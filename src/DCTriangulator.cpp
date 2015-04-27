#include "DCTriangulator.h"

// debug
#ifndef NDEBUG
#include <csignal>
#endif

// std library
#include <fstream>
#include <iostream>
#include <vector>
#include <functional>
#include <set>
#include <stdexcept>
#include <queue>
#include <atomic>
#include <thread>

#include <tbb/parallel_for.h>
#include <tbb/parallel_do.h>
#include <tbb/parallel_sort.h>
#include <tbb/enumerable_thread_specific.h>

// tbb
#include <tbb/spin_mutex.h>

// own
#include "Geometry.h"
#include "CGALTriangulator.h"

#include "utils/Timings.h"
#include "utils/Logger.h"
#include "utils/ASSERT.h"

//**************************

template <uint D, typename Precision>
constexpr Precision DCTriangulator<D, Precision>::SAFETY;

template <uint D, typename Precision>
constexpr Precision DCTriangulator<D, Precision>::ADAPTION_FACTOR;

//**************************

template <uint D, typename Precision>
DCTriangulator<D, Precision>::DCTriangulator(
    const dBox<D, Precision> &_bounds,
    dPoints<D, Precision> &_points,
    const uint _baseThreshold,
    const unsigned char splitter,
    const uint gridOccupancy,
    const bool parallelBaseSolver)
    : Triangulator<D, Precision>(_bounds, _points), baseThreshold(_baseThreshold) {


    // add infinite points to data set
    auto stats = getPointStats(std::size_t(0), this->points.size(), this->points);
    for (uint i = 0; i < pow(2, D); ++i) {
    VLOG("Point stats: " << stats.min << " - " << stats.mid << " - "
         << stats.max << std::endl);

    dPoint<D, Precision> p;
    p.id = dPoint<D, Precision>::cINF | i;
    p.coords = stats.mid;

    for (uint d = 0; d < D; ++d)
    p.coords[d] +=
    (i & (1 << d) ? 1 : -1) * 2 * SAFETY * (stats.max[d] - stats.min[d]);

    this->points.emplace_back(p);
    }

    if(parallelBaseSolver)
        baseTriangulator = std::make_unique<CGALTriangulator<D, Precision, true>>(this->baseBounds, this->points, gridOccupancy);
    else
        baseTriangulator = std::make_unique<CGALTriangulator<D, Precision, false>>(this->baseBounds, this->points, gridOccupancy);

    if(splitter != 0)
        partitioner = Partitioner<D, Precision>::make(splitter);

}

template <uint D, typename Precision>
Ids DCTriangulator<D, Precision>::getEdge(
    const dSimplices<D, Precision> &simplices,
    const Partitioning<D, Precision> &partitioning, const uint &partition) {
  Ids edgeSimplices;
  Ids wqa = simplices.convexHull; // set of already checked simplices
  std::deque<uint> wq(simplices.convexHull.begin(), simplices.convexHull.end());

      /* Walk along the neighbors,
       * testing for each neighbor whether its circumsphere is within the
       * partition or not
       */

      INDENT
      while (!wq.empty()) {
        uint x = wq.front();
        wq.pop_front();

        if (simplices.contains(x)) {
          const auto cs = simplices[x].circumsphere(this->points);
          bool intersectsBounds = false;
          for (uint i = 0; i < partitioning.size(); ++i) {
            if (i != partition && partitioning[i].bounds.intersects(cs))
              intersectsBounds = true;
          }
          if (intersectsBounds) {

            PLOG("Adding " << simplices[x]
                           << " to edge -> circumcircle criterion"
                           << std::endl);
            edgeSimplices.insert(simplices[x].id);

            for (const auto &n : simplices[x].neighbors) {
              if (wqa.insert(n).second) {
                // n was not yet inspected
                wq.push_back(n);
              }
            }
          }
        }
      }
      DEDENT

  return edgeSimplices;
}

template <uint D, typename Precision>
Ids DCTriangulator<D, Precision>::extractPoints(
    const Ids &edgeSimplices, const dSimplices<D, Precision> &simplices,
    bool ignoreInfinite) {
  Ids outPoints;

  for (const auto &id : edgeSimplices) {
    ASSERT(simplices.contains(id));

    for (uint i = 0; i < D + 1; ++i) {
      if (dPoint<D, Precision>::isFinite(simplices[id].vertices[i]))
        outPoints.insert(simplices[id].vertices[i]);
    }
  }

  // add the extreme infinite points to the set
  if (!ignoreInfinite) {
    for (uint k = dPoint<D, Precision>::cINF; k != 0; ++k) {
      outPoints.insert(k);
    }
  }

  return outPoints;
}

template <uint D, typename Precision>
void DCTriangulator<D, Precision>::updateNeighbors(
    dSimplices<D, Precision> &simplices, const Ids &toCheck,
#ifdef NDEBUG
    __attribute__((unused))
#endif
    const std::string &provenance) {

  INDENT
  const uint saveIndent = LOGGER.getIndent();

  tbb::concurrent_unordered_set<uint> wqa(toCheck.size());
  tbb::enumerable_thread_specific<std::set<uint>,
          tbb::cache_aligned_allocator<std::set<uint>>,
          tbb::ets_key_usage_type::ets_key_per_instance> tsNeighborSet;

#ifndef NDEBUG
  std::atomic<uint> checked(0);
  std::atomic<uint> updated(0);
#endif

  tbb::parallel_do(toCheck, [&](const uint &id,
                                tbb::parallel_do_feeder<uint> &feeder) {

    LOGGER.setIndent(saveIndent);

    if (!wqa.insert(id).second) {
      return;
    }

    dSimplex<D, Precision> &simplex = simplices[id];

    PLOG("Checking neighbors of " << simplex << std::endl);
#ifndef NDEBUG
    ++checked;
#endif

    bool requiresUpdate = false;
    for (const auto &n : simplex.neighbors) {
      if (dSimplex<D, Precision>::isFinite(n) &&
          (!simplices.contains(n) || !simplex.isNeighbor(simplices[n]))) {
        requiresUpdate = true;
        break;
      }
    }

    if (!requiresUpdate) {
      INDENT
      PLOG("No update necessary" << std::endl);
      DEDENT

      return;
    }

    PLOG("Updating neighbors of " << simplex << std::endl);
#ifndef NDEBUG
    ++updated;
#endif

#ifndef NDEBUG
    dSimplex<D, Precision> saveSimplex = simplex;
#endif

    INDENT
      uint neighborIdx = 0;
      auto neighborSet = tsNeighborSet.local();
      neighborSet.clear();
    for (uint i = 0; i < D + 1; ++i) {

        uint facetteHash = simplex.vertexFingerprint ^ simplex.vertices[i];
        auto range = simplices.wuFaces.equal_range(facetteHash);
        PLOG("Key: " << facetteHash << " #Results: " << std::distance(range.first, range.second) << std::endl;);
        for (auto it = range.first; it != range.second; ++it) {
            if (it->second != simplex.id && dSimplex<D, Precision>::isFinite(it->second)
                && simplices.contains(it->second) && simplex.isNeighbor(simplices[it->second])) {
                PLOG("Neighbor with " << simplices[it->second] << std::endl);

                if(neighborSet.insert(it->second).second) {
                    simplex.neighbors[neighborIdx++] = it->second;
                    feeder.add(it->second);
                }
            }
        }
    }
    DEDENT

    ASSERT(0 < neighborIdx && neighborIdx <= D + 1);

#ifndef NDEBUG
    if (!(neighborIdx > 0 && neighborIdx <= D + 1)) {
      LOG("Error: wrong number of neighbors for simplex " << simplex
                                                          << std::endl);
    }
#endif

    // ASSERT(simplex.neighbors.size() > 0 && simplex.neighbors.size() <=
    // D+1);

    while(neighborIdx < D + 1) { // if it is a simplex at the border,
      // it might have one infinite
      // simplex as neighbor
      simplex.neighbors[neighborIdx++] = uint(dSimplex<D, Precision>::cINF);
    }
  });

#ifndef NDEBUG
  LOG("Updating neighbors - checked: " << checked << " updated: " << updated
                                       << " simplices: " << simplices.size()
                                       << std::endl);
#endif

  DEDENT
}

template <uint D, typename Precision>
dSimplices<D, Precision> DCTriangulator<D, Precision>::mergeTriangulation(
    std::vector<dSimplices<D, Precision>> &partialDTs, const Ids &edgeSimplices,
    const dSimplices<D, Precision> &edgeDT,
    const Partitioning<D, Precision> &partitioning,
    const std::string &provenance,
#ifdef NDEBUG
    __attribute__((unused))
#endif
    const dSimplices<D, Precision> *realDT) {

  LOG("Merging partial DTs into one triangulation" << std::endl);

  dSimplices<D, Precision> DT;
  // use the first partition as estimator for the size of the triangulation
  DT.reserve(partialDTs.size() * partialDTs[0].size());
  DT.convexHull.reserve(partialDTs.size() * partialDTs[0].convexHull.size());
  DT.wuFaces.reserve(partialDTs.size() * partialDTs[0].wuFaces.size());

  std::vector<dSimplex<D, Precision>> deletedSimplices;
  deletedSimplices.reserve(edgeSimplices.size());

  for (uint i = 0; i < partialDTs.size(); ++i) {

      // copy the simplices not belonging to the edge
      for(auto & s : partialDTs[i]){
          if(!edgeSimplices.count(s.id))
              DT.insert(std::move(s));
          else
              deletedSimplices.push_back(std::move(s));
      }

      // copy the convex hull of the partial DT
      // only copy values not belonging to edgeSimplices
      for(auto & idx : partialDTs[i].convexHull){
          if(!edgeSimplices.count(idx))
              DT.convexHull.insert(std::move(idx));
      }

      // copy the faces where-used list
      // don't copy values belonging to the edge
      for(auto & wu : partialDTs[i].wuFaces){
          if(!edgeSimplices.count(wu.second))
            DT.wuFaces.insert(std::move(wu));
      }
  }

  auto cmpFingerprint =
          [](const dSimplex<D, Precision> &a, const dSimplex<D, Precision> &b) {
              return a.vertexFingerprint < b.vertexFingerprint;
          };
  tbb::parallel_sort(deletedSimplices.begin(), deletedSimplices.end(), cmpFingerprint);

  // merge partial DTs and edge DT
  LOG("Merging triangulations" << std::endl);
  tbb::spin_mutex insertMtx;
  Ids insertedSimplices;

  tbb::parallel_for(std::size_t(0), edgeDT.bucket_count(), [&](const uint i) {

    for (auto it = edgeDT.begin(i); it != edgeDT.end(i); ++it) {

      const dSimplex<D, Precision> &edgeSimplex = *it;

      // check whether edgeSimplex is completely contained in one partition
      uint p0 = partitioning.partition(edgeSimplex.vertices[0]);
      bool inOnePartition = partitioning[p0].contains(edgeSimplex);

      if (!inOnePartition) {
        tbb::spin_mutex::scoped_lock lock(insertMtx);
        DT.insert(edgeSimplex);

          //convex hull treatment
          if(!edgeSimplex.isFinite())
              DT.convexHull.insert(edgeSimplex.id);

        for (uint d = 0; d < D + 1; ++d) {
          DT.wuFaces.emplace((edgeSimplex.vertexFingerprint ^ edgeSimplex.vertices[d]), edgeSimplex.id);
        }

        insertedSimplices.insert(edgeSimplex.id);
      } else {
        // the simplex is completely in one partition -> it must have been found
        // before
        auto range =
            std::equal_range(deletedSimplices.begin(), deletedSimplices.end(),
                             edgeSimplex, cmpFingerprint);

        for (auto it = range.first; it != range.second; ++it) {
          if (edgeSimplex.equalVertices(*it)) {
            tbb::spin_mutex::scoped_lock lock(insertMtx);
            DT.insert(edgeSimplex);

              //convex hull treatment
              if(!edgeSimplex.isFinite())
                  DT.convexHull.insert(edgeSimplex.id);

            for (uint d = 0; d < D + 1; ++d) {
                DT.wuFaces.emplace((edgeSimplex.vertexFingerprint ^ edgeSimplex.vertices[d]), edgeSimplex.id);
            }

            insertedSimplices.insert(edgeSimplex.id);
          }
        }
      }
    }
  });

  ASSERT(DT.countDuplicates() == 0);

  LOG("Updating neighbors" << std::endl);
  updateNeighbors(DT, insertedSimplices, provenance);

  return DT;
}

template <uint D, typename Precision>
void DCTriangulator<D, Precision>::evaluateVerificationReport(
    __attribute__((unused)) const VerificationReport<D, Precision> &vr,
    __attribute__((unused)) const std::string &provenance) const {

}

template <uint D, typename Precision>
void DCTriangulator<D, Precision>::evaluateCrossCheckReport(
    const CrossCheckReport<D, Precision> &ccr, const std::string &provenance,
    const dSimplices<D, Precision> &DT, const dSimplices<D, Precision> &realDT,
    const Partitioning<D, Precision> *partitioning) const {
  if (IS_PROLIX) {
    std::stringstream txt;

    for (const auto &missingSimplex : ccr.missing) {

      auto mySimplices = DT.findSimplices(missingSimplex.vertices);

      // textual description of the missing simplex

      // sort the simplices by descending number of shared vertices
      std::priority_queue<std::pair<uint, uint>> pqSharedVertices;
      for (const auto &s : mySimplices) {
        pqSharedVertices.emplace(s.countSharedVertices(missingSimplex), s.id);
      }

      auto format = [&](const uint vertex) -> std::string {
        std::stringstream s;
        if (!missingSimplex.contains(vertex))
          s << zkr::cc::fore::red;

        s << vertex;

        if (partitioning != nullptr) {
          s << "-" << partitioning->partition(vertex);
        }

        s << zkr::cc::console;

        return s.str();
      };

      txt << "Missing simplex: " << missingSimplex.id << ": V = ["
          << format(missingSimplex.vertices[0]);
      for (uint i = 1; i < D + 1; ++i)
        txt << ", " << format(missingSimplex.vertices[i]);
      txt << "]" << std::endl;

      txt << "Points:" << std::endl;
      for (const auto &p : missingSimplex.vertices) {
        txt << '\t' << this->points[p] << std::endl;
      }

      txt << "Cirumsphere: " << missingSimplex.circumsphere(this->points)
          << std::endl;

      txt << "Possible candidates:" << std::endl;

      while (!pqSharedVertices.empty()) {
        // check whether we have a valid simplex
        const auto &s = mySimplices[pqSharedVertices.top().second];
        auto realS = realDT.findSimplices(s.vertices, true);

        txt << "\t";

        txt << s.id << ": V = [" << format(s.vertices[0]);
        for (uint i = 1; i < D + 1; ++i)
          txt << ", " << format(s.vertices[i]);
        txt << "] - real";

        ASSERT(realS.empty() || realS.size() == 1);

        if (realS.empty()) {
          txt << " INVALID";
        } else {
          for (const auto &r : realS) {
            txt << " " << r;
          }
        }
        txt << std::endl;

        pqSharedVertices.pop();
      }
      txt << std::endl;
    }

    for (const auto &invalidSimplex : ccr.invalid) {
      auto realSimplices = realDT.findSimplices(invalidSimplex.vertices);

      // textual description of the missing simplex
      txt << "Invalid simplex: " << invalidSimplex << std::endl;
      txt << "Possible real candidates:" << std::endl;
      for (const auto &s : realSimplices) {
        txt << "\t" << s << std::endl;
      }
      txt << std::endl;
    }

    if (!txt.str().empty()) {
      std::ofstream o("log/" + provenance + "_CrossCheckReport.log",
                      std::ios::out | std::ios::trunc);
      o << txt.str();
    }
  }
}

template <uint D, typename Precision>
dSimplices<D, Precision>
DCTriangulator<D, Precision>::_triangulateBase(const Ids partitionPoints,
                                            const dBox<D, Precision> &bounds,
                                            const std::string provenance) {

  LOGGER.setIndent(provenance.length());

  LOG("triangulateBASE called on level " << provenance << " with "
                                         << partitionPoints.size() << " points"
                                         << std::endl);

  bool isTOP = this->isTOP(provenance);

  std::unique_ptr<dSimplices<D, Precision>> realDT = nullptr;
  if (isTOP && Triangulator<D, Precision>::VERIFY) {
    LOG("Real triangulation" << std::endl);
    INDENT
    realDT = std::make_unique<dSimplices<D, Precision>>(
        baseTriangulator->_triangulate(partitionPoints, bounds, provenance));
    LOG("Real triangulation contains " << realDT->size() << " tetrahedra"
                                       << std::endl);
    DEDENT
  }

  INDENT
  // if this is the top-most triangulation, ignore infinite vertices
  auto dt = baseTriangulator->_triangulate(partitionPoints, bounds, provenance);
  LOG("Triangulation contains " << dt.size() << " tetrahedra" << std::endl);
  DEDENT

  if (isTOP && Triangulator<D, Precision>::VERIFY) {
    LOG("Verifying CGAL triangulation" << std::endl);
    dt.verify(partitionPoints, this->points);
  }

#ifndef NDEBUG
  auto saveDT = dt;
#endif

  // LOG("Updating neighbors" << std::endl);
  // findNeighbors(dt, provenance);

  ASSERT(saveDT == dt); // only performed if not NDEBUG

  TriangulationReportEntry rep;
  rep.provenance = provenance;
  rep.base_case = true;
  rep.valid = false;
  rep.edge_triangulation = provenance.find('e') != std::string::npos;
  rep.nPoints = partitionPoints.size();
  rep.nSimplices = dt.size();
  rep.nEdgePoints = 0;
  rep.nEdgeSimplices = 0;

  if (isTOP && Triangulator<D, Precision>::VERIFY) {
    LOG("Consistency check of triangulation" << std::endl);
    auto vr = dt.verify(partitionPoints, this->points);
    evaluateVerificationReport(vr, provenance);

    LOG("Cross check with real triangulation" << std::endl);
    auto ccr = dt.crossCheck(*realDT);
    evaluateCrossCheckReport(ccr, provenance, dt, *realDT);

    rep.valid = vr.valid && ccr.valid;
  }

  triangulationReport.push_back(rep);

  return dt;
}

template <uint D, typename Precision>
dSimplices<D, Precision>
DCTriangulator<D, Precision>::_triangulate(const Ids & partitionPoints,
                                          const dBox<D, Precision> &bounds,
                                          const std::string provenance) {

  LOGGER.setIndent(provenance.length());

  LOG("triangulateDAC called on level " << provenance << " with "
                                        << partitionPoints.size() << " points"
                                        << std::endl);

  bool isTOP = this->isTOP(provenance);

  if (partitionPoints.size() > baseThreshold) {
    LOG("Recursive case" << std::endl);

    std::unique_ptr<dSimplices<D, Precision>> realDT = nullptr;
    if (isTOP && Triangulator<D, Precision>::VERIFY) {
      // perform real triangulation
      LOG("Real triangulation" << std::endl);
      INDENT

      realDT = std::make_unique<dSimplices<D, Precision>>(
              baseTriangulator->_triangulate(partitionPoints, bounds, provenance));
      LOG("Real triangulation contains " << realDT->size() << " tetrahedra"
                                         << std::endl);
      DEDENT
    }

    // partition input
    LOG("Partioning" << std::endl);
    INDENT
    const auto partioning =
        partitioner->partition(partitionPoints, this->points, provenance);
    DEDENT

    std::vector<dSimplices<D, Precision>> partialDTs;
    partialDTs.resize(partioning.size());

    tbb::parallel_for(
        std::size_t(0), partioning.size(),
        [&](const uint i)
        {
          LOGGER.setIndent(
              provenance.length()); // new thread, initialize Logger indent
          INDENT
          LOG("Partition " << i << " on thread " << std::this_thread::get_id()
                           << std::endl);
          partialDTs[i] =
              _triangulate(partioning[i].points, partioning[i].bounds,
                          provenance + std::to_string(i));
          LOG("Triangulation " << i << " contains " << partialDTs[i].size()
                               << " tetrahedra" << std::endl);
          DEDENT
        }
        );

    LOG("Extracting edges" << std::endl);
    INDENT

    Ids edgePointIds;
    Ids edgeSimplexIds;

    for (uint i = 0; i < partioning.size(); ++i) {
      // points are in different partitions, there can be no overlap

      auto edge = getEdge(partialDTs[i], partioning, i);
      edgeSimplexIds.insert(edge.begin(), edge.end());

      // ignore infinite vertices if this is the top most triangulation
      auto ep = extractPoints(edge, partialDTs[i]);
      edgePointIds.insert(ep.begin(), ep.end());
    }
    DEDENT

    LOG("Edge has " << edgeSimplexIds.size() << " simplices with "
                    << edgePointIds.size() << " points" << std::endl);

      dSimplices<D, Precision> edgeDT;

      if(edgePointIds.size() < ADAPTION_FACTOR * partitionPoints.size()) {
          LOG("Triangulating edges recursively" << std::endl);
          INDENT
          edgeDT = _triangulate(edgePointIds, bounds, provenance + "e");
      } else {
          LOG("Triangulating edges with basecase" << std::endl);
          INDENT
          edgeDT = _triangulateBase(edgePointIds, bounds, provenance + "e");
      }
    LOG("Edge triangulation contains " << edgeDT.size() << " tetrahedra"
                                       << std::endl
                                       << std::endl);
    DEDENT

    auto mergedDT = mergeTriangulation(partialDTs, edgeSimplexIds, edgeDT,
                                       partioning, provenance, realDT.get());

    TriangulationReportEntry rep;
    rep.provenance = provenance;
    rep.base_case = false;
    rep.valid = false;
    rep.edge_triangulation = provenance.find('e') != std::string::npos;
    rep.nPoints = partitionPoints.size();
    rep.nSimplices = mergedDT.size();
    rep.nEdgePoints = edgePointIds.size();
    rep.nEdgeSimplices = edgeDT.size();

    if (isTOP & Triangulator<D, Precision>::VERIFY) {
      LOG("Consistency check of triangulation" << std::endl);
      auto vr = mergedDT.verify(partitionPoints, this->points);
      evaluateVerificationReport(vr, provenance);

      LOG("Cross check with real triangulation" << std::endl);
      auto ccr = mergedDT.crossCheck(*realDT);
      evaluateCrossCheckReport(ccr, provenance, mergedDT, *realDT, &partioning);

      rep.valid = vr.valid && ccr.valid;
    }

    triangulationReport.push_back(rep);

    return mergedDT;

  } else { // base case
    return _triangulateBase(partitionPoints, bounds, provenance);
  }
}

// specializations

template class DCTriangulator<2, float>;
template class DCTriangulator<3, float>;

template class DCTriangulator<2, double>;
template class DCTriangulator<3, double>;

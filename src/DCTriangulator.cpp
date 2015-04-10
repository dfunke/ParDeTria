#include "DCTriangulator.h"

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

// debug
#ifndef NDEBUG
#include <csignal>
#endif

// own
#include "Geometry.h"
#include "Painter.h"
#include "CGALTriangulator.h"

#include "utils/Timings.h"
#include "utils/Logger.h"
#include "utils/ASSERT.h"

//**************************

template <uint D, typename Precision>
constexpr Precision DCTriangulator<D, Precision>::SAFETY;

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
  Ids wqa; // set of already checked simplices

  // the infinite points are stored as the last ones in the point array
  for (uint infVertex = this->points.size() - dPoint<D, Precision>::nINF;
       infVertex < this->points.size(); ++infVertex) {

    ASSERT(this->points.contains(infVertex));
    ASSERT(!this->points[infVertex].isFinite());

    PLOG("Infinite vertex " << this->points[infVertex] << std::endl);

    uint infId = this->points[infVertex].id;
    LOGGER.logContainer(simplices.whereUsed.at(infId),
                        Logger::Verbosity::PROLIX, "Used in simplices: ");

    INDENT
    for (const auto &s : simplices.whereUsed.at(infId)) {
      if (!dSimplex<D, Precision>::isFinite(s))
        continue;

      ASSERT(simplices.contains(s));

      PLOG("Adding " << simplices[s] << " to edge" << std::endl);
      edgeSimplices.insert(simplices[s].id);

      /* Walk along the neighbors,
       * testing for each neighbor whether its circumsphere is within the
       * partition or not
       */
      wqa.insert(simplices[s].id);
      std::deque<uint> wq;

      // work queue of simplices to check for circum circle criterian
      for (const auto &n : simplices[s].neighbors) {
        if (wqa.insert(n).second) {
          // simplex not yet inspected -> add it to queue
          wq.push_back(n);
        }
      }

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
    }
    DEDENT
  }

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
void DCTriangulator<D, Precision>::findNeighbors(
    dSimplices<D, Precision> &simplices,
#ifdef NDEBUG
    __attribute__((unused))
#endif
    const std::string &provenance) {
#ifndef NDEBUG
  PainterBulkWriter<D, Precision> paintWriter;
#endif

  INDENT
  const uint saveIndent = LOGGER.getIndent();

  tbb::parallel_for(
      std::size_t(0), simplices.bucket_count(), [&](const uint i) {

        LOGGER.setIndent(saveIndent);

        for (auto it = simplices.begin(i); it != simplices.end(i); ++it) {

          dSimplex<D, Precision> &simplex = *it;

          PLOG("Updating neighbors of " << simplex << std::endl);

#ifndef NDEBUG
          dSimplex<D, Precision> saveSimplex = simplex;
#endif

          simplex.neighbors.clear();
          simplex.neighbors.reserve(D + 1);

          std::unordered_map<uint, uint> counters;
          counters.reserve((D + 1) * (D + 1) * (D + 1));

          INDENT
          for (uint v = 0; v < D + 1; ++v) {
            // for every point, look where else its used
            // if the other simplex shares a second point -> it is a neighbor

            const dPoint<D, Precision> &vertex = this->points[simplex.vertices[v]];

            if (IS_PROLIX) {
              LOGGER.logContainer(simplices.whereUsed.at(vertex.id),
                                  Logger::Verbosity::PROLIX,
                                  "Vertex " + to_string(vertex) + " used in");
            }

            INDENT
            for (const uint u : simplices.whereUsed.at(vertex.id)) {
              if (dSimplex<D, Precision>::isFinite(u))
                counters[u] += 1;
            }

            for (const auto &it : counters) {
              if (it.first != simplex.id)
                if (it.second == D)
                  if (simplices.contains(it.first)) {
                    PLOG("Neighbor with " << simplices[it.first] << std::endl);

                    simplex.neighbors.insert(it.first);

                    // LOGGER.logContainer(simplex.neighbors,
                    // Logger::Verbosity::PROLIX);

                    // ASSERT(simplex.neighbors.size() <= D+1);
                    // u will be updated in its own round;
                  }
            }
            DEDENT
          }
          DEDENT

          ASSERT(0 < simplex.neighbors.size() &&
                 simplex.neighbors.size() <= D + 1);

#ifndef NDEBUG
          if (!(simplex.neighbors.size() > 0 &&
                simplex.neighbors.size() <= D + 1)) {
            LOG("Error: wrong number of neighbors for simplex " << simplex
                                                                << std::endl);

            if (IS_PROLIX) {
              paintWriter.add("img/" + provenance + "_neighbors_" +
                                  std::to_string(simplex.id),
                              this->baseBounds);
              paintWriter.top().draw(this->points);
              paintWriter.top().drawPartition(this->points);

              paintWriter.top().setColor(0, 0, 0, 0.4);
              paintWriter.top().draw(simplices, this->points, true);

              paintWriter.top().setColor(1, 0, 0); // simplex in red
              paintWriter.top().draw(simplex, this->points, true);

              paintWriter.top().setColor(1, 1, 0, 0.4); // neighbors in yellow
              paintWriter.top().drawNeighbors(simplex, simplices, this->points, true);

              if (!(simplex.equalVertices(saveSimplex) &&
                    simplex.equalNeighbors(saveSimplex))) {
                LOG("Error: was before " << saveSimplex << std::endl);

                paintWriter.top().setColor(0, 0, 1); // old simplex in blue
                paintWriter.top().draw(saveSimplex, this->points, true);

                paintWriter.top().setColor(
                    0, 1, 1,
                    0.4); // old simplex neighbors in cyan
                paintWriter.top().drawNeighbors(saveSimplex, simplices, this->points,
                                                true);
              }
            }
          }
#endif

          // ASSERT(simplex.neighbors.size() > 0 && simplex.neighbors.size() <=
          // D+1);

          // TODO it might be better NOT to save the infinite simplex as
          // neighbor
          if (simplex.neighbors.size() <
              D + 1) { // if it is a simplex at the border,
            // it might have one infinite
            // simplex as neighbor
            simplex.neighbors.insert(uint(dSimplex<D, Precision>::cINF));
          }
        }
      });
  DEDENT
}

// supplementary class used in updateNeighbors
template <uint D, typename Precision> class tWhereUsedSuppl {
public:
  tWhereUsedSuppl() : whereUsed(nullptr), index(0), size(0), deleted(0) {}
  tWhereUsedSuppl(tbb::concurrent_vector<uint> *_whereUsed)
      : whereUsed(_whereUsed), index(0), size(_whereUsed->size()), deleted(0) {}

public:
  tbb::concurrent_vector<uint> *whereUsed = nullptr;
  uint index = 0;
  uint size = 0;
  uint deleted = 0;
};

template <uint D, typename Precision>
void DCTriangulator<D, Precision>::updateNeighbors(
    dSimplices<D, Precision> &simplices, const Ids &toCheck,
#ifdef NDEBUG
    __attribute__((unused))
#endif
    const std::string &provenance) {
#ifndef NDEBUG
  PainterBulkWriter<D, Precision> paintWriter;
#endif

  INDENT
  const uint saveIndent = LOGGER.getIndent();

  tbb::concurrent_unordered_set<uint> wqa(toCheck.size());

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

    simplex.neighbors.clear();
    simplex.neighbors.reserve(D + 1);

    std::array<tWhereUsedSuppl<D, Precision>, D + 1> suppl;
    for (uint i = 0; i < D + 1; ++i) {
      suppl[i] = std::move(tWhereUsedSuppl<D, Precision>(
          &simplices.whereUsed.at(simplex.vertices[i])));
    }

    uint currentMax = 0;
    uint currentV = 0;
    bool cont = true;

    INDENT
    while (cont) {
      uint currentMin = std::numeric_limits<uint>::max();
      uint countOver = 0;
      uint count = 0;
      uint v = 0;

      for (uint i = 0; i < D + 1; ++i) {

        tWhereUsedSuppl<D, Precision> &supp = suppl[i];

        while (supp.index < supp.size ? (v = supp.whereUsed->at(supp.index),
                                         !dSimplex<D, Precision>::isFinite(v))
                                      : (++countOver, false)) {
          ++supp.index;
          ++supp.deleted;
        }

        if (supp.index < supp.size) {
          currentMin = std::min(v, currentMin);
          currentMax = std::max(v, currentMax);
          count += v == currentV;

          supp.index += v <= currentV;
        }
      }
      if (count == D && simplices.contains(currentV)) {
        PLOG("Neighbor with " << simplices.at(currentV) << std::endl);

        simplex.neighbors.insert(currentV);
        feeder.add(currentV);
      }

      if (currentV == currentMin)
        ++currentV;
      else
        currentV = currentMin;

      cont = countOver <= 2;
    }

    // compactify where used list if necessary
    for (uint i = 0; i < D + 1; ++i) {
      tWhereUsedSuppl<D, Precision> &supp = suppl[i];

      if (supp.deleted / supp.size > .5) {
        tbb::concurrent_vector<uint> v;
        v.reserve(supp.size - supp.deleted);

        for (const auto &s : *supp.whereUsed) {
          if (dSimplex<D, Precision>::isFinite(s))
            v.emplace_back(s);
        }

        *supp.whereUsed = std::move(v);
      }
    }
    DEDENT

    ASSERT(0 < simplex.neighbors.size() && simplex.neighbors.size() <= D + 1);

#ifndef NDEBUG
    if (!(simplex.neighbors.size() > 0 && simplex.neighbors.size() <= D + 1)) {
      LOG("Error: wrong number of neighbors for simplex " << simplex
                                                          << std::endl);

      if (IS_PROLIX) {
        paintWriter.add("img/" + provenance + "_neighbors_" +
                            std::to_string(simplex.id),
                        this->baseBounds);
        paintWriter.top().draw(this->points);
        paintWriter.top().drawPartition(this->points);

        paintWriter.top().setColor(0, 0, 0, 0.4);
        paintWriter.top().draw(simplices, this->points, true);

        paintWriter.top().setColor(1, 0, 0); // simplex in red
        paintWriter.top().draw(simplex, this->points, true);

        paintWriter.top().setColor(1, 1, 0, 0.4); // neighbors in yellow
        paintWriter.top().drawNeighbors(simplex, simplices, this->points, true);

        if (!(simplex.equalVertices(saveSimplex) &&
              simplex.equalNeighbors(saveSimplex))) {
          LOG("Error: was before " << saveSimplex << std::endl);

          paintWriter.top().setColor(0, 0, 1); // old simplex in blue
          paintWriter.top().draw(saveSimplex, this->points, true);

          paintWriter.top().setColor(0, 1, 1,
                                     0.4); // old simplex neighbors in cyan
          paintWriter.top().drawNeighbors(saveSimplex, simplices, this->points, true);
        }
      }
    }
#endif

    // ASSERT(simplex.neighbors.size() > 0 && simplex.neighbors.size() <=
    // D+1);

    // TODO it might be better NOT to save the infinite simplex as neighbor
    if (simplex.neighbors.size() < D + 1) { // if it is a simplex at the border,
      // it might have one infinite
      // simplex as neighbor
      simplex.neighbors.insert(uint(dSimplex<D, Precision>::cINF));
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

  for (uint i = 0; i < partialDTs.size(); ++i) {
    DT.insert(partialDTs[i].begin(), partialDTs[i].end());
    for (const auto &wu : partialDTs[i].whereUsed) {
      DT.whereUsed[wu.first].grow_by(wu.second.begin(), wu.second.end());
    }
  }

  auto edgePointIds = extractPoints(edgeSimplices, DT);

//**********************
#ifndef NDEBUG
  auto paintMerging = [&](const dSimplices<D, Precision> &dt,
                          const std::string &name,
                          bool drawInfinite) -> Painter<D, Precision> {
    Painter<D, Precision> painter(this->baseBounds);

    painter.draw(this->points);
    painter.drawPartition(this->points);

    painter.setColor(0, 1, 0, 0.4);
    painter.draw(dt, this->points, drawInfinite);

    painter.setColor(1, 0, 0);
    painter.draw(this->points.template filter<dPoints<D, Precision>>(edgePointIds));
    painter.save(name);

    if (realDT != nullptr) {
      painter.setColor(0, 0, 0, 0.1);
      painter.draw(*realDT, this->points);
    }

    painter.save(name + "_overlay");

    return painter;
  };
#endif
//**********************

#ifndef NDEBUG
  paintMerging(DT, provenance + "_05a_merging_merged", false);
#endif

  auto cmpFingerprint =
      [](const dSimplex<D, Precision> &a, const dSimplex<D, Precision> &b) {
        return a.vertexFingerprint < b.vertexFingerprint;
      };

  std::vector<dSimplex<D, Precision>> deletedSimplices;
  deletedSimplices.reserve(edgeSimplices.size());
  for (const auto &id : edgeSimplices) {
    deletedSimplices.emplace_back(DT[id]);
  }
  tbb::parallel_sort(deletedSimplices.begin(), deletedSimplices.end(),
                     cmpFingerprint);

#ifndef NDEBUG
// paintMerging(deletedSimplices, provenance + "_05b_merging_deleted", true);
#endif

  // delete all simplices belonging to the edge from DT
  LOG("Striping triangulation from edge" << std::endl);

  tbb::spin_mutex eraseMtx;
  tbb::parallel_for(
      std::size_t(0), edgeSimplices.bucket_count(), [&](const uint i) {

        for (auto it = edgeSimplices.begin(i); it != edgeSimplices.end(i);
             ++it) {

          const uint id = *it;

          ASSERT(DT.contains(id));

          // remove simplex from where used list
          for (uint p = 0; p < D + 1; ++p) {

            auto it = std::find(DT.whereUsed[DT[id].vertices[p]].begin(),
                                DT.whereUsed[DT[id].vertices[p]].end(), id);

            ASSERT(it != DT.whereUsed[DT[id].vertices[p]].end());

            // lock.upgrade_to_writer();
            // no need to upgrade to writer here - no invalidation of iterators
            // possible
            *it = dSimplex<D, Precision>::cINF;
          }

          tbb::spin_mutex::scoped_lock lock(eraseMtx);
          DT.erase(id);
        }
      });

//**********************
#ifndef NDEBUG
  auto painter = paintMerging(DT, provenance + "_05b_merging_stripped", false);

  painter.setColor(0, 0, 1, 0.4);
  painter.setDashed();
  painter.draw(edgeDT, this->points);

  painter.save(provenance + "_05b_merging_stripped+edge_overlay");

  painter.setColor(1, 0, 0, 0.4);
  painter.setDashed();
  // painter.draw(deletedSimplices, this->points);
  painter.setDashed(false);

  painter.save(provenance + "_05b_merging_stripped+edge+deleted_overlay");
#endif
  //**********************

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

        for (uint d = 0; d < D + 1; ++d) {
          DT.whereUsed[edgeSimplex.vertices[d]].emplace_back(edgeSimplex.id);
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

            for (uint d = 0; d < D + 1; ++d) {
              DT.whereUsed[edgeSimplex.vertices[d]].emplace_back(
                  edgeSimplex.id);
            }

            insertedSimplices.insert(edgeSimplex.id);
          }
        }
      }
    }
  });

#ifndef NDEBUG
  paintMerging(DT, provenance + "_05c_merging_edge", false);
#endif

  ASSERT(DT.countDuplicates() == 0);

  // sort the where-used list
  tbb::parallel_for(DT.whereUsed.range(), [&](auto &r) {
    for (auto &it : r)
      tbb::parallel_sort(it.second.begin(), it.second.end());
  });

  LOG("Updating neighbors" << std::endl);
  updateNeighbors(DT, insertedSimplices, provenance);

#ifndef NDEBUG
  paintMerging(DT, provenance + "_05d_merging_finished", false);
#endif

  return DT;
}

template <uint D, typename Precision>
void DCTriangulator<D, Precision>::evaluateVerificationReport(
    const VerificationReport<D, Precision> &vr,
    const std::string &provenance) const {
  if (IS_PROLIX) {
    Painter<D, Precision> basePainter(this->baseBounds);
    basePainter.draw(this->points);
    basePainter.drawPartition(this->points);
    basePainter.setLogging();

    PainterBulkWriter<D, Precision> writer;
    for (const auto &inCircle : vr.inCircle) {
      writer.add("img/" + provenance + "_inCircle_" +
                     std::to_string(inCircle.first.id),
                 basePainter);

      writer.top().setColor(0, 1, 0);
      writer.top().draw(inCircle.first, this->points, true);
      writer.top().drawCircumSphere(inCircle.first, this->points, true);

      writer.top().setColor(1, 0, 0);
      writer.top().draw(
          this->points.template filter<dPoints<D, Precision>>(inCircle.second));
    }

    for (const auto &neighbors : vr.wrongNeighbors) {
      writer.add("img/" + provenance + "_wrongNeighbors_" +
                     std::to_string(neighbors.first.id) + "_" +
                     std::to_string(neighbors.second.id),
                 basePainter);

      writer.top().setColor(0, 1, 0);
      writer.top().draw(neighbors.first, this->points, true);

      writer.top().setColor(1, 0, 0);
      writer.top().draw(neighbors.second, this->points, true);
    }
  }
}

template <uint D, typename Precision>
void DCTriangulator<D, Precision>::evaluateCrossCheckReport(
    const CrossCheckReport<D, Precision> &ccr, const std::string &provenance,
    const dSimplices<D, Precision> &DT, const dSimplices<D, Precision> &realDT,
    const Partitioning<D, Precision> *partitioning) const {
  if (IS_PROLIX) {
    Painter<D, Precision> basePainter(this->baseBounds);
    basePainter.draw(this->points);
    basePainter.drawPartition(this->points);
    basePainter.setLogging();

    PainterBulkWriter<D, Precision> writer;
    std::stringstream txt;

    for (const auto &missingSimplex : ccr.missing) {
      writer.add("img/" + provenance + "_missing_" +
                     std::to_string(missingSimplex.id),
                 basePainter);

      auto mySimplices = DT.findSimplices(missingSimplex.vertices);

      writer.top().setColor(0, 1, 0);
      writer.top().draw(mySimplices, this->points, true);

      writer.top().setColor(1, 0, 0);
      writer.top().draw(missingSimplex, this->points, true);

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
      writer.add("img/" + provenance + "_invalid_" +
                     std::to_string(invalidSimplex.id),
                 basePainter);

      auto realSimplices = realDT.findSimplices(invalidSimplex.vertices);

      writer.top().setColor(0, 1, 0);
      writer.top().draw(realSimplices, this->points, true);

      writer.top().setColor(1, 0, 0);
      writer.top().draw(invalidSimplex, this->points, true);

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

#ifndef NDEBUG
    // setup base painter
    Painter<D, Precision> basePainter(this->baseBounds);
    basePainter.draw(
        this->points.template filter<dPoints<D, Precision>>(partitionPoints));
    basePainter.drawPartition(
        this->points.template filter<dPoints<D, Precision>>(partitionPoints));
    basePainter.save(provenance + "_00_points");
#endif

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

#ifndef NDEBUG
      Painter<D, Precision> paintRealDT = basePainter;

      paintRealDT.setColor(0, 0, 1);
      paintRealDT.draw(*realDT, this->points);
      paintRealDT.setColor(0, 0, 0);

      paintRealDT.save(provenance + "_01_realDT");
#endif
    }

    // partition input
    LOG("Partioning" << std::endl);
    INDENT
    const auto partioning =
        partitioner->partition(partitionPoints, this->points, provenance);
    DEDENT

    std::vector<dSimplices<D, Precision>> partialDTs;
    partialDTs.resize(partioning.size());

#ifdef NDEBUG
    tbb::parallel_for(
        std::size_t(0), partioning.size(),
        [&](const uint i)
#else
    for (uint i = 0; i < partioning.size(); ++i)
#endif
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
#ifdef NDEBUG
        );
#endif

#ifndef NDEBUG
    Painter<D, Precision> paintPartialDTs = basePainter;
    for (uint i = 0; i < partioning.size(); ++i) {
      paintPartialDTs.setColor(Painter<D, Precision>::tetradicColor(i), 0.4);
      paintPartialDTs.draw(partialDTs[i], this->points);

      Painter<D, Precision> paintPartialDT = basePainter;
      paintPartialDT.draw(partialDTs[i], this->points, true);
      paintPartialDT.save(provenance + "_02_partialDT_" + std::to_string(i));
    }

    paintPartialDTs.save(provenance + "_02_partialDTs");
#endif

    LOG("Extracting edges" << std::endl);
    INDENT

    Ids edgePointIds;
    Ids edgeSimplexIds;

#ifndef NDEBUG
    Painter<D, Precision> paintEdges = paintPartialDTs;
    paintEdges.setColor(1, 0, 0);
#endif

    for (uint i = 0; i < partioning.size(); ++i) {
      // points are in different partitions, there can be no overlap

      auto edge = getEdge(partialDTs[i], partioning, i);
      edgeSimplexIds.insert(edge.begin(), edge.end());

      // ignore infinite vertices if this is the top most triangulation
      auto ep = extractPoints(edge, partialDTs[i]);
      edgePointIds.insert(ep.begin(), ep.end());

#ifndef NDEBUG
      paintEdges.draw(partialDTs[i].project(edge), this->points, false);
#endif
    }
    DEDENT

    LOG("Edge has " << edgeSimplexIds.size() << " simplices with "
                    << edgePointIds.size() << " points" << std::endl);

#ifndef NDEBUG
    paintEdges.draw(
        this->points.template filter<dPoints<D, Precision>>(edgePointIds));
    paintEdges.save(provenance + "_03_edgeMarked");
#endif

    LOG("Triangulating edges" << std::endl);
    INDENT
    auto edgeDT = _triangulate(edgePointIds, bounds, provenance + "e");
    LOG("Edge triangulation contains " << edgeDT.size() << " tetrahedra"
                                       << std::endl
                                       << std::endl);
    DEDENT

#ifndef NDEBUG
    Painter<D, Precision> paintEdgeDT = basePainter;

    paintEdgeDT.setColor(0, 1, 0);
    paintEdgeDT.draw(edgeDT, this->points, true);
    paintEdgeDT.setColor(1, 0, 0);
    paintEdgeDT.draw(
        this->points.template filter<dPoints<D, Precision>>(edgePointIds));

    paintEdgeDT.save(provenance + "_04_edgeDT");
#endif

    auto mergedDT = mergeTriangulation(partialDTs, edgeSimplexIds, edgeDT,
                                       partioning, provenance, realDT.get());

#ifndef NDEBUG
    Painter<D, Precision> paintFinal = basePainter;
    paintFinal.setColor(0, 1, 0);
    paintFinal.draw(mergedDT, this->points, true);
    paintFinal.save(provenance + "_06_final");
#endif

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

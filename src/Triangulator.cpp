#include "Triangulator.h"

// std library
#include <fstream>
#include <iostream>
#include <vector>
#include <functional>
#include <set>
#include <stdexcept>
#include <queue>
#include <thread>

#include <tbb/parallel_for.h>

// debug
#ifndef NDEBUG
#include <csignal>
#endif

// own
#include "Geometry.h"
#include "Painter.h"
#include "CGAL_Interface.h"

#include "utils/Timings.h"
#include "utils/Logger.h"
#include "utils/ASSERT.h"

//**************************

template <uint D, typename Precision>
constexpr Precision Triangulator<D, Precision>::SAFETY;

template <uint D, typename Precision>
constexpr char Triangulator<D, Precision>::TOP;

template <uint D, typename Precision>
bool Triangulator<D, Precision>::VERIFY = true;
//**************************

template <uint D, typename Precision>
Triangulator<D, Precision>::Triangulator(
    const dBox<D, Precision> &_bounds, const uint _baseThreshold,
    dPoints<D, Precision> &_points,
    std::unique_ptr<Partitioner<D, Precision>> &&_partitioner)
    : bounds(_bounds), baseThreshold(_baseThreshold), points(_points),
      partitioner(_partitioner) {
  if (!points.contains(dPoint<D, Precision>::cINF)) {
    auto stats = getPointStats(points.begin_keys(), points.end_keys(), points);
    for (uint i = 0; i < pow(2, D); ++i) {
      VLOG("Point stats: " << stats.min << " - " << stats.mid << " - "
                           << stats.max << std::endl);

      dPoint<D, Precision> p;
      p.id = dPoint<D, Precision>::cINF | i;
      p.coords = stats.mid;

      for (uint d = 0; d < D; ++d)
        p.coords[d] += (i & (1 << d) ? 1 : -1) * 2 * SAFETY *
                       (stats.max[d] - stats.min[d]);

      points.insert(p);
    }
  }
}

template <uint D, typename Precision>
Ids Triangulator<D, Precision>::getEdge(
    const dSimplices<D, Precision> &simplices,
    const Partitioning<D, Precision> &partitioning, const uint &partition) {
  Ids edgeSimplices;
  std::set<uint> wqa; // set of already checked simplices

  // we use the overflow of the uint to zero to abort the loop
  for (uint infVertex = dPoint<D, Precision>::cINF; infVertex != 0;
       ++infVertex) {

    ASSERT(points.contains(infVertex));

    PLOG("Infinite vertex " << points[infVertex] << std::endl);
    LOGGER.logContainer(points[infVertex].simplices, Logger::Verbosity::PROLIX,
                        "Used in simplices: ");

    INDENT
    for (const auto &s : points[infVertex].simplices) {
      if (!simplices.contains(s))
        continue;

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
          const auto cs = simplices[x].circumsphere(points);
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
Ids Triangulator<D, Precision>::extractPoints(
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
void
Triangulator<D, Precision>::updateNeighbors(dSimplices<D, Precision> &simplices,
#ifdef NDEBUG
                                            __attribute__((unused))
#endif
                                            const std::string &provenance) {
#ifndef NDEBUG
  PainterBulkWriter<D, Precision> paintWriter;
#endif

  INDENT
  const uint saveIndent = LOGGER.getIndent();

  tbb::parallel_for(std::size_t(0), simplices.bucket_count(),
                    [&](const uint i) {

    LOGGER.setIndent(saveIndent);

    for (auto it = simplices.begin(i); it != simplices.end(i); ++it) {

      dSimplex<D, Precision> &simplex = *it;

      PLOG("Updating neighbors of " << simplex << std::endl);

#ifndef NDEBUG
      dSimplex<D, Precision> saveSimplex = simplex;
#endif

      simplex.neighbors.clear();

      INDENT
      for (uint v = 0; v < D + 1; ++v) {
        // for every point, look where else its used
        // if the other simplex shares a second point -> it is a neighbor

        const dPoint<D, Precision> &vertex = points[simplex.vertices[v]];

        if (IS_PROLIX) {
          LOGGER.logContainer(vertex.simplices, Logger::Verbosity::PROLIX,
                              "Vertex " + to_string(vertex) + " used in");
        }

        INDENT
        for (const uint u : vertex.simplices) {
          if (u != simplex.id && simplices.contains(u) &&
              simplex.isNeighbor(simplices[u])) {
            PLOG("Neighbor with " << simplices[u] << std::endl);

            simplex.neighbors.insert(u);

            // LOGGER.logContainer(simplex.neighbors,
            // Logger::Verbosity::PROLIX);

            // ASSERT(simplex.neighbors.size() <= D+1);
            // u will be updated in its own round;
          }
        }
        DEDENT
      }
      DEDENT

      ASSERT(0 < simplex.neighbors.size() && simplex.neighbors.size() <= D + 1);

#ifndef NDEBUG
      if (!(simplex.neighbors.size() > 0 &&
            simplex.neighbors.size() <= D + 1)) {
        LOG("Error: wrong number of neighbors for simplex " << simplex
                                                            << std::endl);

        if (IS_PROLIX) {
          paintWriter.add("img/" + provenance + "_neighbors_" +
                              std::to_string(simplex.id),
                          bounds);
          paintWriter.top().draw(points);
          paintWriter.top().drawPartition(points);

          paintWriter.top().setColor(0, 0, 0, 0.4);
          paintWriter.top().draw(simplices, points, true);

          paintWriter.top().setColor(1, 0, 0); // simplex in red
          paintWriter.top().draw(simplex, points, true);

          paintWriter.top().setColor(1, 1, 0, 0.4); // neighbors in yellow
          paintWriter.top().drawNeighbors(simplex, simplices, points, true);

          if (!(simplex.equalVertices(saveSimplex) &&
                simplex.equalNeighbors(saveSimplex))) {
            LOG("Error: was before " << saveSimplex << std::endl);

            paintWriter.top().setColor(0, 0, 1); // old simplex in blue
            paintWriter.top().draw(saveSimplex, points, true);

            paintWriter.top().setColor(0, 1, 1,
                                       0.4); // old simplex neighbors in cyan
            paintWriter.top().drawNeighbors(saveSimplex, simplices, points,
                                            true);
          }
        }
      }
#endif

      // ASSERT(simplex.neighbors.size() > 0 && simplex.neighbors.size() <=
      // D+1);

      // TODO it might be better NOT to save the infinite simplex as neighbor
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

template <uint D, typename Precision>
dSimplices<D, Precision> Triangulator<D, Precision>::mergeTriangulation(
    std::vector<dSimplices<D, Precision>> &partialDTs, const Ids &edgeSimplices,
    const dSimplices<D, Precision> &edgeDT,
    const Partitioning<D, Precision> &partitioning,
    const std::string &provenance, const dSimplices<D, Precision> *realDT) {

  LOG("Merging partial DTs into one triangulation" << std::endl);
  dSimplices<D, Precision> DT;
  // use the first partition as estimator for the size of the triangulation
  DT.reserve(partialDTs.size() * partialDTs[0].size());

  for (uint i = 0; i < partialDTs.size(); ++i) {
    DT.insert(partialDTs[i].begin(), partialDTs[i].end());
  }

  auto edgePointIds = extractPoints(edgeSimplices, DT);

  //**********************
  auto paintMerging = [&](const dSimplices<D, Precision> &dt,
                          const std::string &name,
                          bool drawInfinite = false) -> Painter<D, Precision> {
    Painter<D, Precision> painter(bounds);

    painter.draw(points);
    painter.drawPartition(points);

    painter.setColor(0, 1, 0, 0.4);
    painter.draw(dt, points, drawInfinite);

    painter.setColor(1, 0, 0);
    painter.draw(points.project(edgePointIds));
    painter.save(name);

    if (realDT != nullptr) {
      painter.setColor(0, 0, 0, 0.1);
      painter.draw(*realDT, points);
    }

    painter.save(name + "_overlay");

    return painter;
  };
  //**********************

  paintMerging(DT, provenance + "_05a_merging_merged");

  const auto deletedSimplices = DT.project(edgeSimplices);

  paintMerging(deletedSimplices, provenance + "_05b_merging_deleted", true);

  // delete all simplices belonging to the edge from DT
  LOG("Striping triangulation from edge" << std::endl);
  for (const uint id : edgeSimplices) {
    ASSERT(DT.contains(id));

    // remove simplex from where used list
    for (uint p = 0; p < D + 1; ++p) {
      dPoint<D, Precision> &point = points[DT[id].vertices[p]];

      std::lock_guard<SpinMutex> lock(point.mtx);
      point.simplices.erase(id);
    }

    DT.erase(id);
  }

  //**********************
  auto painter = paintMerging(DT, provenance + "_05b_merging_stripped");

  painter.setColor(0, 0, 1, 0.4);
  painter.setDashed();
  painter.draw(edgeDT, points);

  painter.save(provenance + "_05b_merging_stripped+edge_overlay");

  painter.setColor(1, 0, 0, 0.4);
  painter.setDashed();
  painter.draw(deletedSimplices, points);
  painter.setDashed(false);

  painter.save(provenance + "_05b_merging_stripped+edge+deleted_overlay");
  //**********************

  // merge partial DTs and edge DT
  LOG("Merging triangulations" << std::endl);
  SpinMutex insertMtx;
  tbb::parallel_for(std::size_t(0), edgeDT.bucket_count(), [&](const uint i) {

    for (auto it = edgeDT.begin(i); it != edgeDT.end(i); ++it) {

      const dSimplex<D, Precision> &edgeSimplex = *it;

      // check whether edgeSimplex is completely contained in one partition
      uint p0 = partitioning.partition(edgeSimplex.vertices[0]);
      bool inOnePartition = partitioning[p0].contains(edgeSimplex);

      if (!inOnePartition) {
        std::lock_guard<SpinMutex> lock(insertMtx);
        DT.insert(edgeSimplex);
      } else {
        // the simplex is completely in one partition -> it must have been found
        // before
        auto it = std::find_if(deletedSimplices.begin(), deletedSimplices.end(),
                               [&](const dSimplex<D, Precision> &s) {
          return s.equalVertices(edgeSimplex);
        });
        if (it != deletedSimplices.end()) {
          std::lock_guard<SpinMutex> lock(insertMtx);
          DT.insert(edgeSimplex);
        }
      }
    }
  });

  paintMerging(DT, provenance + "_05c_merging_edge");

  ASSERT(DT.countDuplicates() == 0);
  updateNeighbors(DT, provenance);

  paintMerging(DT, provenance + "_05d_merging_finished");

  return DT;
}

template <uint D, typename Precision>
void Triangulator<D, Precision>::evaluateVerificationReport(
    const VerificationReport<D, Precision> &vr,
    const std::string &provenance) const {
  if (IS_PROLIX) {
    Painter<D, Precision> basePainter(bounds);
    basePainter.draw(points);
    basePainter.drawPartition(points);
    basePainter.setLogging();

    PainterBulkWriter<D, Precision> writer;
    for (const auto &inCircle : vr.inCircle) {
      writer.add("img/" + provenance + "_inCircle_" +
                     std::to_string(inCircle.first.id),
                 basePainter);

      writer.top().setColor(0, 1, 0);
      writer.top().draw(inCircle.first, points, true);
      writer.top().drawCircumSphere(inCircle.first, points, true);

      writer.top().setColor(1, 0, 0);
      writer.top().draw(points.project(inCircle.second));
    }
  }
}

template <uint D, typename Precision>
void Triangulator<D, Precision>::evaluateCrossCheckReport(
    const CrossCheckReport<D, Precision> &ccr, const std::string &provenance,
    const dSimplices<D, Precision> &DT, const dSimplices<D, Precision> &realDT,
    const Partitioning<D, Precision> *partitioning) const {
  if (IS_PROLIX) {
    Painter<D, Precision> basePainter(bounds);
    basePainter.draw(points);
    basePainter.drawPartition(points);
    basePainter.setLogging();

    PainterBulkWriter<D, Precision> writer;
    std::stringstream txt;

    for (const auto &missingSimplex : ccr.missing) {
      writer.add("img/" + provenance + "_missing_" +
                     std::to_string(missingSimplex.id),
                 basePainter);

      auto mySimplices = DT.findSimplices(missingSimplex.vertices);

      writer.top().setColor(0, 1, 0);
      writer.top().draw(mySimplices, points, true);

      writer.top().setColor(1, 0, 0);
      writer.top().draw(missingSimplex, points, true);

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
        txt << '\t' << points[p] << std::endl;
      }

      txt << "Cirumsphere: " << missingSimplex.circumsphere(points)
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
      writer.top().draw(realSimplices, points, true);

      writer.top().setColor(1, 0, 0);
      writer.top().draw(invalidSimplex, points, true);

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
Triangulator<D, Precision>::triangulateBase(const Ids partitionPoints,
                                            const std::string provenance) {

  LOGGER.setIndent(provenance.length());

  LOG("triangulateBASE called on level " << provenance << " with "
                                         << partitionPoints.size() << " points"
                                         << std::endl);

  bool isTOP = provenance == std::to_string(TOP);

  std::unique_ptr<dSimplices<D, Precision>> realDT = nullptr;
  if (isTOP && VERIFY) {
    LOG("Real triangulation" << std::endl);
    INDENT
    realDT = std::make_unique<dSimplices<D, Precision>>(
        CGALInterface<D, Precision>::triangulate(points, &partitionPoints));
    LOG("Real triangulation contains " << realDT->size() << " tetrahedra"
                                       << std::endl);
    DEDENT
  }

  INDENT
  // if this is the top-most triangulation, ignore infinite vertices
  auto dt = CGALInterface<D, Precision>::triangulate(points, &partitionPoints);
  LOG("Triangulation contains " << dt.size() << " tetrahedra" << std::endl);
  DEDENT

  if (isTOP && VERIFY) {
    LOG("Verifying CGAL triangulation" << std::endl);
    dt.verify(points.project(partitionPoints));
  }

#ifndef NDEBUG
  auto saveDT = dt;
#endif

  LOG("Updating neighbors" << std::endl);
  updateNeighbors(dt, provenance);

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

  if (isTOP && VERIFY) {
    LOG("Consistency check of triangulation" << std::endl);
    auto vr = dt.verify(points.project(partitionPoints));
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
Triangulator<D, Precision>::triangulateDAC(const Ids partitionPoints,
                                           const std::string provenance) {

  LOGGER.setIndent(provenance.length());

  LOG("triangulateDAC called on level " << provenance << " with "
                                        << partitionPoints.size() << " points"
                                        << std::endl);

  bool isTOP = provenance == std::to_string(TOP);

  if (partitionPoints.size() > baseThreshold) {
    LOG("Recursive case" << std::endl);

    // setup base painter
    Painter<D, Precision> basePainter(bounds);
    basePainter.draw(points.project(partitionPoints));
    basePainter.drawPartition(points.project(partitionPoints));
    basePainter.save(provenance + "_00_points");

    std::unique_ptr<dSimplices<D, Precision>> realDT = nullptr;
    if (isTOP && VERIFY) {
      // perform real triangulation
      LOG("Real triangulation" << std::endl);
      INDENT

      realDT = std::make_unique<dSimplices<D, Precision>>(
          CGALInterface<D, Precision>::triangulate(points, &partitionPoints));
      LOG("Real triangulation contains " << realDT->size() << " tetrahedra"
                                         << std::endl);
      DEDENT

      Painter<D, Precision> paintRealDT = basePainter;

      paintRealDT.setColor(0, 0, 1);
      paintRealDT.draw(*realDT, points);
      paintRealDT.setColor(0, 0, 0);

      paintRealDT.save(provenance + "_01_realDT");
    }

    // partition input
    LOG("Partioning" << std::endl);
    INDENT
    auto partioning =
        partitioner->partition(partitionPoints, points, provenance);
    DEDENT

    std::vector<dSimplices<D, Precision>> partialDTs;
    partialDTs.resize(partioning.size());

    tbb::parallel_for(std::size_t(0), partioning.size(), [&](const uint i) {
      LOGGER.setIndent(
          provenance.length()); // new thread, initialize Logger indent
      INDENT
      LOG("Partition " << i << " on thread " << std::this_thread::get_id()
                       << std::endl);
      partialDTs[i] =
          triangulateDAC(partioning[i].points, provenance + std::to_string(i));
      LOG("Triangulation " << i << " contains " << partialDTs[i].size()
                           << " tetrahedra" << std::endl);
      DEDENT
    });

    Painter<D, Precision> paintPartialDTs = basePainter;
    for (uint i = 0; i < partioning.size(); ++i) {
      paintPartialDTs.setColor(Painter<D, Precision>::tetradicColor(i), 0.4);
      paintPartialDTs.draw(partialDTs[i], points);

      Painter<D, Precision> paintPartialDT = basePainter;
      paintPartialDT.draw(partialDTs[i], points, true);
      paintPartialDT.save(provenance + "_02_partialDT_" + std::to_string(i));
    }

    paintPartialDTs.save(provenance + "_02_partialDTs");

    LOG("Extracting edges" << std::endl);
    INDENT

    Ids edgePointIds;
    Ids edgeSimplexIds;
    Painter<D, Precision> paintEdges = paintPartialDTs;
    paintEdges.setColor(1, 0, 0);

    for (uint i = 0; i < partioning.size(); ++i) {
      // points are in different partitions, there can be no overlap

      auto edge = getEdge(partialDTs[i], partioning, i);
      edgeSimplexIds.insert(edge.begin(), edge.end());

      // ignore infinite vertices if this is the top most triangulation
      auto ep = extractPoints(edge, partialDTs[i]);
      edgePointIds.insert(ep.begin(), ep.end());

      paintEdges.draw(partialDTs[i].project(edge), points, false);
    }
    DEDENT

    LOG("Edge has " << edgeSimplexIds.size() << " simplices with "
                    << edgePointIds.size() << " points" << std::endl);

    paintEdges.draw(points.project(edgePointIds));
    paintEdges.save(provenance + "_03_edgeMarked");

    LOG("Triangulating edges" << std::endl);
    INDENT
    auto edgeDT = triangulateDAC(edgePointIds, provenance + "e");
    LOG("Edge triangulation contains " << edgeDT.size() << " tetrahedra"
                                       << std::endl << std::endl);
    DEDENT

    Painter<D, Precision> paintEdgeDT = basePainter;

    paintEdgeDT.setColor(0, 1, 0);
    paintEdgeDT.draw(edgeDT, points, true);
    paintEdgeDT.setColor(1, 0, 0);
    paintEdgeDT.draw(points.project(edgePointIds));

    paintEdgeDT.save(provenance + "_04_edgeDT");

    auto mergedDT = mergeTriangulation(partialDTs, edgeSimplexIds, edgeDT,
                                       partioning, provenance, realDT.get());

    Painter<D, Precision> paintFinal = basePainter;
    paintFinal.setColor(0, 1, 0);
    paintFinal.draw(mergedDT, points, true);
    paintFinal.save(provenance + "_06_final");

    TriangulationReportEntry rep;
    rep.provenance = provenance;
    rep.base_case = false;
    rep.valid = false;
    rep.edge_triangulation = provenance.find('e') != std::string::npos;
    rep.nPoints = partitionPoints.size();
    rep.nSimplices = mergedDT.size();
    rep.nEdgePoints = edgePointIds.size();
    rep.nEdgeSimplices = edgeDT.size();

    if (isTOP & VERIFY) {
      LOG("Consistency check of triangulation" << std::endl);
      auto vr = mergedDT.verify(points.project(partitionPoints));
      evaluateVerificationReport(vr, provenance);

      LOG("Cross check with real triangulation" << std::endl);
      auto ccr = mergedDT.crossCheck(*realDT);
      evaluateCrossCheckReport(ccr, provenance, mergedDT, *realDT, &partioning);

      rep.valid = vr.valid && ccr.valid;
    }

    triangulationReport.push_back(rep);

    return mergedDT;

  } else { // base case
    return triangulateBase(partitionPoints, provenance);
  }
}

template <uint D, typename Precision>
dSimplices<D, Precision> Triangulator<D, Precision>::triangulate() {
  Ids allPoints(points.begin_keys(), points.end_keys());
  return triangulateDAC(allPoints, std::to_string(TOP));
}

// specializations

template class Triangulator<2, float>;
template class Triangulator<3, float>;

template class Triangulator<2, double>;
template class Triangulator<3, double>;

#include "Triangulator.h"

// std library
#include <fstream>
#include <iostream>
#include <vector>
#include <functional>
#include <set>
#include <stdexcept>
#include <queue>

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

template <uint D> constexpr tCoordinate Triangulator<D>::SAFETY;

template <uint D> constexpr uint Triangulator<D>::BASE_CASE;

template <uint D> constexpr char Triangulator<D>::TOP;
//**************************

template <uint D>
Triangulator<D>::Triangulator(dBox<D> &_bounds, dPoints<D> &_points,
                              std::unique_ptr<Partitioner<D>> &&_partitioner)
    : bounds(_bounds), points(_points), partitioner(_partitioner) {
  if (!points.contains(dPoint<D>::cINF)) {
    auto stats = getPointStats(points.begin_keys(), points.end_keys(), points);
    for (uint i = 0; i < pow(2, D); ++i) {
      VLOG << "Point stats: " << stats.min << " - " << stats.mid << " - "
           << stats.max << std::endl;

      dPoint<D> p = stats.mid;
      p.id = dPoint<D>::cINF | i;

      for (uint d = 0; d < D; ++d)
        p.coords[d] += (i & (1 << d) ? 1 : -1) * 2 * SAFETY *
                       (stats.max.coords[d] - stats.min.coords[d]);

      points.insert(p);
    }
  }
}

template <uint D>
Ids Triangulator<D>::getEdge(const dSimplices<D> &simplices,
                             const Partition<D> &partition) {
  Ids edgeSimplices;
  std::set<uint> wqa; // set of already checked simplices

  // we use the overflow of the uint to zero to abort the loop
  for (uint infVertex = dPoint<D>::cINF; infVertex != 0; ++infVertex) {
    ASSERT(
        std::find(partition.points.begin(), partition.points.end(),
                  infVertex) !=
        partition.points.end()); // the infinite point must be in the partition
    ASSERT(points.contains(infVertex));

    PLOG << "Vertex " << points[infVertex] << std::endl;

    INDENT
    for (const auto &s : points[infVertex].simplices) {
      if (!simplices.contains(s))
        continue;

      PLOG << "Adding " << simplices[s] << " to edge" << std::endl;
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

        if (simplices.contains(x) &&
            !partition.bounds.contains(simplices[x].circumsphere(points))) {
          PLOG << "Adding " << simplices[x]
               << " to edge -> circumcircle criterion" << std::endl;
          edgeSimplices.insert(simplices[x].id);

          for (const auto &n : simplices[x].neighbors) {
            if (wqa.insert(n).second) {
              // n was not yet inspected
              wq.push_back(n);
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

template <uint D>
Ids Triangulator<D>::extractPoints(const Ids &edgeSimplices,
                                   const dSimplices<D> &simplices,
                                   bool ignoreInfinite) {
  Ids outPoints;
  std::set<uint> idx;

  for (const auto &id : edgeSimplices) {
    ASSERT(simplices.contains(id));

    for (uint i = 0; i < D + 1; ++i) {
      if (dPoint<D>::isFinite(simplices[id].vertices[i]) &&
          idx.insert(simplices[id].vertices[i]).second)
        outPoints.insert(simplices[id].vertices[i]);
    }
  }

  // add the extreme infinite points to the set
  if (!ignoreInfinite) {
    for (uint k = dPoint<D>::cINF; k != 0; ++k) {
      outPoints.insert(k);
    }
  }

  return outPoints;
}

template <uint D>
void Triangulator<D>::updateNeighbors(dSimplices<D> &simplices,
                                      const std::string &provenance) {
  PainterBulkWriter<D> paintWriter;

  INDENT
  for (dSimplex<D> &simplex : simplices) {
    PLOG << "Updating neighbors of " << simplex << std::endl;

#ifndef NDEBUG
    dSimplex<D> saveSimplex = simplex;
#endif

    simplex.neighbors.clear();

    INDENT
    for (uint v = 0; v < D + 1; ++v) {
      // for every point, look where else its used
      // if the other simplex shares a second point -> it is a neighbor

      const dPoint<D> &vertex = points[simplex.vertices[v]];

      if (IS_PROLIX) {
        LOGGER.logContainer(vertex.simplices, Logger::Verbosity::PROLIX,
                            "Vertex " + to_string(vertex) + " used in");
      }

      INDENT
      for (const uint u : vertex.simplices) {
        if (u != simplex.id && simplices.contains(u) &&
            simplex.isNeighbor(simplices[u])) {
          PLOG << "Neighbor with " << simplices[u] << std::endl;

          simplex.neighbors.insert(u);

          // LOGGER.logContainer(simplex.neighbors, Logger::Verbosity::PROLIX);

          // ASSERT(simplex.neighbors.size() <= D+1);
          // u will be updated in its own round;
        }
      }
      DEDENT
    }
    DEDENT

    if (!(simplex.neighbors.size() > 0 && simplex.neighbors.size() <= D + 1)) {
      LOG << "Error: wrong number of neighbors for simplex " << simplex
          << std::endl;

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
#ifndef NDEBUG
        if (!(simplex.equalVertices(saveSimplex) &&
              simplex.equalNeighbors(saveSimplex))) {
          LOG << "Error: was before " << saveSimplex << std::endl;

          paintWriter.top().setColor(0, 0, 1); // old simplex in blue
          paintWriter.top().draw(saveSimplex, points, true);

          paintWriter.top().setColor(0, 1, 1,
                                     0.4); // old simplex neighbors in cyan
          paintWriter.top().drawNeighbors(saveSimplex, simplices, points, true);
        }
#endif
      }
    }

    // ASSERT(simplex.neighbors.size() > 0 && simplex.neighbors.size() <= D+1);

    // TODO it might be better NOT to save the infinite simplex as neighbor
    if (simplex.neighbors.size() < D + 1) { // if it is a simplex at the border,
                                            // it might have one infinite
                                            // simplex as neighbor
      simplex.neighbors.insert(uint(dSimplex<D>::cINF));
    }
  }
  DEDENT
}

template <uint D>
dSimplices<D> Triangulator<D>::mergeTriangulation(
    std::vector<dSimplices<D>> &partialDTs, const Ids &edgeSimplices,
    const dSimplices<D> &edgeDT, const Partitioning<D> &partitioning,
    const std::string &provenance, const dSimplices<D> *realDT) {

  LOG << "Merging partial DTs into one triangulation" << std::endl;
  dSimplices<D> DT;
  for (uint i = 0; i < partialDTs.size(); ++i)
    DT.insert(partialDTs[i].begin(), partialDTs[i].end());

  auto edgePointIds = extractPoints(edgeSimplices, DT);

  //**********************
  auto paintMerging =
      [&](const dSimplices<D> &dt, const std::string &name) -> Painter<D> {
    Painter<D> painter(bounds);

    painter.draw(points);
    painter.drawPartition(points);

    painter.setColor(0, 1, 0, 0.4);
    painter.draw(dt, points);

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

  auto deletedSimplices = DT.project(edgeSimplices);

  // delete all simplices belonging to the edge from DT
  LOG << "Striping triangulation from edge" << std::endl;
  for (const uint id : edgeSimplices) {
    ASSERT(DT.contains(id));

    // remove simplex from where used list
    for (uint p = 0; p < D + 1; ++p) {
      points[DT[id].vertices[p]].simplices.erase(id);
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
  LOG << "Merging triangulations" << std::endl;
  for (const auto &edgeSimplex : edgeDT) {
    // check whether edgeSimplex is completely contained in one partition
    uint p0 = partitioning.partition(edgeSimplex.vertices[0]);
    bool inOnePartition = partitioning[p0].contains(edgeSimplex);

    if (!inOnePartition) {
      DT.insert(edgeSimplex);
    } else {
      // the simplex is completely in one partition -> it must have been found
      // before
      auto it = std::find_if(
          deletedSimplices.begin(), deletedSimplices.end(),
          [&](const dSimplex<D> &s) { return s.equalVertices(edgeSimplex); });
      if (it != deletedSimplices.end()) {
        DT.insert(edgeSimplex);
      }
    }
  }

  paintMerging(DT, provenance + "_05c_merging_edge");

  ASSERT(DT.countDuplicates() == 0);
  updateNeighbors(DT, provenance);

  paintMerging(DT, provenance + "_05d_merging_finished");

  return DT;
}

template <uint D>
void Triangulator<D>::evaluateVerificationReport(
    const VerificationReport<D> &vr, const std::string &provenance) const {
  if (IS_PROLIX) {
    Painter<D> basePainter(bounds);
    basePainter.draw(points);
    basePainter.drawPartition(points);
    basePainter.setLogging();

    PainterBulkWriter<D> writer;
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

template <uint D>
void Triangulator<D>::evaluateCrossCheckReport(
    const CrossCheckReport<D> &ccr, const std::string &provenance,
    const dSimplices<D> &DT, const dSimplices<D> &realDT,
    const Partitioning<D> *partitioning) const {
  if (IS_PROLIX) {
    Painter<D> basePainter(bounds);
    basePainter.draw(points);
    basePainter.drawPartition(points);
    basePainter.setLogging();

    PainterBulkWriter<D> writer;
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

template <uint D>
dSimplices<D> Triangulator<D>::triangulateBase(const Ids partitionPoints,
                                               const std::string provenance) {
  LOG << "triangulateBASE called on level " << provenance << " with "
      << partitionPoints.size() << " points" << std::endl;

  LOG << "Real triangulation" << std::endl;
  INDENT
  auto realDT = delaunayCgal(points, &partitionPoints, true);
  LOG << "Real triangulation contains " << realDT.size() << " tetrahedra"
      << std::endl;
  DEDENT

  // if this is the top-most triangulation, ignore infinite vertices
  bool ignoreInfinite = provenance == std::to_string(TOP);

  INDENT
  auto dt = delaunayCgal(points, &partitionPoints, ignoreInfinite);
  LOG << "Triangulation contains " << dt.size() << " tetrahedra" << std::endl;
  DEDENT

  LOG << "Verifying CGAL triangulation" << std::endl;
  dt.verify(points.project(partitionPoints));

#ifndef NDEBUG
  auto saveDT = dt;
#endif

  LOG << "Updating neighbors" << std::endl;
  updateNeighbors(dt, provenance);

  LOG << "Consistency check of triangulation" << std::endl;
  auto vr = dt.verify(points.project(partitionPoints));
  evaluateVerificationReport(vr, provenance);

  LOG << "Cross check with real triangulation" << std::endl;
  auto ccr = dt.crossCheck(realDT);
  evaluateCrossCheckReport(ccr, provenance, dt, realDT);

  ASSERT(saveDT == dt); // only performed if not NDEBUG

  TriangulationReportEntry rep;
  rep.provenance = provenance;
  rep.base_case = true;
  rep.edge_triangulation = provenance.find('e') != std::string::npos;
  rep.valid = vr.valid && ccr.valid;
  rep.nPoints = partitionPoints.size();
  rep.nSimplices = dt.size();
  rep.nEdgePoints = 0;
  rep.nEdgeSimplices = 0;
  triangulationReport.push_back(rep);

  return dt;
}

template <uint D>
dSimplices<D> Triangulator<D>::triangulateDAC(const Ids partitionPoints,
                                              const std::string provenance) {
  LOG << "triangulateDAC called on level " << provenance << " with "
      << partitionPoints.size() << " points" << std::endl;

  if (partitionPoints.size() > BASE_CASE) {
    LOG << "Recursive case" << std::endl;

    // setup base painter
    Painter<D> basePainter(bounds);
    basePainter.draw(points.project(partitionPoints));
    basePainter.drawPartition(points.project(partitionPoints));
    basePainter.save(provenance + "_00_points");

    // perform real triangulation
    LOG << "Real triangulation" << std::endl;
    INDENT
    auto realDT = delaunayCgal(points, &partitionPoints, true);
    LOG << "Real triangulation contains " << realDT.size() << " tetrahedra"
        << std::endl;
    DEDENT

    Painter<D> paintRealDT = basePainter;

    paintRealDT.setColor(0, 0, 1);
    paintRealDT.draw(realDT, points);
    paintRealDT.setColor(0, 0, 0);

    paintRealDT.save(provenance + "_01_realDT");

    // partition input
    LOG << "Partioning" << std::endl;
    INDENT
    auto partioning =
        partitioner->partition(partitionPoints, points, provenance);
    DEDENT

    for (auto &p : partioning)
      PLOG << "Partition " << p << std::endl;

    std::vector<dSimplices<D>> partialDTs;
    partialDTs.resize(partioning.size());
    Painter<D> paintPartialDTs = basePainter;

    INDENT
    for (uint i = 0; i < partioning.size(); ++i) {
      LOG << "Partition " << i << std::endl;
      INDENT
      partialDTs[i] =
          triangulateDAC(partioning[i].points, provenance + std::to_string(i));
      LOG << "Triangulation " << i << " contains " << partialDTs[i].size()
          << " tetrahedra" << std::endl;
      DEDENT

      paintPartialDTs.setColor(Painter<D>::tetradicColor(i), 0.4);
      paintPartialDTs.draw(partialDTs[i], points);

      Painter<D> paintPartialDT = basePainter;
      paintPartialDT.draw(partialDTs[i], points, true);
      paintPartialDT.save(provenance + "_02_partialDT_" + std::to_string(i));
    }
    DEDENT

    paintPartialDTs.save(provenance + "_02_partialDTs");

    LOG << "Extracting edges" << std::endl;
    INDENT

    Ids edgePointIds;
    Ids edgeSimplexIds;
    Painter<D> paintEdges = paintPartialDTs;
    paintEdges.setColor(1, 0, 0);

    // ignore infinite vertices if this is the top most triangulation
    bool ignoreInfinite = provenance == std::to_string(TOP);

    for (uint i = 0; i < partioning.size(); ++i) {
      // points are in different partitions, there can be no overlap

      auto edge = getEdge(partialDTs[i], partioning[i]);
      edgeSimplexIds.insert(edge.begin(), edge.end());

      auto ep = extractPoints(edge, partialDTs[i], ignoreInfinite);
      edgePointIds.insert(ep.begin(), ep.end());

      paintEdges.draw(partialDTs[i].project(edge), points, false);
    }
    DEDENT

    LOG << "Edge has " << edgeSimplexIds.size() << " simplices with "
        << edgePointIds.size() << " points" << std::endl;

    paintEdges.draw(points.project(edgePointIds));
    paintEdges.save(provenance + "_03_edgeMarked");

    LOG << "Triangulating edges" << std::endl;
    INDENT
    auto edgeDT = triangulateBase(edgePointIds, provenance + "e");
    LOG << "Edge triangulation contains " << edgeDT.size() << " tetrahedra"
        << std::endl << std::endl;
    DEDENT

    Painter<D> paintEdgeDT = basePainter;

    paintEdgeDT.setColor(0, 1, 0);
    paintEdgeDT.draw(edgeDT, points, true);
    paintEdgeDT.setColor(1, 0, 0);
    paintEdgeDT.draw(points.project(edgePointIds));

    paintEdgeDT.save(provenance + "_04_edgeDT");

    auto mergedDT = mergeTriangulation(partialDTs, edgeSimplexIds, edgeDT,
                                       partioning, provenance, &realDT);

    Painter<D> paintFinal = basePainter;
    paintFinal.setColor(0, 1, 0);
    paintFinal.draw(mergedDT, points, true);
    paintFinal.save(provenance + "_06_final");

    LOG << "Consistency check of triangulation" << std::endl;
    auto vr = mergedDT.verify(points.project(partitionPoints));
    evaluateVerificationReport(vr, provenance);

    LOG << "Cross check with real triangulation" << std::endl;
    auto ccr = mergedDT.crossCheck(realDT);
    evaluateCrossCheckReport(ccr, provenance, mergedDT, realDT, &partioning);

    TriangulationReportEntry rep;
    rep.provenance = provenance;
    rep.base_case = false;
    rep.edge_triangulation = provenance.find('e') != std::string::npos;
    rep.valid = vr.valid && ccr.valid;
    rep.nPoints = partitionPoints.size();
    rep.nSimplices = mergedDT.size();
    rep.nEdgePoints = edgePointIds.size();
    rep.nEdgeSimplices = edgeDT.size();
    triangulationReport.push_back(rep);

    return mergedDT;

  } else { // base case
    return triangulateBase(partitionPoints, provenance);
  }
}

template <uint D> dSimplices<D> Triangulator<D>::triangulate() {
  Ids allPoints(points.begin_keys(), points.end_keys());
  return triangulateDAC(allPoints, std::to_string(TOP));
}

// specializations

template class Triangulator<2>;
template class Triangulator<3>;

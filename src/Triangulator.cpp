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
#include "Timings.h"
#include "Geometry.h"
#include "Painter.h"
#include "Logger.h"
#include "CGAL_Interface.h"

//**************************

constexpr tCoordinate Triangulator::SAFETY;
constexpr uint Triangulator::BASE_CASE;
constexpr char Triangulator::TOP;
//**************************

Triangulator::Triangulator(dBox &_bounds, dPoints &_points,
                           std::unique_ptr<Partitioner> &&_partitioner)
    : bounds(_bounds), points(_points), partitioner(_partitioner) {
  if (!points.contains(dPoint::cINF)) {
    auto stats = getPointStats(points.begin_keys(), points.end_keys(), points);
    for (uint i = 0; i < pow(2, D); ++i) {
      VLOG << "Point stats: " << stats.min << " - " << stats.mid << " - "
           << stats.max << std::endl;

      dPoint p = stats.mid;
      p.id = dPoint::cINF | i;

      for (uint d = 0; d < D; ++d)
        p.coords[d] += (i & (1 << d) ? 1 : -1) * 2 * SAFETY *
                       (stats.max.coords[d] - stats.min.coords[d]);

      points.insert(p);
    }
  }
}

Ids Triangulator::getEdge(const dSimplices &simplices,
                          const Partition &partition) {
  Ids edgeSimplices;
  std::set<uint> wqa; // set of already checked simplices

  // we use the overflow of the uint to zero to abort the loop
  for (uint infVertex = dPoint::cINF; infVertex != 0; ++infVertex) {
    assert(
        std::find(partition.points.begin(), partition.points.end(),
                  infVertex) !=
        partition.points.end()); // the infinite point must be in the partition
    assert(points.contains(infVertex));

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
            !partition.bounds.contains(simplices[x].circumcircle(points))) {
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

Ids Triangulator::extractPoints(const Ids &edgeSimplices,
                                const dSimplices &simplices,
                                bool ignoreInfinite) {
  Ids outPoints;
  std::set<uint> idx;

  for (const auto &id : edgeSimplices) {
    assert(simplices.contains(id));

    for (uint i = 0; i < D + 1; ++i) {
      if (dPoint::isFinite(simplices[id].vertices[i]) &&
          idx.insert(simplices[id].vertices[i]).second)
        outPoints.insert(simplices[id].vertices[i]);
    }
  }

  // add the extreme infinite points to the set
  if (!ignoreInfinite) {
    for (uint k = dPoint::cINF; k != 0; ++k) {
      outPoints.insert(k);
    }
  }

  return outPoints;
}

void Triangulator::updateNeighbors(dSimplices &simplices,
                                   const std::string &provenance) {
  PainterBulkWriter paintWriter;

  INDENT
  for (dSimplex &simplex : simplices) {
    PLOG << "Updating neighbors of " << simplex << std::endl;

#ifndef NDEBUG
    dSimplex saveSimplex = simplex;
#endif

    simplex.neighbors.clear();

    INDENT
    for (uint v = 0; v < D + 1; ++v) {
      // for every point, look where else its used
      // if the other simplex shares a second point -> it is a neighbor

      const dPoint &vertex = points[simplex.vertices[v]];

      if (IS_PROLIX) {
        PLOG << "Vertex " << vertex << " used in ";
        LOGGER.logContainer(vertex.simplices, Logger::Verbosity::PROLIX);
      }

      INDENT
      for (const uint u : vertex.simplices) {
        if (u != simplex.id && simplices.contains(u) &&
            simplex.isNeighbor(simplices[u])) {
          PLOG << "Neighbor with " << simplices[u] << std::endl;

          simplex.neighbors.insert(u);

          LOGGER.logContainer(simplex.neighbors, Logger::Verbosity::PROLIX);

          // assert(simplex.neighbors.size() <= D+1);
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
                            std::to_string(simplex.id) + ".png",
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

    // assert(simplex.neighbors.size() > 0 && simplex.neighbors.size() <= D+1);

    // TODO it might be better NOT to save the infinite simplex as neighbor
    if (simplex.neighbors.size() < D + 1) { // if it is a simplex at the border,
                                            // it might have one infinite
                                            // simplex as neighbor
      simplex.neighbors.insert(uint(dSimplex::cINF));
    }
  }
  DEDENT
}

dSimplices Triangulator::mergeTriangulation(std::vector<dSimplices> &partialDTs,
                                            const Ids &edgeSimplices,
                                            const dSimplices &edgeDT,
                                            const Partitioning &partitioning,
                                            const std::string &provenance,
                                            const dSimplices *realDT) {
  auto partitionPoint = [&](const uint &point) -> uint {

    for (uint p = 0; p < partitioning.size(); ++p) {
      if (std::find(partitioning[p].points.begin(),
                    partitioning[p].points.end(),
                    point) != partitioning[p].points.end())
        return p;
    }

    throw std::out_of_range(std::to_string(point) + "not found");

  };

  auto partitionContains = [&](const dSimplex &simplex) -> bool {
    uint p0 = partitionPoint(simplex.vertices[0]);
    uint p1 = partitionPoint(simplex.vertices[1]);
    uint p2 = partitionPoint(simplex.vertices[2]);

    return p0 == p1 && p1 == p2;
  };

  LOG << "Merging partial DTs into one triangulation" << std::endl;
  dSimplices DT;
  for (uint i = 0; i < partialDTs.size(); ++i)
    DT.insert(partialDTs[i].begin(), partialDTs[i].end());

  auto edgePointIds = extractPoints(edgeSimplices, DT);

  auto paintMerging =
      [&](const dSimplices &dt, const std::string &name) -> Painter {
    Painter painter(bounds);

    painter.draw(points);
    painter.drawPartition(points);

    painter.setColor(0, 1, 0, 0.4);
    painter.draw(dt, points);

    painter.setColor(1, 0, 0);
    painter.draw(points.project(edgePointIds));
    painter.savePNG(name + ".png");

    if (realDT != nullptr) {
      painter.setColor(0, 0, 0, 0.1);
      painter.draw(*realDT, points);
    }

    painter.savePNG(name + "_overlay.png");

    return painter;
  };

  paintMerging(DT, provenance + "_05a_merging_merged");

  auto deletedSimplices = DT.project(edgeSimplices);

  // delete all simplices belonging to the edge from DT
  LOG << "Striping triangulation from edge" << std::endl;
  for (const uint id : edgeSimplices) {
    assert(DT.contains(id));

    // remove simplex from where used list
    for (uint p = 0; p < D + 1; ++p) {
      points[DT[id].vertices[p]].simplices.erase(id);
    }

    DT.erase(id);
  }

  auto painter = paintMerging(DT, provenance + "_05b_merging_stripped");

  painter.setColor(0, 0, 1, 0.4);
  painter.setLineDash({2, 4});
  painter.draw(edgeDT, points);

  painter.savePNG(provenance + "_05b_merging_stripped+edge_overlay.png");

  painter.setColor(1, 0, 0, 0.4);
  painter.setLineDash({2, 4});
  painter.draw(deletedSimplices, points);
  painter.unsetLineDash();

  painter.savePNG(provenance +
                  "_05b_merging_stripped+edge+deleted_overlay.png");

  // merge partial DTs and edge DT
  LOG << "Merging triangulations" << std::endl;
  for (const auto &edgeSimplex : edgeDT) {
    if (!partitionContains(edgeSimplex)) {
      DT.insert(edgeSimplex);
    } else {
      auto it = std::find_if(
          deletedSimplices.begin(), deletedSimplices.end(),
          [&](const dSimplex &s) { return s.equalVertices(edgeSimplex); });
      if (it != deletedSimplices.end())
        DT.insert(edgeSimplex);
    }
  }

  paintMerging(DT, provenance + "_05c_merging_edge");

  assert(DT.countDuplicates() == 0);
  updateNeighbors(DT, provenance);

  paintMerging(DT, provenance + "_05d_merging_finished");

  return DT;
}

void
Triangulator::evaluateVerificationReport(const VerificationReport &vr,
                                         const std::string &provenance) const {
  if (IS_PROLIX) {
    Painter basePainter(bounds);
    basePainter.draw(points);
    basePainter.drawPartition(points);
    basePainter.setLogging();

    PainterBulkWriter writer;
    for (const auto &inCircle : vr.inCircle) {
      writer.add("img/" + provenance + "_inCircle_" +
                     std::to_string(inCircle.first.id) + ".png",
                 basePainter);

      writer.top().setColor(0, 1, 0);
      writer.top().draw(inCircle.first, points, true);
      writer.top().drawCircumCircle(inCircle.first, points, true);

      writer.top().setColor(1, 0, 0);
      writer.top().draw(points.project(inCircle.second));
    }
  }
}

void Triangulator::evaluateCrossCheckReport(const CrossCheckReport &ccr,
                                            const std::string &provenance,
                                            const dSimplices &DT,
                                            const dSimplices &realDT) const {
  if (IS_PROLIX) {
    Painter basePainter(bounds);
    basePainter.draw(points);
    basePainter.drawPartition(points);
    basePainter.setLogging();

    PainterBulkWriter writer;
    for (const auto &missingSimplex : ccr.missing) {
      writer.add("img/" + provenance + "_missing_" +
                     std::to_string(missingSimplex.id) + ".png",
                 basePainter);

      auto mySimplices = DT.findSimplices(missingSimplex.vertices);

      writer.top().setColor(0, 1, 0);
      writer.top().draw(mySimplices, points, true);

      writer.top().setColor(1, 0, 0);
      writer.top().draw(missingSimplex, points, true);
    }

    for (const auto &invalidSimplex : ccr.invalid) {
      writer.add("img/" + provenance + "_invalid_" +
                     std::to_string(invalidSimplex.id) + ".png",
                 basePainter);

      auto realSimplices = realDT.findSimplices(invalidSimplex.vertices);

      writer.top().setColor(0, 1, 0);
      writer.top().draw(realSimplices, points, true);

      writer.top().setColor(1, 0, 0);
      writer.top().draw(invalidSimplex, points, true);
    }
  }
}

dSimplices Triangulator::triangulateBase(const Ids partitionPoints,
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

  assert(saveDT == dt); // only performed if not NDEBUG

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

dSimplices Triangulator::triangulateDAC(const Ids partitionPoints,
                                        const std::string provenance) {
  LOG << "triangulateDAC called on level " << provenance << " with "
      << partitionPoints.size() << " points" << std::endl;

  if (partitionPoints.size() > BASE_CASE) {
    LOG << "Recursive case" << std::endl;

    LOG << "Partioning" << std::endl;
    INDENT
    auto partioning =
        partitioner->partition(partitionPoints, points, provenance);
    DEDENT

    Painter basePainter(bounds);
    basePainter.draw(points.project(partitionPoints));
    basePainter.drawPartition(points.project(partitionPoints));
    basePainter.savePNG(provenance + "_00_points.png");

    for (auto &p : partioning)
      PLOG << "Partition " << p << std::endl;

    LOG << "Real triangulation" << std::endl;
    INDENT
    auto realDT = delaunayCgal(points, &partitionPoints, true);
    LOG << "Real triangulation contains " << realDT.size() << " tetrahedra"
        << std::endl;
    DEDENT

    Painter paintRealDT = basePainter;

    paintRealDT.setColor(0, 0, 1);
    paintRealDT.draw(realDT, points);
    paintRealDT.setColor(0, 0, 0);

    paintRealDT.savePNG(provenance + "_01_realDT.png");

    std::vector<dSimplices> partialDTs;
    partialDTs.resize(partioning.size());
    Painter paintPartialDTs = basePainter;

    INDENT
    for (uint i = 0; i < partioning.size(); ++i) {
      LOG << "Partition " << i << std::endl;
      INDENT
      partialDTs[i] =
          triangulateDAC(partioning[i].points, provenance + std::to_string(i));
      LOG << "Triangulation " << i << " contains " << partialDTs[i].size()
          << " tetrahedra" << std::endl;
      DEDENT

      paintPartialDTs.setColor(Painter::tetradicColor(i), 0.4);
      paintPartialDTs.draw(partialDTs[i], points);

      Painter paintPartialDT = basePainter;
      paintPartialDT.draw(partialDTs[i], points, true);
      paintPartialDT.savePNG(provenance + "_02_partialDT_" + std::to_string(i) +
                             ".png");
    }
    DEDENT

    paintPartialDTs.savePNG(provenance + "_02_partialDTs.png");

    LOG << "Extracting edges" << std::endl;
    INDENT

    Ids edgePointIds;
    Ids edgeSimplexIds;
    Painter paintEdges = paintPartialDTs;
    paintEdges.setColor(1, 0, 0);

    // ignore infinite vertices of this is the top most triangulation
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
    paintEdges.savePNG(provenance + "_03_edgeMarked.png");

    LOG << "Triangulating edges" << std::endl;
    INDENT
    auto edgeDT = triangulateBase(edgePointIds, provenance + "e");
    LOG << "Edge triangulation contains " << edgeDT.size() << " tetrahedra"
        << std::endl << std::endl;
    DEDENT

    Painter paintEdgeDT = basePainter;

    paintEdgeDT.setColor(0, 1, 0);
    paintEdgeDT.draw(edgeDT, points, true);
    paintEdgeDT.setColor(1, 0, 0);
    paintEdgeDT.draw(points.project(edgePointIds));

    paintEdgeDT.savePNG(provenance + "_04_edgeDT.png");

    auto mergedDT = mergeTriangulation(partialDTs, edgeSimplexIds, edgeDT,
                                       partioning, provenance, &realDT);

    Painter paintFinal = basePainter;
    paintFinal.setColor(0, 1, 0);
    paintFinal.draw(mergedDT, points, true);
    paintFinal.savePNG(provenance + "_06_final.png");

    LOG << "Consistency check of triangulation" << std::endl;
    auto vr = mergedDT.verify(points.project(partitionPoints));
    evaluateVerificationReport(vr, provenance);

    LOG << "Cross check with real triangulation" << std::endl;
    auto ccr = mergedDT.crossCheck(realDT);
    evaluateCrossCheckReport(ccr, provenance, mergedDT, realDT);

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

dSimplices Triangulator::triangulate() {
  Ids allPoints(points.begin_keys(), points.end_keys());
  return triangulateDAC(allPoints, std::to_string(TOP));
}

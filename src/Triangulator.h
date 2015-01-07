#pragma once

// stl
#include <functional>
#include <vector>

// own
#include "Geometry.h"
#include "Partitioner.h"

struct TriangulationReportEntry {
  std::string provenance;
  bool base_case;
  bool edge_triangulation;
  bool valid;

  uint nPoints;
  uint nSimplices;
  uint nEdgePoints;
  uint nEdgeSimplices;
};

typedef std::vector<TriangulationReportEntry> TriangulationReport;

class Triangulator {
public:
  Triangulator(dBox &_bounds, dPoints &_points,
               std::unique_ptr<Partitioner> &&_partitioner);

  dSimplices triangulate();

  const TriangulationReport &getTriangulationReport() const {
    return triangulationReport;
  }

private:
  dSimplices triangulateBase(const Ids partitionPoints,
                             const std::string provenance);

  dSimplices triangulateDAC(const Ids partitionPoints,
                            const std::string provenance);

  Ids getEdge(const dSimplices &simplices, const Partition &partition);

  Ids extractPoints(const Ids &edgeSimplices, const dSimplices &simplices,
                    bool ignoreInfinite = false);

  void updateNeighbors(dSimplices &simplices, const std::string &provenance);

  dSimplices mergeTriangulation(std::vector<dSimplices> &partialDTs,
                                const Ids &edgeSimplices,
                                const dSimplices &edgeDT,
                                const Partitioning &partitioning,
                                const std::string &provenance,
                                const dSimplices *realDT = nullptr);

  void evaluateVerificationReport(const VerificationReport &vr,
                                  const std::string &provenance) const;

  void evaluateCrossCheckReport(const CrossCheckReport &ccr,
                                const std::string &provenance,
                                const dSimplices &DT,
                                const dSimplices &realDT) const;

private:
  dBox bounds;
  dPoints points;

  TriangulationReport triangulationReport;

  std::unique_ptr<Partitioner> &partitioner;

  static constexpr tCoordinate SAFETY = 100;
  static constexpr uint BASE_CASE = 100;
  static constexpr char TOP = 0;
};

#pragma once

// stl
#include <functional>

// own
#include "Geometry.h"
#include "Partitioner.h"

class Triangulator {

public:
  Triangulator(dBox &_bounds, dPoints &_points,
               const Partitioner &_partitioner);

  dSimplices triangulate();

private:
  dSimplices triangulateBase(const Ids partitionPoints,
                             const std::string provenance);

  dSimplices triangulateDAC(const Ids partitionPoints,
                            const std::string provenance);

  Ids getEdge(const dSimplices &simplices, const Partition &partition);

  Ids extractPoints(const Ids &edgeSimplices, const dSimplices &simplices);

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

  const Partitioner &partitioner;

  static constexpr tCoordinate SAFETY = 100;
  static const uint BASE_CASE = 100;
};

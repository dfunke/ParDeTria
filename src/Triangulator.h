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

template <uint D> class Triangulator {
public:
  Triangulator(dBox<D> &_bounds, dPoints<D> &_points,
               std::unique_ptr<Partitioner<D>> &&_partitioner);

  dSimplices<D> triangulate();

  const TriangulationReport &getTriangulationReport() const {
    return triangulationReport;
  }

private:
  dSimplices<D> triangulateBase(const Ids partitionPoints,
                                const std::string provenance);

  dSimplices<D> triangulateDAC(const Ids partitionPoints,
                               const std::string provenance);

  Ids getEdge(const dSimplices<D> &simplices, const Partition<D> &partition);

  Ids extractPoints(const Ids &edgeSimplices, const dSimplices<D> &simplices,
                    bool ignoreInfinite = false);

  void updateNeighbors(dSimplices<D> &simplices, const std::string &provenance);

  dSimplices<D> mergeTriangulation(std::vector<dSimplices<D>> &partialDTs,
                                   const Ids &edgeSimplices,
                                   const dSimplices<D> &edgeDT,
                                   const Partitioning<D> &partitioning,
                                   const std::string &provenance,
                                   const dSimplices<D> *realDT = nullptr);

  void evaluateVerificationReport(const VerificationReport<D> &vr,
                                  const std::string &provenance) const;

  void
  evaluateCrossCheckReport(const CrossCheckReport<D> &ccr,
                           const std::string &provenance,
                           const dSimplices<D> &DT, const dSimplices<D> &realDT,
                           const Partitioning<D> *partitioning = nullptr) const;

private:
  bool verify = false;
  dBox<D> bounds;
  dPoints<D> points;

  TriangulationReport triangulationReport;

  std::unique_ptr<Partitioner<D>> &partitioner;

  static constexpr tCoordinate SAFETY = 100;
  static constexpr uint BASE_CASE = 100;
  static constexpr char TOP = 0;
};

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

template <uint D, typename Precision> class Triangulator {
public:
  Triangulator(const dBox<D, Precision> &_bounds,
               dPoints<D, Precision> &_points,
               std::unique_ptr<Partitioner<D, Precision>> &&_partitioner);

  dSimplices<D, Precision> triangulate();

  const TriangulationReport &getTriangulationReport() const {
    return triangulationReport;
  }

protected:
  dSimplices<D, Precision> triangulateBase(const Ids partitionPoints,
                                           const std::string provenance);

  dSimplices<D, Precision> triangulateDAC(const Ids partitionPoints,
                                          const std::string provenance);

  Ids getEdge(const dSimplices<D, Precision> &simplices,
              const dBox<D, Precision> &bounds);

  Ids extractPoints(const Ids &edgeSimplices,
                    const dSimplices<D, Precision> &simplices,
                    bool ignoreInfinite = false);

  void updateNeighbors(dSimplices<D, Precision> &simplices,
                       const std::string &provenance);

  dSimplices<D, Precision>
  mergeTriangulation(std::vector<dSimplices<D, Precision>> &partialDTs,
                     const Ids &edgeSimplices,
                     const dSimplices<D, Precision> &edgeDT,
                     const Partitioning<D, Precision> &partitioning,
                     const std::string &provenance,
                     const dSimplices<D, Precision> *realDT = nullptr);

  void evaluateVerificationReport(const VerificationReport<D, Precision> &vr,
                                  const std::string &provenance) const;

  void evaluateCrossCheckReport(
      const CrossCheckReport<D, Precision> &ccr, const std::string &provenance,
      const dSimplices<D, Precision> &DT,
      const dSimplices<D, Precision> &realDT,
      const Partitioning<D, Precision> *partitioning = nullptr) const;

protected:
  const dBox<D, Precision> bounds;
  dPoints<D, Precision> points;
  TriangulationReport triangulationReport;
  std::unique_ptr<Partitioner<D, Precision>> &partitioner;

public:
  static bool isTOP(const std::string &provenance) {
    return provenance == std::to_string(TOP);
  }

  static bool VERIFY;

protected:
  static constexpr Precision SAFETY = 100;
  static constexpr uint BASE_CASE = 100;
  static constexpr char TOP = 0;
};

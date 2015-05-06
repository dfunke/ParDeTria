#pragma once

// stl
#include <functional>
#include <vector>
#include <memory>

// own
#include "Geometry.h"
#include "Partitioner.h"
#include "Triangulator.h"
#include "CGALTriangulator.h"
#include "utils/TBB_Containers.h"

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

typedef tbb::concurrent_vector<TriangulationReportEntry> TriangulationReport;

template<uint D, typename Precision>
class DCTriangulator : public Triangulator<D, Precision> {
public:


    DCTriangulator(const dBox<D, Precision> &_bounds, dPoints<D, Precision> &_points,
                   const uint _recursionDepth,
                   const unsigned char splitter,
                   const uint gridOccupancy = 1,
                   const bool parallelBaseSolver = false);

protected:
    dSimplices<D, Precision> _triangulateBase(const Ids partitionPoints,
                                              const dBox<D, Precision> &bounds,
                                              const std::string provenance);

    dSimplices<D, Precision> _triangulate(const Ids &partitionPoints,
                                          const dBox<D, Precision> &bounds,
                                          const std::string provenance);

    Ids getEdge(dSimplices<D, Precision> &simplices,
                const Partitioning<D, Precision> &partitioning,
                const uint &partition);

    Ids extractPoints(const Ids &edgeSimplices,
                      const dSimplices<D, Precision> &simplices,
                      bool ignoreInfinite = false);

    void updateNeighbors(dSimplices<D, Precision> &simplices, const Ids &toCheck,
                         const std::string &provenance);

    dSimplices<D, Precision>
            mergeTriangulation(std::vector<dSimplices<D, Precision>> &partialDTs,
                               const Ids &edgeSimplices,
                               const dSimplices<D, Precision> &edgeDT,
                               const Partitioning<D, Precision> &partitioning,
                               const std::string &provenance);


protected:
    const uint recursionDepth;

    std::unique_ptr<Partitioner<D, Precision>> partitioner;
    std::unique_ptr<Triangulator<D, Precision>> baseTriangulator;

public:
    static constexpr Precision SAFETY = 100;
    static constexpr uint BASE_CUTOFF = 500;
};

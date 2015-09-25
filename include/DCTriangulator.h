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

    tIdType nPoints;
    tIdType nSimplices;
    tIdType nEdgePoints;
    tIdType nEdgeSimplices;
};

typedef tbb::concurrent_vector<TriangulationReportEntry> TriangulationReport;

template<uint D, typename Precision>
class DCTriangulator : public Triangulator<D, Precision> {
public:


    DCTriangulator(const dBox<D, Precision> &_bounds, dPoints<D, Precision> &_points,
                   const uint _recursionDepth,
                   const unsigned char splitter,
                   const uint gridOccupancy = 1,
                   const bool parallelBaseSolver = false,
                   const bool parallelEdgeTria = true,
                   const bool addInfinitePoints = true);

protected:
    dSimplices<D, Precision> _triangulateBase(const Point_Ids &partitionPoints,
                                              const dBox<D, Precision> &bounds,
                                              const std::string provenance);

    dSimplices<D, Precision> _triangulate(const Point_Ids &partitionPoints,
                                          const dBox<D, Precision> &bounds,
                                          const std::string provenance,
                                          const unsigned char _splitter
    );

    dSimplices<D, Precision> _triangulate(const Point_Ids &partitionPoints,
                                          const dBox<D, Precision> &bounds,
                                          const std::string provenance) {
        return _triangulate(partitionPoints, bounds, provenance, (unsigned char) 0);
    }

    void getEdge(const dSimplices<D, Precision> &simplices,
                 const Partitioning<D, Precision> &partitioning,
                 const uint &partition,
                 Concurrent_Growing_Point_Ids &edgePoints, Concurrent_Simplex_Ids &edgeSimplices);

    cWuFaces buildWhereUsed(const dSimplices<D, Precision> &simplices,
                            const Simplex_Ids &edgeSimplices);

    void updateNeighbors(dSimplices<D, Precision> &simplices,
                         const Simplex_Ids &toCheck,
                         const cWuFaces &wuFaces,
                         const std::string &provenance);

    dSimplices<D, Precision> mergeTriangulation(std::vector<dSimplices<D, Precision>> &&partialDTs,
                                                const Simplex_Ids &edgeSimplices,
                                                dSimplices<D, Precision> &&edgeDT,
                                                const Partitioning<D, Precision> &partitioning,
                                                const std::string &provenance
    );


protected:
    const uint recursionDepth;
    const bool parallelEdgeTria;
    const unsigned char splitter;

    std::unique_ptr<Triangulator<D, Precision>> baseTriangulator;

public:
    static constexpr Precision SAFETY = 100;
    static constexpr uint BASE_CUTOFF = 500;
};

//
// Created by dfunke on 4/8/15.
//

#pragma once

#include "Geometry.h"

#include "utils/VTuneAdapter.h"
#include "datastructures/LP_Set.hxx"
#include "datastructures/Growing_LP.hxx"

struct PartialTriangulation {
    GrowingHashTable<Concurrent_LP_Set> simplices;
    GrowingHashTable<Concurrent_LP_Set> convexHull;

    PartialTriangulation(const std::size_t nSimplices, const std::size_t nConvexHull)
            : simplices(nSimplices), convexHull(nConvexHull) { }

    PartialTriangulation()
            : simplices(1), convexHull(1) { }
};

template<uint D, typename Precision>
class Triangulator {

    template<uint D2, typename Precision2>
    friend
    class DCTriangulator;

public:
    virtual ~Triangulator() = default;

protected:
    Triangulator(const dBox<D, Precision> &_bounds, dPoints<D, Precision> &_points)
            : baseBounds(_bounds), points(_points) { }

public:

    std::pair<dSimplices<D, Precision>, PartialTriangulation> triangulate() {
        dSimplices<D, Precision> DT;
        DT.reserve(points.size());

        PartialTriangulation pt = this->_triangulate(DT, allPoints(), baseBounds, std::to_string(TOP));

        return std::make_pair(std::move(DT), std::move(pt));
    };

protected:
    virtual PartialTriangulation _triangulate(dSimplices<D, Precision> &DT,
                                              const Ids &partitionPoints,
                                              const dBox<D, Precision> &bounds,
                                              const std::string provenance) = 0;

    Ids allPoints() const;

protected:
    const dBox<D, Precision> baseBounds;
    dPoints<D, Precision> &points;

public:
    static bool isTOP(const std::string &provenance) {
        return provenance == std::to_string(TOP);
    }

protected:
    static constexpr char TOP = 0;
};

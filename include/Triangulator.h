//
// Created by dfunke on 4/8/15.
//

#pragma once

#include <datastructures/Bit_Set.hxx>
#include "Geometry.h"

#include "utils/VTuneAdapter.h"
#include "datastructures/LP_Set.hxx"
#include "datastructures/Growing_LP.hxx"

struct PartialTriangulation {
    Simplex_Ids simplices;
    Simplex_Ids convexHull;

    PartialTriangulation(const std::size_t lowerBound, const std::size_t upperBound)
            : simplices(lowerBound, upperBound),
              convexHull(lowerBound, upperBound) { }

    PartialTriangulation()
            : simplices(0, 1), convexHull(0, 1) { }
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
                                              const Point_Ids &partitionPoints,
                                              const dBox<D, Precision> &bounds,
                                              const std::string provenance) = 0;

    Point_Ids allPoints() const;

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

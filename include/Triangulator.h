//
// Created by dfunke on 4/8/15.
//

#pragma once

#include <datastructures/Bit_Set.hxx>
#include "Geometry.h"

#include "utils/VTuneAdapter.h"
#include "datastructures/LP_Set.hxx"
#include "datastructures/Growing_LP.hxx"

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

    dSimplices<D, Precision> triangulate() {

        return this->_triangulate(allPoints(), baseBounds, std::to_string(TOP));
    };

protected:
    virtual dSimplices<D, Precision> _triangulate(const Point_Ids &partitionPoints,
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

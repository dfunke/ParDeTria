//
// Created by dfunke on 4/8/15.
//

#pragma once

#include "Geometry.h"

template<uint D, typename Precision>
class Triangulator {

    template<uint D2, typename Precision2>
    friend
    class DCTriangulator;

protected:
    Triangulator(const dBox<D, Precision> &_bounds, dPoints<D, Precision> &_points)
            : baseBounds(_bounds), points(_points) { }

public:

    dSimplices<D, Precision> triangulate() {
        return this->_triangulate(allPoints(), baseBounds, std::to_string(TOP));
    };

protected:
    virtual dSimplices<D, Precision> _triangulate(const Ids &partitionPoints,
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
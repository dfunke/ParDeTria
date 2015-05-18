//
// Created by dfunke on 4/8/15.
//

#pragma once

#include "Geometry.h"

struct PartialTriangulation {
    Ids simplices;
    Ids convexHull;
};

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
        dSimplices<D, Precision> tmpDT;
        PartialTriangulation pt = this->_triangulate(tmpDT, allPoints(), baseBounds, std::to_string(TOP));

        dSimplices<D, Precision> DT = tmpDT.project(pt.simplices);
        DT.convexHull.insert(pt.convexHull.begin(), pt.convexHull.end());

        return DT;
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
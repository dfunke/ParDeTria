//
// Created by dfunke on 4/8/15.
//

#pragma once

#include "Geometry.h"

#include "utils/VTuneAdapter.h"
namespace LoadBalancing
{
    template<uint D, typename Precision, typename Monitor>
    class DCTriangulator;
}

template<uint D, typename Precision>
class Triangulator {

    template<uint D2, typename Precision2>
    friend
    class DCTriangulator;
    
    template<uint D2, typename Precision2, typename Monitor>
    friend
    class LoadBalancing::DCTriangulator;


public:
    virtual ~Triangulator() = default;

protected:
    Triangulator(const dBox<D, Precision> &_bounds, dPoints<D, Precision> &_points)
            : baseBounds(_bounds), points(_points) { }

public:

    dSimplices<D, Precision> triangulate() {

        return this->_triangulate(allPoints(), baseBounds, std::to_string(TOP));
    }

    template <typename CONTAINER>
    dSimplices<D, Precision> triangulate(const CONTAINER & a){
        Point_Ids idSet(a.size());
        for(const auto & id : a){
            idSet.insert(id);
        }

        return _triangulate(idSet, baseBounds, std::to_string(TOP));
    }

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

    static constexpr char TOP = 0;
};

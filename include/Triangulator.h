//
// Created by dfunke on 4/8/15.
//

#pragma once

#include "Geometry.h"

template <uint D, typename Precision>
class Triangulator {

protected:
    Triangulator(const dBox<D, Precision> &_bounds, dPoints<D, Precision> &_points)
            : baseBounds(_bounds), points(_points) { }

public:

    dSimplices<D, Precision> triangulate() {
        return this->_triangulate(allPoints(), baseBounds, std::to_string(TOP));
    };

    template<typename... Args>
    static std::unique_ptr<Triangulator<D, Precision>> make(const unsigned char type, Args &... args);

protected:
    virtual dSimplices<D, Precision> _triangulate(const Ids & partitionPoints,
                                         const dBox<D, Precision> &bounds,
                                         const std::string provenance) = 0;

    Ids allPoints() const;

protected:
    const dBox<D, Precision> baseBounds;
    dPoints<D, Precision> & points;

public:
    static bool isTOP(const std::string &provenance) {
        return provenance == std::to_string(TOP);
    }

    static bool VERIFY;

protected:
    static constexpr char TOP = 0;
};

// forward declare triangulator types
template <uint D, typename Precision> class DCTriangulator;
template <uint D, typename Precision, bool Parallel> class CGALTriangulator;

template<uint D, typename Precision> template <typename... Args>
std::unique_ptr<Triangulator<D, Precision>> Triangulator<D, Precision>::make
        (const unsigned char type, Args &... args) {

    std::unique_ptr<Triangulator<D, Precision>> triangulator_ptr;
    switch (type) {
        case 'd':
            triangulator_ptr = std::make_unique<DCTriangulator<D, Precision>>(args...);
            break;
        case 'c':
            triangulator_ptr = std::make_unique<CGALTriangulator<D, Precision, false>>(args...);
            break;
        case 'm':
            triangulator_ptr = std::make_unique<CGALTriangulator<D, Precision, true>>(args...);
            break;
    }

    ASSERT(triangulator_ptr);
    return triangulator_ptr;
}
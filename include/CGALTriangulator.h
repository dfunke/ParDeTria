#pragma once

#include "Geometry.h"
#include "Triangulator.h"
#include "utils/IndexHandler.h"

template<uint D, typename Precision,
        bool Parallel = false>
class CGALTriangulator : public Triangulator<D, Precision> {

public:
    CGALTriangulator(const dBox<D, Precision> &_bounds, dPoints<D, Precision> &_points,
                     const uint _gridOccupancy = 1);

    dSimplices<D, Precision> triangulateSome(const Point_Ids& ids, const dBox<D, Precision>& bounds) {
        return this->_triangulate(ids, bounds, std::to_string(Triangulator<D, Precision>::TOP));
    }

protected:
    dSimplices<D, Precision> _triangulate(const Point_Ids &ids,
                                          const dBox<D, Precision> &bounds,
                                          const std::string provenance
            /*, bool filterInfinite = false */);

protected:
    const uint gridOccupancy;
};

template<uint D, typename Precision,
        bool Parallel = false>
class PureCGALTriangulator : public Triangulator<D, Precision> {

public:
    PureCGALTriangulator(const dBox<D, Precision> &_bounds, dPoints<D, Precision> &_points,
                         const uint _gridOccupancy = 1);

protected:
    dSimplices<D, Precision> _triangulate(const Point_Ids &ids,
                                          const dBox<D, Precision> &bounds,
                                          const std::string provenance
            /*, bool filterInfinite = false */);

protected:
    const uint gridOccupancy;
};

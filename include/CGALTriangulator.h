#pragma once

#include "Geometry.h"
#include "Triangulator.h"

template<uint D, typename Precision,
        bool Parallel = false>
class CGALTriangulator : public Triangulator<D, Precision> {

public:
    CGALTriangulator(const dBox<D, Precision> &_bounds, dPoints<D, Precision> &_points,
                     const uint _gridOccupancy = 1);

protected:
    PartialTriangulation _triangulate(dSimplices<D, Precision> &DT,
                                      const Ids &ids,
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
    PartialTriangulation _triangulate(dSimplices<D, Precision> &DT,
                                      const Ids &ids,
                                      const dBox<D, Precision> &bounds,
                                      const std::string provenance
            /*, bool filterInfinite = false */);

protected:
    const uint gridOccupancy;
};
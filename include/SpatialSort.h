#pragma once

#include "Geometry.h"

template<uint D, typename Precision>
class SpatialSorter {

public:
    void sort(dPoints<D, Precision>, const Ids &ids
            /*, bool filterInfinite = false */);
};

template<uint D, typename Precision>
class CGALSpatialSorter : public SpatialSorter<D, Precision> {

public:
    void sort(dPoints<D, Precision> & points
            /*, bool filterInfinite = false */);
};
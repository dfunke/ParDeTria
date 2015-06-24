#pragma once

#include "Geometry.h"

template<uint D, typename Precision>
class SpatialSorter {

public:
    virtual void sort(dPoints<D, Precision> &points
            /*, bool filterInfinite = false */) = 0;
};

template<uint D, typename Precision>
class CGALSpatialSorter : public SpatialSorter<D, Precision> {

public:
    void sort(dPoints<D, Precision> &points
            /*, bool filterInfinite = false */);
};
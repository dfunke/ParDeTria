#pragma once

#include "Geometry.h"

template <uint D>
dSimplices<D> delaunayCgal(dPoints<D> &points, const Ids *ids = nullptr,
                           bool filterInfinite = false);

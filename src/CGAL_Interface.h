#pragma once

#include "Geometry.h"

template <uint D, typename Precision> class CGALInterface {
public:
  static dSimplices<D, Precision> triangulate(const Ids &ids,
                                              dPoints<D, Precision> &points,
                                              bool filterInfinite = false);
};

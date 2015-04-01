#pragma once

#include "Geometry.h"

template <uint D, typename Precision> class CGALInterface {
public:
  static dSimplices<D, Precision> triangulate(const Ids &ids,
                                              dPoints<D, Precision> &points,
                                              const dBox<D, Precision> &bounds,
                                              bool filterInfinite = false);
};

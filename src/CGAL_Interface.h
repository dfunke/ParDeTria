#pragma once

#include "Geometry.h"

template <uint D, typename Precision> class CGALInterface {
public:
  static dSimplices<D, Precision> triangulate(dPoints<D, Precision> &points,
                                              const Ids *ids = nullptr,
                                              bool filterInfinite = false);
};

#include "Painter.h"
#include "Partitioner.h"

#include <valgrind/callgrind.h>

// static variables

constexpr uint DataHolder<3>::PADDING;
constexpr uint DataHolder<3>::FONT_SIZE;

// TODO remove
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"

template <> tCoordinate Painter<3>::imgDim(const uint dim) const {
  return data.img.dim[dim] - data.img.coords[dim];
}

template <>
tCoordinate Painter<3>::translatePoint(tCoordinate in, uint dim) const {
  return ((data.offset[dim] + in - data.bounds.coords[dim]) /
          (3 * data.bounds.dim[dim])) *
         imgDim(dim);
}

template <>
tCoordinate Painter<3>::translateLength(tCoordinate in, uint dim) const {
  return (in / (3 * data.bounds.dim[dim])) * imgDim(dim);
}

template <> template <typename Object> void Painter<3>::_log(const Object &o) {
  if (logging) {
    // TODO logging
  }
}

template <>
void Painter<3>::setColor(tCoordinate r, tCoordinate g, tCoordinate b,
                          tCoordinate alpha) {
  // TODO set color
}

//
template <> void Painter<3>::draw(const dPoint<3> &point, bool drawInfinite) {

  // TODO draw point

  _log(point);
}

template <>
void Painter<3>::draw(const dSimplex<3> &simplex, const dPoints<3> &points,
                      bool drawInfinite) {

  // TODO draw simplex

  _log(simplex);
}

template <>
void Painter<3>::drawCircumSphere(const dSimplex<3> &s,
                                  const dPoints<3> &points, bool drawInfinite) {

  // TODO draw circumsphere
}

template <> void Painter<3>::drawPartition(const dPoints<3> &points) {
  // TODO draw partition
}

template <>
void Painter<3>::drawNeighbors(const dSimplex<3> &simplex,
                               const dSimplices<3> &neighbors,
                               const dPoints<3> &points, bool drawInfinite) {
  for (uint n : simplex.neighbors) {
    if (dSimplex<3>::isFinite(n) && neighbors.contains(n)) {
      draw(neighbors[n], points, drawInfinite);
    }
  }
}

template <>
void Painter<3>::save(
#ifdef NO_OUTPUT
    __attribute((unused))
#endif
    const std::string &file) const {
// stop callgrind instrumentation for saving png

#ifndef NO_OUTPUT // disable png output
  CALLGRIND_STOP_INSTRUMENTATION;

  // TODO save painting

  CALLGRIND_START_INSTRUMENTATION;
#endif
}

template <> void Painter<3>::_init(const dBox<3> &_bounds, uint _resolution) {
  data.bounds = _bounds;
  logging = false;
  dashed = false;
  data.logLine = 1;

  for (uint d = 0; d < 3; ++d) {
    data.offset[d] = data.bounds.dim[d]; // offset bounds in middle of image

    data.img.coords[d] = 0;                                 // image starts at 0
    data.img.dim[d] = 3 * data.bounds.dim[d] * _resolution; // 9 quadrants
  }

  // TODO init painter
}

template <> void Painter<3>::_copy(const Painter<3> &a) {
  data.bounds = a.data.bounds;
  data.img = a.data.img;
  data.offset = a.data.offset;
  logging = a.logging;
  data.logLine = a.data.logLine;

  // TODO copy painting
}

// TODO remove
#pragma clang diagnostic pop

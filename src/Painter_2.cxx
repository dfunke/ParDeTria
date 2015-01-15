#include "Painter.h"
#include "Partitioner.h"

#include <valgrind/callgrind.h>

// static variables

constexpr uint DataHolder<2>::PADDING;
constexpr uint DataHolder<2>::FONT_SIZE;

const int D = 2;

template <> tCoordinate Painter<2>::imgDim(const uint dim) const {
  return data.img.dim[dim] - data.img.coords[dim];
}

template <>
tCoordinate Painter<2>::translatePoint(tCoordinate in, uint dim) const {
  return ((data.offset[dim] + in - data.bounds.coords[dim]) /
          (3 * data.bounds.dim[dim])) *
         imgDim(dim);
}

template <>
tCoordinate Painter<2>::translateLength(tCoordinate in, uint dim) const {
  return (in / (3 * data.bounds.dim[dim])) * imgDim(dim);
}

template <> template <typename Object> void Painter<2>::_log(const Object &o) {
  if (logging) {
    data.cr->move_to(translatePoint(data.bounds.dim[0] + data.PADDING, 0),
                     translatePoint(data.bounds.coords[1], 1) +
                         1.5 * data.FONT_SIZE * data.logLine++);
    data.cr->show_text(to_string(o));
  }
}

template <>
void Painter<2>::setColor(tCoordinate r, tCoordinate g, tCoordinate b,
                          tCoordinate alpha) {
  data.cr->set_source_rgba(r, g, b, alpha);
}

//
template <> void Painter<2>::draw(const dPoint<2> &point) {

  if (!point.isFinite()) {
    data.cr->save();
    data.cr->set_source_rgb(1, 1, 0);
  }

  data.cr->arc(translatePoint(point.coords[0], 0),
               translatePoint(point.coords[1], 1), 5, 0, 2 * M_PI);
  data.cr->fill();

  data.cr->move_to(translatePoint(point.coords[0], 0) + 7,
                   translatePoint(point.coords[1], 1) + 7);
  data.cr->set_font_size(10);
  data.cr->show_text(std::to_string(point.id));
  data.cr->set_font_size(data.FONT_SIZE);

  if (!point.isFinite())
    data.cr->restore();

  _log(point);
}

template <> void Painter<2>::draw(const dPoints<2> &points) {

  for (const auto &p : points)
    draw(p);
}

template <>
void Painter<2>::draw(const dSimplex<2> &simplex, const dPoints<2> &points,
                      bool drawInfinite) {

  auto line = [&](uint a, uint b) {

    auto A = points[a];
    auto B = points[b];

    data.cr->save();

    if (!A.isFinite() || !B.isFinite()) {
      if (!drawInfinite)
        return;

      data.cr->set_dash(std::vector<double>({2, 2}), 0);
    }

    if (dashed)
      data.cr->set_dash(std::vector<double>({2, 2}), 0);
    else
      data.cr->unset_dash();

    data.cr->move_to(translatePoint(A.coords[0], 0),
                     translatePoint(A.coords[1], 1));
    data.cr->line_to(translatePoint(B.coords[0], 0),
                     translatePoint(B.coords[1], 1));
    data.cr->stroke();

    data.cr->restore();
  };

  for (uint d = 0; d < D; ++d) {
    line(simplex.vertices[d], simplex.vertices[d + 1]);
  }
  // close the loop
  line(simplex.vertices[D], simplex.vertices[0]);

  _log(simplex);
}

template <>
void Painter<2>::draw(const dSimplices<2> &simplices, const dPoints<2> &points,
                      bool drawInfinite) {

  for (const auto &s : simplices)
    draw(s, points, drawInfinite);
}

template <>
void Painter<2>::drawCircumSphere(const dSimplex<2> &s,
                                  const dPoints<2> &points, bool drawInfinite) {

  if (!drawInfinite && !s.isFinite())
    return;

  auto circumcircle = s.circumsphere(points);

  data.cr->save();

  // draw circumcenter
  data.cr->arc(translatePoint(circumcircle.center[0], 0),
               translatePoint(circumcircle.center[1], 1), 5, 0, 2 * M_PI);
  data.cr->fill();

  if (dashed)
    data.cr->set_dash(std::vector<double>({2, 2}), 0);
  else
    data.cr->unset_dash();

  data.cr->arc(translatePoint(circumcircle.center[0], 0),
               translatePoint(circumcircle.center[1], 1),
               translateLength(circumcircle.radius, 0), 0, 2 * M_PI);
  data.cr->stroke();

  data.cr->restore();
}

template <>
void Painter<2>::drawCircumSphere(const dSimplices<2> &s,
                                  const dPoints<2> &points, bool drawInfinite) {
  for (const auto &x : s)
    drawCircumSphere(x, points, drawInfinite);
}

template <> void Painter<2>::drawPartition(const dPoints<2> &points) {
  auto stats = getPointStats(points.begin_keys(), points.end_keys(), points);

  // draw partition borders
  data.cr->move_to(translatePoint(stats.mid.coords[0], 0), 0);
  data.cr->line_to(translatePoint(stats.mid.coords[0], 0), imgDim(1));

  data.cr->move_to(0, translatePoint(stats.mid.coords[1], 1));
  data.cr->line_to(imgDim(0), translatePoint(stats.mid.coords[1], 1));

  data.cr->stroke();
}

template <>
void Painter<2>::drawNeighbors(const dSimplex<2> &simplex,
                               const dSimplices<2> &neighbors,
                               const dPoints<2> &points, bool drawInfinite) {
  for (uint n : simplex.neighbors) {
    if (dSimplex<2>::isFinite(n) && neighbors.contains(n)) {
      draw(neighbors[n], points, drawInfinite);
    }
  }
}

template <>
void Painter<2>::drawNeighbors(const dSimplices<2> &simplices,
                               const dSimplices<2> &neighbors,
                               const dPoints<2> &points, bool drawInfinite) {
  for (const auto &x : simplices)
    drawNeighbors(x, neighbors, points, drawInfinite);
}

template <>
void Painter<2>::save(
#ifdef NO_OUTPUT
    __attribute((unused))
#endif
    const std::string &file) const {
// stop callgrind instrumentation for saving png

#ifndef NO_OUTPUT // disable png output
  CALLGRIND_STOP_INSTRUMENTATION;

  data.cs->write_to_png(file);

  CALLGRIND_START_INSTRUMENTATION;
#endif
}

template <> void Painter<2>::_init(const dBox<2> &_bounds, uint _resolution) {
  data.bounds = _bounds;
  logging = false;
  dashed = false;
  data.logLine = 1;

  for (uint d = 0; d < D; ++d) {
    data.offset[d] = data.bounds.dim[d]; // offset bounds in middle of image

    data.img.coords[d] = 0;                                 // image starts at 0
    data.img.dim[d] = 3 * data.bounds.dim[d] * _resolution; // 9 quadrants
  }

  data.cs =
      Cairo::ImageSurface::create(Cairo::FORMAT_ARGB32, imgDim(0), imgDim(1));
  data.cr = Cairo::Context::create(data.cs);

  // draw background white
  data.cr->set_source_rgb(1, 1, 1);
  data.cr->paint();

  data.cr->set_line_width(1.0);
  data.cr->set_source_rgb(0, 0, 0);

  // set font options
  data.cr->select_font_face("serif", Cairo::FONT_SLANT_NORMAL,
                            Cairo::FONT_WEIGHT_NORMAL);
  data.cr->set_font_size(data.FONT_SIZE);

  // draw bounds
  data.cr->rectangle(translatePoint(data.bounds.coords[0], 0),
                     translatePoint(data.bounds.coords[1], 1),
                     translateLength(data.bounds.dim[0], 0),
                     translateLength(data.bounds.dim[1], 1));
  data.cr->stroke();
}

template <> void Painter<2>::_copy(const Painter<2> &a) {
  data.bounds = a.data.bounds;
  data.img = a.data.img;
  data.offset = a.data.offset;
  logging = a.logging;
  data.logLine = a.data.logLine;

  data.cs =
      Cairo::ImageSurface::create(Cairo::FORMAT_ARGB32, imgDim(0), imgDim(1));
  data.cr = Cairo::Context::create(data.cs);

  // set font options
  data.cr->select_font_face("serif", Cairo::FONT_SLANT_NORMAL,
                            Cairo::FONT_WEIGHT_NORMAL);
  data.cr->set_font_size(data.FONT_SIZE);

  // copy a's surface
  data.cr->save();
  data.cr->set_source(a.data.cs, 0, 0);
  data.cr->paint();
  data.cr->restore();
}

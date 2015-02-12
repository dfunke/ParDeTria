#include "Painter.h"
#include "Partitioner.h"

#include <valgrind/callgrind.h>

// static variables

template <typename Precision>
constexpr uint PainterImplementation<2, Precision>::PADDING;
template <typename Precision>
constexpr uint PainterImplementation<2, Precision>::FONT_SIZE;

template <typename Precision>
float PainterImplementation<2, Precision>::imgDim(const uint dim) {
  return img.high[dim] - img.low[dim];
}

template <typename Precision>
float PainterImplementation<2, Precision>::translatePoint(float in, uint dim) {
  return ((offset[dim] + in - bounds.low[dim]) /
          (3 * (bounds.high[dim] - bounds.low[dim]))) *
         imgDim(dim);
}

template <typename Precision>
float PainterImplementation<2, Precision>::translateLength(float in, uint dim) {
  return (in / (3 * (bounds.high[dim] - bounds.low[dim]))) * imgDim(dim);
}

template <typename Precision>
template <typename Object>
void PainterImplementation<2, Precision>::_log(const Object &o) {
  if (logging) {
    cr->move_to(translatePoint(bounds.high[0] + PADDING, 0),
                translatePoint(bounds.low[1], 1) + 1.5 * FONT_SIZE * logLine++);
    cr->show_text(to_string(o));
  }
}

template <typename Precision>
void PainterImplementation<2, Precision>::setColor(float r, float g, float b,
                                                   float alpha) {
  cr->set_source_rgba(r, g, b, alpha);
}

//
template <typename Precision>
void
PainterImplementation<2, Precision>::draw(const dPoint<2, Precision> &point,
                                          bool drawInfinite) {

  if (!point.isFinite()) {
    if (!drawInfinite)
      return;

    cr->save();
    cr->set_source_rgb(1, 1, 0);
  }

  cr->arc(translatePoint(point.coords[0], 0),
          translatePoint(point.coords[1], 1), 5, 0, 2 * M_PI);
  cr->fill();

  cr->move_to(translatePoint(point.coords[0], 0) + 7,
              translatePoint(point.coords[1], 1) + 7);
  cr->set_font_size(10);
  cr->show_text(std::to_string(point.id));
  cr->set_font_size(FONT_SIZE);

  if (!point.isFinite())
    cr->restore();

  _log(point);
}

template <typename Precision>
void
PainterImplementation<2, Precision>::draw(const dSimplex<2, Precision> &simplex,
                                          const dPoints<2, Precision> &points,
                                          bool drawInfinite) {

  auto line = [&](uint a, uint b) {

    auto A = points[a];
    auto B = points[b];

    cr->save();

    if (!A.isFinite() || !B.isFinite()) {
      if (!drawInfinite)
        return;

      cr->set_dash(std::vector<double>({2, 2}), 0);
    }

    if (dashed)
      cr->set_dash(std::vector<double>({2, 2}), 0);
    else
      cr->unset_dash();

    cr->move_to(translatePoint(A.coords[0], 0), translatePoint(A.coords[1], 1));
    cr->line_to(translatePoint(B.coords[0], 0), translatePoint(B.coords[1], 1));
    cr->stroke();

    cr->restore();
  };

  for (uint d = 0; d < 2; ++d) {
    line(simplex.vertices[d], simplex.vertices[d + 1]);
  }
  // close the loop
  line(simplex.vertices[2], simplex.vertices[0]);

  _log(simplex);
}

template <typename Precision>
void PainterImplementation<2, Precision>::drawCircumSphere(
    const dSimplex<2, Precision> &s, const dPoints<2, Precision> &points,
    bool drawInfinite) {

  if (!drawInfinite && !s.isFinite())
    return;

  auto circumcircle = s.circumsphere(points);

  cr->save();

  // draw circumcenter
  cr->arc(translatePoint(circumcircle.center[0], 0),
          translatePoint(circumcircle.center[1], 1), 5, 0, 2 * M_PI);
  cr->fill();

  if (dashed)
    cr->set_dash(std::vector<double>({2, 2}), 0);
  else
    cr->unset_dash();

  cr->arc(translatePoint(circumcircle.center[0], 0),
          translatePoint(circumcircle.center[1], 1),
          translateLength(circumcircle.radius, 0), 0, 2 * M_PI);
  cr->stroke();

  cr->restore();
}

template <typename Precision>
void PainterImplementation<2, Precision>::drawPartition(
    const dPoints<2, Precision> &points) {
  auto stats = getPointStats(points.begin_keys(), points.end_keys(), points);

  // draw partition borders
  cr->move_to(translatePoint(stats.mid.coords[0], 0), 0);
  cr->line_to(translatePoint(stats.mid.coords[0], 0), imgDim(1));

  cr->move_to(0, translatePoint(stats.mid.coords[1], 1));
  cr->line_to(imgDim(0), translatePoint(stats.mid.coords[1], 1));

  cr->stroke();
}

template <typename Precision>
void PainterImplementation<2, Precision>::save(const std::string &file) const {
  // stop callgrind instrumentation for saving png

  CALLGRIND_STOP_INSTRUMENTATION;

  cs->write_to_png(file + ".png");

  CALLGRIND_START_INSTRUMENTATION;
}

template <typename Precision>
void
PainterImplementation<2, Precision>::_init(const dBox<2, Precision> &_bounds,
                                           uint _resolution) {
  bounds = _bounds;
  logging = false;
  dashed = false;
  logLine = 1;

  for (uint d = 0; d < 2; ++d) {
    offset[d] =
        bounds.high[d] - bounds.low[d]; // offset bounds in middle of image

    img.low[d] = 0; // image starts at 0
    img.high[d] =
        img.low[d] +
        3 * (bounds.high[d] - bounds.low[d]) * _resolution; // 9 quadrants
  }

  cs = Cairo::ImageSurface::create(Cairo::FORMAT_ARGB32, imgDim(0), imgDim(1));
  cr = Cairo::Context::create(cs);

  // draw background white
  cr->set_source_rgb(1, 1, 1);
  cr->paint();

  cr->set_line_width(1.0);
  cr->set_source_rgb(0, 0, 0);

  // set font options
  cr->select_font_face("serif", Cairo::FONT_SLANT_NORMAL,
                       Cairo::FONT_WEIGHT_NORMAL);
  cr->set_font_size(FONT_SIZE);

  // draw bounds
  cr->rectangle(translatePoint(bounds.low[0], 0),
                translatePoint(bounds.low[1], 1),
                translateLength(bounds.high[0] - bounds.low[0], 0),
                translateLength(bounds.high[1] - bounds.low[1], 1));
  cr->stroke();
}

template <typename Precision>
void
PainterImplementation<2, Precision>::_copy(const Painter<2, Precision> &a) {
  bounds = a.impl.bounds;
  img = a.impl.img;
  offset = a.impl.offset;
  logging = a.impl.logging;
  logLine = a.impl.logLine;

  cs = Cairo::ImageSurface::create(Cairo::FORMAT_ARGB32, imgDim(0), imgDim(1));
  cr = Cairo::Context::create(cs);

  // set font options
  cr->select_font_face("serif", Cairo::FONT_SLANT_NORMAL,
                       Cairo::FONT_WEIGHT_NORMAL);
  cr->set_font_size(FONT_SIZE);

  // copy a's surface
  cr->save();
  cr->set_source(a.impl.cs, 0, 0);
  cr->paint();
  cr->restore();
}
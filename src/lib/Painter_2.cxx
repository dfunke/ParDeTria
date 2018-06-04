#include "Painter.h"
#include "Painter_impl.hxx"
#include "Partitioner.h"

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
void PainterImplementation<2, Precision>::stroke() {
  while (strokeMtx.test_and_set(std::memory_order_acquire)) // acquire lock
    ;                                                       // spin
  cr->stroke();

  strokeMtx.clear(std::memory_order_release); // release lock
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

template <typename Precision>
void PainterImplementation<2, Precision>::drawLine(const dVector<2, Precision>& a, const dVector<2, Precision>& b) {
  cr->move_to(translatePoint(a[0], 0),
	      translatePoint(a[1], 1));
  cr->line_to(translatePoint(b[0], 0),
	      translatePoint(b[1], 1));
  stroke();
}

template <typename Precision>
void PainterImplementation<2, Precision>::drawBox(const dBox<2, Precision>& bounds) {
  cr->rectangle(translatePoint(bounds.low[0], 0),
                translatePoint(bounds.low[1], 1),
                translateLength(bounds.high[0] - bounds.low[0], 0),
                translateLength(bounds.high[1] - bounds.low[1], 1));
  stroke();
}

template <typename Precision>
void PainterImplementation<2, Precision>::drawText(const std::string& text,
                                                   dVector<2, Precision> coords) {
  cr->move_to(translatePoint(coords[0], 0),
              translatePoint(coords[1], 1));
  cr->set_font_size(20);
  cr->show_text(text);
  cr->set_font_size(FONT_SIZE);
}

//
template <typename Precision>
void PainterImplementation<2, Precision>::draw(
    const dPoint<2, Precision> &point, bool /*drawInfinite*/) {

  cr->arc(translatePoint(point.coords[0], 0),
          translatePoint(point.coords[1], 1), 5, 0, 2 * M_PI);
  cr->fill();

  cr->move_to(translatePoint(point.coords[0], 0) + 7,
              translatePoint(point.coords[1], 1) + 7);
  cr->set_font_size(10);
  //cr->show_text(std::to_string(point.id));
  cr->set_font_size(FONT_SIZE);

  _log(point);
}

template <typename Precision>
void PainterImplementation<2, Precision>::draw(
    const dSimplex<2, Precision> &simplex, const dPoints<2, Precision> &points,
    bool /*drawInfinite*/) {
	
  using sizeT = typename dPoints<2, Precision>::vector::size_type;
  auto line = [&](sizeT a, sizeT b) {

    const auto &A = points[a];
    const auto &B = points[b];

    cr->save();

    if (dashed)
      cr->set_dash(std::vector<double>({2, 2}), 0);

    cr->move_to(translatePoint(A.coords[0], 0), translatePoint(A.coords[1], 1));
    cr->line_to(translatePoint(B.coords[0], 0), translatePoint(B.coords[1], 1));
    stroke();

    cr->unset_dash();
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

  cr->arc(translatePoint(circumcircle.center[0], 0),
          translatePoint(circumcircle.center[1], 1),
          translateLength(circumcircle.radius, 0), 0, 2 * M_PI);
  stroke();

  cr->unset_dash();
  cr->restore();
}

template <typename Precision>
void PainterImplementation<2, Precision>::drawPartition(
    const dPoints<2, Precision> &/*points*/) {
    assert(false);
  /*auto stats = getPointStats(std::size_t(0), points.size(), points);

  // draw partition borders
  cr->move_to(translatePoint(stats.mid[0], 0), 0);
  cr->line_to(translatePoint(stats.mid[0], 0), imgDim(1));

  //cr->move_to(0, translatePoint(stats.mid[1], 1));
  //cr->line_to(imgDim(0), translatePoint(stats.mid[1], 1));

  stroke();*/
}

template <typename Precision>
void PainterImplementation<2, Precision>::save(const std::string &file) const {
  // stop callgrind instrumentation for saving png

  //CALLGRIND_STOP_INSTRUMENTATION;

  cs->write_to_png(file + ".png");

  //CALLGRIND_START_INSTRUMENTATION;
}

template <typename Precision>
void PainterImplementation<2, Precision>::_init(
    const dBox<2, Precision> &_bounds, uint _resolution) {
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
  stroke();
}

template <typename Precision>
void PainterImplementation<2, Precision>::_copy(
#ifdef NDEBUG
    __attribute__((unused))
#endif
    const Painter<2, Precision> &a) {

#ifndef NDEBUG
  bounds = a.impl.bounds;
  img = a.impl.img;
  offset = a.impl.offset;
  logging = a.impl.logging;
  dashed = a.impl.dashed;
  logLine = a.impl.logLine;
#endif

  cs = Cairo::ImageSurface::create(Cairo::FORMAT_ARGB32, imgDim(0), imgDim(1));
  cr = Cairo::Context::create(cs);

  // set font options
  cr->select_font_face("serif", Cairo::FONT_SLANT_NORMAL,
                       Cairo::FONT_WEIGHT_NORMAL);
  cr->set_font_size(FONT_SIZE);

  // copy a's surface
  cr->save();

#ifndef NDEBUG
  cr->set_source(a.impl.cs, 0, 0);
#endif

  while (strokeMtx.test_and_set(std::memory_order_acquire)) // acquire lock
    ;                                                       // spin
  cr->paint();

  strokeMtx.clear(std::memory_order_release); // release lock

  cr->restore();
}

// specializations

template struct PainterImplementation<2, float>;
template struct PainterImplementation<2, double>;

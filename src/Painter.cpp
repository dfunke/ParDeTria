#include "Painter.h"
#include "Partitioner.h"

// static variables

constexpr uint Painter::PADDING;
constexpr uint Painter::FONT_SIZE;

//

void Painter::draw(const dPoint &point) {

  if (!point.isFinite()) {
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

void Painter::draw(const dPoints &points) {

  for (const auto &p : points)
    draw(p);
}

void Painter::draw(const dSimplex &simplex, const dPoints &points,
                   bool drawInfinite) {

  for (uint d = 0; d < D; ++d) {
    line(points, simplex.vertices[d], simplex.vertices[d + 1], drawInfinite);
  }
  // close the loop
  line(points, simplex.vertices[D], simplex.vertices[0], drawInfinite);

  cr->stroke();

  _log(simplex);
}

void Painter::draw(const dSimplices &simplices, const dPoints &points,
                   bool drawInfinite) {

  for (const auto &s : simplices)
    draw(s, points, drawInfinite);
}

void Painter::drawCircumCircle(const dSimplex &s, const dPoints &points,
                               bool drawInfinite) {

  if (!drawInfinite && !s.isFinite())
    return;

  auto circumcircle = s.circumcircle(points);

  // draw circumcenter
  cr->arc(translatePoint(circumcircle.center[0], 0),
          translatePoint(circumcircle.center[1], 1), 5, 0, 2 * M_PI);
  cr->fill();

  cr->arc(translatePoint(circumcircle.center[0], 0),
          translatePoint(circumcircle.center[1], 1),
          translateLength(circumcircle.radius, 0), 0, 2 * M_PI);
  cr->stroke();
}

void Painter::drawCircumCircle(const dSimplices &s, const dPoints &points,
                               bool drawInfinite) {
  for (const auto &x : s)
    drawCircumCircle(x, points, drawInfinite);
}

void Painter::drawPartition(const dPoints &points) {
  auto stats = getPointStats(points.begin_keys(), points.end_keys(), points);

  // draw partition borders
  cr->move_to(translatePoint(stats.mid.coords[0], 0), 0);
  cr->line_to(translatePoint(stats.mid.coords[0], 0), height());

  cr->move_to(0, translatePoint(stats.mid.coords[1], 1));
  cr->line_to(width(), translatePoint(stats.mid.coords[1], 1));

  cr->stroke();
}

void Painter::drawNeighbors(const dSimplex &simplex,
                            const dSimplices &neighbors, const dPoints &points,
                            bool drawInfinite) {
  for (uint n : simplex.neighbors) {
    if (dSimplex::isFinite(n) && neighbors.contains(n)) {
      draw(neighbors[n], points, drawInfinite);
    }
  }
}

void Painter::drawNeighbors(const dSimplices &simplices,
                            const dSimplices &neighbors, const dPoints &points,
                            bool drawInfinite) {
  for (const auto &x : simplices)
    drawNeighbors(x, neighbors, points, drawInfinite);
}

void Painter::_init(const dBox &_bounds, uint _resolution) {
  bounds = _bounds;
  logging = false;
  logLine = 1;

  for (uint d = 0; d < D; ++d) {
    offset.coords[d] = bounds.dim[d]; // offset bounds in middle of image

    img.coords[d] = 0;                            // image starts at 0
    img.dim[d] = 3 * bounds.dim[d] * _resolution; // 9 quadrants
  }

  cs = Cairo::ImageSurface::create(Cairo::FORMAT_ARGB32, width(), height());
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
  cr->rectangle(
      translatePoint(bounds.coords[0], 0), translatePoint(bounds.coords[1], 1),
      translateLength(bounds.dim[0], 0), translateLength(bounds.dim[1], 1));
  cr->stroke();
}

void Painter::_copy(const Painter &a) {
  bounds = a.bounds;
  img = a.img;
  offset = a.offset;
  logging = a.logging;
  logLine = a.logLine;

  cs = Cairo::ImageSurface::create(Cairo::FORMAT_ARGB32, width(), height());
  cr = Cairo::Context::create(cs);

  // set font options
  cr->select_font_face("serif", Cairo::FONT_SLANT_NORMAL,
                       Cairo::FONT_WEIGHT_NORMAL);
  cr->set_font_size(FONT_SIZE);

  // copy a's surface
  cr->save();
  cr->set_source(a.cs, 0, 0);
  cr->paint();
  cr->restore();
}

void Painter::line(const dPoints &points, uint a, uint b, bool drawInfinite) {

  auto A = points[a];
  auto B = points[b];

  cr->save();

  if (!A.isFinite() || !B.isFinite()) {
    if (!drawInfinite)
      return;

    setLineDash();
  }

  cr->move_to(translatePoint(A.coords[0], 0), translatePoint(A.coords[1], 1));
  cr->line_to(translatePoint(B.coords[0], 0), translatePoint(B.coords[1], 1));
  cr->stroke();

  cr->restore();
}

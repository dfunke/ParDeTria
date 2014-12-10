#pragma once

#include <cairo/cairo.h>
#include <cairomm/cairomm.h>

#include <algorithm>

#include "Geometry.h"
#include "Logger.h"

class Painter {

public:
  typedef std::tuple<tCoordinate, tCoordinate, tCoordinate> tRGB;

public:
  Painter() {
    // non working version;
  }

  Painter(const dBox &_bounds, uint _resolution = 10) {

    _init(_bounds, _resolution);
  }

  Painter(const Painter &a) { _copy(a); }

  Painter &operator=(const Painter &a) {

    _copy(a);

    return *this;
  }

  void draw(const dPoint &point);
  void draw(const dPoints &points);

  void draw(const dSimplex &simplex, const dPoints &points,
            bool drawInfinite = false);
  void draw(const dSimplices &simplices, const dPoints &points,
            bool drawInfinite = false);

  void drawCircumCircle(const dSimplex &simplex, const dPoints &points,
                        bool drawInfinite = false);
  void drawCircumCircle(const dSimplices &simplices, const dPoints &points,
                        bool drawInfinite = false);

  void drawPartition(const dPoints &points);

  void drawNeighbors(const dSimplex &simplex, const dSimplices &neighbors,
                     const dPoints &points, bool drawInfinite = false);
  void drawNeighbors(const dSimplices &simplices, const dSimplices &neighbors,
                     const dPoints &points, bool drawInfinite = false);

  void setColor(const tRGB &rgb, tCoordinate alpha = 1.0) {
    setColor(std::get<0>(rgb), std::get<1>(rgb), std::get<2>(rgb), alpha);
  }

  void setColor(tCoordinate r, tCoordinate g, tCoordinate b,
                tCoordinate alpha = 1.0) {
    cr->set_source_rgba(r, g, b, alpha);
  }

  void setLineWidth(tCoordinate w) { cr->set_line_width(w); }

  void setLineDash(const std::vector<double> &dash = {2, 2},
                   const uint offset = 0) {
    cr->set_dash(dash, offset);
  }

  void unsetLineDash() { cr->unset_dash(); }

  void savePNG(const std::string &file) const { cs->write_to_png(file); }

public:
  static tRGB tetradicColor(uint i) { return tRGB(!(i & 2), i & 2, i & 1); }

private:
  void _init(const dBox &_bounds, uint _resolution);
  void _copy(const Painter &a);

  void line(const dPoints &points, uint a, uint b, bool drawInfinite);

  inline tCoordinate width() const { return imgDim(0); }

  inline tCoordinate height() const { return imgDim(1); }

  inline tCoordinate imgDim(const uint dim) const {
    return img.dim[dim] - img.coords[dim];
  }

  inline tCoordinate translatePoint(tCoordinate in, uint dim) const {
    return ((offset.coords[dim] + in - bounds.coords[dim]) /
            (3 * bounds.dim[dim])) *
           imgDim(dim);
  }

  inline tCoordinate translateLength(tCoordinate in, uint dim) const {
    return (in / (3 * bounds.dim[dim])) * imgDim(dim);
  }

private:
  dBox bounds;
  dBox img;
  dPoint offset;

  Cairo::RefPtr<Cairo::ImageSurface> cs;
  Cairo::RefPtr<Cairo::Context> cr;
};

#include <boost/progress.hpp>
#include <deque>

class PainterBulkWriter {

public:
  const uint FLUSH_THRESHOLD = 100;

public:
  ~PainterBulkWriter() { flush(); }

  template <class... Types>
  void add(const std::string &file, Types &&... args) {
    if (m_items.size() > FLUSH_THRESHOLD)
      flush();

    m_items.emplace_front(std::piecewise_construct, std::forward_as_tuple(file),
                          std::forward_as_tuple(args...));
  }

  Painter &top() { return m_items.front().second; }

private:
  void flush() {

    if (!m_items.empty()) {
      LOG << "Writing " << m_items.size() << " items" << std::endl;

      boost::progress_display progress(m_items.size(), LOG);

      while (!m_items.empty()) {
        m_items.back().second.savePNG(m_items.back().first);
        m_items.pop_back();
        ++progress;
      }
      // clear the memory
      m_items.shrink_to_fit();
    }
  }

private:
  std::deque<std::pair<std::string, Painter>> m_items;
};

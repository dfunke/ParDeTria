#pragma once

#include <algorithm>
#include <vector>
#include <tuple>

#include <cairo/cairo.h>
#include <cairomm/cairomm.h>

#include "Geometry.h"
#include "utils/Logger.h"

// define tuple for color
typedef std::tuple<tCoordinate, tCoordinate, tCoordinate> tRGB;
typedef std::tuple<tCoordinate, tCoordinate, tCoordinate, tCoordinate> tRGBa;

template <uint D> struct DataHolder;

template <> struct DataHolder<2> {
  dBox<2> bounds;
  dBox<2> img;
  dVector<2> offset;

  uint logLine;

  static constexpr uint PADDING = 10;
  static constexpr uint FONT_SIZE = 15;

  Cairo::RefPtr<Cairo::ImageSurface> cs;
  Cairo::RefPtr<Cairo::Context> cr;
};

template <> struct DataHolder<3> {

  template <class Object> class ColoredObject : public Object {

  public:
    ColoredObject() : Object() {}
    ColoredObject(const Object &o) : Object(o) {}
    ColoredObject(const Object &o, const tRGBa &c) : Object(o), color(c) {}

    ColoredObject &operator=(const Object &a) {

      Object::operator=(a);

      return *this;
    }

  public:
    tRGBa color;
  };

  IndexedVector<ColoredObject<dSimplex<3>>> pSimplices;
  IndexedVector<ColoredObject<dPoint<3>>> pPoints;
  std::vector<ColoredObject<dSphere<3>>> pSpheres;
  std::vector<ColoredObject<std::string>> pLog;
  tRGBa cColor;
};

template <uint D> class Painter {

public:
  Painter() : logging(false), dashed(false) {
    // non working version;
  }

  Painter(const dBox<D> &_bounds, uint _resolution = 10) {

    _init(_bounds, _resolution);
  }

  Painter(const Painter<D> &a) { _copy(a); }

  Painter &operator=(const Painter<D> &a) {

    _copy(a);

    return *this;
  }

  void draw(const dPoint<D> &point, bool drawInfinite = false);
  void draw(const dPoints<D> &points, bool drawInfinite = false) {
    for (const auto &p : points)
      draw(p, drawInfinite);
  }

  void draw(const dSimplex<D> &simplex, const dPoints<D> &points,
            bool drawInfinite = false);
  void draw(const dSimplices<D> &simplices, const dPoints<D> &points,
            bool drawInfinite = false) {
    for (const auto &s : simplices)
      draw(s, points, drawInfinite);
  }

  void drawCircumSphere(const dSimplex<D> &simplex, const dPoints<D> &points,
                        bool drawInfinite = false);
  void drawCircumSphere(const dSimplices<D> &simplices,
                        const dPoints<D> &points, bool drawInfinite = false) {
    for (const auto &x : simplices)
      drawCircumSphere(x, points, drawInfinite);
  }

  void drawNeighbors(const dSimplex<D> &simplex, const dSimplices<D> &neighbors,
                     const dPoints<D> &points, bool drawInfinite = false) {
    for (uint n : simplex.neighbors) {
      if (dSimplex<3>::isFinite(n) && neighbors.contains(n)) {
        draw(neighbors[n], points, drawInfinite);
      }
    }
  }

  void drawNeighbors(const dSimplices<D> &simplices,
                     const dSimplices<D> &neighbors, const dPoints<D> &points,
                     bool drawInfinite = false) {
    for (const auto &x : simplices)
      drawNeighbors(x, neighbors, points, drawInfinite);
  }

  void drawPartition(const dPoints<D> &points);

  void setColor(const tRGBa &rgba) {
    setColor(std::get<0>(rgba), std::get<1>(rgba), std::get<2>(rgba),
             std::get<3>(rgba));
  }

  void setColor(const tRGB &rgb, tCoordinate alpha = 1.0) {
    setColor(std::get<0>(rgb), std::get<1>(rgb), std::get<2>(rgb), alpha);
  }

  void setColor(tCoordinate r, tCoordinate g, tCoordinate b,
                tCoordinate alpha = 1.0);

  void save(const std::string &file) const;

  void setLogging(bool _logging = true) { logging = _logging; }
  void setDashed(bool _dashed = true) { dashed = _dashed; }

public:
  static tRGB tetradicColor(uint i) {
    return tRGB((bool)!(i & 2), (bool)i & 2, (bool)i & 1);
  }

private:
  void _init(const dBox<D> &_bounds, uint _resolution);
  void _copy(const Painter<D> &a);

  template <typename Object> void _log(const Object &o);

private:
  DataHolder<D> data;
  bool logging;
  bool dashed;
};

#include <boost/progress.hpp>
#include <deque>

template <uint D> class PainterBulkWriter {

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

  Painter<D> &top() { return m_items.front().second; }

private:
  void flush() {

    if (!m_items.empty()) {
      LOG << "Writing " << m_items.size() << " items" << std::endl;

      boost::progress_display progress(m_items.size(), LOG);

      while (!m_items.empty()) {
        m_items.back().second.save(m_items.back().first);
        m_items.pop_back();
        ++progress;
      }
      // clear the memory
      m_items.shrink_to_fit();
    }
  }

private:
  std::deque<std::pair<std::string, Painter<D>>> m_items;
};

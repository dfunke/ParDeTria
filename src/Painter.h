#pragma once

#include <algorithm>
#include <vector>
#include <tuple>

#include <cairo/cairo.h>
#include <cairomm/cairomm.h>

#include "Geometry.h"
#include "utils/Logger.h"

// define tuple for color
typedef std::tuple<float, float, float> tRGB;
typedef std::tuple<float, float, float, float> tRGBa;

// forward declarations
template <uint D, typename Precision> class Painter;
template <uint D, typename Precision> struct PainterImplementation;

#include "Painter_impl.hxx"

template <uint D, typename Precision> class Painter {

public:
  friend struct PainterImplementation<D, Precision>;

  Painter() {
    // non working version;
  }

  Painter(const dBox<D, Precision> &_bounds, uint _resolution = 10) {

    impl._init(_bounds, _resolution);
  }

  Painter(const Painter<D, Precision> &a) { impl._copy(a); }

  Painter &operator=(const Painter<D, Precision> &a) {

    impl._copy(a);

    return *this;
  }

  void draw(const dPoint<D, Precision> &point, bool drawInfinite = false) {
    impl.draw(point, drawInfinite);
  }
  void draw(const dPoints<D, Precision> &points, bool drawInfinite = false) {
    for (const auto &p : points)
      draw(p, drawInfinite);
  }

  void draw(const dSimplex<D, Precision> &simplex,
            const dPoints<D, Precision> &points, bool drawInfinite = false) {
    impl.draw(simplex, points, drawInfinite);
  }
  void draw(const dSimplices<D, Precision> &simplices,
            const dPoints<D, Precision> &points, bool drawInfinite = false) {
    for (const auto &s : simplices)
      draw(s, points, drawInfinite);
  }

  void drawCircumSphere(const dSimplex<D, Precision> &simplex,
                        const dPoints<D, Precision> &points,
                        bool drawInfinite = false) {
    impl.drawCircumSphere(simplex, points, drawInfinite);
  }
  void drawCircumSphere(const dSimplices<D, Precision> &simplices,
                        const dPoints<D, Precision> &points,
                        bool drawInfinite = false) {
    for (const auto &x : simplices)
      drawCircumSphere(x, points, drawInfinite);
  }

  void drawNeighbors(const dSimplex<D, Precision> &simplex,
                     const dSimplices<D, Precision> &neighbors,
                     const dPoints<D, Precision> &points,
                     bool drawInfinite = false) {
    for (uint n : simplex.neighbors) {
      if (dSimplex<D, Precision>::isFinite(n) && neighbors.contains(n)) {
        draw(neighbors[n], points, drawInfinite);
      }
    }
  }

  void drawNeighbors(const dSimplices<D, Precision> &simplices,
                     const dSimplices<D, Precision> &neighbors,
                     const dPoints<D, Precision> &points,
                     bool drawInfinite = false) {
    for (const auto &x : simplices)
      drawNeighbors(x, neighbors, points, drawInfinite);
  }

  void drawPartition(const dPoints<D, Precision> &points) {
    impl.drawPartition(points);
  }

  void setColor(const tRGBa &rgba) {
    setColor(std::get<0>(rgba), std::get<1>(rgba), std::get<2>(rgba),
             std::get<3>(rgba));
  }

  void setColor(const tRGB &rgb, float alpha = 1.0) {
    setColor(std::get<0>(rgb), std::get<1>(rgb), std::get<2>(rgb), alpha);
  }

  void setColor(float r, float g, float b, float alpha = 1.0) {
    impl.setColor(r, g, b, alpha);
  }

  void save(const std::string &file) const {
    if (ENABLED)
      impl.save(file);
  }

  void setLogging(bool _logging = true) { impl.logging = _logging; }
  void setDashed(bool _dashed = true) { impl.dashed = _dashed; }

private:
  PainterImplementation<D, Precision> impl;

public:
  static tRGB tetradicColor(uint i) {
    return tRGB((bool)!(i & 2), (bool)i & 2, (bool)i & 1);
  }

  static bool ENABLED;
};

#include <boost/progress.hpp>
#include <deque>

template <uint D, typename Precision> class PainterBulkWriter {

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

  Painter<D, Precision> &top() { return m_items.front().second; }

private:
  void flush() {

    if (!m_items.empty()) {
      LOG("Writing " << m_items.size() << " items" << std::endl);

      boost::progress_display progress(
          m_items.size(), LOGGER.addLogEntry(Logger::Verbosity::NORMAL));

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
  std::deque<std::pair<std::string, Painter<D, Precision>>> m_items;
};

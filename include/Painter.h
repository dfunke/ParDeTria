#pragma once

#include <algorithm>
#include <vector>
#include <tuple>

#include "Geometry.h"
#include "utils/Logger.h"

// define tuple for color
typedef std::tuple<float, float, float> tRGB;
typedef std::tuple<float, float, float, float> tRGBa;

// forward declarations
template <uint D, typename Precision> class Painter;
template <uint D, typename Precision> struct PainterImplementation;

#ifndef NDEBUG
#include "Painter_impl.hxx"
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

template <uint D, typename Precision> class Painter {

public:
  friend struct PainterImplementation<D, Precision>;

  Painter() {
    // non working version;
  }

  Painter(const dBox<D, Precision> &_bounds, uint _resolution = 10) {
#ifndef NDEBUG
    if (ENABLED)
      impl._init(_bounds, _resolution);
#endif
  }

  Painter(const Painter<D, Precision> &a) {
#ifndef NDEBUG
    if (ENABLED)
      impl._copy(a);
#endif
  }

  Painter &operator=(const Painter<D, Precision> &a) {

#ifndef NDEBUG
    if (ENABLED)
      impl._copy(a);
#endif

    return *this;
  }
  
  void drawBox(const dBox<2, Precision>& bounds) {
    #ifndef NDEBUG
        if (ENABLED)
        impl.drawBox(bounds);
    #endif
  }
  
  void drawText(const std::string& text, dVector<2, Precision> coords) {
    #ifndef NDEBUG
        if (ENABLED)
        impl.drawText(text, coords);
    #endif
  }

  void draw(const dPoint<D, Precision> &point, bool drawInfinite = false) {
#ifndef NDEBUG
    if (ENABLED)
      impl.draw(point, drawInfinite);
#endif
  }
  void draw(const dPoints<D, Precision> &points, bool drawInfinite = false) {
    if (ENABLED)
      for (const auto &p : points)
        draw(p, drawInfinite);
  }

  void draw(const dSimplex<D, Precision> &simplex,
            const dPoints<D, Precision> &points, bool drawInfinite = false) {
#ifndef NDEBUG
    if (ENABLED)
      impl.draw(simplex, points, drawInfinite);
#endif
  }
  void draw(const dSimplices<D, Precision> &simplices,
            const dPoints<D, Precision> &points, bool drawInfinite = false) {
    if (ENABLED)
      for (const auto &s : simplices)
        draw(s, points, drawInfinite);
  }

  void drawCircumSphere(const dSimplex<D, Precision> &simplex,
                        const dPoints<D, Precision> &points,
                        bool drawInfinite = false) {
#ifndef NDEBUG
    if (ENABLED)
      impl.drawCircumSphere(simplex, points, drawInfinite);
#endif
  }
  void drawCircumSphere(const dSimplices<D, Precision> &simplices,
                        const dPoints<D, Precision> &points,
                        bool drawInfinite = false) {
    if (ENABLED)
      for (const auto &x : simplices)
        drawCircumSphere(x, points, drawInfinite);
  }

  void drawNeighbors(const dSimplex<D, Precision> &simplex,
                     const dSimplices<D, Precision> &neighbors,
                     const dPoints<D, Precision> &points,
                     bool drawInfinite = false) {
      assert(false);
    /*if (ENABLED)
      for (uint n : simplex.neighbors) {
        if (dSimplex<D, Precision>::isFinite(n) && neighbors.contains(n)) {
          draw(neighbors[n], points, drawInfinite);
        }
      }*/
  }

  void drawNeighbors(const dSimplices<D, Precision> &simplices,
                     const dSimplices<D, Precision> &neighbors,
                     const dPoints<D, Precision> &points,
                     bool drawInfinite = false) {
    if (ENABLED)
      for (const auto &x : simplices)
        drawNeighbors(x, neighbors, points, drawInfinite);
  }

  void drawPartition(const dPoints<D, Precision> &points) {
#ifndef NDEBUG
    if (ENABLED)
      impl.drawPartition(points);
#endif
  }

  void setColor(const tRGBa &rgba) {
    setColor(std::get<0>(rgba), std::get<1>(rgba), std::get<2>(rgba),
             std::get<3>(rgba));
  }

  void setColor(const tRGB &rgb, float alpha = 1.0) {
    setColor(std::get<0>(rgb), std::get<1>(rgb), std::get<2>(rgb), alpha);
  }

  void setColor(float r, float g, float b, float alpha = 1.0) {
#ifndef NDEBUG
    if (ENABLED)
      impl.setColor(r, g, b, alpha);
#endif
  }

  void save(const std::string &file) const {
#ifndef NDEBUG
    if (ENABLED)
      impl.save(file);
#endif
  }

  void setLogging(bool _logging = true) {
#ifndef NDEBUG
    impl.logging = _logging;
#endif
  }
  void setDashed(bool _dashed = true) {
#ifndef NDEBUG
    impl.dashed = _dashed;
#endif
  }

private:
#ifndef NDEBUG
  PainterImplementation<D, Precision> impl;
#endif

public:
  static tRGB tetradicColor(uint i) {
    return tRGB((bool)(i & 2), (bool)!(i & 2), (bool)i & 1);
  }

  static bool ENABLED;
};

#pragma GCC diagnostic pop

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

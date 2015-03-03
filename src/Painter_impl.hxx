#pragma once

#include <algorithm>
#include <vector>
#include <tuple>

#include <cairo/cairo.h>
#include <cairomm/cairomm.h>
#include <bits/atomic_base.h>

#include "Geometry.h"
#include "utils/Logger.h"

// define tuple for color
typedef std::tuple<float, float, float> tRGB;
typedef std::tuple<float, float, float, float> tRGBa;

// forward declarations
template <uint D, typename Precision> class Painter;
template <uint D, typename Precision> struct PainterImplementation;

template <typename Precision> struct PainterImplementation<2, Precision> {
public:
  __attribute__((noinline)) void draw(const dPoint<2, Precision> &point,
                                      bool drawInfinite = false);

  __attribute__((noinline)) void draw(const dSimplex<2, Precision> &simplex,
                                      const dPoints<2, Precision> &points,
                                      bool drawInfinite = false);

  __attribute__((noinline)) void drawCircumSphere(
      const dSimplex<2, Precision> &simplex,
      const dPoints<2, Precision> &points, bool drawInfinite = false);

  __attribute__((noinline)) void drawPartition(
      const dPoints<2, Precision> &points);

  __attribute__((noinline)) void setColor(float r, float g, float b,
                                          float alpha = 1.0);

  __attribute__((noinline)) void save(const std::string &file) const;

  __attribute__((noinline)) void _init(const dBox<2, Precision> &_bounds,
                                       uint _resolution);

  __attribute__((noinline)) void _copy(const Painter<2, Precision> &a);

  template <typename Object>
  __attribute__((noinline)) void _log(const Object &o);

private:
  float imgDim(const uint dim);
  float translatePoint(float in, uint dim);
  float translateLength(float in, uint dim);
  void stroke();

private:
  dBox<2, Precision> bounds;
  dBox<2, Precision> img;
  dVector<2, Precision> offset;

  uint logLine;
  std::atomic_flag strokeMtx = ATOMIC_FLAG_INIT;

public:
  bool logging;
  bool dashed;

private:
  static constexpr uint PADDING = 10;
  static constexpr uint FONT_SIZE = 15;

  Cairo::RefPtr<Cairo::ImageSurface> cs;
  Cairo::RefPtr<Cairo::Context> cr;
};

template <typename Precision> struct PainterImplementation<3, Precision> {

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

public:
  __attribute__((noinline)) void draw(const dPoint<3, Precision> &point,
                                      bool drawInfinite = false);
  __attribute__((noinline)) void draw(const dSimplex<3, Precision> &simplex,
                                      const dPoints<3, Precision> &points,
                                      bool drawInfinite = false);
  __attribute__((noinline)) void drawCircumSphere(
      const dSimplex<3, Precision> &simplex,
      const dPoints<3, Precision> &points, bool drawInfinite = false);
  __attribute__((noinline)) void drawPartition(
      const dPoints<3, Precision> &points);

  __attribute__((noinline)) void setColor(float r, float g, float b,
                                          float alpha = 1.0);

  __attribute__((noinline)) void save(const std::string &file) const;

  __attribute__((noinline)) void _init(const dBox<3, Precision> &_bounds,
                                       uint _resolution);
  __attribute__((noinline)) void _copy(const Painter<3, Precision> &a);

  template <typename Object>
  __attribute__((noinline)) void _log(const Object &o);

public:
  bool logging;
  bool dashed;

private:
  std::vector<ColoredObject<dSimplex<3, Precision>>> pSimplices;
  std::vector<ColoredObject<dPoint<3, Precision>>> pPoints;
  std::vector<ColoredObject<dSphere<3, Precision>>> pSpheres;
  std::vector<ColoredObject<std::string>> pLog;
  tRGBa cColor;
};

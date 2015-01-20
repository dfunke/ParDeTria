#include "Painter.h"
#include "Partitioner.h"

// disable warnings in VTK library
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wextra-semi"

#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPoints.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkTetra.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>

#pragma clang diagnostic pop

#include <valgrind/callgrind.h>

// static variables

constexpr uint DataHolder<3>::PADDING;
constexpr uint DataHolder<3>::FONT_SIZE;

// TODO remove
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"

template <> template <typename Object> void Painter<3>::_log(const Object &o) {
  if (logging) {
    // TODO logging
  }
}

template <>
void Painter<3>::setColor(tCoordinate r, tCoordinate g, tCoordinate b,
                          tCoordinate alpha) {
  data.cColor = std::make_tuple(r, g, b, alpha);
}

//
template <> void Painter<3>::draw(const dPoint<3> &point, bool drawInfinite) {

  if (point.isFinite() || drawInfinite) {

    if (data.pPoints.contains(point)) {
      data.pPoints[point.id].color = data.cColor; // update color of object
      return;
    }

    DataHolder<3>::ColoredObject<dPoint<3>> p = point;
    p.color = data.cColor;

    data.pPoints[p.id] = p;
  }

  _log(point);
}

template <>
void Painter<3>::draw(const dSimplex<3> &simplex, const dPoints<3> &points,
                      bool drawInfinite) {

  if (simplex.isFinite() || drawInfinite) {

    if (data.pSimplices.contains(simplex)) {
      data.pSimplices[simplex.id].color = data.cColor; // update color of object
      return;
    }

    DataHolder<3>::ColoredObject<dSimplex<3>> p = simplex;
    p.color = data.cColor;

    data.pSimplices[p.id] = p;
  }

  _log(simplex);
}

template <>
void Painter<3>::drawCircumSphere(const dSimplex<3> &simplex,
                                  const dPoints<3> &points, bool drawInfinite) {

  if (simplex.isFinite() || drawInfinite) {

    DataHolder<3>::ColoredObject<dSphere<3>> s = simplex.circumsphere(points);
    s.color = data.cColor;

    data.pSpheres.push_back(s);
  }
}

template <> void Painter<3>::drawPartition(const dPoints<3> &points) {
  // TODO draw partition
}

template <>
void Painter<3>::save(
#ifdef NO_OUTPUT
    __attribute((unused))
#endif
    const std::string &file) const {
// stop callgrind instrumentation for saving VTU

#ifndef NO_OUTPUT // disable png output
  CALLGRIND_STOP_INSTRUMENTATION;

  auto points = vtkSmartPointer<vtkPoints>::New();

  std::map<uint, uint> pointIdMap;

  for (const auto &o : data.pPoints) {
    pointIdMap[o.id] = points->InsertNextPoint(o.coords.data());
  }

  auto cellArray = vtkSmartPointer<vtkCellArray>::New();

  auto colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName("Colors");

  for (const auto &o : data.pSimplices) {
    auto tetra = vtkSmartPointer<vtkTetra>::New();

    for (uint d = 0; d < 4; ++d) {
      tetra->GetPointIds()->SetId(d, pointIdMap[o.vertices[d]]);
    }

    cellArray->InsertNextCell(tetra);
    colors->InsertNextTuple3(255 * std::get<0>(o.color),
                             255 * std::get<1>(o.color),
                             255 * std::get<2>(o.color));
  }

  auto polyData = vtkSmartPointer<vtkPolyData>::New();
  polyData->SetPoints(points);
  polyData->SetPolys(cellArray);
  polyData->GetCellData()->SetScalars(colors);

  // Write file
  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName((file + ".vtu").c_str());
  writer->SetInputData(polyData);
  writer->Write();

  CALLGRIND_START_INSTRUMENTATION;
#endif
}

template <> void Painter<3>::_init(const dBox<3> &_bounds, uint _resolution) {
  logging = false;
  dashed = false;
  data.logLine = 1;
}

template <> void Painter<3>::_copy(const Painter<3> &a) {
  logging = a.logging;
  data.logLine = a.data.logLine;

  // copy data
  data.cColor = a.data.cColor;
  data.pPoints.insert(a.data.pPoints.begin(), a.data.pPoints.end());
  data.pSimplices.insert(a.data.pSimplices.begin(), a.data.pSimplices.end());
  data.pSpheres.insert(data.pSpheres.end(), a.data.pSpheres.begin(),
                       a.data.pSpheres.end());
}

// TODO remove
#pragma clang diagnostic pop

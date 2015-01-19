#include "Painter.h"
#include "Partitioner.h"

// disable warnings in VTK library
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wextra-semi"

#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

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
  // TODO set color
}

//
template <> void Painter<3>::draw(const dPoint<3> &point, bool drawInfinite) {

  if (point.isFinite() && point.id > data.points->GetNumberOfPoints()) {
    data.points->InsertPoint(point.id, point.coords.data());
  }

  _log(point);
}

template <>
void Painter<3>::draw(const dSimplex<3> &simplex, const dPoints<3> &points,
                      bool drawInfinite) {

  if (simplex.isFinite() && simplex.id > data.cells->GetNumberOfCells()) {
    auto tetra = vtkSmartPointer<vtkTetra>::New();

    for (uint d = 0; d < 3; ++d) {
      draw(points[simplex.vertices[d]]); // ensure point is in the list of
                                         // points
      tetra->GetPointIds()->SetId(d, simplex.vertices[d]);
    }

    data.cells->InsertNextCell(tetra);
  }

  _log(simplex);
}

template <>
void Painter<3>::drawCircumSphere(const dSimplex<3> &s,
                                  const dPoints<3> &points, bool drawInfinite) {

  // TODO draw circumsphere
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

  auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  grid->SetPoints(data.points);
  grid->SetCells(VTK_TETRA, data.cells);

  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
      vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName((file + ".vtu").c_str());
  writer->SetInputData(grid);
  writer->Write();

  CALLGRIND_START_INSTRUMENTATION;
#endif
}

template <> void Painter<3>::_init(const dBox<3> &_bounds, uint _resolution) {
  logging = false;
  dashed = false;
  data.logLine = 1;

  data.points = vtkSmartPointer<vtkPoints>::New();
  data.cells = vtkSmartPointer<vtkCellArray>::New();
}

template <> void Painter<3>::_copy(const Painter<3> &a) {
  logging = a.logging;
  data.logLine = a.data.logLine;

  data.points = vtkSmartPointer<vtkPoints>::New();
  data.cells = vtkSmartPointer<vtkCellArray>::New();

  data.points->DeepCopy(a.data.points);
  data.cells->DeepCopy(a.data.cells);
}

// TODO remove
#pragma clang diagnostic pop

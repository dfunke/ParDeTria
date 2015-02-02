#include "Painter.h"
#include "Partitioner.h"

#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkTetra.h>
#include <vtkPolyData.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPolyDataWriter.h>

#include <valgrind/callgrind.h>

template <class ColoredObject>
void updateColor(ColoredObject &o, const tRGBa &color) {

  tRGBa res;
  // alpha_0 = alpha_a + alpha_b (1 - alpha_a)
  std::get<3>(res) =
      std::get<3>(o.color) + std::get<3>(color) * (1 - std::get<3>(o.color));

  tCoordinate invAlpha = 1 / std::get<3>(res);

  // color_0 = 1/alpha_0 (color_a * alpha_a + color_b * alpha_b (1 - alpha_a))
  std::get<0>(res) = invAlpha * (std::get<0>(o.color) * std::get<3>(o.color) +
                                 std::get<0>(color) * std::get<3>(color) *
                                     (1 - std::get<3>(o.color)));

  std::get<1>(res) = invAlpha * (std::get<1>(o.color) * std::get<3>(o.color) +
                                 std::get<1>(color) * std::get<3>(color) *
                                     (1 - std::get<3>(o.color)));

  std::get<2>(res) = invAlpha * (std::get<2>(o.color) * std::get<3>(o.color) +
                                 std::get<2>(color) * std::get<3>(color) *
                                     (1 - std::get<3>(o.color)));

  o.color = res;
}

template <> template <typename Object> void Painter<3>::_log(const Object &o) {
  if (logging) {
    data.pLog.emplace_back(to_string(o), data.cColor);
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
    data.pPoints.emplace_back(point, data.cColor);
  }

  _log(point);
}

template <>
void Painter<3>::draw(const dSimplex<3> &simplex,
                      __attribute((unused)) const dPoints<3> &points,
                      bool drawInfinite) {

  if (simplex.isFinite() || drawInfinite) {
    data.pSimplices.emplace_back(simplex, data.cColor);
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

template <>
void Painter<3>::drawPartition(__attribute((unused)) const dPoints<3> &points) {
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

  PLOG << "Writing file " << (file + ".vtu") << std::endl;

  auto points = vtkSmartPointer<vtkPoints>::New();

  auto colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(4);
  colors->SetName("Colors");

  std::map<uint, uint> pointIdMap;

  for (const auto &o : data.pPoints) {
    pointIdMap[o.id] = points->InsertNextPoint(o.coords.data());
    colors->InsertNextTuple4(
        255 * std::get<0>(o.color), 255 * std::get<1>(o.color),
        255 * std::get<2>(o.color), 255 * std::get<3>(o.color));
  }

  auto pointsPolyData = vtkSmartPointer<vtkPolyData>::New();
  pointsPolyData->SetPoints(points);

  auto vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
  vertexFilter->SetInputData(pointsPolyData);
  vertexFilter->Update();

  auto cellArray = vtkSmartPointer<vtkCellArray>::New();

  for (const auto &o : data.pSimplices) {
    auto tetra = vtkSmartPointer<vtkTetra>::New();

    for (uint d = 0; d < 4; ++d) {
      tetra->GetPointIds()->SetId(d, pointIdMap[o.vertices[d]]);
    }

    cellArray->InsertNextCell(tetra);
    colors->InsertNextTuple4(
        255 * std::get<0>(o.color), 255 * std::get<1>(o.color),
        255 * std::get<2>(o.color), 255 * std::get<3>(o.color));

    PLOG << o.id << ": " << 255 * std::get<0>(o.color) << " "
         << 255 * std::get<1>(o.color) << " " << 255 * std::get<2>(o.color)
         << " " << 255 * std::get<3>(o.color) << std::endl;
  }

  auto polyData = vtkSmartPointer<vtkPolyData>::New();
  polyData->ShallowCopy(vertexFilter->GetOutput());
  // polyData->SetPoints(points);
  polyData->SetPolys(cellArray);
  // polyData->GetPointData()->SetScalars(pColors);
  polyData->GetCellData()->SetScalars(colors);

  // Write file
  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName((file + ".vtu").c_str());
  writer->SetInputData(polyData);
  writer->Write();

  CALLGRIND_START_INSTRUMENTATION;
#endif
}

template <>
void Painter<3>::_init(__attribute((unused)) const dBox<3> &_bounds,
                       __attribute((unused)) uint _resolution) {
  logging = false;
  dashed = false;

  data.cColor = std::make_tuple(0, 0, 0, 1);
}

template <> void Painter<3>::_copy(const Painter<3> &a) {
  logging = a.logging;

  // copy data
  data.cColor = a.data.cColor;
  data.pPoints.insert(data.pPoints.end(), a.data.pPoints.begin(),
                      a.data.pPoints.end());
  data.pSimplices.insert(data.pSimplices.end(), a.data.pSimplices.begin(),
                         a.data.pSimplices.end());
  data.pSpheres.insert(data.pSpheres.end(), a.data.pSpheres.begin(),
                       a.data.pSpheres.end());
  data.pLog.insert(data.pLog.end(), a.data.pLog.begin(), a.data.pLog.end());
}

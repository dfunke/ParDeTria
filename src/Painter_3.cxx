#include "Painter.h"
#include "Painter_impl.hxx"
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

//

template <class ColoredObject>
void updateColor(ColoredObject &o, const tRGBa &color) {

  tRGBa res;
  // alpha_0 = alpha_a + alpha_b (1 - alpha_a)
  std::get<3>(res) =
      std::get<3>(o.color) + std::get<3>(color) * (1 - std::get<3>(o.color));

  float invAlpha = 1 / std::get<3>(res);

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

template <typename Precision>
template <typename Object>
void PainterImplementation<3, Precision>::_log(const Object &o) {
  if (logging) {
    pLog.emplace_back(to_string(o), cColor);
  }
}

template <typename Precision>
void PainterImplementation<3, Precision>::setColor(float r, float g, float b,
                                                   float alpha) {
  cColor = std::make_tuple(r, g, b, alpha);
}

//
template <typename Precision>
void PainterImplementation<3, Precision>::draw(
    const dPoint<3, Precision> &point, bool drawInfinite) {

  if (point.isFinite() || drawInfinite) {
    pPoints.emplace_back(point, cColor);
  }

  _log(point);
}

template <typename Precision>
void PainterImplementation<3, Precision>::draw(
    const dSimplex<3, Precision> &simplex,
    __attribute((unused)) const dPoints<3, Precision> &points,
    bool drawInfinite) {

  if (simplex.isFinite() || drawInfinite) {
    pSimplices.emplace_back(simplex, cColor);
  }

  _log(simplex);
}

template <typename Precision>
void PainterImplementation<3, Precision>::drawCircumSphere(
    const dSimplex<3, Precision> &simplex, const dPoints<3, Precision> &points,
    bool drawInfinite) {

  if (simplex.isFinite() || drawInfinite) {

    PainterImplementation<3, Precision>::ColoredObject<dSphere<3, Precision>>
        s = simplex.circumsphere(points);
    s.color = cColor;

    pSpheres.push_back(s);
  }
}

template <typename Precision>
void PainterImplementation<3, Precision>::drawPartition(
    const dPoints<3, Precision> &points) {

  auto stats = getPointStats(std::size_t(0), points.size(), points);
  pPoints.emplace_back(stats.mid, std::make_tuple(1, 0, 0, 1));
}

template <typename Precision>
void PainterImplementation<3, Precision>::save(const std::string &file) const {
  // stop callgrind instrumentation for saving VTU

  CALLGRIND_STOP_INSTRUMENTATION;

  PLOG("Writing file " << (file + ".vtu") << std::endl);

  auto points = vtkSmartPointer<vtkPoints>::New();

  auto colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(4);
  colors->SetName("Colors");

  std::map<uint, uint> pointIdMap;

  for (const auto &o : pPoints) {
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

  for (const auto &o : pSimplices) {
    auto tetra = vtkSmartPointer<vtkTetra>::New();

    for (uint d = 0; d < 4; ++d) {
      tetra->GetPointIds()->SetId(d, pointIdMap[o.vertices[d]]);
    }

    cellArray->InsertNextCell(tetra);
    colors->InsertNextTuple4(
        255 * std::get<0>(o.color), 255 * std::get<1>(o.color),
        255 * std::get<2>(o.color), 255 * std::get<3>(o.color));

    /*PLOG(o.id << ": " << 255 * std::get<0>(o.color) << " "
         << 255 * std::get<1>(o.color) << " " << 255 * std::get<2>(o.color)
         << " " << 255 * std::get<3, Precision>(o.color) << std::endl);*/
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
}

template <typename Precision>
void PainterImplementation<3, Precision>::_init(
    const dBox<3, Precision> &_bounds, __attribute((unused)) uint _resolution) {

  logging = false;
  dashed = false;

  pPoints.emplace_back(_bounds.low, std::make_tuple(1, 0, 0, 1));
  pPoints.emplace_back(_bounds.high, std::make_tuple(1, 0, 0, 1));

  cColor = std::make_tuple(0, 0, 0, 1);
}

template <typename Precision>
void PainterImplementation<3, Precision>::_copy(
#ifdef NDEBUG
    __attribute__((unused))
#endif
    const Painter<3, Precision> &a) {

#ifndef NDEBUG
  logging = a.impl.logging;

  // copy data
  cColor = a.impl.cColor;
  pPoints.insert(pPoints.end(), a.impl.pPoints.begin(), a.impl.pPoints.end());
  pSimplices.insert(pSimplices.end(), a.impl.pSimplices.begin(),
                    a.impl.pSimplices.end());
  pSpheres.insert(pSpheres.end(), a.impl.pSpheres.begin(),
                  a.impl.pSpheres.end());
  pLog.insert(pLog.end(), a.impl.pLog.begin(), a.impl.pLog.end());
#endif
}

// specializations

template struct PainterImplementation<3, float>;
template struct PainterImplementation<3, double>;

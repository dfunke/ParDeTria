// disable warnings in VTK library
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wextra-semi"

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkCellArray.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPointData.h>
#include <vtkVertexGlyphFilter.h>

#pragma clang diagnostic pop

int main(int argc, char *argv[]) {
  // Parse command line arguments
  if (argc != 2) {
    std::cout << "Required arguments: Input Filename" << std::endl;
    return EXIT_FAILURE;
  }

  std::string filename = argv[1];

  // Read and display file for verification that it was written correclty
  auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(reader->GetOutput());

  auto actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  auto renderer = vtkSmartPointer<vtkRenderer>::New();
  auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);

  auto renderWindowInteractor =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  renderer->AddActor(actor);
  renderer->SetBackground(255, 255, 255);

  renderWindow->Render();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}

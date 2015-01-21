#include <vtkObjectFactory.h>
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
#include <vtkProperty.h>
#include <vtkInteractorStyleTrackballCamera.h>

// boost
#include <boost/filesystem.hpp>

//******************************************************************************

vtkSmartPointer<vtkXMLPolyDataReader> reader;
vtkSmartPointer<vtkRenderWindow> renderWindow;

std::vector<boost::filesystem::path> files;

//******************************************************************************

// Define interaction style
class CustomInteractorStyle : public vtkInteractorStyleTrackballCamera {

public:
  static CustomInteractorStyle *New();

  vtkTypeMacro(CustomInteractorStyle, vtkInteractorStyleTrackballCamera);

  virtual void OnKeyPress() {
    // Get the keypress
    vtkRenderWindowInteractor *rwi = this->Interactor;
    std::string key = rwi->GetKeySym();

    // Handle an arrow key
    if (key == "Next") {
      ++currentFile;
    }

    if (key == "Prior") {
      --currentFile;
    }

    if (key == "Next" || key == "Prior") {
      currentFile %= files.size();

      reader->SetFileName(files[currentFile].native().c_str());

      reader->Modified();
      reader->Update();

      renderWindow->SetWindowName(files[currentFile].filename().c_str());
      renderWindow->Render();
    }

    // Forward events
    vtkInteractorStyleTrackballCamera::OnKeyPress();
  }

private:
  uint currentFile = 0;
};

vtkStandardNewMacro(CustomInteractorStyle);

//******************************************************************************

void displayFiles() {
  // Read and display file for verification that it was written correclty
  reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(files[0].native().c_str());
  reader->Update();

  auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(reader->GetOutput());

  auto actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetRepresentationToWireframe();

  auto renderer = vtkSmartPointer<vtkRenderer>::New();
  renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);

  auto renderWindowInteractor =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  auto interactor = vtkSmartPointer<CustomInteractorStyle>::New();
  renderWindowInteractor->SetInteractorStyle(interactor);

  renderer->AddActor(actor);
  renderer->SetBackground(255, 255, 255);

  renderWindow->SetWindowName(files[0].filename().c_str());
  renderWindow->Render();
  renderWindowInteractor->Start();
}

int main(int argc, char *argv[]) {
  // Parse command line arguments
  if (argc != 2) {
    std::cout << "required argument: input filename or directory" << std::endl;
    return EXIT_FAILURE;
  }

  std::string p = argv[1];

  try {
    if (boost::filesystem::exists(p)) // does p actually exist?
    {
      if (boost::filesystem::is_regular_file(p)) // is p a regular file?
        files.emplace_back(p);
      else if (boost::filesystem::is_directory(p)) { // is p a directory?
        for (auto it = boost::filesystem::directory_iterator(p);
             it != boost::filesystem::directory_iterator(); ++it) {

          if (it->path().extension() == ".vtu") { // we have a vtu file
            files.push_back(it->path());
          }
        }

        std::sort(files.begin(), files.end());

      }

      else
        cout << p << " exists, but is neither a regular file nor a directory\n";
    } else
      cout << p << " does not exist\n";
  }

  catch (const boost::filesystem::filesystem_error &ex) {
    cout << ex.what() << '\n';
  }

  for (const auto &f : files) {
    std::cout << f << std::endl;
  }

  displayFiles();

  return EXIT_SUCCESS;
}

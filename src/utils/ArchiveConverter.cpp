
// own
#include "Geometry.h"
#include "Partitioner.h"
#include "utils/Serialization.hxx"

#include <boost/program_options.hpp>
#include <boost/progress.hpp>

//**************************

int main(int argc, char *argv[]) {

  // define program options

  namespace po = boost::program_options;

  std::string inFile;
  std::string outFile;

  po::options_description cCommandLine("Command Line Options");
  cCommandLine.add_options()("in", po::value<std::string>(&inFile),
                             "load float points from file");
  cCommandLine.add_options()("out", po::value<std::string>(&outFile),
                             "store double points in file");
  cCommandLine.add_options()("help", "produce help message");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cCommandLine), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << cCommandLine << std::endl;
    return EXIT_SUCCESS;
  }

  if (!vm.count("in") || !vm.count("out")) {
    std::cerr << "Specify in- and output" << std::endl;
    return EXIT_FAILURE;
  }

  //***************************************************************************

  dPoints<3, float> inPoints;
  inPoints = loadObject<dPoints<3, float>>(inFile);

  dPoints<3, double> outPoints;
  for (const auto &i : inPoints) {
    outPoints[i.id].id = i.id;

    for (uint k = 0; k < i.coords.size(); ++k)
      outPoints[i.id].coords[k] = i.coords[k];

    // outPoints[i.id].simplices = i.simplices;
  }

  storeObject(outPoints, outFile);

  return EXIT_SUCCESS;
}


// own
#include "Geometry.h"
#include "Triangulator.h"
#include "Partitioner.h"
#include "Painter.h"

#include "utils/Random.h"
#include "utils/Logger.h"
#include "utils/CSV.h"
#include "utils/Timings.h"
#include "utils/RSS.h"
#include "utils/ASSERT.h"
#include "utils/Serialization.hxx"
#include "utils/ProgressDisplay.h"

#include <boost/program_options.hpp>
#include <boost/utility/in_place_factory.hpp>

#include <mutex>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

//**************************

#define D 2
#define Precision double

//**************************

struct TriangulateReturn {
  TriangulationReport tr;
  bool exception;
  tDuration time;
};

//**************************

TriangulateReturn triangulate(const dBox<D, Precision> &bounds,
                              const uint baseCase,
                              dPoints<D, Precision> &points,
                              const unsigned char splitter) {
  std::unique_ptr<Partitioner<D, Precision>> partitioner_ptr;
  switch (splitter) {
  case 'd':
    partitioner_ptr = std::make_unique<dPartitioner<D, Precision>>();
    break;
  case 'c':
    partitioner_ptr = std::make_unique<CyclePartitioner<D, Precision>>();
    break;
  default:
    // p must be a dimension - subtract '0' to get integer value
    int d = splitter - '0';
    ASSERT(0 <= d && d < D);
    partitioner_ptr = std::make_unique<kPartitioner<D, Precision>>(d);
    break;
  }

  Triangulator<D, Precision> triangulator(bounds, baseCase, points,
                                          std::move(partitioner_ptr));

  TriangulateReturn ret;

  try {

    auto t1 = Clock::now();
    auto dt = triangulator.triangulate();
    auto t2 = Clock::now();

    ret.tr = triangulator.getTriangulationReport();
    ret.exception = false;
    ret.time = std::chrono::duration_cast<tDuration>(t2 - t1);

  } catch (AssertionException &e) {
    std::cerr << "Assertion failed - ABORTING - n =  " << points.size()
              << " p = " << splitter << std::endl;
    std::cerr << e.what() << std::endl;

    // output points
    storeObject(points, "assertionFailed_" + std::to_string(points.size()) +
                            "_" + (char)splitter + ".dat");

    ret.exception = true;
  }

  return ret;
}

//**************************

int doStudy(const uint N, const dBox<D, Precision> &bounds, const uint baseCase,
            std::function<Precision()> &&dice) {
  // define splitter test cases
  std::vector<unsigned char> splitters = {'d', 'c', '0', '1'};

  LOGGER.setLogLevel(Logger::Verbosity::SILENT);

  Painter<D, Precision>::ENABLED = false;

  uint nIterations = 9 * (log10(N) - 1) + 1;
  ProgressDisplay progress(splitters.size() * nIterations, std::cout);

  std::ofstream f("triangulation_report.csv", std::ios::out | std::ios::trunc);
  f << CSV::csv("n", "splitter", "provenance", "base_case", "edge_tria",
                "valid", "nPoints", "nSimplices", "nEdgePoints",
                "nEdgeSimplices", "time", "mem", "max_mem") << std::endl;

  // generate test points deterministically
  std::vector<dPoints<D, Precision>> pointss;
  for (uint i = 0; i < nIterations; ++i) {
    uint n = ((i % 9) + 1) * std::pow(10, std::floor(i / 9 + 1));
    pointss.push_back(genPoints(n, bounds, dice));
  }

  std::mutex mtx;
  tbb::parallel_for(std::size_t(0), pointss.size(), [&](const uint i) {
    tbb::parallel_for(std::size_t(0), splitters.size(), [&](const uint j) {

      auto &points = pointss[i];
      unsigned char p = splitters[j];

      auto ret = triangulate(bounds, baseCase, points, p);

      // evaluate the triangulation report

      for (const auto &tr : ret.tr) {
        {
          std::lock_guard<std::mutex> lock(mtx);
          f << CSV::csv(points.size(), p, tr.provenance, tr.base_case,
                        tr.edge_triangulation, tr.valid, tr.nPoints,
                        tr.nSimplices, tr.nEdgePoints, tr.nEdgeSimplices,
                        ret.time.count(), getCurrentRSS(), getPeakRSS())
            << std::endl;
        }
        if (Triangulator<D, Precision>::VERIFY &&
            Triangulator<D, Precision>::isTOP(tr.provenance) && !tr.valid) {

          // output points
          storeObject(points, "invalid_" + std::to_string(points.size()) + "_" +
                                  (char)p + ".dat");
        }
      }

      ++progress;
    });
  });

  std::cout << LOGGER << std::endl;

  LOG("Finished" << std::endl);

  return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {

  // define program options

  namespace po = boost::program_options;

  int verbosity = -1;

  unsigned char p;
  uint N;
  uint baseCase;
  uint threads = tbb::task_scheduler_init::default_num_threads();

  bool study = false;

  std::string pointFile;
  bool loadPoints = false;

  po::options_description cCommandLine("Command Line Options");
  cCommandLine.add_options()("n", po::value<uint>(&N), "number of points");
  cCommandLine.add_options()("basecase", po::value<uint>(&baseCase),
                             "threshold for base case");
  cCommandLine.add_options()(
      "splitter", po::value<unsigned char>(&p),
      "splitter - _c_ycle, _d_-dimensional, _[0-d-1]_ fixed dimension");
  cCommandLine.add_options()("study", "study mode");
  cCommandLine.add_options()("no-verify", "don't verify triangulation");
  cCommandLine.add_options()("no-output", "don't output triangulation");
  cCommandLine.add_options()("verbosity", po::value<int>(&verbosity),
                             "verbosity");
  cCommandLine.add_options()("points", po::value<std::string>(&pointFile),
                             "load points from file");
  cCommandLine.add_options()("threads", po::value<uint>(&threads),
                             "specify number of threads");
  cCommandLine.add_options()("help", "produce help message");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cCommandLine), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << cCommandLine << std::endl;
    return EXIT_SUCCESS;
  }

  // evaluate program options
  LOGGER.setLogLevel(static_cast<Logger::Verbosity>(verbosity));

  // plausability checks
  bool valid = true;
  if ((!vm.count("basecase"))) {
    std::cout << "Please specify basecase threshold" << std::endl;
    valid = false;
  }

  if ((!vm.count("study") && !vm.count("splitter"))) {
    std::cout << "Please specify splitter" << std::endl;
    valid = false;
  }

  if ((!vm.count("study") && !(vm.count("n") || vm.count("points")))) {
    std::cout << "Please specify number of points or point file" << std::endl;
    valid = false;
  }

  if (!valid)
    return EXIT_FAILURE;

  if (vm.count("points")) {
    loadPoints = true; // filename will be in pointFile
  }

  if (vm.count("study")) {
    study = true;
  }

  if (vm.count("no-verify")) {
    Triangulator<D, Precision>::VERIFY = false;
  }

  if (vm.count("no-output")) {
    Painter<D, Precision>::ENABLED = false;
  }

  // load scheduler with specified number of threads
  tbb::task_scheduler_init init(threads);

  //***************************************************************************

  dBox<D, Precision> bounds;
  for (uint i = 0; i < D; ++i) {
    bounds.low[i] = 0;
    bounds.high[i] = 100;
  }

  std::uniform_real_distribution<Precision> distribution(0, 1);
  std::function<Precision()> dice = std::bind(distribution, generator);

  /*
   * File format for triangulation report
   * input_n splitter provenance base_case edge_triangulation nPoints
   * nSimplices
   * nEdgePoints nEdgeSimplices
   */

  uint returnCode = 0;
  if (study) {
    returnCode = doStudy(N, bounds, baseCase, std::move(dice));
  } else {

    dPoints<D, Precision> points;
    if (loadPoints) {
      points = loadObject<dPoints<D, Precision>>(pointFile);
    } else {
      points = genPoints(N, bounds, dice);
    }

    auto ret = triangulate(bounds, baseCase, points, p);

    LOG("Triangulating "
        << points.size() << " points took "
        << std::chrono::duration_cast<std::chrono::milliseconds>(ret.time)
               .count() << " ms" << std::endl);

    VLOG(std::endl
         << "Triangulation report: " << std::endl);

    VLOG(CSV::csv("n", "splitter", "provenance", "base_case", "edge_tria",
                  "valid", "nPoints", "nSimplices", "nEdgePoints",
                  "nEdgeSimplices", "time", "mem", "max_mem"));

    for (const auto &tr : ret.tr) {
      VLOG(CSV::csv(points.size(), p, tr.provenance, tr.base_case,
                    tr.edge_triangulation, tr.valid, tr.nPoints, tr.nSimplices,
                    tr.nEdgePoints, tr.nEdgeSimplices, ret.time.count(),
                    getCurrentRSS(), getPeakRSS()));
    }

    returnCode = EXIT_SUCCESS;
  }

  if (LOGGER.getLogLevel() >= Logger::Verbosity::NORMAL) {
    LOGGER.printLog(std::cout);
  }

  return returnCode;
}

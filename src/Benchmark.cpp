
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

#include <boost/program_options.hpp>

//**************************

#define D 3
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

int benchmark(const uint N, const dBox<D, Precision> &bounds,
              const uint baseCase, const unsigned char splitter,
              std::function<Precision()> &&dice) {

  LOGGER.setLogLevel(Logger::Verbosity::SILENT);
  Triangulator<D, Precision>::VERIFY = false;
  Painter<D, Precision>::ENABLED = false;

  uint nIterations = 9 * (log10(N) - 1) + 1;

  for (uint i = 0; i < nIterations; ++i) {
    uint n = ((i % 9) + 1) * std::pow(10, std::floor(i / 9 + 1));

    auto points = genPoints(n, bounds, dice);

    auto ret = triangulate(bounds, baseCase, points, splitter);

    std::cout << "RESULT"
              << " n=" << n << " splitter=" << splitter
              << " basecase=" << baseCase << " time=" << ret.time.count()
              << std::endl;
  };

  return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {

  // define program options

  namespace po = boost::program_options;

  unsigned char p;
  uint N;
  uint baseCase;

  po::options_description cCommandLine("Command Line Options");
  cCommandLine.add_options()("n", po::value<uint>(&N),
                             "maximum number of points");
  cCommandLine.add_options()("basecase", po::value<uint>(&baseCase),
                             "threshold for base case");
  cCommandLine.add_options()(
      "splitter", po::value<unsigned char>(&p),
      "splitter - _c_ycle, _d_-dimensional, _[0-d-1]_ fixed dimension");
  cCommandLine.add_options()("help", "produce help message");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cCommandLine), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << cCommandLine << std::endl;
    return EXIT_SUCCESS;
  }

  // plausability checks
  bool valid = true;
  if ((!vm.count("basecase"))) {
    std::cout << "Please specify basecase threshold" << std::endl;
    valid = false;
  }

  if ((!vm.count("splitter"))) {
    std::cout << "Please specify splitter" << std::endl;
    valid = false;
  }

  if (!valid)
    return EXIT_FAILURE;

  //***************************************************************************

  dBox<D, Precision> bounds;
  for (uint i = 0; i < D; ++i) {
    bounds.low[i] = 0;
    bounds.high[i] = 100;
  }

  std::uniform_real_distribution<Precision> distribution(0, 1);
  std::function<Precision()> dice = std::bind(distribution, generator);

  return benchmark(N, bounds, baseCase, p, std::move(dice));
}

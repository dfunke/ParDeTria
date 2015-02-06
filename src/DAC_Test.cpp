
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
#include <boost/progress.hpp>

/* Test the output of orient3D and inSphere3D
 */

void testGeometry3D() {

  dPoints<3> points;

  const uint A = 0;
  const uint B = 1;
  const uint C = 2;

  const uint Du = 3;
  const uint Dd = 4;

  const uint I = 5;
  const uint O = 6;

  // A (0,0,0)
  points[A].id = A;
  points[A].coords = {{0, 0, 0}};

  // B (2,0,0)
  points[B].id = B;
  points[B].coords = {{2, 0, 0}};

  // C (1,2,0)
  points[C].id = C;
  points[C].coords = {{1, 2, 0}};

  // D_u (1,1,1)
  points[Du].id = Du;
  points[Du].coords = {{1, 1, 1}};

  // D_d (1,1,-1)
  points[Dd].id = Dd;
  points[Dd].coords = {{1, 1, -1}};

  // I (1,1,0)
  points[I].id = I;
  points[I].coords = {{1, 1, 0}};

  // O (5,5,5)
  points[O].id = O;
  points[O].coords = {{5, 5, 5}};

  auto test = [&](const std::array<uint, 4> &v) {
    dSimplex<3> s;
    s.vertices = v;
    std::cout << s << std::endl;
    std::cout << "orientation: " << s.orientation(points) << std::endl;
    std::cout << "inSphere(I): " << s.inSphere(points[I], points) << std::endl;
    std::cout << "inSphere(O): " << s.inSphere(points[O], points) << std::endl;

    dSphere<3> cs = s.circumsphere(points);
    std::cout << "circumsphere: " << cs << std::endl;
  };

  // triangle ABC is counter-clockwise
  // D_d is below -> orientation > 0
  test({{A, B, C, Dd}});

  // triangle ABC is counter-clockwise
  // D_u is above -> orientation < 0
  test({{A, B, C, Du}});

  // triangle ACB is clockwise
  // D_d is below -> orientation < 0
  test({{A, C, B, Dd}});

  // triangle ACB is clockwise
  // D_u is above -> orientation > 0
  test({{A, C, B, Du}});
}

//**************************

#define D 3

//**************************

struct TriangulateReturn {
  TriangulationReport tr;
  bool exception;
  tDuration time;
};

//**************************

TriangulateReturn triangulate(const dBox<D> &bounds, dPoints<D> &points,
                              const unsigned char splitter) {
  std::unique_ptr<Partitioner<D>> partitioner_ptr;
  switch (splitter) {
  case 'd':
    partitioner_ptr = std::make_unique<dPartitioner<D>>();
    break;
  case 'c':
    partitioner_ptr = std::make_unique<CyclePartitioner<D>>();
    break;
  default:
    // p must be a dimension - subtract '0' to get integer value
    uint d = splitter - '0';
    ASSERT(0 <= d && d < D);
    partitioner_ptr = std::make_unique<kPartitioner<D>>(d);
    break;
  }

  Triangulator<D> triangulator(bounds, points, std::move(partitioner_ptr));

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

int main(int argc, char *argv[]) {

  // define program options

  namespace po = boost::program_options;

  int verbosity = -1;

  unsigned char p = 'd';
  uint N = 1e3;

  bool study = false;

  std::string pointFile;
  bool loadPoints = false;

  std::string simplexFile;
  bool loadSimplices = false;

  po::options_description cCommandLine("Command Line Options");
  cCommandLine.add_options()("n", po::value<uint>(&N), "number of points");
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
  cCommandLine.add_options()("simplices", po::value<std::string>(&simplexFile),
                             "load simplices from file");
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

  if (vm.count("points")) {
    loadPoints = true; // filename will be in pointFile
  }

  if (vm.count("simplices")) {
    loadSimplices = true; // filename will be in simplexFile
  }

  if (vm.count("study")) {
    study = true;
  }

  if (vm.count("no-verify")) {
    Triangulator<D>::VERIFY = false;
  }

  if (vm.count("no-output")) {
    Painter<D>::ENABLED = false;
  }

  //***************************************************************************

  dBox<D> bounds;
  for (uint i = 0; i < D; ++i) {
    bounds.low[i] = 0;
    bounds.high[i] = 100;
  }

  std::uniform_real_distribution<tCoordinate> distribution(0, 1);
  std::function<tCoordinate()> dice = std::bind(distribution, generator);

  /*
   * File format for triangulation report
   * input_n splitter provenance base_case edge_triangulation nPoints
   * nSimplices
   * nEdgePoints nEdgeSimplices
   */

  if (study) {
    // define splitter test cases
    std::vector<unsigned char> splitters = {'d', 'c', '0', '1'};

    if (vm.count("verbosity")) {
      // silence std::cout until STUDY is done - preserve log level
      LOGGER.setLogLevel(Logger::abs(LOGGER.getLogLevel()));
    } else {
      LOGGER.setLogLevel(Logger::Verbosity::SILENT);
    }

    Painter<D>::ENABLED = false;

    // abort flag if an exception occurred
    bool abort = false;

    boost::progress_display progress(
        splitters.size() * (9 * (log10(N) - 1) + 1), std::cout);

    std::ofstream f("triangulation_report.csv",
                    std::ios::out | std::ios::trunc);
    f << CSV::csv("n", "splitter", "provenance", "base_case", "edge_tria",
                  "valid", "nPoints", "nSimplices", "nEdgePoints",
                  "nEdgeSimplices", "time", "mem", "max_mem") << std::endl;

    for (uint n = 10; n <= N && !abort; n += pow(10, floor(log10(n)))) {
      for (uint i = 0; i < splitters.size() && !abort; ++i) {

        unsigned char p = splitters[i];

        dPoints<D> points;
        if (loadPoints) {
          points = loadObject<dPoints<D>>(pointFile);
        } else {
          points = genPoints(n, bounds, dice);
        }

        auto ret = triangulate(bounds, points, p);

        if (ret.exception) {
          abort = true;
        }

        // evaluate the triangulation report

        for (const auto &tr : ret.tr) {
          f << CSV::csv(points.size(), p, tr.provenance, tr.base_case,
                        tr.edge_triangulation, tr.valid, tr.nPoints,
                        tr.nSimplices, tr.nEdgePoints, tr.nEdgeSimplices,
                        ret.time.count(), getCurrentRSS(), getPeakRSS())
            << std::endl;

          if (Triangulator<D>::VERIFY &&
              Triangulator<D>::isTOP(tr.provenance) && !tr.valid) {

            // output points
            storeObject(points, "invalid_" + std::to_string(n) + "_" + (char)p +
                                    ".dat");
          }
        }

        ++progress;
      }
    }

    std::cout << LOGGER << std::endl;

    LOG << "Finished" << std::endl;

    return abort ? EXIT_FAILURE : EXIT_SUCCESS;

  } else {

    dPoints<D> points;
    if (loadPoints) {
      points = loadObject<dPoints<D>>(pointFile);
    } else {
      points = genPoints(N, bounds, dice);
    }

    auto ret = triangulate(bounds, points, p);

    LOG << "Triangulating " << points.size() << " points took "
        << std::chrono::duration_cast<std::chrono::seconds>(ret.time).count()
        << " s" << std::endl;

    VLOG << std::endl << "Triangulation report: " << std::endl;

    VLOG << CSV::csv("n", "splitter", "provenance", "base_case", "edge_tria",
                     "valid", "nPoints", "nSimplices", "nEdgePoints",
                     "nEdgeSimplices", "time", "mem", "max_mem") << std::endl;
    for (const auto &tr : ret.tr) {
      VLOG << CSV::csv(points.size(), p, tr.provenance, tr.base_case,
                       tr.edge_triangulation, tr.valid, tr.nPoints,
                       tr.nSimplices, tr.nEdgePoints, tr.nEdgeSimplices,
                       ret.time.count(), getCurrentRSS(), getPeakRSS())
           << std::endl;
    }

    return EXIT_SUCCESS;
  }
}

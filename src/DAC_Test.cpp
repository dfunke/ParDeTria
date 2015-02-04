
// own
#include "Geometry.h"
#include "Triangulator.h"
#include "Partitioner.h"

#include "utils/Random.h"
#include "utils/Logger.h"
#include "utils/CSV.h"
#include "utils/Timings.h"
#include "utils/RSS.h"
#include "utils/ASSERT.h"

#ifdef STUDY
#include "utils/Serialization.h"
#include <boost/progress.hpp>
#endif

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
const uint N = 1e3;

int main(int argc, char *argv[]) {

  if (argc == 2) {
    LOGGER.setLogLevel(static_cast<Logger::Verbosity>(std::stoi(argv[1])));
    std::cout << "Output level set to " << argv[1] << std::endl;
  } else
    LOGGER.setLogLevel(Logger::Verbosity::LIVE);

  dBox<D> bounds;
  for (uint i = 0; i < D; ++i) {
    bounds.coords[i] = 0;
    bounds.dim[i] = 100;
  }

  std::uniform_real_distribution<tCoordinate> distribution(0, 1);
  std::function<tCoordinate()> dice = std::bind(distribution, generator);

/*
 * File format for triangulation report
 * input_n splitter provenance base_case edge_triangulation nPoints
 * nSimplices
 * nEdgePoints nEdgeSimplices
 */

#ifdef STUDY
  // define splitter test cases
  std::vector<unsigned char> splitters = {'d', 'c', '0', '1'};
  // silence std::cout until STUDY is done - preserve log level
  LOGGER.setLogLevel(Logger::abs(LOGGER.getLogLevel()));

  // abort flag if an exception occurred
  bool abort = false;

  boost::progress_display progress(splitters.size() * (9 * (log10(N) - 1) + 1),
                                   std::cout);

  std::ofstream f("triangulation_report.csv", std::ios::out | std::ios::trunc);
  f << CSV::csv("n", "splitter", "provenance", "base_case", "edge_tria",
                "valid", "nPoints", "nSimplices", "nEdgePoints",
                "nEdgeSimplices", "time", "mem", "max_mem") << std::endl;

  for (uint n = 10; n <= N && !abort; n += pow(10, floor(log10(n)))) {
    for (uint i = 0; i < splitters.size(); ++i) {

      unsigned char p = splitters[i];
#else  // STUDY
  unsigned char p = 'c';
  uint n = N;
#endif // STUDY

      auto points = genPoints(n, bounds, dice);

      std::unique_ptr<Partitioner<D>> partitioner_ptr;
      switch (p) {
      case 'd':
        partitioner_ptr = std::make_unique<dPartitioner<D>>();
        break;
      case 'c':
        partitioner_ptr = std::make_unique<CyclePartitioner<D>>();
        break;
      default:
        // p must be a dimension - subtract '0' to get integer value
        uint d = p - '0';
        ASSERT(0 <= d && d < D);
        partitioner_ptr = std::make_unique<kPartitioner<D>>(d);
        break;
      }

      Triangulator<D> triangulator(bounds, points, std::move(partitioner_ptr));

      Clock::time_point t1, t2;
      INDENT
#ifdef STUDY
      try {
#endif

        t1 = Clock::now();
        auto dt = triangulator.triangulate();
        t2 = Clock::now();

#ifdef STUDY
      } catch (AssertionException e) {
        std::cerr << "Assertion failed - ABORTING - n =  " << n << " p = " << p
                  << std::endl;
        std::cerr << e.what() << std::endl;

        // output points
        storeObject(points, "testPoints_" + std::to_string(n) + "_" +
                                std::to_string(p) + ".dat");

        // set abort flag
        abort = true;
      }
#endif

      LOG << "Triangulating " << n << " points took "
          << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
          << " s" << std::endl;
      DEDENT

#ifdef STUDY
      // evaluate the triangulation report
      const auto &trReport = triangulator.getTriangulationReport();

      for (const auto &tr : trReport) {
        f << CSV::csv(points.size(), p, tr.provenance, tr.base_case,
                      tr.edge_triangulation, tr.valid, tr.nPoints,
                      tr.nSimplices, tr.nEdgePoints, tr.nEdgeSimplices,
                      std::chrono::duration_cast<tDuration>(t2 - t1).count(),
                      getCurrentRSS(), getPeakRSS()) << std::endl;
      }

      ++progress;
    }
  }

  std::cout << LOGGER << std::endl;
#endif // STUDY

  LOG << "Finished" << std::endl;

#ifdef STUDY
  return abort ? EXIT_FAILURE : EXIT_SUCCESS;
#else
  return EXIT_SUCCESS;
#endif
}


// own
#include "Random.h"
#include "Geometry.h"
#include "Logger.h"
#include "Triangulator.h"
#include "Partitioner.h"
#include "CSV.h"

// boost
#include <boost/progress.hpp>

const uint N = 1e3;
//**************************

int main(int argc, char *argv[]) {
  if (argc == 2)
    LOGGER.setLogLevel(static_cast<Logger::Verbosity>(std::stoi(argv[1])));
  else
    LOGGER.setLogLevel(Logger::Verbosity::LIVE);

  dBox bounds;
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
  std::vector<unsigned char> splitters = {'d', 'c', 0, 1};
  // silence std::cout until STUDY is done - preserve log level
  LOGGER.setLogLevel(Logger::abs(LOGGER.getLogLevel()));

  boost::progress_display progress(splitters.size() * (9 * (log10(N) - 1)),
                                   std::cout);

  std::ofstream f("triangulation_report.csv", std::ios::out | std::ios::trunc);
  f << CSV::csv("n", "splitter", "provenance", "base_case", "edge_tria",
                "valid", "nPoints", "nSimplices", "nEdgePoints",
                "nEdgeSimplices") << std::endl;

  for (uint n = 10; n < N; n += pow(10, floor(log10(n)))) {
    for (auto p : splitters) {
#else  // STUDY
  unsigned char p = 'd';
  uint n = N;
#endif // STUDY

      auto points = genPoints(n, bounds, dice);

      std::unique_ptr<Partitioner> partitioner_ptr;
      switch (p) {
      case 'd':
        partitioner_ptr = std::make_unique<dPartitioner>();
        break;
      case 'c':
        partitioner_ptr = std::make_unique<CyclePartitioner>();
        break;
      default:
        assert(0 <= p && p < D);
        partitioner_ptr = std::make_unique<kPartitioner>(p);
        break;
      }

      Triangulator triangulator(bounds, points, std::move(partitioner_ptr));

      INDENT
      auto dt = triangulator.triangulate();
      DEDENT

#ifdef STUDY
      // evaluate the triangulation report
      const auto &trReport = triangulator.getTriangulationReport();

      for (const auto &tr : trReport) {
        f << CSV::csv(points.size(), p, tr.provenance, tr.base_case,
                      tr.edge_triangulation, tr.valid, tr.nPoints,
                      tr.nSimplices, tr.nEdgePoints, tr.nEdgeSimplices)
          << std::endl;
      }

      ++progress;
    }
  }

  std::cout << LOGGER << std::endl;
#endif // STUDY

  LOG << "Finished" << std::endl;
}

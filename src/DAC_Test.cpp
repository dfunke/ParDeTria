
// own
#include "Random.h"
#include "Geometry.h"
#include "Logger.h"
#include "Triangulator.h"
#include "Partitioner.h"

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

  auto points = genPoints(N, bounds, dice);

  Triangulator triangulator(bounds, points, CyclePartitioner());

  INDENT
  auto dt = triangulator.triangulate();
  DEDENT

  LOG << "Finished" << std::endl;
}

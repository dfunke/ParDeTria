#include <fstream>
#include <utils/CSV.h>
#include <utils/Random.h>

#include "Triangulator.h"
#include "utils/CSV.h"

template <uint D, typename Precision>
class EdgeStudy : public Triangulator<D, Precision> {

public:
  EdgeStudy(std::ofstream &_out, const dBox<D, Precision> &_bounds,
            dPoints<D, Precision> &_points,
            std::unique_ptr<Partitioner<D, Precision>> &&_partitioner = nullptr)
      : Triangulator<D, Precision>(_bounds, _points, std::move(_partitioner)),
        out(_out) {

    Ids allPoints(this->points.begin_keys(), this->points.end_keys());
    auto DT = this->triangulateBase(allPoints, "e");

    auto edgeSimplices = this->getEdge(DT, this->bounds);
    auto edgePoints = this->extractPoints(edgeSimplices, DT);

    auto s =
        CSV::csv(this->points.size(), edgeSimplices.size(), edgePoints.size());

    out << s << std::endl;
    LOG << s << std::endl;
  }

private:
  std::ostream &out;
};

#define D 3
#define Precision double

#define N 1e5

int main() {

  LOGGER.setLogLevel(Logger::Verbosity::LIVE);

  std::ofstream f("edge_study.csv", std::ios::out | std::ios::trunc);
  f << CSV::csv("n", "eS", "eP") << std::endl;

  dBox<D, Precision> bounds;
  for (uint i = 0; i < D; ++i) {
    bounds.low[i] = 0;
    bounds.high[i] = 100;
  }

  std::uniform_real_distribution<Precision> distribution(0, 1);
  std::function<Precision()> dice = std::bind(distribution, generator);

  Triangulator<D, Precision>::VERIFY = false;

  for (uint n = 10; n <= N; n += pow(10, floor(log10(n)))) {
    auto points = genPoints(n, bounds, dice);
    EdgeStudy<D, Precision> es(f, bounds, points);
  }

  return EXIT_SUCCESS;
}

#include <fstream>
#include <utils/CSV.h>
#include <utils/Random.h>

#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <future>
#include <bits/stl_deque.h>

#include "Triangulator.h"
#include "utils/CSV.h"
#include "Painter.h"

#define LLOG(msg) std::cout << msg;

template <uint D, typename Precision>
class WhereUsedListSizeStudy : Triangulator<D, Precision> {

public:
  WhereUsedListSizeStudy(
      std::ofstream &_out, std::mutex &_outMtx,
      const dBox<D, Precision> &_bounds, dPoints<D, Precision> &_points,
      std::unique_ptr<Partitioner<D, Precision>> &&_partitioner)
      : Triangulator<D, Precision>(_bounds, 500, _points,
                                   std::move(_partitioner)),
        out(_out), outMtx(_outMtx) {

    auto DT = this->triangulate();

    std::stringstream ss;

    for (const auto &p : this->points) {
      uint countInfiniteS = 0;
      for (const auto &s : p.simplices) {
        countInfiniteS += !dSimplex<D, Precision>::isFinite(s);
      }

      ss << CSV::csv(this->points.size(), p.simplices.size(), countInfiniteS)
         << '\n';
    }

    {
      std::lock_guard<std::mutex> lock(outMtx);
      out << ss.str() << std::endl;
    }
  }

private:
  std::ostream &out;
  std::mutex &outMtx;
};

template <uint D, typename Precision>
void studyWhereUsedListSize(const uint N, const char splitter) {

  std::ofstream f("whereUsed_study_set_" + std::to_string(D) + "_" + splitter +
                      ".csv",
                  std::ios::out | std::ios::trunc);
  f << CSV::csv("nP", "nWU", "nIWU") << std::endl;

  dBox<D, Precision> bounds;
  for (uint i = 0; i < D; ++i) {
    bounds.low[i] = 0;
    bounds.high[i] = 100;
  }

  std::uniform_real_distribution<Precision> distribution(0, 1);
  std::function<Precision()> dice = std::bind(distribution, generator);

  Triangulator<D, Precision>::VERIFY = false;
  Painter<D, Precision>::ENABLED = false;

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
    uint d = splitter - '0';
    ASSERT(0 <= d && d < D);
    partitioner_ptr = std::make_unique<kPartitioner<D, Precision>>(d);
    break;
  }

  std::mutex mtx;
  tbb::parallel_for(
      std::size_t(0), std::size_t((9 * (log10(N) - 1) + 1)), [&](const uint i) {
        uint n = ((i % 9) + 1) * std::pow(10, std::floor(i / 9 + 1));
        LLOG("Testing with " << n << " points" << std::endl);
        auto points = genPoints(n, bounds, dice);
        WhereUsedListSizeStudy<D, Precision> es(f, mtx, bounds, points,
                                                std::move(partitioner_ptr));
      });
}

#define Precision double
#define N 1e4

int main() {

  LOGGER.setLogLevel(Logger::Verbosity::SILENT);
  const unsigned char splitter = 'c';

  std::async(std::launch::async, [&]() {
    LLOG("2D" << std::endl);
    INDENT
    studyWhereUsedListSize<2, Precision>(N, splitter);
    DEDENT
  });

  std::async(std::launch::async, [&]() {
    LLOG("3D" << std::endl);
    INDENT
    studyWhereUsedListSize<3, Precision>(N, splitter);
    DEDENT
  });

  return EXIT_SUCCESS;
}

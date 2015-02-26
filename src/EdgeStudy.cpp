#include <fstream>
#include <utils/CSV.h>
#include <utils/Random.h>

#include <tbb/parallel_for.h>
#include <future>
#include <bits/stl_deque.h>

#include "Triangulator.h"
#include "utils/CSV.h"

template <uint D, typename Precision>
class EdgeSizeStudy : public Triangulator<D, Precision> {

public:
  EdgeSizeStudy(
      std::ofstream &_out, std::mutex &_outMtx,
      const dBox<D, Precision> &_bounds, dPoints<D, Precision> &_points,
      std::unique_ptr<Partitioner<D, Precision>> &&_partitioner = nullptr)
      : Triangulator<D, Precision>(_bounds, _points, std::move(_partitioner)),
        out(_out), outMtx(_outMtx) {

    Ids allPoints(this->points.begin_keys(), this->points.end_keys());

    INDENT
    auto DT = this->triangulateBase(allPoints, "e");
    DEDENT

    auto edgeSimplices = this->getEdge(DT, this->bounds);
    auto edgePoints = this->extractPoints(edgeSimplices, DT);

    auto s = CSV::csv(this->points.size(), DT.size(), edgeSimplices.size(),
                      edgePoints.size());

    LOG << s << std::endl;
    {
      std::lock_guard<std::mutex> lock(outMtx);
      out << s << std::endl;
    }
  }

private:
  std::ostream &out;
  std::mutex &outMtx;
};

template <uint D, typename Precision>
class PartitionEdgeSizeStudy : public Triangulator<D, Precision> {

public:
  PartitionEdgeSizeStudy(
      std::ofstream &_out, std::mutex &_outMtx,
      const dBox<D, Precision> &_bounds, dPoints<D, Precision> &_points,
      std::unique_ptr<Partitioner<D, Precision>> &&_partitioner)
      : Triangulator<D, Precision>(_bounds, _points, std::move(_partitioner)),
        out(_out), outMtx(_outMtx) {

    Ids allPoints(this->points.begin_keys(), this->points.end_keys());

    auto partitioning =
        this->partitioner->partition(allPoints, this->points, "0");

    std::vector<dSimplices<D, Precision>> partialDTs;
    partialDTs.resize(partitioning.size());

    Ids edgePoints;
    Ids edgeSimplices;

    INDENT
    for (uint i = 0; i < partitioning.size(); ++i) {
      partialDTs[i] = this->triangulateBase(partitioning[i].points,
                                            "0" + std::to_string(i));

      std::vector<dBox<D, Precision>> pBounds;
      for (uint j = 0; j < partitioning.size(); ++j) {
        if (i != j)
          pBounds.push_back(partitioning[j].bounds);
      }

      auto eS = getEdge(partialDTs[i], pBounds);
      edgeSimplices.insert(eS.begin(), eS.end());

      auto eP = this->extractPoints(eS, partialDTs[i]);
      edgePoints.insert(eP.begin(), eP.end());
    }
    DEDENT

    auto s =
        CSV::csv(this->points.size(), edgeSimplices.size(), edgePoints.size());

    LOG << s << std::endl;
    {
      std::lock_guard<std::mutex> lock(outMtx);
      out << s << std::endl;
    }
  }

private:
  Ids getEdge(const dSimplices<D, Precision> &simplices,
              const std::vector<dBox<D, Precision>> &pBounds) {
    Ids edgeSimplices;
    std::set<uint> wqa; // set of already checked simplices

    // we use the overflow of the uint to zero to abort the loop
    for (uint infVertex = dPoint<D, Precision>::cINF; infVertex != 0;
         ++infVertex) {

      ASSERT(this->points.contains(infVertex));

      PLOG << "Infinite vertex " << this->points[infVertex] << std::endl;
      LOGGER.logContainer(this->points[infVertex].simplices,
                          Logger::Verbosity::PROLIX, "Used in simplices: ");

      INDENT
      for (const auto &s : this->points[infVertex].simplices) {
        if (!simplices.contains(s))
          continue;

        PLOG << "Adding " << simplices[s] << " to edge" << std::endl;
        edgeSimplices.insert(simplices[s].id);

        /* Walk along the neighbors,
         * testing for each neighbor whether its circumsphere is within the
         * partition or not
         */
        wqa.insert(simplices[s].id);
        std::deque<uint> wq;

        // work queue of simplices to check for circum circle criterian
        for (const auto &n : simplices[s].neighbors) {
          if (wqa.insert(n).second) {
            // simplex not yet inspected -> add it to queue
            wq.push_back(n);
          }
        }

        INDENT
        while (!wq.empty()) {
          uint x = wq.front();
          wq.pop_front();

          if (simplices.contains(x)) {
            const auto cs = simplices[x].circumsphere(this->points);
            bool intersectsBounds = false;
            for (const auto &p : pBounds) {
              if (p.intersects(cs))
                intersectsBounds = true;
            }
            if (intersectsBounds) {
              PLOG << "Adding " << simplices[x]
                   << " to edge -> circumcircle criterion" << std::endl;
              edgeSimplices.insert(simplices[x].id);

              for (const auto &n : simplices[x].neighbors) {
                if (wqa.insert(n).second) {
                  // n was not yet inspected
                  wq.push_back(n);
                }
              }
            }
          }
        }
        DEDENT
      }
      DEDENT
    }

    return edgeSimplices;
  }

private:
  std::ostream &out;
  std::mutex &outMtx;
};

template <uint D, typename Precision> void studyEdgeSize(const uint N) {

  std::ofstream f("edge_study_" + std::to_string(D) + ".csv",
                  std::ios::out | std::ios::trunc);
  f << CSV::csv("nP", "nS", "eS", "eP") << std::endl;

  dBox<D, Precision> bounds;
  for (uint i = 0; i < D; ++i) {
    bounds.low[i] = 0;
    bounds.high[i] = 100;
  }

  std::uniform_real_distribution<Precision> distribution(0, 1);
  std::function<Precision()> dice = std::bind(distribution, generator);

  Triangulator<D, Precision>::VERIFY = false;

  std::mutex mtx;
  tbb::parallel_for(std::size_t(0), std::size_t((9 * (log10(N) - 1) + 1)),
                    [&](const uint i) {
    uint n = ((i % 9) + 1) * std::pow(10, std::floor(i / 9 + 1));
    auto points = genPoints(n, bounds, dice);
    EdgeSizeStudy<D, Precision> es(f, mtx, bounds, points);
  });
}

template <uint D, typename Precision>
void studyPartitionEdgeSize(const uint N, const char splitter) {

  std::ofstream f("partition_study_" + std::to_string(D) + "_" + splitter +
                      ".csv",
                  std::ios::out | std::ios::trunc);
  f << CSV::csv("nP", "eS", "eP") << std::endl;

  dBox<D, Precision> bounds;
  for (uint i = 0; i < D; ++i) {
    bounds.low[i] = 0;
    bounds.high[i] = 100;
  }

  std::uniform_real_distribution<Precision> distribution(0, 1);
  std::function<Precision()> dice = std::bind(distribution, generator);

  Triangulator<D, Precision>::VERIFY = false;

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
  tbb::parallel_for(std::size_t(0), std::size_t((9 * (log10(N) - 1) + 1)),
                    [&](const uint i) {
    uint n = ((i % 9) + 1) * std::pow(10, std::floor(i / 9 + 1));
    auto points = genPoints(n, bounds, dice);
    PartitionEdgeSizeStudy<D, Precision> es(f, mtx, bounds, points,
                                            std::move(partitioner_ptr));
  });
}

#define Precision double
#define N 1e4

int main() {

  LOGGER.setLogLevel(Logger::Verbosity::LIVE);
  std::vector<unsigned char> splitters = {'d', 'c'};

  std::async(std::launch::async, [&]() {
    LOG << "2D" << std::endl;
    INDENT
    // studyEdgeSize<2, Precision>(N);
    tbb::parallel_for(std::size_t(0), splitters.size(), [&](const uint i) {
      const auto s = splitters[i];
      LOG << "Splitter " << std::to_string(s) << std::endl;
      studyPartitionEdgeSize<2, Precision>(N, s);
    });
    DEDENT
  });

  std::async(std::launch::async, [&]() {
    LOG << "3D" << std::endl;
    INDENT
    // studyEdgeSize<3, Precision>(N);
    tbb::parallel_for(std::size_t(0), splitters.size(), [&](const uint i) {
      const auto s = splitters[i];
      LOG << "Splitter " << std::to_string(s) << std::endl;
      studyPartitionEdgeSize<3, Precision>(N, s);
    });
    DEDENT
  });

  return EXIT_SUCCESS;
}

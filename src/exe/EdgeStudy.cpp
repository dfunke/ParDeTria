#include <fstream>
#include <utils/CSV.h>
#include <utils/Random.h>

#include <tbb/parallel_for.h>
#include <future>
#include <bits/stl_deque.h>

#include "DCTriangulator.h"
#include "utils/CSV.h"

template<uint D, typename Precision>
class PartitionEdgeSizeStudy : public DCTriangulator<D, Precision> {

public:
    PartitionEdgeSizeStudy(
            std::ofstream &_out, std::mutex &_outMtx,
            const dBox<D, Precision> &_bounds, dPoints<D, Precision> &_points,
            const unsigned char splitter)
            : DCTriangulator<D, Precision>(_bounds, _points, 0, splitter),
              out(_out), outMtx(_outMtx) {

        auto partitioning =
                this->partitioner->partition(this->allPoints(), this->points, "0");

        std::vector<dSimplices<D, Precision>> partialDTs;
        partialDTs.resize(partitioning.size());

        Ids edgePoints;
        Ids edgeSimplices;

        INDENT
        for (uint i = 0; i < partitioning.size(); ++i) {
            partialDTs[i] =
                    this->_triangulateBase(partitioning[i].points, partitioning[i].bounds,
                                           "0" + std::to_string(i));

            auto eS = this->getEdge(partialDTs[i], partitioning, i);
            edgeSimplices.insert(eS.begin(), eS.end());

            auto eP = this->extractPoints(eS, partialDTs[i]);
            edgePoints.insert(eP.begin(), eP.end());
        }
        DEDENT

        auto s =
                CSV::csv(this->points.size(), edgeSimplices.size(), edgePoints.size());

        LOG(s << std::endl);
        {
            std::lock_guard<std::mutex> lock(outMtx);
            out << s << std::endl;
        }
    }

private:
    std::ostream &out;
    std::mutex &outMtx;
};

template<uint D, typename Precision>
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
    std::function<Precision()> dice = std::bind(distribution, startGen);

    DCTriangulator<D, Precision>::VERIFY = false;

    std::mutex mtx;
    tbb::parallel_for(
            std::size_t(0), std::size_t((9 * (log10(N) - 1) + 1)), [&](const uint i) {
        uint n = ((i % 9) + 1) * std::pow(10, std::floor(i / 9 + 1));
        LOG("Testing with " << n << " points" << std::endl);
        auto points = genPoints(n, bounds, dice);
        PartitionEdgeSizeStudy<D, Precision> es(f, mtx, bounds, points, splitter);
    });
}

#define Precision double
#define N 1e4

int main() {

    LOGGER.setLogLevel(Logger::Verbosity::LIVE);
    std::vector<unsigned char> splitters = {'d', 'c'};

    std::async(std::launch::async, [&]() {
        LOG("2D" << std::endl);
        INDENT
        tbb::parallel_for(std::size_t(0), splitters.size(), [&](const uint i) {
            const auto s = splitters[i];
            LOG("Splitter " << s << std::endl);
            studyPartitionEdgeSize<2, Precision>(N, s);
        });
        DEDENT
    });

    std::async(std::launch::async, [&]() {
        LOG("3D" << std::endl);
        INDENT
        tbb::parallel_for(std::size_t(0), splitters.size(), [&](const uint i) {
            const auto s = splitters[i];
            LOG("Splitter " << s << std::endl);
            studyPartitionEdgeSize<3, Precision>(N, s);
        });
        DEDENT
    });

    return EXIT_SUCCESS;
}

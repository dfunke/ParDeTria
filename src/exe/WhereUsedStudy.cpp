#include <fstream>
#include <utils/CSV.h>
#include <utils/Random.h>

#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <future>
#include <bits/stl_deque.h>

#include "DCTriangulator.h"
#include "utils/CSV.h"

#define LLOG(msg) std::cout << msg;

template<uint D, typename Precision>
class WhereUsedListSizeStudy : DCTriangulator<D, Precision> {

public:
    WhereUsedListSizeStudy(
            std::ofstream &_out, std::mutex &_outMtx,
            const dBox<D, Precision> &_bounds, dPoints<D, Precision> &_points,
            const unsigned char splitter)
            : DCTriangulator<D, Precision>(_bounds, _points, 500, splitter),
              out(_out), outMtx(_outMtx) {

        auto DT = this->triangulate();

        std::stringstream ss;

        for (auto it = DT.wuFaces.begin();
             it != DT.wuFaces.end();) {

            auto range = DT.wuFaces.equal_range(it->first);
            ss << CSV::csv(this->points.size(), std::distance(range.first, range.second)) << '\n';
            it = range.second;

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

template<uint D, typename Precision>
void studyWhereUsedListSize(const uint N, const char splitter) {

    std::ofstream f("whereUsed_collisions_" + std::to_string(D) + "_" +
                    splitter + ".csv",
                    std::ios::out | std::ios::trunc);
    f << CSV::csv("n", "collisions") << std::endl;

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
        LLOG("Testing with " << n << " points" << std::endl);
        auto points = genPoints(n, bounds, dice);
        WhereUsedListSizeStudy<D, Precision> es(f, mtx, bounds, points, splitter);
    });
}

#define Precision double

int main(int argc, char *argv[]) {

    LOGGER.setLogLevel(Logger::Verbosity::SILENT);
    const unsigned char splitter = 'c';

    int N = 1e4;
    if(argc == 2)
        N = std::stoi(argv[1]);

    /*std::async(std::launch::async, [&]() {
        LLOG("2D" << std::endl);
        INDENT
        studyWhereUsedListSize<2, Precision>(N, splitter);
        DEDENT
    });*/

    std::async(std::launch::async, [&]() {
        LLOG("3D" << std::endl);
        INDENT
        studyWhereUsedListSize<3, Precision>(N, splitter);
        DEDENT
    });

    return EXIT_SUCCESS;
}

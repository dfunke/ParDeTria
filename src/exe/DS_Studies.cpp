// own
#include "Geometry.h"
#include "DCTriangulator.h"
#include "Partitioner.h"

#include "utils/Random.h"
#include "utils/Generator.h"
#include "utils/Logger.h"
#include "utils/CSV.h"
#include "utils/Timings.h"
#include "utils/System.h"
#include "utils/ASSERT.h"
#include "utils/Serialization.hxx"
#include "utils/ProgressDisplay.h"
#include "utils/version.h"

#include <boost/program_options.hpp>
#include <boost/utility/in_place_factory.hpp>

#include <mutex>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

//**************************

#define D 3
#define Precision double

//**************************

struct TriangulateReturn {
    bool exception;
    bool valid;
    tDuration time;
    std::size_t currRSS;
    std::size_t peakRSS;
    ExperimentRun run;

    //
    std::size_t countSimplices;
    std::size_t countDeletedSimplices;
    std::size_t cvSize;
    std::size_t cvCapacity;

    double dtMB;
    double cvMB;
    double pointMB;
};

TriangulateReturn triangulate(const dBox<D, Precision> &bounds,
                              const uint recursionDepth,
                              dPoints<D, Precision> &points,
                              const unsigned char splitter,
                              const bool parallelBase,
                              const bool verify) {

    DCTriangulator<D, Precision> triangulator(bounds, points, recursionDepth, splitter, 100, parallelBase, true, true);

    TriangulateReturn ret;
    PROFILER.setRun(&ret.run);

    //try {

    auto t1 = Clock::now();
    auto dt = triangulator.triangulate();
    auto t2 = Clock::now();

    ret.currRSS = getCurrentRSS();
    ret.peakRSS = getPeakRSS();

    ret.exception = false;
    ret.time = std::chrono::duration_cast<tDuration>(t2 - t1);

    if (verify) {
        LOGGER.setIndent(0);
        LOG("Generate Reference Triangulation" << std::endl);

        INDENT;
        CGALTriangulator<D, Precision, false> cgal(bounds, points);
        auto realDT = cgal.triangulate();
        DEDENT;

        auto vr = dt.verify(points);
        auto ccr = dt.crossCheck(realDT);

        ret.valid = vr.valid && ccr.valid;
    } else {
        ret.valid = true;
    }

    std::atomic<std::size_t> aInf(0);
    std::atomic<std::size_t> aFin(0);

    tbb::parallel_for(tbb::blocked_range<std::size_t>(dt.lowerBound(), dt.upperBound(), 1e3), [&](const auto &range) {

        uint hint = 0;

        for (auto it = range.begin(); it != range.end(); ++it) {
            bool inf = dt.at(it, hint).id == dSimplex<D, Precision>::cINF;

            aInf += inf;
            aFin += !inf;
        }

    });

    ret.countSimplices = aFin.load();
    ret.countDeletedSimplices = aInf.load();

    ret.cvSize = dt.convexHull.size();
    ret.cvCapacity = dt.convexHull.capacity();

    ret.dtMB = ((dt.upperBound() - dt.lowerBound()) * sizeof(dSimplex<D, Precision>)) / 1e6;
    ret.pointMB = ((points.size()) * sizeof(dPoint<D, Precision>)) / 1e6;
    ret.cvMB = ((dt.convexHull.capacity()) * sizeof(tIdType)) / 1e6;

    return ret;
}

//**************************

int main(int argc, char *argv[]) {

    // define program options

    namespace po = boost::program_options;

    unsigned char p = 'c';
    tIdType minN;
    tIdType maxN;

    uint minThreads = 1;
    uint maxThreads = tbb::task_scheduler_init::default_num_threads();

    uint recursionDepth;
    uint threads = tbb::task_scheduler_init::default_num_threads();
    tIdType N;
    bool parallelBase = false;
    bool verify = false;

    po::options_description cCommandLine("Command Line Options");
    cCommandLine.add_options()("N", po::value<tIdType>(&N), "number of points");
    cCommandLine.add_options()("minN", po::value<tIdType>(&minN), "minimum number of points");
    cCommandLine.add_options()("maxN", po::value<tIdType>(&maxN), "maximum number of points");
    cCommandLine.add_options()("minThreads", po::value<uint>(&minThreads), "minimum number of threads = 1");
    cCommandLine.add_options()("maxThreads", po::value<uint>(&maxThreads), "maximum number of threads = #cores");
    cCommandLine.add_options()("recDepth", po::value<uint>(&recursionDepth),
                               "maximum levels of recursion");
    cCommandLine.add_options()(
            "splitter", po::value<unsigned char>(&p),
            "splitter - _c_ycle, _d_-dimensional, _[0-d-1]_ fixed dimension");
    cCommandLine.add_options()("threads", po::value<uint>(&threads),
                               "specify number of threads");
    cCommandLine.add_options()("par-base", po::value<bool>(&parallelBase),
                               "use parallel base solver");
    cCommandLine.add_options()("verify", po::value<bool>(&verify),
                               "verify triangulation");
    cCommandLine.add_options()("help", "produce help message");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cCommandLine), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << cCommandLine << std::endl;
        return EXIT_SUCCESS;
    }

    // evaluate program options
    LOGGER.setLogLevel(Logger::Verbosity::SILENT);

    // plausability checks
    bool valid = true;
    if (!vm.count("recDepth")) {
        recursionDepth = log2(threads);
    }

    if (!valid)
        return EXIT_FAILURE;


    //***************************************************************************

    dBox<D, Precision> bounds;
    for (uint i = 0; i < D; ++i) {
        bounds.low[i] = 0;
        bounds.high[i] = 100;
    }

    std::vector<unsigned char> dists = {'n', 'b'};
    if (vm.count("minN") && vm.count("maxN")) {
        std::cout << "Point study from " << minN << " to " << maxN << " points" << std::endl;

        // load scheduler with specified number of threads
        tbb::task_scheduler_init init(threads);
        
        std::ofstream f("ds_study_points_" + std::string(GIT_COMMIT) + ".csv");
        f << "dist n simplices delSimplices cvSize cvCap pointMB dtMB cvMB currRSS peakRSS time" << std::endl;

        for (auto dist : dists) {
            for (std::size_t n = minN; n <= maxN; n += pow(10, floor(log10(n)))) {
                auto pg = GeneratorFactory<D, Precision>::make(dist);
                auto points = pg->generate(n, bounds, startGen);

                TriangulateReturn ret = triangulate(bounds, recursionDepth, points, p, parallelBase, verify);

                f << dist << " " << n << " " << ret.countSimplices << " " << ret.countDeletedSimplices << " "
                << ret.cvSize << " " << ret.cvCapacity << " "
                << ret.pointMB << " " << ret.dtMB << " " << ret.cvMB << " "
                << ret.currRSS / 1e6 << " " << ret.peakRSS / 1e6 << " " << ret.time.count()
                << std::endl;

                std::cout << "\tTriangulating "
                << points.size() << " points took "
                << std::chrono::duration_cast<std::chrono::milliseconds>(ret.time).count() << " ms and "
                << ret.currRSS / 1e6 << "/" << ret.peakRSS / 1e6 << " MB "
                << "with " << threads << " threads  valid: " << ret.valid << std::endl;
            }
        }
    }

    if (vm.count("minThreads") && vm.count("maxThreads")) {
        std::cout << "Thread study from " << minThreads << " to " << maxThreads << " threads" << std::endl;

        std::ofstream f("ds_study_threads_" + std::string(GIT_COMMIT) + ".csv");
        f << "threads simplices delSimplices cvSize cvCap pointMB dtMB cvMB currRSS peakRSS time" << std::endl;

        if (vm.count("maxN")) {
            N = maxN / 2;
        } else {
            if (!vm.count("N")) {
                std::cout << "no number of points specified" << std::endl;
                return EXIT_FAILURE;
            }
        }

        for (uint t = minThreads; t <= maxThreads; t <<= 1) {
            auto pg = GeneratorFactory<D, Precision>::make('u');
            auto points = pg->generate(N, bounds, startGen);

            tbb::task_scheduler_init init(t);
            uint recD = log2(t);

            TriangulateReturn ret = triangulate(bounds, recD, points, p, parallelBase, verify);

            f << t << " " << ret.countSimplices << " " << ret.countDeletedSimplices << " "
            << ret.cvSize << " " << ret.cvCapacity << " "
            << ret.pointMB << " " << ret.dtMB << " " << ret.cvMB << " "
            << ret.currRSS / 1e6 << " " << ret.peakRSS / 1e6 << " " << ret.time.count()
            << std::endl;

            std::cout << "\tTriangulating "
            << points.size() << " points took "
            << std::chrono::duration_cast<std::chrono::milliseconds>(ret.time).count() << " ms and "
            << ret.currRSS / 1e6 << "/" << ret.peakRSS / 1e6 << " MB "
            << "with " << t << " threads  valid: " << ret.valid << std::endl;
        }
    }

    return EXIT_SUCCESS;
}

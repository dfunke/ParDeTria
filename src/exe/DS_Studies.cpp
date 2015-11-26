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
    std::size_t countSimplices = 0;
    std::size_t countDeletedSimplices = 0;
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
                              const bool parallelBase) {

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
    ret.valid = true;

    for (std::size_t i = dt.lowerBound(); i < dt.upperBound(); ++i) {
        ret.countSimplices += dt.unsafe_at(i).id != dSimplex<D, Precision>::cINF;
        ret.countDeletedSimplices += dt.unsafe_at(i).id == dSimplex<D, Precision>::cINF;
    }

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
    uint recursionDepth;
    uint threads = tbb::task_scheduler_init::default_num_threads();
    bool parallelBase = false;

    po::options_description cCommandLine("Command Line Options");
    cCommandLine.add_options()("minN", po::value<tIdType>(&minN), "minimum number of points");
    cCommandLine.add_options()("maxN", po::value<tIdType>(&maxN), "maximum number of points");
    cCommandLine.add_options()("recDepth", po::value<uint>(&recursionDepth),
                               "maximum levels of recursion");
    cCommandLine.add_options()(
            "splitter", po::value<unsigned char>(&p),
            "splitter - _c_ycle, _d_-dimensional, _[0-d-1]_ fixed dimension");
    cCommandLine.add_options()("threads", po::value<uint>(&threads),
                               "specify number of threads");
    cCommandLine.add_options()("par-base", po::value<bool>(&parallelBase),
                               "use parallel base solver");
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

    if (!vm.count("minN") || !vm.count("maxN")) {
        std::cout << "Please specify number of points or point file" << std::endl;
        valid = false;
    }

    if (!valid)
        return EXIT_FAILURE;

    // load scheduler with specified number of threads
    tbb::task_scheduler_init init(threads);

    //***************************************************************************

    dBox<D, Precision> bounds;
    for (uint i = 0; i < D; ++i) {
        bounds.low[i] = 0;
        bounds.high[i] = 100;
    }

    std::ofstream f("ds_study_" + std::string(GIT_COMMIT) + ".csv");
    f << "n simplices delSimplices cvSize cvCap pointMB dtMB cvMB currRSS peakRSS time" << std::endl;

    for (std::size_t n = minN; n <= maxN; n += pow(10, floor(log10(n)))) {
        auto pg = GeneratorFactory<D, Precision>::make('u');
        auto points = pg->generate(n, bounds, startGen);

        TriangulateReturn ret = triangulate(bounds, recursionDepth, points, p, parallelBase);

        f << n << " " << ret.countSimplices << " " << ret.countDeletedSimplices << " "
          << ret.cvSize << " " << ret.cvCapacity << " "
          << ret.pointMB << " " << ret.dtMB << " " << ret.cvMB << " "
          << ret.currRSS / 1e6 << " " << ret.peakRSS / 1e6 << " " << ret.time.count()
          << std::endl;

        std::cout << "Triangulating "
        << points.size() << " points took "
        << std::chrono::duration_cast<std::chrono::milliseconds>(ret.time)
                .count() << " ms and " << ret.currRSS / 1e6 << "/" << ret.peakRSS / 1e6 << " MB" << std::endl;
    }

    return 0;
}

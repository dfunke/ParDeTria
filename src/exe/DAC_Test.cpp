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
    ExperimentRun run;
};

TriangulateReturn triangulate(const dBox<D, Precision> &bounds,
                              const uint recursionDepth,
                              dPoints<D, Precision> &points,
                              const unsigned char splitter,
                              const unsigned char alg = 'd',
                              const bool verify = true) {

    std::unique_ptr<Triangulator<D, Precision>> triangulator_ptr;
    if (alg == 'c') {
        triangulator_ptr =
                std::make_unique<PureCGALTriangulator<D, Precision, false>>(bounds, points);
    } else {
        if (alg == 'm') {
            triangulator_ptr =
                    std::make_unique<PureCGALTriangulator<D, Precision, true>>(bounds, points, 100);
        } else {
            triangulator_ptr =
                    std::make_unique<DCTriangulator<D, Precision>>(bounds, points, recursionDepth, splitter, 100,
                                                                   false);
        }
    }

    TriangulateReturn ret;
    PROFILER.setRun(&ret.run);

    //try {

    auto t1 = Clock::now();
    auto dt = triangulator_ptr->triangulate();
    auto t2 = Clock::now();

    ret.exception = false;
    ret.time = std::chrono::duration_cast<tDuration>(t2 - t1);

    if (verify) {
        LOGGER.setIndent(0);
        LOG("Generate Reference Triangulation" << std::endl);

        INDENT;
        CGALTriangulator<D, Precision, false> cgal(bounds, points);
        auto realDT = cgal.triangulate();
        DEDENT;

        auto vr = dt.first.verify(dt.second, points);
        auto ccr = dt.first.crossCheck(dt.second, realDT.first, realDT.second);

        ret.valid = vr.valid && ccr.valid;
    } else {
        ret.valid = true;
    }

    /*} catch (AssertionException &e) {
      std::cerr << "Assertion failed - ABORTING - n =  " << points.size()
                << " p = " << splitter << std::endl;
      std::cerr << e.what() << std::endl;

      // output points
      storeObject(points, "assertionFailed_" + std::to_string(points.size()) +
                              "_" + (char)splitter + ".dat");

      ret.exception = true;
    }*/

    return ret;
}

//**************************

int main(int argc, char *argv[]) {

    // define program options

    namespace po = boost::program_options;

    int verbosity = -1;

    unsigned char p;
    uint N;
    uint recursionDepth;
    uint threads = tbb::task_scheduler_init::default_num_threads();
    unsigned char alg = 'd';

    bool verify = true;

    std::string pointFile;
    bool loadPoints = false;

    po::options_description cCommandLine("Command Line Options");
    cCommandLine.add_options()("n", po::value<uint>(&N), "number of points");
    cCommandLine.add_options()("recDepth", po::value<uint>(&recursionDepth),
                               "maximum levels of recursion");
    cCommandLine.add_options()(
            "splitter", po::value<unsigned char>(&p),
            "splitter - _c_ycle, _d_-dimensional, _[0-d-1]_ fixed dimension");
    cCommandLine.add_options()("no-verify", "don't verify triangulation");
    cCommandLine.add_options()("no-output", "don't write images");
    cCommandLine.add_options()("seq-fault", "run triangulation until seq fault occurs");
    cCommandLine.add_options()("verbosity", po::value<int>(&verbosity),
                               "verbosity");
    cCommandLine.add_options()("points", po::value<std::string>(&pointFile),
                               "load points from file");
    cCommandLine.add_options()("threads", po::value<uint>(&threads),
                               "specify number of threads");
    cCommandLine.add_options()("alg", po::value<unsigned char>(&alg),
                               "specify algorithm to use");
    cCommandLine.add_options()("help", "produce help message");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cCommandLine), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << cCommandLine << std::endl;
        return EXIT_SUCCESS;
    }

    // evaluate program options
    LOGGER.setLogLevel(static_cast<Logger::Verbosity>(verbosity));

    // plausability checks
    bool valid = true;
    if ((!vm.count("recDepth"))) {
        std::cout << "Please specify the maximum depth of recursion" << std::endl;
        valid = false;
    }

    if (!vm.count("splitter")) {
        std::cout << "Please specify splitter" << std::endl;
        valid = false;
    }

    if (!(vm.count("n") || vm.count("points"))) {
        std::cout << "Please specify number of points or point file" << std::endl;
        valid = false;
    }

    if (!valid)
        return EXIT_FAILURE;

    if (vm.count("points")) {
        loadPoints = true; // filename will be in pointFile
    }

    verify = !vm.count("no-verify");

    // load scheduler with specified number of threads
    tbb::task_scheduler_init init(threads);

    //***************************************************************************

    dBox<D, Precision> bounds;
    for (uint i = 0; i < D; ++i) {
        bounds.low[i] = 0;
        bounds.high[i] = 100;
    }

    std::uniform_real_distribution<Precision> distribution(0, 1);
    std::function<Precision()> dice = std::bind(distribution, startGen);

    /*
     * File format for triangulation report
     * input_n splitter provenance base_case edge_triangulation nPoints
     * nSimplices
     * nEdgePoints nEdgeSimplices
     */

    uint returnCode = 0;

    dPoints<D, Precision> points;
    if (loadPoints) {
        points = loadObject<dPoints<D, Precision>>(pointFile);
    } else {
        points = genPoints(N, bounds, dice);
    }

    TriangulateReturn ret;
    if (vm.count("seq-fault")) {
        ulong i = 0;
        std::random_device rd;
        tGenerator gen(rd());
        std::function<Precision()> dice = std::bind(distribution, gen);
        for (; ;) {
            std::cout << "." << std::flush;
            if (++i % 80 == 0)
                std::cout << std::endl;

            points = genPoints(N, bounds, dice);
            triangulate(bounds, recursionDepth, points, p, alg, verify);
        }
    } else
        ret = triangulate(bounds, recursionDepth, points, p, alg, verify);

    LOG("Triangulating "
        << points.size() << " points took "
        << std::chrono::duration_cast<std::chrono::milliseconds>(ret.time)
                .count() << " ms" << std::endl);

    returnCode = ret.valid ? EXIT_SUCCESS : EXIT_FAILURE;

#ifdef ENABLE_PROFILING
        VLOG(std::endl << "Profiling Counters:" << std::endl);
        for(const auto &c : ret.run.counters()){
            VLOG(c.first << ": " << c.second << std::endl);
        }
#endif

    if (LOGGER.getLogLevel() >= Logger::Verbosity::NORMAL) {
        LOGGER.printLog(std::cout);
    }

#ifdef ENABLE_PROFILING
        std::cout << "Profiling enabled!" << std::endl;
#endif

    return returnCode;
}

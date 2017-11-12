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

#include "LoadBalancedDCTriangulator.h"
#include "load_balancing/OldPartitionerPartitioner.h"
#include "load_balancing/SimplePartitioner.h"
#include "load_balancing/sample_partitioner/BinaryBoxEstimatingSamplePartitioner.h"

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
    std::size_t nSimplices;
    ExperimentRun run;
};

std::unique_ptr<LoadBalancing::Partitioner<D, Precision>> getPartitioner(unsigned char type, size_t threads) {
    std::unique_ptr<LoadBalancing::Partitioner<D, Precision>> result;
    switch (type) {
        case 'D':
            result = std::make_unique<LoadBalancing::SimplePartitioner<D, Precision>>();
            break;
        case 's': {
            std::random_device rand;
            result = std::make_unique<LoadBalancing::BinaryBoxEstimatingSamplePartitioner<D, Precision>>(100, rand(), 2500);
            break;
        }
        case 'd':
        case 'c':
        case 'e':
        default:
            auto oldPartitioner = Partitioner<D, Precision>::make(type);
            auto maxRecursions = oldPartitioner->getRecursionDepth(threads);
            auto baseCutoff = DCTriangulator<D, Precision>::BASE_CUTOFF;
            result = std::make_unique<LoadBalancing::OldPartitionerPartitioner<D, Precision>>(
                std::move(oldPartitioner), maxRecursions, baseCutoff);
            break;
    }
    return result;
}

TriangulateReturn triangulate(const dBox<D, Precision> &bounds,
                              const uint threads,
                              dPoints<D, Precision> &points,
                              const unsigned char splitter,
                              const unsigned char alg,
                              const bool parallelBase,
                              const bool verify) {

    std::unique_ptr<Triangulator<D, Precision>> triangulator_ptr;
    switch(alg) {
        case 'c':
            triangulator_ptr = 
            std::make_unique<PureCGALTriangulator<D, Precision, false>>(bounds, points);
            break;
            
        case 'm':
            triangulator_ptr =
            std::make_unique<PureCGALTriangulator<D, Precision, true>>(bounds, points, 100);
            break;
            
        case 'l': {
            auto partitioner = getPartitioner(splitter, threads);
            triangulator_ptr =
            std::make_unique<LoadBalancing::DCTriangulator<D, Precision>>(bounds, points, std::move(partitioner),
                                                            100, parallelBase, false);
            break;
        }
            
        case 'd': { // with same parameters as 'l'
            auto splitter_ptr = Partitioner<D, Precision>::make(splitter);
            triangulator_ptr =
            std::make_unique<DCTriangulator<D, Precision>>(bounds, points, threads,
                                                           std::move(splitter_ptr), 100, parallelBase,
                                                           false);
            break;
        }
            
        default:
            auto splitter_ptr = Partitioner<D, Precision>::make(splitter);
            triangulator_ptr =
            std::make_unique<DCTriangulator<D, Precision>>(bounds, points, threads,
                                                           std::move(splitter_ptr), 100, parallelBase);
            break;
    }

    TriangulateReturn ret;
    PROFILER.setRun(&ret.run);

    //try {

    auto t1 = Clock::now();
    auto dt = triangulator_ptr->triangulate();
    auto t2 = Clock::now();

    ret.currRSS = getCurrentRSS();
    ret.peakRSS = getPeakRSS();
    ret.nSimplices = dt.exact_size();

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

    unsigned char p = 'c';
    unsigned char dist = 'u';
    tIdType N;
    uint threads = tbb::task_scheduler_init::default_num_threads();
    unsigned char alg = 'd';
    bool parallelBase = false;

    bool verify = true;

    std::string pointFile;
    bool loadPoints = false;

    po::options_description cCommandLine("Command Line Options");
    cCommandLine.add_options()("n", po::value<tIdType>(&N), "number of points");
    cCommandLine.add_options()(
            "splitter", po::value<unsigned char>(&p),
            "splitter - _c_ycle, _d_-dimensional, _[0-d-1]_ fixed dimension");
    cCommandLine.add_options()(
            "dist", po::value<unsigned char>(&dist),
            "distribution: _u_niform, _e_llipsoid, _l_ines");
    cCommandLine.add_options()("no-verify", "don't verify triangulation");
    cCommandLine.add_options()("no-output", "don't write images");
    cCommandLine.add_options()("stats", "write run statistics even with silent");
    cCommandLine.add_options()("seq-fault", "run triangulation until seq fault occurs");
    cCommandLine.add_options()("verbosity", po::value<int>(&verbosity),
                               "verbosity");
    cCommandLine.add_options()("points", po::value<std::string>(&pointFile),
                               "load points from file");
    cCommandLine.add_options()("threads", po::value<uint>(&threads),
                               "specify number of threads");
    cCommandLine.add_options()("alg", po::value<unsigned char>(&alg),
                               "specify algorithm to use");
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
    LOGGER.setLogLevel(static_cast<Logger::Verbosity>(verbosity));

    // plausability checks
    bool valid = true;

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
        auto pg = GeneratorFactory<D, Precision>::make(dist);
        points = pg->generate(N, bounds, startGen);
    }

    TriangulateReturn ret;
    if (vm.count("seq-fault")) {
        ulong i = 0;
        std::random_device rd;
        tGenerator gen(rd());
        auto pg = GeneratorFactory<D, Precision>::make(dist);
        for (; ;) {
            std::cout << "." << std::flush;
            if (++i % 80 == 0)
                std::cout << std::endl;

            points = pg->generate(N, bounds, gen);
            triangulate(bounds, threads, points, p, alg, parallelBase, verify);
        }
    } else
        ret = triangulate(bounds, threads, points, p, alg, parallelBase, verify);

    LOG("Triangulating "
        << points.size() << " points to " << ret.nSimplices  << " simplices took "
        << std::chrono::duration_cast<std::chrono::milliseconds>(ret.time)
                .count() << " ms and " << ret.currRSS / 1e6 << "/" << ret.peakRSS / 1e6 << " MB" << std::endl);

    if(IS_SILENT && vm.count("stats")){
      std::cout << "Triangulating "
        << points.size() << " points to " << ret.nSimplices  << " simplices took "
        << std::chrono::duration_cast<std::chrono::milliseconds>(ret.time)
                .count() << " ms and " << ret.currRSS / 1e6 << "/" << ret.peakRSS / 1e6 << " MB" << std::endl;
    }

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

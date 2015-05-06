// own
#include "Geometry.h"
#include "DCTriangulator.h"
#include "Partitioner.h"

#include "utils/Random.h"
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
    tDuration time;
};

//**************************

TriangulateReturn triangulate(const dBox<D, Precision> &bounds,
                              const uint baseCase,
                              dPoints<D, Precision> &points,
                              const unsigned char splitter,
                              const unsigned char alg = 'd') {

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
                    std::make_unique<DCTriangulator<D, Precision>>(bounds, points, baseCase, splitter, 100,
                                                                   false);
        }
    }

    TriangulateReturn ret;

    //try {

    auto t1 = Clock::now();
    auto dt = triangulator_ptr->triangulate();
    auto t2 = Clock::now();

    ret.exception = false;
    ret.time = std::chrono::duration_cast<tDuration>(t2 - t1);

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

int doStudy(const uint N, const dBox<D, Precision> &bounds, const uint baseCase,
            std::function<Precision()> &&dice) {
    // define splitter test cases
    std::vector<unsigned char> splitters = {'d', 'c', '0', '1'};

    LOGGER.setLogLevel(Logger::Verbosity::SILENT);

    uint nIterations = 9 * (log10(N) - 1) + 1;
    ProgressDisplay progress(splitters.size() * nIterations, std::cout);

    // generate test points deterministically
    std::vector<dPoints<D, Precision>> pointss;
    for (uint i = 0; i < nIterations; ++i) {
        uint n = ((i % 9) + 1) * std::pow(10, std::floor(i / 9 + 1));
        pointss.push_back(genPoints(n, bounds, dice));
    }

    tbb::parallel_for(std::size_t(0), pointss.size(), [&](const uint i) {
        tbb::parallel_for(std::size_t(0), splitters.size(), [&](const uint j) {

            auto &points = pointss[i];
            unsigned char p = splitters[j];

            triangulate(bounds, baseCase, points, p);

            ++progress;
        });
    });

    std::cout << LOGGER << std::endl;

    LOG("Finished" << std::endl);

    return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {

    // define program options

    namespace po = boost::program_options;

    int verbosity = -1;

    unsigned char p;
    uint N;
    uint baseCase;
    uint threads = tbb::task_scheduler_init::default_num_threads();
    unsigned char alg = 'd';

    bool study = false;

    std::string pointFile;
    bool loadPoints = false;

    po::options_description cCommandLine("Command Line Options");
    cCommandLine.add_options()("n", po::value<uint>(&N), "number of points");
    cCommandLine.add_options()("basecase", po::value<uint>(&baseCase),
                               "threshold for base case");
    cCommandLine.add_options()(
            "splitter", po::value<unsigned char>(&p),
            "splitter - _c_ycle, _d_-dimensional, _[0-d-1]_ fixed dimension");
    cCommandLine.add_options()("study", "study mode");
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
    if ((!vm.count("basecase"))) {
        std::cout << "Please specify basecase threshold" << std::endl;
        valid = false;
    }

    if ((!vm.count("study") && !vm.count("splitter"))) {
        std::cout << "Please specify splitter" << std::endl;
        valid = false;
    }

    if ((!vm.count("study") && !(vm.count("n") || vm.count("points")))) {
        std::cout << "Please specify number of points or point file" << std::endl;
        valid = false;
    }

    if (!valid)
        return EXIT_FAILURE;

    if (vm.count("points")) {
        loadPoints = true; // filename will be in pointFile
    }

    if (vm.count("study")) {
        study = true;
    }

    if (vm.count("no-verify")) {
        DCTriangulator<D, Precision>::VERIFY = false;
    }

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
    if (study) {
        returnCode = doStudy(N, bounds, baseCase, std::move(dice));
    } else {

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
                triangulate(bounds, baseCase, points, p, alg);
            }
        } else
            ret = triangulate(bounds, baseCase, points, p, alg);

        LOG("Triangulating "
            << points.size() << " points took "
            << std::chrono::duration_cast<std::chrono::milliseconds>(ret.time)
                    .count() << " ms" << std::endl);

        returnCode = EXIT_SUCCESS;
    }

    if (LOGGER.getLogLevel() >= Logger::Verbosity::NORMAL) {
        LOGGER.printLog(std::cout);
    }

    return returnCode;
}

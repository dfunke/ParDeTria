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
#include "utils/DBConnection.h"
#include "utils/Serialization.hxx"

#include <boost/program_options.hpp>
#include <tbb/task_scheduler_init.h>

//**************************

#define D 3
#define Precision double

std::vector<unsigned char> splitters = {'c'};
std::vector<unsigned char> triangulators = {'c', 'm', 'd'};
std::vector<unsigned char> distributions = {'u', 'e', 'l'};
std::vector<uint> occupancies = {10, 50, 100, 1000};

enum class TriState : char { INDEF = -1, FALSE = 0, TRUE = 1};

DBConnection db("db_" + getHostname() + ".dat",
#ifndef ENABLE_PROFILING
                "benchmarks"
#else
                "profiling"
#endif
);

dBox<D, Precision> bounds(dVector<D, Precision>( {{ 0, 0, 0 }}

), dVector<D, Precision>({{ 100,100,100}}));

void runExperiment(ExperimentRun &run, const uint reps = 10) {

    //quite all output and unnecessary computations
    LOGGER.setLogLevel(Logger::Verbosity::SILENT);

    //get point set
    unsigned char dist = run.getTrait<unsigned char>("dist");
    uint nPoints = run.getTrait<uint>("nP");

    //use same start seed for all experiment runs
    tGenerator gen(START_SEED);
    auto pg = GeneratorFactory<D, Precision>::make(dist);
    auto points = pg->generate(nPoints, bounds, gen);

    std::unique_ptr<Triangulator<D, Precision>> triangulator_ptr;
    unsigned char alg = run.getTrait<unsigned char>("alg");

    uint threads = run.hasTrait("threads") ? run.getTrait<uint>("threads") : 1;
    // load scheduler with specified number of threads
    tbb::task_scheduler_init init(threads);

    if (alg == 'c') {
        triangulator_ptr =
                std::make_unique<PureCGALTriangulator<D, Precision, false>>(bounds, points);
    } else {
        uint gridOccupancy = run.getTrait<uint>("occupancy");

        if (alg == 'm') {
            triangulator_ptr =
                    std::make_unique<PureCGALTriangulator<D, Precision, true>>(bounds, points, gridOccupancy);
        } else {
            uint recursionDepth = run.getTrait<uint>("recursion-depth");
            unsigned char splitter = run.getTrait<unsigned char>("splitter");
            bool parBase = run.getTrait<bool>("parallel-base");
            bool parEdge = run.getTrait<bool>("parallel-edge");

            triangulator_ptr =
                    std::make_unique<DCTriangulator<D, Precision>>(bounds, points, recursionDepth, splitter, gridOccupancy,
                                                                   parBase, parEdge);
        }
    }

    run.addTrait("start-time", getDatetime());

    PROFILER.setRun(&run);

    try {

        for (uint i = 0; i < reps; ++i) {
            auto t1 = Clock::now();
            auto dt = triangulator_ptr->triangulate();
            auto t2 = Clock::now();

            run.addMeasurement("memory", getCurrentRSS());
            run.addMeasurement("peakMem", getPeakRSS());
            run.addMeasurement("times", std::chrono::duration_cast<tDuration>(t2 - t1).count());
            
        }

    } catch (std::exception &e) {
        std::cerr << "Experiment: " << run.str() << std::endl;
        std::cerr << "\tException raised: " << e.what() << std::endl;

        // output points
        std::ofstream errFile("failedExperiments", std::ios::out | std::ios::app);
        errFile << run.str() << std::endl;
    }

    run.addTrait("end-time", getDatetime());

    db.save(run);
}

//**************************

void runExperiments(std::vector<ExperimentRun> &runs, const uint reps = 10, bool reverse = false) {

    signed long start = !reverse ? 0 : runs.size() - 1;
    signed long end = !reverse ? runs.size() : -1;

    for (auto i = start; i != end; i += (!reverse ? 1 : -1)) {
        std::cout << i+1 << "/" << runs.size() << ": " << runs[i].str() << std::endl;

        runExperiment(runs[i], reps);

        std::cout << "\tAverage time: "
        << runs[i].avgMeasurement("times") / 1e6
        << " ms\tAverage mem: "
        << runs[i].avgMeasurement("memory") / 1e6 << " MB" << std::endl;
    }

}

//**************************

std::vector<ExperimentRun> generateExperimentRuns(const uint maxN, const uint minN = 10,
                                                  int maxThreads = -1, int minThreads = 1,
                                                  int minRecDepth = 0, uint maxRecDepth = std::numeric_limits<uint>::max(),
                                                  bool parallel_base = false,
                                                  TriState parallel_edge = TriState::INDEF,
                                                  int runNumber = -1) {

    std::vector<ExperimentRun> runs;

    //determine maximum number of threads
    if(maxThreads == -1)
        maxThreads = tbb::task_scheduler_init::default_num_threads();
    if(minThreads == -1)
        minThreads = maxThreads;

    //determine the latest run number
    if(runNumber == -1)
        runNumber = db.getMaximum<uint>("run-number") + 1;

    bool optimalRecDepth = (minRecDepth == 0 && maxRecDepth == std::numeric_limits<uint>::max());

    //loop over distributions
    for (const unsigned char dist : distributions) {

        //create ExperimentRuns
        //loop over number of points
        for (uint nPoints = minN; nPoints <= maxN; nPoints += pow(10, floor(log10(nPoints)))) {

            //loop over triangulators
            for (const unsigned char alg : triangulators) {

                if (alg == 'c') {

                    ExperimentRun run;
                    run.addTrait("run-number", runNumber);
                    run.addTrait("dist", dist);
                    run.addTrait("nP", nPoints);
                    run.addTrait("alg", alg);
                    run.addTrait("threads", 1);

                    runs.emplace_back(run);

                } else {

                    //loop over number of threads, if algo is multi-threaded
                    for (uint threads = minThreads; threads <= (unsigned) maxThreads; threads <<= 1) {

                        //loop over gridOccupancy
                        bool firstOccupancy = true;
                        for (const uint occ : occupancies) {

                            if (alg == 'm') {
                                ExperimentRun run;
                                run.addTrait("run-number", runNumber);
                                run.addTrait("dist", dist);
                                run.addTrait("nP", nPoints);
                                run.addTrait("alg", alg);
                                run.addTrait("threads", threads);
                                run.addTrait("occupancy", occ);

                                runs.emplace_back(run);

                            } else {

                                //loop over splitters
                                for (const unsigned char splitter : splitters) {

                                    if(optimalRecDepth){
                                        minRecDepth = maxRecDepth = log2(threads);
                                    }

                                    //loop over recursion depth
                                    for (uint recursionDepth = minRecDepth;
                                         recursionDepth <= std::min(maxRecDepth, (uint) std::ceil(std::log2(nPoints / DCTriangulator<D, Precision>::BASE_CUTOFF)));
                                         ++recursionDepth) {

                                        if (firstOccupancy) { // we do not need to generate runs for the following occupancies
                                            ExperimentRun runSeq;
                                            runSeq.addTrait("run-number", runNumber);
                                            runSeq.addTrait("dist", dist);
                                            runSeq.addTrait("nP", nPoints);
                                            runSeq.addTrait("alg", alg);
                                            runSeq.addTrait("threads", threads);
                                            runSeq.addTrait("occupancy", 1);
                                            runSeq.addTrait("splitter", splitter);
                                            runSeq.addTrait("recursion-depth", recursionDepth);
                                            runSeq.addTrait("parallel-base", false);

                                            if(parallel_edge == TriState::INDEF){
                                                runSeq.addTrait("parallel-edge", false);
                                            runs.emplace_back(runSeq);

                                                runSeq.updateTrait("parallel-edge", true);
                                                runs.emplace_back(runSeq);
                                            } else {
                                                runSeq.addTrait("parallel-edge", parallel_edge == TriState::TRUE ? true : false);
                                                runs.emplace_back(runSeq);
                                            }
                                        }

                                        if (parallel_base) {
                                            ExperimentRun runPar;
                                            runPar.addTrait("run-number", runNumber);
                                            runPar.addTrait("dist", dist);
                                            runPar.addTrait("nP", nPoints);
                                            runPar.addTrait("alg", alg);
                                            runPar.addTrait("threads", threads);
                                            runPar.addTrait("occupancy", occ);
                                            runPar.addTrait("splitter", splitter);
                                            runPar.addTrait("recursion-depth", recursionDepth);
                                            runPar.addTrait("parallel-base", true);

                                            if(parallel_edge == TriState::INDEF){
                                                runPar.addTrait("parallel-edge", false);
                                            runs.emplace_back(runPar);

                                                runPar.updateTrait("parallel-edge", true);
                                                runs.emplace_back(runPar);
                                            } else {
                                                runPar.addTrait("parallel-edge", parallel_edge == TriState::TRUE ? true : false);
                                                runs.emplace_back(runPar);
                                            }
                                        }

                                    }
                                }
                            }

                            firstOccupancy = false;
                        }
                    }
                }
            }
        }
    }

    return runs;
}

int main(int argc, char *argv[]) {

    // define program options

    namespace po = boost::program_options;

    uint maxN, minN = 10;
    int maxThreads = -1;
    int minThreads = 1;
    uint minRecDepth = 0;
    uint maxRecDepth = std::numeric_limits<uint>::max();
    uint occupancy = 1;
    unsigned char alg;
    unsigned char dist;
    bool parallelBase;
    TriState parallelEdge = TriState::TRUE;

    uint reps = 10;
    std::string runFile;
    std::string run;
    int runNumber = -1;

    po::options_description cCommandLine("Command Line Options");
    // point options
    cCommandLine.add_options()("maxN", po::value<uint>(&maxN),
                               "maximum number of points");
    cCommandLine.add_options()("minN", po::value<uint>(&minN),
                               "minimum number of points, default 10");

    // thread options
    cCommandLine.add_options()("minThreads", po::value<int>(&minThreads),
                               "minimum number of threads, default 1");
    cCommandLine.add_options()("maxThreads", po::value<int>(&maxThreads),
                               "maximum number of threads, default -1 = automatic");

    // rec depth options
    // thread options
    cCommandLine.add_options()("minRecDepth", po::value<uint>(&minRecDepth),
                               "minimum recursion depth, default 0");
    cCommandLine.add_options()("maxRecDepth", po::value<uint>(&maxRecDepth),
                               "maximum recursion depth, default: expected size of min base case");

    // algorithm options
    cCommandLine.add_options()("algorithm", po::value<unsigned char>(&alg),
                               "algorithm to use");

    // distribution options
    cCommandLine.add_options()("dist", po::value<unsigned char>(&dist),
                               "point distribution");

    // occupancy options
    cCommandLine.add_options()("occupancy", po::value<uint>(&occupancy),
                               "specify occupancy of grid lock data structure");

    // parallel base options
    cCommandLine.add_options()("no-parallel-base", "don't use parallel base solver");

    // parallel edge
    cCommandLine.add_options()("no-parallel-edge", "don't use parallel edge solver");
    cCommandLine.add_options()("parallel-edge-study", "use both parallel and sequential edge solver");

    // operative options
    cCommandLine.add_options()("reps", po::value<uint>(&reps),
                               "repetitions of experiments");
    cCommandLine.add_options()("runs", po::value<std::string>(&runFile),
                               "file containing experiments to run");
    cCommandLine.add_options()("run-string", po::value<std::string>(&run),
                               "string describing an experiment to run");
    cCommandLine.add_options()("run-number", po::value<int>(&runNumber),
                               "specify run-number, default -1 = automatic");
    cCommandLine.add_options()("gen-only", "just generate test-cases");
    cCommandLine.add_options()("reverse", "reverse test-case execution");

#ifdef ENABLE_PROFILING
    cCommandLine.add_options()("profiling",
                               "perform fine-grained profiling");
#endif

    cCommandLine.add_options()("help", "produce help message");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cCommandLine), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << cCommandLine << std::endl;
        return EXIT_SUCCESS;
    }

    //***************************************************************************

    std::vector<ExperimentRun> runs;
    if (vm.count("runs")) {
        std::ifstream input(runFile);

        for (std::string line; getline(input, line);) {
            runs.emplace_back(line);
        }
    } else if (vm.count("run-string")) {
        ExperimentRun exRun(run);
        if(runNumber != -1){
            // a run-number was specified on the command line
            exRun.addTrait("run-number", runNumber);
        }

        runs.push_back(std::move(exRun));
    } else {
        if ((!vm.count("maxN"))) {
            std::cout << "Please specify number of points" << std::endl;
            return EXIT_FAILURE;
        }

        if (vm.count("occupancy")) {
            occupancies = {occupancy};
        }

        if (vm.count("algorithm")) {
            triangulators = {alg};
        }

        if (vm.count("dist")) {
            distributions = {dist};
        }

        parallelBase = !vm.count("no-parallel-base");

        if(vm.count("no-parallel-edge"))
            parallelEdge = TriState::FALSE;

        if(vm.count("parallel-edge-study"))
            parallelEdge = TriState::INDEF;

        if (vm.count("profiling")) {
            // only our algorithm is instrumented
            triangulators = {'d'};
            occupancies = {100};
            parallelBase = false;

            //operations don't change with number of threads
            minThreads = -1;
            maxThreads = -1;
        }

        runs = generateExperimentRuns(maxN, minN, maxThreads, minThreads, minRecDepth, maxRecDepth, parallelBase, parallelEdge, runNumber);
    }

#ifdef ENABLE_PROFILING
    if (!vm.count("profiling")) {
        std::cout << "Compiled with profiling, but not specified on command line. Exiting" << std::endl;
        return EXIT_FAILURE;
    } else {
        // only 1 reptition, profiling counters don't change
        reps = 1;
    }
#endif

    if (vm.count("gen-only")) {
        for (const auto &r : runs)
            std::cout << r.str() << std::endl;
    } else {
        runExperiments(runs, reps, vm.count("reverse"));
    }

    return EXIT_SUCCESS;
}

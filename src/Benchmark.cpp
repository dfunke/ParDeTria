
// own
#include "Geometry.h"
#include "DCTriangulator.h"
#include "Partitioner.h"
#include "Painter.h"

#include "utils/Random.h"
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
#define REPS 10

std::vector<unsigned char> splitters = { 'c' };
std::vector<unsigned char> triangulators = { 'c', 'm', 'd' };
std::vector<unsigned char> distributions = { 'u' };
std::vector<uint> occupancies = { 10, 50, 100, 1000 };

DBConnection db("db_" + getHostname() + ".dat", "benchmarks");
RandomPointHolder<D, Precision> randomPoints;
dBox<D, Precision> bounds(dVector<D, Precision>({{0,0,0}}), dVector<D, Precision>({{100,100,100}}));

void runExperiment(ExperimentRun & run) {

  //quite all output and unnecessary computations
  LOGGER.setLogLevel(Logger::Verbosity::SILENT);
  DCTriangulator<D, Precision>::VERIFY = false;
  Painter<D, Precision>::ENABLED = false;

  //get point set
  unsigned char dist = run.getTrait<unsigned char>("dist");
  uint nPoints = run.getTrait<uint>("nP");
  auto points = randomPoints.get(dist, nPoints);

  std::unique_ptr<Triangulator<D, Precision>> triangulator_ptr;
  unsigned char alg = run.traits().at("alg")[0];

  uint threads = 1;

  if(alg == 'c'){
      triangulator_ptr =
              std::make_unique<CGALTriangulator<D, Precision, false>>(bounds, points);
  } else {
    uint gridOccupancy = run.getTrait<uint>("occupancy");
    threads = run.getTrait<uint>("threads");

    if(alg == 'm'){
      triangulator_ptr =
              std::make_unique<CGALTriangulator<D, Precision, true>>(bounds, points, gridOccupancy);
    } else {
      uint basecase = run.getTrait<uint>("basecase");
      unsigned char splitter = run.getTrait<unsigned char>("splitter");
      bool parBase = run.getTrait<bool>("parallel-base");

      triangulator_ptr =
              std::make_unique<DCTriangulator<D, Precision>>(bounds, points, basecase, splitter, gridOccupancy, parBase);
    }
  }

  // load scheduler with specified number of threads
  tbb::task_scheduler_init init(threads);

  try {

    for(uint i = 0; i < REPS; ++i) {
      auto t1 = Clock::now();
      auto dt = triangulator_ptr->triangulate();
      auto t2 = Clock::now();

      run.addMemory(getCurrentRSS());
      run.addTime(std::chrono::duration_cast<tDuration>(t2 - t1));
    }

  } catch (AssertionException &e) {
    std::cerr << "Experiment: " << run.str() << std::endl;
    std::cerr << "\tAssertion failed: " << e.what() << std::endl;

    // output points
    storeObject(points, "assertionFailed_" + run.str("_") + ".dat");
  }

  db.save(run);
}

//**************************

void runExperiments(std::vector<ExperimentRun> & runs){

  for(uint i = 0; i < runs.size(); ++i){
    std::cout << i << "/" << runs.size() << ": " << runs[i].str() << std::endl;

    runExperiment(runs[i]);

    std::cout << "\tAverage time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(runs[i].avgTime()).count()
              << " ms" << std::endl;
  }

}

//**************************

std::vector<ExperimentRun> generateExperimentRuns(const uint N) {

  std::vector<ExperimentRun> runs;

  //determine maximum number of threads
  uint maxThreads = tbb::task_scheduler_init::default_num_threads();

  //loop over distributions
  for(const unsigned char dist : distributions){

    //fill random point holder
    randomPoints.fill(dist, N, bounds);

    //create ExperimentRuns
    //loop over number of points
    uint nPointsIterations = 9 * (log10(N) - 1) + 1;
    for (uint i = 0; i < nPointsIterations; ++i) {
      uint nPoints = ((i % 9) + 1) * std::pow(10, std::floor(i / 9 + 1));

      //loop over triangulators
      for(const unsigned char alg : triangulators){

        if(alg == 'c'){

          ExperimentRun run;
          run.addTrait("dist", dist);
          run.addTrait("nP", nPoints);
          run.addTrait("alg", alg);

          runs.emplace_back(run);

        } else {

          //loop over number of threads, if algo is multi-threaded
          for (uint threads = 1; threads <= maxThreads; threads <<= 1) {

            //loop over gridOccupancy
            bool firstOccupancy = true;
            for(const uint occ : occupancies){

              if(alg == 'm'){
                ExperimentRun run;
                run.addTrait("dist", dist);
                run.addTrait("nP", nPoints);
                run.addTrait("alg", alg);
                run.addTrait("threads", threads);
                run.addTrait("occupancy", occ);

                runs.emplace_back(run);

              } else {

                //loop over splitters
                for(const unsigned char splitter : splitters){

                  //loop over base cases
                  for(uint basecase = 500; basecase <= std::max(500u, nPoints / 2); basecase <<= 2){

                    if(firstOccupancy) { // we do not need to generate runs for the following occupancies
                      ExperimentRun runSeq;
                      runSeq.addTrait("dist", dist);
                      runSeq.addTrait("nP", nPoints);
                      runSeq.addTrait("alg", alg);
                      runSeq.addTrait("threads", threads);
                      runSeq.addTrait("occupancy", 1);
                      runSeq.addTrait("splitter", splitter);
                      runSeq.addTrait("basecase", basecase);
                      runSeq.addTrait("parallel-base", false);

                      runs.emplace_back(runSeq);
                    }

                    ExperimentRun runPar;
                    runPar.addTrait("dist", dist);
                    runPar.addTrait("nP", nPoints);
                    runPar.addTrait("alg", alg);
                    runPar.addTrait("threads", threads);
                    runPar.addTrait("occupancy", occ);
                    runPar.addTrait("splitter", splitter);
                    runPar.addTrait("basecase", basecase);
                    runPar.addTrait("parallel-base", true);

                    runs.emplace_back(runPar);

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

  uint N;

  po::options_description cCommandLine("Command Line Options");
  cCommandLine.add_options()("n", po::value<uint>(&N),
                             "maximum number of points");
  cCommandLine.add_options()("help", "produce help message");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cCommandLine), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << cCommandLine << std::endl;
    return EXIT_SUCCESS;
  }

  // plausability checks
  if ((!vm.count("n"))) {
    std::cout << "Please specify number of points" << std::endl;
    return EXIT_FAILURE;
  }

  //***************************************************************************

  auto runs = generateExperimentRuns(N);
  runExperiments(runs);

  return EXIT_SUCCESS;
}

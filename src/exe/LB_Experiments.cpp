#include <iostream>
#include <unordered_map>
#include <memory>
#include <string>
#include <vector>

#include <tbb/task_scheduler_init.h>

#include "load_balancing/CommandLineInterface.h"
#include "load_balancing/Experiment.h"
#include "load_balancing/monitors/ComprehensiveMonitor.h"


template <uint D, typename Precision>
void execute(const po::variables_map& vm, uint threads, const std::string& argString) {
    dBox<D, Precision> bounds;
    bounds.low.fill(0);
    bounds.high.fill(100);

    dPoints<D, Precision> points;
    std::string distribution;
    if (vm.count("points")) {
        std::string filename = vm["points"].as<std::string>();
        bounds = points.loadFromFile(filename);
        distribution = filename;
    } else {
        auto pg = createGenerator<D, Precision>(vm);
        auto N = vm["n"].as<tIdType>();
        points = pg->generate(N, bounds, startGen);
        distribution = vm["distribution"].as<std::string>();
    }
    
    auto partitioner = createPartitioner<D, Precision, lb::ComprehensiveMonitor>(vm,
																				 threads,
																				 startGen);
    assert(partitioner);
    
    typename lb::Experiment<D, Precision>::Setup setup {
        argString,
        bounds,
        points,
        vm.count("validate") > 0,
        distribution,
        threads,
        vm.count("seed") ? vm["seed"].as<uint>() : START_SEED
    };
        
    lb::Experiment<D, Precision> exp(std::move(partitioner), std::move(setup), std::cout);
    bool valid = exp.runOnce();
    if(!valid) {
	    LOGGER.printLog(std::cerr);
	}
}

int main(int argc, char *argv[]) {
    uint threads = tbb::task_scheduler_init::default_num_threads();
    std::string pointFile;
    std::string partitionerName;
    uint dim;
    int verbosity = (short)Logger::Verbosity::NORMAL;

    auto cCommandLine = defaultOptions<double>();
    cCommandLine.add_options()("dim", po::value(&dim), "dimension");
    cCommandLine.add_options()("threads", po::value(&threads), "specify number of threads");
    cCommandLine.add_options()("points", po::value<std::string>(&pointFile),
                               "load points from file");
    cCommandLine.add_options()("verbosity", po::value<int>(&verbosity),
                               "verbosity");
    LOGGER.setLogLevel(static_cast<Logger::Verbosity>(verbosity));

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cCommandLine), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << cCommandLine << std::endl;
        return EXIT_SUCCESS;
    }

	if(vm.count("seed")) {
		startGen.seed(vm["seed"].as<uint>());
	}
    
    if (!(vm.count("n") || vm.count("points"))) {
        std::cout << "Please specify number of points or point file" << std::endl;
        return EXIT_FAILURE;
    }
    
    if(!vm.count("partitioner")) {
        std::cout << "Please specify a partitioner" << std::endl;
        return EXIT_FAILURE;
    }
    
    if(!(vm.count("distribution") || vm.count("points"))) {
        std::cout << "Please specify a distribution" << std::endl;
        return EXIT_FAILURE;
    }

    tbb::task_scheduler_init init(threads);
    
    std::string argString = std::accumulate(argv + 2, argv + argc, std::string(argc > 1 ? argv[1] : ""),
                                        [](std::string acc, const std::string& element) {
                                            return acc += " " + element;
                                        });

    switch(dim) {
    case 2:
        execute<2, double>(vm, threads, argString);
        break;
    case 3:
        execute<3, double>(vm, threads, argString);
        break;
    default:
        std::cerr << "unsupported dimension\n";
    }
    
    return EXIT_SUCCESS;
}


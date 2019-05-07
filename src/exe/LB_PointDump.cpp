#include <cstdio>
#include "load_balancing/CommandLineInterface.h"

constexpr uint D = 3;
using Precision = double;

int main(int argc, char** argv) {
    auto cCommandLine = defaultOptions<double>();
    LOGGER.setLogLevel(static_cast<Logger::Verbosity>(Logger::Verbosity::SILENT));

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
    
    if (vm.count("n") == 0) {
        std::cout << "Please specify number of points or point file" << std::endl;
        return EXIT_FAILURE;
    }
    
    if(vm.count("distribution") == 0) {
        std::cout << "Please specify a distribution" << std::endl;
        return EXIT_FAILURE;
    }
    

    dBox<D, Precision> bounds;
    bounds.low.fill(0);
    bounds.high.fill(100);
    
    dPoints<D, Precision> points;
    auto pg = createGenerator<D, Precision>(vm);
    auto N = vm["n"].as<tIdType>();
    points = pg->generate(N, bounds, startGen);
    
    for(const auto& point : points) {
        printf("%f %f %f\n", point.coords[0], point.coords[1], point.coords[2]);
    }
}

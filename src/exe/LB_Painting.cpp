#include <string>
#include <tbb/task_scheduler_init.h>
#include "Painter.h"
#include "utils/MakeIds.h"
#include "load_balancing/CommandLineInterface.h"

using Precision = double;
constexpr uint D = 2;


int main(int argc, char *argv[]) {
    uint threads = tbb::task_scheduler_init::default_num_threads();
    //std::string pointFile;    
    std::string partitionerName;

    po::options_description cCommandLine("Command Line Options");
    cCommandLine.add_options()("n", po::value<tIdType>(), "number of points");
    cCommandLine.add_options()("partitioner", po::value<std::string>());
    cCommandLine.add_options()("distribution", po::value<std::string>());
    //cCommandLine.add_options()("points", po::value(&pointFile), "load points from file");
    cCommandLine.add_options()("threads", po::value(&threads), "specify number of threads");
    //cCommandLine.add_options()("sample-size", po::value<uint>());
    cCommandLine.add_options()("split-dimension", po::value<uint>());
    cCommandLine.add_options()("out", po::value<std::string>());
    cCommandLine.add_options()("help", "produce help message");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cCommandLine), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << cCommandLine << std::endl;
        return EXIT_SUCCESS;
    }

    LOGGER.setLogLevel(static_cast<Logger::Verbosity>(0));
    
    if (!(vm.count("n") || vm.count("points"))) {
        std::cout << "Please specify number of points or point file" << std::endl;
        return EXIT_FAILURE;
    }
    
    if(!vm.count("partitioner")) {
        std::cout << "Please specify a partitioner" << std::endl;
        return EXIT_FAILURE;
    }
    
    if(!vm.count("distribution")) {
        std::cout << "Please specify a distribution" << std::endl;
        return EXIT_FAILURE;
    }
    
    auto filename = vm["out"].as<std::string>();
    
    dBox<D, Precision> bounds;
    bounds.low.fill(0);
    bounds.high.fill(100);

    dPoints<D, Precision> points;
    /*if (vm.count("points")) {
        points = loadObject<dPoints<D, Precision>>(pointFile);
    } else {*/
        auto pg = createGenerator<D, Precision>(vm["distribution"].as<std::string>());
        auto N = vm["n"].as<tIdType>();
        points = pg->generate(N, bounds, startGen);
    //}
    
    std::random_device rand;
    auto partitioner = createPartitioner<D, Precision>(vm, threads, rand);
    assert(partitioner);

    tbb::task_scheduler_init init(threads);
    
    auto pointIds = makePointIds(points);
    auto partitioning = partitioner->partition(bounds, points, pointIds);
    
    Painter<2, double> painter(bounds);
    for(const auto& point : points) {
        painter.draw(point);
    }
    painter.save(filename + "_points");
}

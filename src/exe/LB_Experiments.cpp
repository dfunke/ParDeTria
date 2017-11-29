#include <iostream>
#include <unordered_map>
#include <memory>
#include <string>
#include <vector>

#include <tbb/task_scheduler_init.h>

//#include "Geometry.h"

#include "utils/Generator.h"
#include "utils/MakeIds.h"
//#include "utils/Serialization.hxx"

#include "LoadBalancedDCTriangulator.h"
#include "load_balancing/Partitioner.h"
#include "load_balancing/OldPartitionerPartitioner.h"
#include "load_balancing/SimplePartitioner.h"
#include "load_balancing/sample_partitioner/BinaryBoxEstimatingSamplePartitioner.h"
#include "load_balancing/sample_partitioner/CenterDistancePointAssigningSamplePartitioner.h"
#include "load_balancing/sample_partitioner/BoundsDistancePointAssigningSamplePartitioner.h"
#include "load_balancing/Experiment.h"
#include <boost/program_options.hpp>

constexpr auto D = 3;
using Precision = double;

namespace po = boost::program_options;
namespace lb = LoadBalancing;

template <uint D, typename Precision>
std::unique_ptr<PointGenerator<D, Precision>> createGenerator(const std::string& distributionName) {
    std::unique_ptr<PointGenerator<D, Precision>> result = nullptr;
    if("ellipsoid" == distributionName)
        result = std::make_unique<EllipsoidPointGenerator<D, Precision>>();
    else if("skew-line" == distributionName)
        result = std::make_unique<SkewLinePointGenerator<D, Precision>>();
    else if("normal" == distributionName)
        result = std::make_unique<NormalPointGenerator<D, Precision>>();
    else if("bubble" == distributionName)
        result = std::make_unique<BubblePointGenerator<D, Precision>>();
    else if("uniform" == distributionName)
        result = std::make_unique<UniformPointGenerator<D, Precision>>();
    
    return result;
}

int main(int argc, char *argv[]) {
    std::string distribution;
    tIdType N;
    uint threads = tbb::task_scheduler_init::default_num_threads();
    std::string pointFile;    
    std::string partitionerName;
    uint dim;

    po::options_description cCommandLine("Command Line Options");
    cCommandLine.add_options()("n", po::value(&N), "number of points");
    cCommandLine.add_options()("dim", po::value(&dim), "dimension");
    cCommandLine.add_options()("partitioner", po::value(&partitionerName));
    cCommandLine.add_options()("distribution", po::value<std::string>());
    cCommandLine.add_options()("points", po::value(&pointFile), "load points from file");
    cCommandLine.add_options()("threads", po::value(&threads), "specify number of threads");
    cCommandLine.add_options()("sample-size", po::value<uint>());
    cCommandLine.add_options()("split-dimension", po::value<uint>());
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

    tbb::task_scheduler_init init(threads);

    dBox<D, Precision> bounds;
    bounds.low.fill(0);
    bounds.high.fill(100);

    dPoints<D, Precision> points;
    /*if (vm.count("points")) {
        points = loadObject<dPoints<D, Precision>>(pointFile);
    } else {*/
        auto pg = createGenerator<D, Precision>(vm["distribution"].as<std::string>());
        points = pg->generate(N, bounds, startGen);
    //}
    
    std::random_device rand;
    
    std::unique_ptr<lb::Partitioner<D, Precision>> partitioner = nullptr;
    if("binary-besp" == partitionerName){
        uint sampleSize = vm["sample-size"].as<uint>();
        auto baseCutoff = lb::DCTriangulator<D, Precision>::BASE_CUTOFF;
        partitioner = std::make_unique<lb::BinaryBoxEstimatingSamplePartitioner<D, Precision>>(sampleSize, rand(), baseCutoff);
    } else if("center-distance-pasp" == partitionerName){
        uint sampleSize = vm["sample-size"].as<uint>();
        partitioner = std::make_unique<lb::CenterDistancePointAssigningSamplePartitioner<D, Precision>>(sampleSize, rand(), threads);
    } else if("bounds-distance-pasp" == partitionerName){
        uint sampleSize = vm["sample-size"].as<uint>();
        partitioner = std::make_unique<lb::BoundsDistancePointAssigningSamplePartitioner<D, Precision>>(sampleSize, rand(), threads);
    } else {
        std::unique_ptr<Partitioner<D, Precision>> oldPartitioner = nullptr;
        if("dWay" == partitionerName) {
            oldPartitioner = std::make_unique<dWayPartitioner<D, Precision>>();
        } else if("cycle" == partitionerName) {
            oldPartitioner = std::make_unique<CyclePartitioner<D, Precision>>();
        } else if("extend" == partitionerName) {
            oldPartitioner = std::make_unique<ExtendPartitioner<D, Precision>>();
        } else if("one-dim" == partitionerName) {
            uint d = vm["split-dimension"].as<uint>();
            ASSERT(d < D);
            oldPartitioner = std::make_unique<OneDimPartitioner<D, Precision>>(d);
        }
        auto maxRecursions = oldPartitioner->getRecursionDepth(threads);
        auto baseCutoff = lb::DCTriangulator<D, Precision>::BASE_CUTOFF;
        partitioner = std::make_unique<LoadBalancing::OldPartitionerPartitioner<D, Precision>>(
            std::move(oldPartitioner), maxRecursions, baseCutoff);
    }
    assert(partitioner);
    
    lb::Experiment<D, Precision>::Setup setup {
        bounds,
        points,
        true
    };
        
    lb::Experiment<D, Precision> exp(std::move(partitioner), std::move(setup), std::cout);
    exp.runOnce();

    return EXIT_SUCCESS;
}

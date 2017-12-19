#pragma once

#include "utils/Generator.h"
#include "load_balancing/sample_partitioner/Sampler.h"
#include "load_balancing/UnboundBubbleGenerator.h"
#include "load_balancing/PatchBubbleGenerator.h"

#include "LoadBalancedDCTriangulator.h"
#include "load_balancing/Partitioner.h"
#include "load_balancing/OldPartitionerPartitioner.h"
#include "load_balancing/SimplePartitioner.h"
#include "load_balancing/sample_partitioner/BinaryBoxEstimatingSamplePartitioner.h"
#include "load_balancing/sample_partitioner/CenterDistancePointAssigningSamplePartitioner.h"
#include "load_balancing/sample_partitioner/BoundsDistancePointAssigningSamplePartitioner.h"
#include "load_balancing/sample_partitioner/NearestSamplePointAssigningSamplePartitioner.h"

#include <boost/program_options.hpp>

namespace po = boost::program_options;
namespace lb = LoadBalancing;

template <typename T>
T extract(const std::string& name, const po::variables_map& vm, T defaultValue) {
    T result;
    if(vm.count(name))
        result = vm[name].as<T>();
    else
        result = defaultValue;
    return result;
}

template <uint D, typename Precision>
std::unique_ptr<PointGenerator<D, Precision>> createGenerator(const po::variables_map& vm) {
    std::string distributionName = vm["distribution"].as<std::string>();
    std::unique_ptr<PointGenerator<D, Precision>> result = nullptr;
    if("ellipsoid" == distributionName)
        result = std::make_unique<EllipsoidPointGenerator<D, Precision>>();
    else if("skew-line" == distributionName)
        result = std::make_unique<SkewLinePointGenerator<D, Precision>>();
    else if("normal" == distributionName)
        result = std::make_unique<NormalPointGenerator<D, Precision>>();
    else if("uniform" == distributionName)
        result = std::make_unique<UniformPointGenerator<D, Precision>>();
    else if("bubble" == distributionName){
        result = std::make_unique<BubblePointGenerator<D, Precision>>();
    } else if("unbound-bubble" == distributionName){
        uint numBubbles = extract<uint>("num-bubbles", vm, std::pow(4, D));
        Precision radius = extract<Precision>("bubble-radius", vm, 10);
        result = std::make_unique<lb::UnboundBubbleGenerator<D, Precision>>(numBubbles, radius);
    } else if("patch-bubble" == distributionName){
        uint numBubbles = extract<uint>("num-bubbles", vm, std::pow(4, D));
        result = std::make_unique<lb::PatchBubbleGenerator<D, Precision>>(numBubbles);
    }
    
    return result;
}

template <uint D, typename Precision, typename RandomGenerator>
std::unique_ptr<lb::Partitioner<D, Precision>> createPartitioner(const po::variables_map& vm, uint threads, RandomGenerator& rand) {
    std::string partitionerName = vm["partitioner"].as<std::string>();
    lb::Sampler<D, Precision> sampler(rand(), [](size_t n) -> size_t { return std::sqrt(n); });
    std::unique_ptr<lb::Partitioner<D, Precision>> partitioner = nullptr;
    if("binary-besp" == partitionerName){
        auto baseCutoff = lb::DCTriangulator<D, Precision>::BASE_CUTOFF;
        partitioner = std::make_unique<lb::BinaryBoxEstimatingSamplePartitioner<D, Precision>>(baseCutoff, std::move(sampler));
    } else if("center-distance-pasp" == partitionerName){
        partitioner = std::make_unique<lb::CenterDistancePointAssigningSamplePartitioner<D, Precision>>(threads, std::move(sampler));
    } else if("bounds-distance-pasp" == partitionerName){
        partitioner = std::make_unique<lb::BoundsDistancePointAssigningSamplePartitioner<D, Precision>>(threads, std::move(sampler));
    } else if("nearest-sample-pasp" == partitionerName) {
        partitioner = std::make_unique<lb::NearestSamplePointAssigningSamplePartitioner<D, Precision>>(threads, std::move(sampler));
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
    return partitioner;
}

template <typename Precision>
po::options_description defaultOptions() {
    po::options_description cCommandLine("Command Line Options");
    cCommandLine.add_options()("n", po::value<tIdType>(), "number of points");
    cCommandLine.add_options()("partitioner", po::value<std::string>());
    cCommandLine.add_options()("distribution", po::value<std::string>());
    //cCommandLine.add_options()("points", po::value(&pointFile), "load points from file");
    //cCommandLine.add_options()("sample-size", po::value<uint>());
    cCommandLine.add_options()("num-bubbles", po::value<uint>());
    cCommandLine.add_options()("bubble-radius", po::value<Precision>());
    cCommandLine.add_options()("split-dimension", po::value<uint>());
    cCommandLine.add_options()("help", "produce help message");
    return cCommandLine;
}

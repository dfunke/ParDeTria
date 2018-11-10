#pragma once
#include <kaHIP_interface.h>

#include "utils/Generator.h"
#include "load_balancing/sample_partitioner/Sampler.h"
#include "load_balancing/UnboundBubbleGenerator.h"
#include "load_balancing/PatchBubbleGenerator.h"
#include "load_balancing/MaliciousBubbleGenerator.h"

#include "LoadBalancedDCTriangulator.h"
#include "load_balancing/Partitioner.h"
#include "load_balancing/OldPartitionerPartitioner.h"
#include "load_balancing/SimplePartitioner.h"
#include "load_balancing/sample_partitioner/BinaryBoxEstimatingSamplePartitioner.h"
#include "load_balancing/sample_partitioner/CenterDistancePointAssigningSamplePartitioner.h"
#include "load_balancing/sample_partitioner/BoundsDistancePointAssigningSamplePartitioner.h"
#include "load_balancing/sample_partitioner/NearestSamplePointAssigningSamplePartitioner.h"
#include "load_balancing/sample_partitioner/NearestSamplePointAssigningSampleBipartitioner.h"
#include "load_balancing/sample_partitioner/HyperplaneBipartitioner.h"
#include "load_balancing/BoundsIntersectionChecker.h"
#include "load_balancing/GridIntersectionChecker.h"

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
    } else if("malicious-bubble" == distributionName){
        uint numBubbles = extract<uint>("num-bubbles", vm, std::pow(4, D));
        result = std::make_unique<lb::MaliciousBubbleGenerator<D, Precision>>(numBubbles);
    }
    
    return result;
}

template <uint D, typename Precision>
auto createIntersectionPartitionMakerFunction(const po::variables_map& vm) {
	lb::IntersectionPartitionMakerFunction<D, Precision> pmf;
    if(vm.count("bounds") <= 0 || vm["bounds"].as<std::string>() == "grid") {
	    Precision cellWidth = vm["cell-width"].as<Precision>();
		//pm = std::make_unique<lb::GridIntersectionPartitionMaker<D, Precision>>
		//	(lb::Grid<D, Precision>{cellWidth});
		pmf = lb::GridIntersectionPartitionMaker<D, Precision>(lb::Grid<D, Precision>{cellWidth});
	} else if(vm["bounds"].as<std::string>() == "box" || vm["bounds"].as<std::string>() == "boxes") {
		pmf = lb::BoundsIntersectionPartitionMaker<D, Precision>();
	} else {
		throw std::runtime_error("unknown option for argument 'bounds'");
	}
    return pmf;
}
template <uint D, typename Precision, typename MonitorT, typename RandomGenerator>
std::unique_ptr<lb::Partitioner<D, Precision, MonitorT>> createPartitioner(const po::variables_map& vm,
																		   uint threads,
																		   RandomGenerator& rand) {
    std::string partitionerName = vm["partitioner"].as<std::string>();

	std::function<size_t(size_t)> f;
	if(vm.count("sample-size") == 0 && vm.count("sampling") == 0) {
		f = [](size_t n) -> size_t { return std::sqrt(n); };
	} else if(vm.count("sample-size") == 0) {
		auto sampling = vm["sampling"].as<std::string>();
		if("sqrt" == sampling) {
			f = [](size_t n) -> size_t { return std::sqrt(n); };
		} else if("log" == sampling) {
			f = [](size_t n) -> size_t { return std::log2(n); };
		} else {
			throw std::runtime_error("unknown option for argument 'sampling'");
		}
	} else if(vm.count("sampling") == 0) {
		f = [c = vm["sample-size"].as<size_t>()](size_t /*n*/) -> size_t { return c; };
	} else {
		throw std::runtime_error("illegaly set mutually exclusive"
		                         "options 'sample-size' and 'sampling'");
	}

	bool uniformEdges = vm.count("edge-weights") == 0
		|| "uniform" == vm["edge-weights"].as<std::string>();

	int mode = FAST;
	if(vm.count("kaffpa-mode") > 0) {
		auto modeString = vm["kaffpa-mode"].as<std::string>();
		if("ECO" == modeString) {
			mode = ECO;
		} else if("STRONG" == modeString) {
			mode = STRONG;
		} else if("PAR" == modeString) {
            mode = FASTSOCIALMULTITRY_PARALLEL;
        } else if("FAST" != modeString) {
			std::cerr << "unsupported kaffpa-mode '" << modeString << "'\n";
		}
	}

	typename lb::Sampler<D, Precision>::DistanceToEdgeWeightFunction weight;
	if(!uniformEdges) {
		if("squared-inverted" == vm["edge-weights"].as<std::string>()) {
			weight = [](Precision d) -> Precision { return 1/(d + 0.0000001); };
		} else if("inverted" == vm["edge-weights"].as<std::string>()) {
			weight = [](Precision d) -> Precision { return 1/std::sqrt(d + 0.0000001); };
		} else if("linear" == vm["edge-weights"].as<std::string>()) {
			weight = [](Precision d) -> Precision { return 1.0 - d; };
		} else if("log" == vm["edge-weights"].as<std::string>()) {
			weight = [](Precision d) -> Precision { return -std::log(d + 0.0000001); };
		}
	}
	
	double imbalance = 0.05;
	if(vm.count("imbalance") > 0) {
		imbalance = vm["imbalance"].as<double>();
	}

    auto sampler = uniformEdges ? lb::Sampler<D, Precision>(rand(), f, mode, imbalance)
	                            : lb::Sampler<D, Precision>(rand(), f, mode, imbalance, weight);
	//std::unique_ptr<lb::GridIntersectionPartitionMaker<D, Precision>> pm = nullptr;
	auto createIPMF = [](const auto& vm) -> auto {
		return createIntersectionPartitionMakerFunction<D, Precision>(vm);
	};
	
    auto baseCutoff = lb::DCTriangulator<D, Precision>::BASE_CUTOFF;
    std::unique_ptr<lb::Partitioner<D, Precision, MonitorT>> partitioner = nullptr;
    if("binary-besp" == partitionerName){
        auto baseCutoff = lb::DCTriangulator<D, Precision>::BASE_CUTOFF;
        partitioner =
	        std::make_unique<lb::BinaryBoxEstimatingSamplePartitioner<D, Precision, MonitorT>>
	        (baseCutoff, std::move(sampler));
    } else if("center-distance-pasp" == partitionerName){
        partitioner = std::make_unique<
	        lb::CenterDistancePointAssigningSamplePartitioner<D, Precision, MonitorT>>
	        (threads, std::move(sampler), createIPMF(vm));
    } else if("bounds-distance-pasp" == partitionerName){
        partitioner = std::make_unique<
	        lb::BoundsDistancePointAssigningSamplePartitioner<D, Precision, MonitorT>>
	        (threads, std::move(sampler), createIPMF(vm));
    } else if("nearest-sample-pasp" == partitionerName) {
        partitioner = std::make_unique<
	        lb::NearestSamplePointAssigningSamplePartitioner<D, Precision, MonitorT>>
	        (threads, std::move(sampler), createIPMF(vm));
    } else if("nearest-sample-pasb" == partitionerName) {
        partitioner = std::make_unique<
	        lb::NearestSamplePointAssigningSampleBipartitioner<D, Precision, MonitorT>>
	        (threads, std::move(sampler), baseCutoff, createIPMF(vm));
    } else if("hyperplane-sample-pasb" == partitionerName) {
        partitioner = std::make_unique<
	        lb::HyperplaneSamplingBipartitioner<D, Precision, MonitorT>>
	        (threads, std::move(sampler));
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
        partitioner =
	        std::make_unique<LoadBalancing::OldPartitionerPartitioner<D, Precision, MonitorT>>(
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
    cCommandLine.add_options()("sample-size", po::value<size_t>());
    cCommandLine.add_options()("sampling", po::value<std::string>());
    cCommandLine.add_options()("edge-weights", po::value<std::string>());
    cCommandLine.add_options()("imbalance", po::value<double>());
    cCommandLine.add_options()("num-bubbles", po::value<uint>());
    cCommandLine.add_options()("bubble-radius", po::value<Precision>());
    cCommandLine.add_options()("cell-width", po::value<Precision>());
    cCommandLine.add_options()("bounds", po::value<std::string>());
    cCommandLine.add_options()("kaffpa-mode", po::value<std::string>());
    cCommandLine.add_options()("split-dimension", po::value<uint>());
    cCommandLine.add_options()("validate", po::value<uint>());
    cCommandLine.add_options()("seed", po::value<uint>());
    cCommandLine.add_options()("help", "produce help message");
    return cCommandLine;
}

#pragma once

#include <chrono>
#include "load_balancing/PartitionTree.h"
#include "load_balancing/Partitioner.h"

namespace LoadBalancing
{
    struct Monitor
    {
        void registerPartitionStart() {}
        void registerPartitionEnd() {}
        template <uint D, typename Precision, typename MonitorT>
	    void registerPartition(const PartitionTree<D, Precision>& /*tree*/,
	                           const Partitioner<D, Precision, MonitorT>& /*partitioner*/) {}
        void registerSampleTriangulation(size_t /*numPoints*/, const std::string& /*provenance*/) {}
        template <uint D, typename Precision>
        void registerSampling(const Sampling<D, Precision>& /*sampling*/,
                              const std::string& /*provenance*/) {}
        template <uint D, typename Precision>
        void registerPartialTriangulation
	        (const std::vector<dSimplices<D, Precision>>& /*partialTriangulations*/,
	         const dSimplices<D, Precision>& /*borderTriangulation*/,
	         const dPoints<D, Precision>& /*points*/,
	         const std::string& /*provenance*/) {}
        void registerTriangulationStart() {}
        void registerTriangulationEnd() {}
        void registerBaseTriangulation(size_t /*numPoints*/, const std::string& /*provenance*/) {}
    };
}

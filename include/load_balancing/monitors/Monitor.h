#pragma once

#include <chrono>
#include "load_balancing/PartitionTree.h"
#include "load_balancing/Partitioner.h"

namespace LoadBalancing
{
    struct Monitor
    {
        void registerPartitionStart() {};
        void registerPartitionEnd() {};
        template <uint D, typename Precision>
	    void registerPartition(const PartitionTree<D, Precision>& /*tree*/, const Partitioner<D, Precision>& /*partitioner*/) {}
        void registerTriangulationStart() {};
        void registerTriangulationEnd() {};
        void registerBaseTriangulation(size_t /*numPoints*/, const std::string& /*provenance*/) {};
    };
}

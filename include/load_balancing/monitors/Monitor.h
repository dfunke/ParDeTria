#pragma once

#include <chrono>
#include "load_balancing/PartitionTree.h"

namespace LoadBalancing
{
    struct Monitor
    {
        void registerPartitionStart() {};
        void registerPartitionEnd() {};
	template <uint D, typename Precision>
	void registerPartition(const PartitionTree<D, Precision>& /*tree*/) {}
        void registerTriangulationStart() {};
        void registerTriangulationEnd() {};
        void registerBaseTriangulation(size_t /*numPoints*/, const std::string& /*provenance*/) {};
    };
}

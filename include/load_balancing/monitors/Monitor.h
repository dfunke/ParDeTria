#pragma once

#include <chrono>

namespace LoadBalancing
{
    struct Monitor
    {
        void registerPartitionStart() {};
        void registerPartitionEnd() {};
        void registerTriangulationStart() {};
        void registerTriangulationEnd() {};
        void registerBaseTriangulation(size_t /*numPoints*/, const std::string& /*provenance*/) {};
    };
}

#pragma once
#include "TimingMonitor.h"
#include <chrono>

namespace LoadBalancing
{
    struct ComprehensiveAccumulator : public TimingAccumulator {
        size_t numTriangulatedPoints = 0;
    };
    
    struct ComprehensiveMonitor : public TimingMonitor
    {
        ComprehensiveMonitor(ComprehensiveAccumulator& accumulator)
          : TimingMonitor(accumulator),
            compAcc(&accumulator) {
        }
        
        void registerBaseTriangulation(size_t numPoints, const std::string& /*provenance*/) {
            compAcc->numTriangulatedPoints += numPoints;
        };
        
    private:
        ComprehensiveAccumulator* compAcc;
    };
}

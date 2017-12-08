#pragma once
#include "Monitor.h"
#include <chrono>

namespace LoadBalancing
{
    struct TimingAccumulator {
        using Duration = std::chrono::milliseconds;
        Duration partitionTime;
        Duration triangulationTime;
    };
    
    struct TimingMonitor : public Monitor
    {
        TimingMonitor(TimingAccumulator& accumulator)
          :  accumulator(&accumulator) {
              this->accumulator->partitionTime = std::chrono::milliseconds(0);
              this->accumulator->triangulationTime = std::chrono::milliseconds(0);
        }
        
        void registerPartitionStart() {
            partitionStart = Clock::now();
        };
        
        void registerPartitionEnd() {
            auto now = Clock::now();
            accumulator->partitionTime = std::chrono::duration_cast<TimingAccumulator::Duration>(now - partitionStart);
        };
        
        void registerTriangulationStart() {
            triangulationStart = Clock::now();
        };
            
        void registerTriangulationEnd() {
            auto now = Clock::now();
            accumulator->triangulationTime = std::chrono::duration_cast<TimingAccumulator::Duration>(now - triangulationStart);
        };
        
        static const std::string unitString;
        
    private:
        using Clock = std::chrono::high_resolution_clock;
        TimingAccumulator* accumulator;
        std::chrono::time_point<Clock> partitionStart;
        std::chrono::time_point<Clock> triangulationStart;
    };
}

#pragma once
#include <ostream>
#include <chrono>
#include "Partitioner.h"
#include "../LoadBalancedDCTriangulator.h"
#include "../utils/MakeIds.h"
#include "monitors/TimingMonitor.h"

namespace LoadBalancing
{
    template <uint D, typename Precision>
    struct Experiment
    {
        struct Setup
        {
            dBox<D, Precision> bounds;
            dPoints<D, Precision> points;
            bool verify;
        };
        
        Experiment(std::unique_ptr<Partitioner<D, Precision>> partitioner, Setup setup, std::ostream& out);
        ~Experiment();
        void runOnce();
        
    private:
        Setup setup;
        std::ostream& out;
        DCTriangulator<D, Precision, TimingMonitor> triangulator;
        TimingMonitor::Accumulator timings;
        
        void finish();
        bool verify(const dSimplices<D, Precision>& dt);
    };
    
    template<uint D, typename Precision>
    Experiment<D, Precision>::Experiment(std::unique_ptr<Partitioner<D, Precision>> partitioner, Setup setup, std::ostream& out)
        : setup(std::move(setup)),
          out(out),
          triangulator(this->setup.bounds, this->setup.points, std::move(partitioner), 100, false, false, true, TimingMonitor(timings)) {
        out << "{\n"
        << "    'results' = [\n";
    }
    
    template<uint D, typename Precision>
    Experiment<D, Precision>::~Experiment() {
        finish();
    }
    
    template<uint D, typename Precision>
    void Experiment<D, Precision>::runOnce() {
        auto dt = triangulator.triangulate();
            
        out
        << "        {\n"
        << "            'partitiontime' = " << timings.partitionTime.count() << ",\n"
        << "            'triangulationtime' = " << timings.triangulationTime.count() << ",\n";
        if(setup.verify) {
            out
            << "            'valid' = " << (verify(dt) ? "true" : "false") << "\n";
        }out
        << "        },\n";
    }
    
    template<uint D, typename Precision>
    void Experiment<D, Precision>::finish() {
        out
        << "    ]\n"
        << "}\n";
    }
    
    template<uint D, typename Precision>
    bool Experiment<D, Precision>::verify(const dSimplices<D, Precision>& dt) {
        CGALTriangulator<D, Precision, false> cgal(setup.bounds, setup.points);
        auto realDT = cgal.triangulate();

        auto vr = dt.verify(setup.points);
        auto ccr = dt.crossCheck(realDT);

        return vr.valid && ccr.valid;
    }
    
}

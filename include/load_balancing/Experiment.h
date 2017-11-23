#pragma once
#include <ostream>
#include <chrono>
#include "Partitioner.h"
#include "../LoadBalancedDCTriangulator.h"
#include "../utils/MakeIds.h"

namespace LoadBalancing
{
    template <uint D, typename Precision>
    struct Experiment
    {
        struct Setup
        {
            DCTriangulator<D, Precision>* triangulator;
            std::unique_ptr<Partitioner<D, Precision>> partitioner;
            dBox<D, Precision> bounds;
            dPoints<D, Precision> points;
        };
        
        Experiment(Setup setup, std::ostream& out);
        ~Experiment();
        void runOnce();
        
    private:
        Setup setup;
        std::ostream& out;
        
        void finish();
    };
    
    template<uint D, typename Precision>
    Experiment<D, Precision>::Experiment(Setup setup, std::ostream& out)
        : setup(std::move(setup)), out(out) {
        out << "{\n"
        << "    'results' = [\n";
    }
    
    template<uint D, typename Precision>
    Experiment<D, Precision>::~Experiment() {
        finish();
    }
    
    template<uint D, typename Precision>
    void Experiment<D, Precision>::runOnce() {
        auto pointIds = makePointIds(setup.points);
        auto t0 = std::chrono::high_resolution_clock::now();
        auto partitioning = setup.partitioner->partition(setup.bounds, setup.points, pointIds);
        auto t1 = std::chrono::high_resolution_clock::now();
        auto dt = setup.triangulator->triangulateTree(partitioning);
        auto t2 = std::chrono::high_resolution_clock::now();
        
        namespace chrono = std::chrono;
        using ms = chrono::milliseconds;
    
        out
        << "        {\n"
        << "            'partitiontime' = " << chrono::duration_cast<ms>(t1 - t0).count() << ",\n"
        << "            'triangulationtime' = " << chrono::duration_cast<ms>(t2 - t1).count() << "\n"
        << "        },\n";
    }
    
    template<uint D, typename Precision>
    void Experiment<D, Precision>::finish() {
        out
        << "    ]\n"
        << "}\n";
    }
    
}

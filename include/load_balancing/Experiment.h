#pragma once
#include <ctime>
#include <ostream>
#include <chrono>
#include <sstream>
#include "utils/version.h"
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
            std::string distribution;
        };
        
        Experiment(std::unique_ptr<Partitioner<D, Precision>> partitioner, Setup setup, std::ostream& out);
        ~Experiment();
        void runOnce();
        
    private:
        Setup setup;
        std::ostream& out;
        std::string partitionerName;
        DCTriangulator<D, Precision, TimingMonitor> triangulator;
        TimingMonitor::Accumulator timings;
        size_t runs;
        
        void start();
        void finish();
        bool verify(const dSimplices<D, Precision>& dt);
        std::string makeDateTimeString() const;
    };
    
    template<uint D, typename Precision>
    Experiment<D, Precision>::Experiment(std::unique_ptr<Partitioner<D, Precision>> partitioner, Setup setup, std::ostream& out)
        : setup(std::move(setup)),
          out(out),
          partitionerName(partitioner->info()),
          triangulator(this->setup.bounds, this->setup.points, std::move(partitioner), 100, false, false, true, TimingMonitor(timings)),
          runs(0) {
        start();
    }
    
    template<uint D, typename Precision>
    Experiment<D, Precision>::~Experiment() {
        finish();
    }
    
    template<uint D, typename Precision>
    void Experiment<D, Precision>::runOnce() {
        auto dt = triangulator.triangulate();
            
        out
            << (runs == 0 ? "" : ",") << " {\n"
            << "            \"partitiontime\": " << timings.partitionTime.count() << ",\n"
            << "            \"triangulationtime\": " << timings.triangulationTime.count() << ",\n";
            if(setup.verify) {
                out
                << "            \"valid\": " << (verify(dt) ? "true" : "false") << "\n";
            }
            out
            << "        }\n";
            
            ++runs;
    }
    
    template<uint D, typename Precision>
    void Experiment<D, Precision>::start() {
        out << "{\n"
            << "    \"commit\": \"" << GIT_COMMIT << "\",\n"
            << "    \"date\": \"" << makeDateTimeString() << "\",\n"
            << "    \"triangulator\": \"" << "load balanced DC Triangulator\",\n"
            << "    \"partitioner\": \"" << partitionerName << "\",\n"
            << "    \"distribution\": \"" << setup.distribution << "\",\n"
            << "    \"numpoints\": " << setup.points.finite_size() << ",\n"
            << "    \"results\": [\n";
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
    
    template<uint D, typename Precision>
    std::string Experiment<D, Precision>::makeDateTimeString() const {
        using SystemClock = std::chrono::system_clock;
        auto now = SystemClock::now();
        std::time_t time = SystemClock::to_time_t(now);
        std::tm tm = *std::localtime(&time);
        std::stringstream stream;
        stream << std::put_time(&tm, "%FT%T%z");
        return stream.str();
    }
}

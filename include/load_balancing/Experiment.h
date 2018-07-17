#pragma once
#include <ctime>
#include <ostream>
#include <chrono>
#include <sstream>
#include "utils/version.h"
#include "Partitioner.h"
#include "../LoadBalancedDCTriangulator.h"
#include "../utils/MakeIds.h"
#include "monitors/ComprehensiveMonitor.h"

namespace LoadBalancing
{
    template <uint D, typename Precision>
    struct Experiment
    {
        struct Setup
        {
            std::string argumentString;
            dBox<D, Precision> bounds;
            dPoints<D, Precision> points;
            bool verify;
            std::string distribution;
            size_t threads;
        };
        
        Experiment(std::unique_ptr<Partitioner<D, Precision>> partitioner, Setup setup, std::ostream& out);
        ~Experiment();
        void runOnce();
        
    private:
        Setup setup;
        std::ostream& out;
        std::string partitionerName;
        DCTriangulator<D, Precision, ComprehensiveMonitor> triangulator;
        ComprehensiveAccumulator acc;
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
          //triangulator(this->setup.bounds, this->setup.points, std::move(partitioner), 100, false, (partitioner->info().find("bipartitioner") == std::string::npos), true, ComprehensiveMonitor(acc)),
          triangulator(this->setup.bounds, this->setup.points, std::move(partitioner), 100, false, false, true, ComprehensiveMonitor(acc)),
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
            << "            \"partitiontime\": " << acc.partitionTime.count() << ",\n"
            << "            \"triangulationtime\": " << acc.triangulationTime.count() << ",\n"
            << "            \"triangulatedPoints\": " << acc.numTriangulatedPoints << ",\n"
            << "            \"sampleSize\": " << acc.sampleSize << ",\n";
            if(setup.verify) {
                out
                << "            \"valid\": " << (verify(dt) ? "true" : "false") << ",\n";
            }
	    size_t partitionCount = 0;
            out
	    << "            \"partitions\": [";
	    for(size_t size : acc.partitionSizes) {
		out
		<< (partitionCount == 0 ? "" : ",") << " {\n"
	        << "                      \"size\": " << size << "\n"
		<< "                }";
		++partitionCount;
	    }
	    out
	    << "\n"
	    << "            ],\n";

		out
	    << "            \"basetriangulations\": [";
		size_t baseCount = 0;
		for(const auto& base : acc.baseTriangulations) {
			out
			<<  (baseCount == 0 ? "" : ",") << " {\n"
			<< "                        \"points\": " << base.first << ",\n"
			<< "                        \"provenance\": \"" << base.second << "\"\n"
			<< "                }";
			++baseCount;
		}
	    out
	    << "\n"
	    << "            ]\n"
        << "        }\n";
            
	    ++runs;
    }
    
    template<uint D, typename Precision>
    void Experiment<D, Precision>::start() {
        out << "{\n"
            << "    \"commit\": \"" << GIT_COMMIT << "\",\n"
            << "    \"date\": \"" << makeDateTimeString() << "\",\n"
            << "    \"args\": \"" << setup.argumentString <<"\",\n"
            << "    \"triangulator\": \"" << "load balanced DC Triangulator\",\n"
            << "    \"partitioner\": \"" << partitionerName << "\",\n"
            << "    \"distribution\": \"" << setup.distribution << "\",\n"
            << "    \"numpoints\": " << setup.points.finite_size() << ",\n"
            << "    \"threads\": " << setup.threads << ",\n"
            << "    \"results\": [";
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

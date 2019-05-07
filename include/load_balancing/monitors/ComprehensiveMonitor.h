#pragma once
#include <cstdlib>
#include <iostream>
#include <vector>
#include <tuple>
#include <mutex>
#include "TimingMonitor.h"

namespace LoadBalancing
{
	struct ComprehensiveAccumulator : public TimingAccumulator {
		size_t numTriangulatedPoints = 0;
		std::vector<size_t> partitionSizes;
		std::vector<std::pair<size_t, std::string>> baseTriangulations;
		size_t sampleSize = 0;
		
		std::mutex mutex;
    };
    
    struct ComprehensiveMonitor : public TimingMonitor
    {
        ComprehensiveMonitor(ComprehensiveAccumulator& accumulator)
          : TimingMonitor(accumulator),
            compAcc(&accumulator) {
        }
        
        void registerBaseTriangulation(size_t numPoints, const std::string& provenance) {
	        std::lock_guard<std::mutex> lock(compAcc->mutex);
	        
            compAcc->numTriangulatedPoints += numPoints;
            compAcc->baseTriangulations.emplace_back(numPoints, provenance);
        }
        
        void registerSampleTriangulation(size_t numPoints, const std::string& provenance) {
	        registerBaseTriangulation(numPoints, provenance);
        }

        template <uint D, typename Precision, typename MonitorT>
	    void registerPartition(const PartitionTree<D, Precision>& tree,
	                           const Partitioner<D, Precision, MonitorT>& partitioner) {
	        registerPartitionRecursively(tree);
	        auto sp = dynamic_cast<const SamplePartitioner<D, Precision, MonitorT>*>(&partitioner);
	        if(sp) {
		        compAcc->sampleSize = sp->sampling().partition.size();
			}
        }
		
    private:
        ComprehensiveAccumulator* compAcc;

        template <uint D, typename Precision>
		void registerPartitionRecursively(const PartitionTree<D, Precision>& tree) {
				if(tree.isLeaf()) {
					const auto& pointIds = std::get<Point_Ids>(tree.attachment);
				compAcc->partitionSizes.push_back(pointIds.size());
				} else {
					const auto& children = std::get<typename lb::PartitionTree<D, Precision>::ChildContainer>(tree.attachment);
					for(const auto& child : children) {
						registerPartitionRecursively(child);
					}
			}
        }
    };
}

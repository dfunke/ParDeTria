#pragma once
#include <cstdlib>
#include <iostream>
#include <vector>
#include <tuple>
#include "TimingMonitor.h"

namespace LoadBalancing
{
	struct ComprehensiveAccumulator : public TimingAccumulator {
		size_t numTriangulatedPoints = 0;
		std::vector<size_t> partitionSizes;
		std::vector<std::pair<size_t, std::string>> baseTriangulations;
    };
    
    struct ComprehensiveMonitor : public TimingMonitor
    {
        ComprehensiveMonitor(ComprehensiveAccumulator& accumulator)
          : TimingMonitor(accumulator),
            compAcc(&accumulator) {
        }
        
        void registerBaseTriangulation(size_t numPoints, const std::string& provenance) {
            compAcc->numTriangulatedPoints += numPoints;
            compAcc->baseTriangulations.emplace_back(numPoints, provenance);
        };

	template <uint D, typename Precision>
	void registerPartition(const PartitionTree<D, Precision>& tree) {
	    registerPartitionRecursively(tree);
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

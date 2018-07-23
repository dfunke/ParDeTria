#pragma once

#include <cassert>
#include <random>
#include <algorithm>
#include <set>
#include "SamplePartitioner.h"
#include "Sampler.h"
#include "load_balancing/IntersectionChecker.h"

namespace LoadBalancing
{
    template <uint D, typename Precision>
    std::vector<dVector<D, Precision>> findPartitionCenters(const std::vector<int>& partitioning,
                                                            size_t partitions,
                                                            const std::vector<dVector<D, Precision>>& samplePoints) {
        std::vector<size_t> numPoints(partitions, 0);
        std::vector<dVector<D, Precision>> result(partitions, dVector<D, Precision>());
        for(size_t i = 0; i < partitioning.size(); ++i) {
            int partition = partitioning[i];
            ++numPoints[partition];
            result[partition] = result[partition] * ((Precision)numPoints[partition]/(1 + numPoints[partition]))
                              + samplePoints[i] * ((Precision)1/(1 + numPoints[partition]));
        }
        return result;
    }
    
    template <uint D, typename Precision>
    struct CenterDistancePointAssigningSamplePartitioner : public SamplePartitioner<D, Precision>
    {
        CenterDistancePointAssigningSamplePartitioner(size_t partitionSize, Sampler<D, Precision> sampler, IntersectionPartitionMakerFunction<D, Precision> intersectionCheckerMaker)
            : mPartitionSize(partitionSize), mSampler(std::move(sampler)), mMakePartition(std::move(intersectionCheckerMaker))			  
        { }
        PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds) override
        {
            mSampling = mSampler(bounds, points, pointIds, mPartitionSize);
            mPartitionCenters = findPartitionCenters<D, Precision>(mSampling.partition, mPartitionSize, mSampling.points);
            auto partitioning = makePartitioning(mPartitionCenters, points, pointIds, mPartitionSize);
            
            typename PartitionTree<D, Precision>::ChildContainer subtrees;
            std::transform(partitioning.begin(), partitioning.end(), std::back_inserter(subtrees), [](auto& partition) {
                PartitionTree<D, Precision> subtree;
                subtree.attachment = std::move(partition.pointIds);
            	subtree.intersectionChecker = std::move(partition.intersectionChecker);
                return subtree;
            });
            
            PartitionTree<D, Precision> tree;
            tree.intersectionChecker = std::make_unique<BoundsIntersectionChecker<D, Precision>>(bounds);
            tree.attachment = std::move(subtrees);
            return tree;
        }
        
        std::string info() const override
        {
            return "center-distance-point-assigning sample partitioner";
        }
        
        const Sampling<D, Precision>& sampling() const {
            return mSampling;
        }

	const std::vector<dVector<D, Precision>>& partitionCenters() const {
	    return mPartitionCenters;
	}
        
    private:
        Sampling<D, Precision> mSampling;
        size_t mPartitionSize;
        Sampler<D, Precision> mSampler;
		std::vector<dVector<D, Precision>> mPartitionCenters;
		IntersectionPartitionMakerFunction<D, Precision> mMakePartition;
		
		std::vector<IntersectionPartition<D, Precision>>
		makePartitioning(const std::vector<dVector<D, Precision>>& partitionPoints,
		                 const dPoints<D, Precision>& points,
		                 const Point_Ids& pointIds,
		                 size_t partitions) {
			return mMakePartition(points, pointIds, partitions, [&] (auto id) -> size_t {
                const auto& coords = points[id].coords;
                auto it = std::min_element(partitionPoints.begin(), partitionPoints.end(),
                                        [&coords](const dVector<D, Precision>& left, const dVector<D, Precision>& right) {
                                            return lenSquared(left - coords) < lenSquared(right - coords);
                                        });
                                            
                auto partition = std::distance(partitionPoints.begin(), it);
				return partition;
			});
		}   
    };
}

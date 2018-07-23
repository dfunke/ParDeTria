#pragma once

#include <cassert>
#include <random>
#include <algorithm>
#include <set>
#include "SamplePartitioner.h"
#include "load_balancing/BoxUtils.h"
#include "Sampler.h"
#include "load_balancing/Grid.h"
#include "load_balancing/IntersectionChecker.h"

namespace LoadBalancing
{
    template <uint D, typename Precision>
    std::vector<dBox<D, Precision>> makeBoundingBoxes(const std::vector<int>& partitioning,
                                                      size_t partitions,
                                                      const std::vector<dVector<D, Precision>>& samplePoints) {
        
        std::vector<std::vector<dVector<D, Precision>>> samplePartitions(partitions);
        for(size_t i = 0; i < partitioning.size(); ++i) {
            int partition = partitioning[i];
            samplePartitions[partition].push_back(samplePoints[i]);
        }
        
        std::vector<dBox<D, Precision>> result(partitions);
        for(size_t i = 0; i < partitions; ++i) {
            result[i] = makeBoundingBox(samplePartitions[i].begin(), samplePartitions[i].end());            
        }
        return result;
    }
    
    template <uint D, typename Precision>
    struct BoundsDistancePointAssigningSamplePartitioner : public SamplePartitioner<D, Precision>
    {
        BoundsDistancePointAssigningSamplePartitioner(size_t partitionSize, Sampler<D, Precision> sampler, IntersectionPartitionMakerFunction<D, Precision> intersectionCheckerMaker)
            : mPartitionSize(partitionSize), mSampler(std::move(sampler)), mMakePartition(std::move(intersectionCheckerMaker))

        { }
        PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds) override
        {
            mSampling = mSampler(bounds, points, pointIds, mPartitionSize);
            auto boxes = makeBoundingBoxes<D, Precision>(mSampling.partition, mPartitionSize, mSampling.points);
            auto partitioning = makePartitioning(boxes, points, mPartitionSize, pointIds);
            
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
            return "bounds-distance-point-assigning sample partitioner";
        }
        
        const Sampling<D, Precision>& sampling() const {
            return mSampling;
        }
        
    private:
        Sampling<D, Precision> mSampling;
        size_t mPartitionSize;
        Sampler<D, Precision> mSampler;
		IntersectionPartitionMakerFunction<D, Precision> mMakePartition;
    
		std::vector<IntersectionPartition<D, Precision>>
		makePartitioning(const std::vector<dBox<D, Precision>>& boxes,
		                 const dPoints<D, Precision>& points,
		                 size_t partitions,
		                 const Point_Ids& pointIds) {
			return mMakePartition(points, pointIds, partitions, [&] (auto id) -> size_t {
				const auto& coords = points[id].coords;
				std::vector<Precision> boxDistances(boxes.size());
				std::transform(boxes.begin(), boxes.end(), boxDistances.begin(),
							   [&coords](const dBox<D, Precision>& box) {
									return lenSquared(vecToBox<D, Precision>(coords, box));
								});
				auto it = std::min_element(begin(boxDistances), end(boxDistances));
				auto partition = std::distance(boxDistances.begin(), it);
				return partition;
			});
		}   
    };
}

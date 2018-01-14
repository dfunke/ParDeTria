#pragma once

#include <cassert>
#include <random>
#include <algorithm>
#include <set>
#include "SamplePartitioner.h"
#include "Sampler.h"

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
    auto makePartitioning(const std::vector<dVector<D, Precision>>& partitionPoints,
                                            const dPoints<D, Precision>& points, const Point_Ids& pointIds) {
        struct ResultType {
            Point_Ids pointIds;
            dBox<D, Precision> bounds;
        };
        
        std::vector<ResultType> result(partitionPoints.size());
        std::vector<bool> initialized(result.size(), false);
        for(auto id : pointIds) {
            if(dPoint<D, Precision>::isFinite(id)) {
                const auto& coords = points[id].coords;
                auto it = std::min_element(partitionPoints.begin(), partitionPoints.end(),
                                        [&coords](const dVector<D, Precision>& left, const dVector<D, Precision>& right) {
                                            return lenSquared(left - coords) < lenSquared(right - coords);
                                        });
                                            
                auto partition = std::distance(partitionPoints.begin(), it);
                result[partition].pointIds.insert(id);
                if(initialized[partition]) {
                    enlargeBoxAroundVector<D, Precision>(result[partition].bounds, coords);
                } else {
                    result[partition].bounds.low = coords;
                    result[partition].bounds.high = result[partition].bounds.low;
                    initialized[partition] = true;
                }
            }
        }
        
        for (tIdType k = dPoint<D, Precision>::cINF; k != 0; ++k) {
            for(auto& e : result) {
                e.pointIds.insert(k);
            }
        }
        
        return result;
    }    
    
    template <uint D, typename Precision>
    struct CenterDistancePointAssigningSamplePartitioner : public SamplePartitioner<D, Precision>
    {
        CenterDistancePointAssigningSamplePartitioner(size_t partitionSize, Sampler<D, Precision> sampler)
            : mPartitionSize(partitionSize), mSampler(std::move(sampler))
        { }
        PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        const dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds) override
        {
            mSampling = mSampler(bounds, points, pointIds, mPartitionSize);
            mPartitionCenters = findPartitionCenters<D, Precision>(mSampling.partition, mPartitionSize, mSampling.points);
            auto partitioning = makePartitioning<D, Precision>(mPartitionCenters, points, pointIds);
            
            typename PartitionTree<D, Precision>::ChildContainer subtrees;
            std::transform(partitioning.begin(), partitioning.end(), std::back_inserter(subtrees), [](auto partition) {
                PartitionTree<D, Precision> subtree;
                subtree.attachment = std::move(partition.pointIds);
                subtree.bounds = std::move(partition.bounds);
                return subtree;
            });
            
            PartitionTree<D, Precision> tree;
            tree.bounds = bounds;
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
    };
}

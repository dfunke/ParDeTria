#pragma once

#include <cassert>
#include <random>
#include <algorithm>
#include <set>
#include "../Partitioner.h"
#include "Sampler.h"

namespace LoadBalancing
{
    template <uint D, typename Precision>
    std::vector<dVector<D, Precision>> findPartitionCenters(const std::vector<int>& partitioning,
                                                            size_t partitions,
                                                            const dPoints<D, Precision>& samplePoints) {
        std::vector<size_t> numPoints(partitions);
        std::vector<dVector<D, Precision>> result(partitions);        
        for(size_t i = 0; i < partitioning.size(); ++i) {
            int partition = partitioning[i];
            ++numPoints[partition];
            result[partition] = result[partition] * ((Precision)numPoints[partition]/(1 + numPoints[partition]))
                              + samplePoints[i].coords * ((Precision)1/(1 + numPoints[partition]));
        }
        
        std::transform(result.begin(), result.end(), numPoints.begin(), result.begin(), [](const dVector<D, Precision>& v, size_t points) {
            return v * (1.0/points);
        });
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
    struct CenterDistancePointAssigningSamplePartitioner : public Partitioner<D, Precision>
    {
        CenterDistancePointAssigningSamplePartitioner(size_t partitionSize, Sampler<D, Precision> sampler)
            : mPartitionSize(partitionSize), mSampler(std::move(sampler))
        { }
        PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        const dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds) override
        {
            auto sampling = mSampler(bounds, points, pointIds, mPartitionSize);
            auto centers = findPartitionCenters(sampling.partition, mPartitionSize, sampling.points);
            auto partitioning = makePartitioning<D, Precision>(centers, points, pointIds);
            
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
        
    private:
        size_t mPartitionSize;
        Sampler<D, Precision> mSampler;
    };
}

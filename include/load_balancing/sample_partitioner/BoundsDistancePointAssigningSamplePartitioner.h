#pragma once

#include <cassert>
#include <random>
#include <algorithm>
#include <set>
#include "SamplePartitioner.h"
#include "load_balancing/BoxUtils.h"
#include "Sampler.h"

namespace LoadBalancing
{
    template <uint D, typename Precision>
    std::vector<dBox<D, Precision>> makeBoundingBoxes(const std::vector<int>& partitioning,
                                                      size_t partitions,
                                                      const dPoints<D, Precision>& samplePoints) {
        
        std::vector<std::vector<dVector<D, Precision>>> samplePartitions(partitions);
        for(size_t i = 0; i < partitioning.size(); ++i) {
            int partition = partitioning[i];
            samplePartitions[partition].push_back(samplePoints[i].coords);
        }
        
        std::vector<dBox<D, Precision>> result(partitions);
        for(size_t i = 0; i < partitions; ++i) {
            result[i] = makeBoundingBox(samplePartitions[i].begin(), samplePartitions[i].end());            
        }
        return result;
    }
    
    template <uint D, typename Precision>
    auto makePartitioningFromBoxes(const std::vector<dBox<D, Precision>>& boxes,
                                            const dPoints<D, Precision>& points, const Point_Ids& pointIds) {
        struct ResultType {
            Point_Ids pointIds;
            dBox<D, Precision> bounds;
        };
        
        std::vector<ResultType> result(boxes.size());
        std::vector<bool> initialized(result.size(), false);
        for(auto id : pointIds) {
            if(dPoint<D, Precision>::isFinite(id)) {
                const auto& coords = points[id].coords;
                /*auto it = std::min_element(boxes.begin(), boxes.end(),
                                        [&coords](const dBox<D, Precision>& left, const dBox<D, Precision>& right) {
                                            return lenSquared(vecToBox<D, Precision>(coords, left))
                                            < lenSquared(vecToBox<D, Precision>(coords, right));
                                        });*/
                std::vector<Precision> boxDistances(boxes.size());
                std::transform(boxes.begin(), boxes.end(), boxDistances.begin(),
                               [&coords](const dBox<D, Precision>& box) {
                                   return lenSquared(vecToBox<D, Precision>(coords, box));
                                   
                            });
                auto it = std::min_element(begin(boxDistances), end(boxDistances));
                                            
                auto partition = std::distance(boxDistances.begin(), it);
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
    struct BoundsDistancePointAssigningSamplePartitioner : public SamplePartitioner<D, Precision>
    {
        BoundsDistancePointAssigningSamplePartitioner(size_t partitionSize, Sampler<D, Precision> sampler)
            : mPartitionSize(partitionSize), mSampler(std::move(sampler))
        { }
        PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        const dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds) override
        {
            mSampling = mSampler(bounds, points, pointIds, mPartitionSize);
            auto boxes = makeBoundingBoxes(mSampling.partition, mPartitionSize, mSampling.points);
            auto partitioning = makePartitioningFromBoxes<D, Precision>(boxes, points, pointIds);
            
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
            return "bounds-distance-point-assigning sample partitioner";
        }
        
        const Sampling<D, Precision>& sampling() const {
            return mSampling;
        }
        
    private:
        Sampling<D, Precision> mSampling;
        size_t mPartitionSize;
        Sampler<D, Precision> mSampler;
    };
}

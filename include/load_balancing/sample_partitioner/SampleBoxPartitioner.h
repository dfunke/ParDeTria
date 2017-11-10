#pragma once

#pragma once

#include <cassert>
#include <random>
#include <algorithm>
#include <set>
#include "../Partitioner.h"
#include "steps.h"

namespace LoadBalancing
{
    template <typename ForwardIt>
    bool hasDuplicate(ForwardIt begin, ForwardIt end) {
        std::unordered_set<typename ForwardIt::value_type> set;
        bool foundDuplicate = false;
        while(begin != end) {
            if(set.count(*begin) > 0)
                foundDuplicate = true;
            set.insert(*begin++);
        }
        return foundDuplicate;
    }
    
    template <uint D, typename Precision>
    struct SampleBoxPartitioner : public Partitioner<D, Precision>
    {
        SampleBoxPartitioner(size_t sampleSize, size_t sampleSeed, size_t minNumOfPoints)
            : mSampleSize(sampleSize), mRand(sampleSeed), mPointsCutoff(minNumOfPoints)
        { }
        PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        const dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds) override
        {
            PartitionTree<D, Precision> tree;
                
            if(pointIds.size() >= std::max(mPointsCutoff, mSampleSize)) {
                auto sample = generateSample<D, Precision>(mSampleSize, pointIds, mRand);
                dPoints<D, Precision> samplePoints;
                for(auto id : sample) {
                    assert((dPoint<D, Precision>::isFinite(id)));
                    samplePoints.emplace_back(points[id]);
                }
                auto simplices = triangulateSample(bounds, samplePoints);
                auto graph = makeGraph(simplices);
                auto partitioning = partitionGraph(graph, mRand);
                auto centerEdges = findPartitionCenterEdges(graph, partitioning);
                auto centerPoints = makePartitionCenterPoints(centerEdges, samplePoints);
                auto boundingBoxes = estimateBoundingBoxes(centerPoints, bounds);
                auto pointIdsPair = seperatePointIds(points, pointIds, std::get<0>(boundingBoxes));

                // todo: alternative:
                // estimate seperation plain
                // separate point pointCloudSize

                // recurse
                
                /*assert((std::all_of(pointIds.begin(), pointIds.end(), [&pointIdsPair](tIdType id) -> bool {
                    return std::get<0>(pointIdsPair).contains(id) || std::get<1>(pointIdsPair).contains(id);
                })));*/
                
                PartitionTree<D, Precision> leftSubtree = partition(std::get<0>(boundingBoxes), points, std::get<0>(pointIdsPair));
                PartitionTree<D, Precision> rightSubtree = partition(std::get<1>(boundingBoxes), points, std::get<1>(pointIdsPair));
            
                tree.bounds = bounds;
                typename PartitionTree<D, Precision>::ChildContainer children{std::move(leftSubtree), std::move(rightSubtree)};
                tree.attachment = std::move(children);
            } else {
                tree.bounds = bounds;
                tree.attachment = pointIds;
            }
            
            return tree;
        }
        
        std::string info() const override
        {
            return "sample partitioner";
        }
        
    private:
        size_t mSampleSize;
        std::mt19937 mRand;
        size_t mPointsCutoff;
    };
}

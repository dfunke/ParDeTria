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
        SampleBoxPartitioner(size_t sampleSize, size_t sampleSeed)
            : mSampleSize(sampleSize), mSampleSeed(sampleSeed)
        { }
        PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        const dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds) const override
        {
            // assert(pointIds == {i + points.offset() | i elem N, i < number of points});
            /*assert(*std::min_element(pointIds.begin(), pointIds.end()) == points.offset());
            assert(*std::max_element(pointIds.begin(), pointIds.end()) == points.offset() + points.size());
            assert(!hasDuplicate(pointIds.begin(), pointIds.end()));*/
            assert(points.size() >= mSampleSize);
            
            std::mt19937 rand(mSampleSeed);
            
            auto sample = generateSample<D, Precision>(mSampleSize, pointIds, rand);
            dPoints<D, Precision> samplePoints;
            for(auto id : sample) {
                assert((dPoint<D, Precision>::isFinite(id)));
                samplePoints.emplace_back(points[id]);
            }
            auto simplices = triangulateSample(bounds, samplePoints);
            auto graph = makeGraph(simplices);
            auto partition = partitionGraph(graph, rand);
            auto centerEdges = findPartitionCenterEdges(graph, partition);
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
            
            PartitionTree<D, Precision> leftSubtree;
            leftSubtree.bounds = std::move(std::get<0>(boundingBoxes));
            leftSubtree.attachment = std::move(std::get<0>(pointIdsPair));

            PartitionTree<D, Precision> rightSubtree;
            rightSubtree.bounds = std::move(std::get<1>(boundingBoxes));
            rightSubtree.attachment = std::move(std::get<1>(pointIdsPair));
            
            PartitionTree<D, Precision> tree;
            tree.bounds = bounds;
            typename PartitionTree<D, Precision>::ChildContainer children{std::move(leftSubtree), std::move(rightSubtree)};
            tree.attachment = std::move(children);
            
            return tree;
        }
        
        std::string info() const override
        {
            return "sample partitioner";
        }
        
    private:
        size_t mSampleSize;
        size_t mSampleSeed;
    };
}

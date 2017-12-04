#pragma once

#include <cassert>
#include <random>
#include <algorithm>
#include <set>
#include <kdtree++/kdtree.hpp>
#include "../Partitioner.h"
#include "steps.h"

namespace LoadBalancing
{ 
    
    template <uint D, typename Precision>
    struct NearestSamplePointAssigningSamplePartitioner : public Partitioner<D, Precision>
    {
        NearestSamplePointAssigningSamplePartitioner(size_t sampleSize, size_t sampleSeed, size_t partitionSize)
            : mSampleSize(sampleSize), mRand(sampleSeed), mPartitionSize(partitionSize)
        { }
        PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        const dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds) override
        {
            auto sample = generateSample<D, Precision>(mSampleSize, pointIds, mRand);
            dPoints<D, Precision> samplePoints;
            for(auto id : sample) {
                assert((dPoint<D, Precision>::isFinite(id)));
                samplePoints.emplace_back(points[id]);
            }
            auto simplices = triangulateSample(bounds, samplePoints);
            auto graph = makeGraph(simplices);
            auto graphPartitioning = partitionGraph(graph, mPartitionSize, mRand);
            auto kdtree = buildKdTree(graphPartitioning, samplePoints);
            auto partitioning = makePartitioning(kdtree, mPartitionSize, points, pointIds);
            
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
            return "nearest-sample-point-assigning sample partitioner";
        }
        
    private:
        struct AssignedVector
        {
            typedef Precision value_type;
            dVector<D, Precision> position;
            int partitionId;
            Precision operator[](size_t d) const {
                return position[d];
            }
        };
        
        struct Partition {
            Point_Ids pointIds;
            dBox<D, Precision> bounds;
        };
        using Tree = KDTree::KDTree<D, AssignedVector>;
        Tree buildKdTree(const std::vector<int>& partitioning, const dPoints<D, Precision>& samplePoints);
        std::vector<Partition> makePartitioning(const Tree& tree, size_t numPartitions,
                              const dPoints<D, Precision>& points, const Point_Ids& pointIds);
         
        size_t mSampleSize;
        std::mt19937 mRand;
        size_t mPartitionSize;
    };
    
    template <uint D, typename Precision>
    typename NearestSamplePointAssigningSamplePartitioner<D, Precision>::Tree
    NearestSamplePointAssigningSamplePartitioner<D, Precision>::buildKdTree(
        const std::vector<int>& partitioning, const dPoints<D, Precision>& samplePoints) {
        
        Tree tree;
        for(size_t i = 0; i < partitioning.size(); ++i) {
            int partition = partitioning[i];
            tree.insert({samplePoints[i].coords, partition});
        }
        return tree;
    }
    
    template <uint D, typename Precision>
    std::vector<typename NearestSamplePointAssigningSamplePartitioner<D, Precision>::Partition>
    NearestSamplePointAssigningSamplePartitioner<D, Precision>::makePartitioning(
        const Tree& tree, size_t numPartitions, const dPoints<D, Precision>& points, const Point_Ids& pointIds) {
        std::vector<Partition> result(numPartitions);
        std::vector<bool> initialized(result.size(), false);
        for(auto id : pointIds) {
            if(dPoint<D, Precision>::isFinite(id)) {
                const auto& coords = points[id].coords;
                auto it_distance = tree.find_nearest(AssignedVector{coords, -1});
                assert(it_distance.first != tree.end());
                const auto& assignedVector = *(it_distance.first);
                                            
                auto partition = assignedVector.partitionId;
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
}

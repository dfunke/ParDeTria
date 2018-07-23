#pragma once

#include <cassert>
#include <random>
#include <algorithm>
#include <set>
#include <kdtree++/kdtree.hpp>
#include "SamplePartitioner.h"
#include "Sampler.h"
#include "load_balancing/IntersectionChecker.h"

namespace LoadBalancing
{ 
    template <uint D, typename Precision>
    struct NearestSamplePointAssigningSamplePartitioner : public SamplePartitioner<D, Precision>
    {
        NearestSamplePointAssigningSamplePartitioner(size_t partitionSize, Sampler<D, Precision> sampler, IntersectionPartitionMakerFunction<D, Precision> intersectionCheckerMaker)
            : mPartitionSize(partitionSize), mSampler(std::move(sampler)), mMakePartition(std::move(intersectionCheckerMaker))
        { }
        PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds) override
        {
            mSampling = mSampler(bounds, points, pointIds, mPartitionSize);
            auto kdtree = buildKdTree(mSampling.partition, mSampling.points);
            auto partitioning = makePartitioning(kdtree, points, mPartitionSize, pointIds);
            
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
            return "nearest-sample-point-assigning sample partitioner";
        }
        
        const Sampling<D, Precision>& sampling() const {
            return mSampling;
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
        
        using Tree = KDTree::KDTree<D, AssignedVector>;
        Tree buildKdTree(const std::vector<int>& partitioning, const std::vector<dVector<D, Precision>>& samplePoints);
        std::vector<IntersectionPartition<D, Precision>> makePartitioning(const Tree& tree,
															  const dPoints<D, Precision>& points,
															  size_t partitions,
															  const Point_Ids& pointIds);
         
        Sampling<D, Precision> mSampling;
        size_t mPartitionSize;
        Sampler<D, Precision> mSampler;
		IntersectionPartitionMakerFunction<D, Precision> mMakePartition;
    };
    
    template <uint D, typename Precision>
    typename NearestSamplePointAssigningSamplePartitioner<D, Precision>::Tree
   NearestSamplePointAssigningSamplePartitioner<D, Precision>::buildKdTree(
        const std::vector<int>& partitioning, const std::vector<dVector<D, Precision>>& samplePoints) {

        VTUNE_TASK(BuildTree);
        Tree tree;
        for(size_t i = 0; i < partitioning.size(); ++i) {
            int partition = partitioning[i];
            tree.insert({samplePoints[i], partition});
        }
        return tree;
    }
    
    template <uint D, typename Precision>
    std::vector<IntersectionPartition<D, Precision>>
    NearestSamplePointAssigningSamplePartitioner<D, Precision>::makePartitioning(const Tree& tree,
																				 const dPoints<D, Precision>& points,
																				 size_t partitions,
																				 const Point_Ids& pointIds) {
        VTUNE_TASK(AssignPointsByTree);

		return mMakePartition(points, pointIds, partitions, [&] (auto id) -> size_t {
			const auto& coords = points[id].coords;
			auto it_distance = tree.find_nearest(AssignedVector{coords, -1});
			assert(it_distance.first != tree.end());
			const auto& assignedVector = *(it_distance.first);
			return assignedVector.partitionId;
		});
    }   
}

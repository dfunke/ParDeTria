#pragma once

#include <cassert>
#include <random>
#include <algorithm>
#include <set>
#include <nanoflann.hpp>
#include "load_balancing/KDTreeVectorOfVectorsAdaptor.h"
#include "SamplePartitioner.h"
#include "Sampler.h"
#include "load_balancing/IntersectionChecker.h"

namespace LoadBalancing
{ 
    template <uint D, typename Precision, typename MonitorT>
    struct NearestSamplePointAssigningSamplePartitioner : SamplePartitioner<D, Precision, MonitorT>
    {
        NearestSamplePointAssigningSamplePartitioner(size_t partitionSize, Sampler<D, Precision> sampler, IntersectionPartitionMakerFunction<D, Precision> intersectionCheckerMaker)
            : mPartitionSize(partitionSize), mSampler(std::move(sampler)), mMakePartition(std::move(intersectionCheckerMaker))
        { }
        PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds,
                                        MonitorT& monitor) override
        {
            mSampling = mSampler(bounds, points, pointIds, mPartitionSize);
			monitor.registerSampleTriangulation(mSampling.size(), "s0");

            auto kdtree = buildKdTree(mSampling.partition, mSampling.points);
            auto partitioning = makePartitioning(kdtree, points, mPartitionSize, pointIds,
                                                 mSampling.partition);
            
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
        
        using Tree = KDTreeVectorOfVectorsAdaptor<std::vector<dVector<D, Precision>>, Precision, D>;
        Tree buildKdTree(const std::vector<int>& partitioning, const std::vector<dVector<D, Precision>>& samplePoints);
        std::vector<IntersectionPartition<D, Precision>> makePartitioning(const Tree& tree,
															  const dPoints<D, Precision>& points,
															  size_t partitions,
															  const Point_Ids& pointIds,
															  const std::vector<int>& partitioning);
         
        Sampling<D, Precision> mSampling;
        size_t mPartitionSize;
        Sampler<D, Precision> mSampler;
		IntersectionPartitionMakerFunction<D, Precision> mMakePartition;
    };
    
    template <uint D, typename Precision, typename MonitorT>
    typename NearestSamplePointAssigningSamplePartitioner<D, Precision, MonitorT>::Tree
   NearestSamplePointAssigningSamplePartitioner<D, Precision, MonitorT>::buildKdTree(
        const std::vector<int>& /*partitioning*/,
        const std::vector<dVector<D, Precision>>& samplePoints) {

        VTUNE_TASK(PartitionerBuildTree);
        Tree tree(D, samplePoints);
        return tree;
    }
    
    template <uint D, typename Precision, typename MonitorT>
    std::vector<IntersectionPartition<D, Precision>>
    NearestSamplePointAssigningSamplePartitioner<D, Precision, MonitorT>::
    makePartitioning(const Tree& tree,
                     const dPoints<D, Precision>& points,
                     size_t partitions,
                     const Point_Ids& pointIds,
                     const std::vector<int>& partitioning) {
        VTUNE_TASK(PartitionerAssignPointsByTree);

		return mMakePartition(points, pointIds, partitions, [&] (auto id) -> size_t {
			const auto& coords = points[id].coords;
			tIdType index;
			Precision distSquared;
			nanoflann::KNNResultSet<Precision> resultSet(1);
			resultSet.init(&index, &distSquared);
			tree.index->findNeighbors(resultSet, coords.data(), nanoflann::SearchParams());
			return partitioning[index];
		});
    }   
}

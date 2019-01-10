#pragma once

#include <cassert>
#include <random>
#include <algorithm>
#include <set>
#include <kdtree++/kdtree.hpp>
#include "SamplePartitioner.h"
#include "Sampler.h"
#include "load_balancing/IntersectionChecker.h"
#include "NearestSamplePointAssigningSamplePartitioner.h"

namespace LoadBalancing
{ 
    template <uint D, typename Precision, typename MonitorT>
    struct NearestSamplePointAssigningSampleBipartitioner : SamplePartitioner<D, Precision, MonitorT>
    {
        NearestSamplePointAssigningSampleBipartitioner(size_t partitionSize, Sampler<D, Precision> sampler, size_t baseCutoff, IntersectionPartitionMakerFunction<D, Precision> intersectionCheckerMaker)
            : mPartitionSize(partitionSize), mSampler(std::move(sampler)),
			mBaseCutoff(baseCutoff),
			mMakePartition(std::move(intersectionCheckerMaker))
        {
	        assert(partitionSize > 0);
	        assert((partitionSize & (partitionSize - 1)) == 0); // power of two
        }

        PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds,
                                        MonitorT& monitor) override
        {
            return makePartitioning(bounds, points, pointIds, mPartitionSize, "s0", monitor);
        }
        
        std::string info() const override
        {
            return "nearest-sample-point-assigning sample bipartitioner";
        }
        
        const Sampling<D, Precision>& sampling() const {
            return mSampling;
        }
        
    private:
		PartitionTree<D, Precision> makePartitioning(const dBox<D, Precision>& bounds,
		                                             dPoints<D, Precision>& points,
		                                             const Point_Ids& pointIds,
		                                             size_t numPartitions,
		                                             std::string provenance,
		                                             MonitorT& monitor) {
			PartitionTree<D, Precision> result;
			if(numPartitions > 1 && pointIds.size() > mBaseCutoff) {
				using BasePartitioner
					=	NearestSamplePointAssigningSamplePartitioner<D, Precision, UnpluggedMonitor>;
				BasePartitioner basePartitioner(2, mSampler, mMakePartition);

				UnpluggedMonitor unpluggedMonitor;
				auto tree = basePartitioner.partition(bounds, points, pointIds, unpluggedMonitor);
				mSampling = basePartitioner.sampling().clone();
				monitor.registerSampleTriangulation(mSampling.size(), provenance);
				monitor.registerSampling(mSampling, provenance);
				
				assert(!tree.isLeaf());
				auto& children =
					std::get<typename decltype(tree)::ChildContainer>(tree.attachment);
				assert(children.size() == 2);
				assert(children[0].isLeaf() && children[1].isLeaf());
				auto leftTree = makePartitioning(children[0].intersectionChecker->bounds(),
				                                 points,
				                                 std::move(std::get<Point_Ids>(children[0].attachment
																			   )),
				                                 numPartitions / 2,
				                                 provenance + "0",
				                                 monitor);
				auto rightTree = makePartitioning(children[1].intersectionChecker->bounds(),
				                                 points,
				                                 std::move(std::get<Point_Ids>(children[1].attachment
																			   )),
				                                 numPartitions / 2,
				                                 provenance + "1",
				                                 monitor);

				leftTree.intersectionChecker = std::move(children[0].intersectionChecker);
				rightTree.intersectionChecker = std::move(children[1].intersectionChecker);
				result.attachment = typename decltype(tree)::ChildContainer {
					std::move(leftTree), std::move(rightTree)};
			} else {
				result.intersectionChecker =
					std::make_unique<BoundsIntersectionChecker<D, Precision>>(bounds);
				result.attachment = pointIds;
			}
			return result;
		}
         
        size_t mPartitionSize;
        Sampler<D, Precision> mSampler;
        size_t mBaseCutoff;
		IntersectionPartitionMakerFunction<D, Precision> mMakePartition;
        Sampling<D, Precision> mSampling;
    };
}


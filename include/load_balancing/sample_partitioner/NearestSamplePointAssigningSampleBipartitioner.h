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
    template <uint D, typename Precision>
    struct NearestSamplePointAssigningSampleBipartitioner : public SamplePartitioner<D, Precision>
    {
        NearestSamplePointAssigningSampleBipartitioner(size_t partitionSize, Sampler<D, Precision> sampler, IntersectionPartitionMakerFunction<D, Precision> intersectionCheckerMaker)
            : mPartitionSize(partitionSize),
			mBasePartitioner(2, std::move(sampler), std::move(intersectionCheckerMaker))
        {
	        assert(partitionSize > 0);
	        assert((partitionSize & (partitionSize - 1)) == 0); // power of two
        }
        PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        const dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds) override
        {
            return makePartitioning(bounds, points, pointIds, mPartitionSize);
        }
        
        std::string info() const override
        {
            return "nearest-sample-point-assigning sample bipartitioner";
        }
        
        const Sampling<D, Precision>& sampling() const {
            return mBasePartitioner.sampling();
        }
        
    private:
		PartitionTree<D, Precision> makePartitioning(const dBox<D, Precision>& bounds,
		                                             const dPoints<D, Precision>& points,
		                                             const Point_Ids& pointIds,
		                                             size_t numPartitions) {
			PartitionTree<D, Precision> result;
			if(numPartitions > 1) {
				auto tree = mBasePartitioner.partition(bounds, points, pointIds);
				assert(!tree.isLeaf());
				const auto& children =
					std::get<typename decltype(tree)::ChildContainer>(tree.attachment);
				assert(children.size() == 2);
				assert(children[0].isLeaf() && children[1].isLeaf());
				auto leftTree = makePartitioning(children[0].intersectionChecker->bounds(),
				                                 points,
				                                 std::move(std::get<Point_Ids>(children[0].attachment
																			   )),
				                                 numPartitions / 2);
				auto rightTree = makePartitioning(children[1].intersectionChecker->bounds(),
				                                 points,
				                                 std::move(std::get<Point_Ids>(children[1].attachment
																			   )),
				                                 numPartitions / 2);
				result.intersectionChecker =
					std::make_unique<BoundsIntersectionChecker<D, Precision>>(bounds);
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
		NearestSamplePointAssigningSamplePartitioner<D, Precision> mBasePartitioner;
    };
}


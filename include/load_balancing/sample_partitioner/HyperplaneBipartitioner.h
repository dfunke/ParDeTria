#pragma once

#include <cassert>
#include <random>
#include <algorithm>
#include <set>
#include <Eigen/Eigen>
#include <kdtree++/kdtree.hpp>
#include "SamplePartitioner.h"
#include "Sampler.h"
#include "load_balancing/IntersectionChecker.h"
#include "NearestSamplePointAssigningSamplePartitioner.h"
#include "load_balancing/VectorOperations.h"
#include "load_balancing/HyperplaneIntersectionChecker.h"

namespace LoadBalancing
{ 
    template <uint D, typename Precision>
    struct HyperplaneSamplingBipartitioner : public SamplePartitioner<D, Precision>
    {
        HyperplaneSamplingBipartitioner(size_t partitionSize, Sampler<D, Precision> sampler)
            : mPartitionSize(partitionSize), mSampler(std::move(sampler))
        {
	        assert(partitionSize > 0);
	        assert((partitionSize & (partitionSize - 1)) == 0); // power of two
        }

        PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        const dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds) override
        {
			auto is = std::make_unique<BoundsIntersectionChecker<D, Precision>>(bounds);
            return makePartitioning(bounds, points, pointIds, std::move(is), mPartitionSize);
        }
        
        std::string info() const override
        {
            return "hyperplane sampling bipartitioner";
        }
        
        const Sampling<D, Precision>& sampling() const {
            return mSampling;
        }
        
    private:
		PartitionTree<D, Precision>
			makePartitioning(const dBox<D, Precision>& bounds,
			                 const dPoints<D, Precision>& points,
			                 const Point_Ids& pointIds,
			                 std::unique_ptr<IntersectionChecker<D, Precision>> intersectionChecker,
			                 size_t numPartitions) {

			PartitionTree<D, Precision> result;
			if(numPartitions > 1) {
           		mSampling = mSampler(bounds, points, pointIds, mPartitionSize);

				// calc split points ...
				std::vector<dVector<D, Precision>> splitPoints;
				const auto& records = mSampling.graph.nodeRecords;
				const auto& adjacency = mSampling.graph.adjacency;
				for(size_t a = 0; a + 1 < records.size(); ++a) {
					for(int i = records[a]; i < records[a+1]; ++i) {
						size_t b = adjacency[i];
						if(mSampling.partition[a] != mSampling.partition[b]) {
							splitPoints.push_back(0.5f * (mSampling.points[a] + mSampling.points[b]));
						}
					}
				}

				// make hyperplane
				const auto m = splitPoints.size();
				constexpr auto n = D + 1;

				using Vector = Eigen::Matrix<Precision, Eigen::Dynamic, 1>;
				using Matrix = Eigen::Matrix<Precision, Eigen::Dynamic, Eigen::Dynamic>;
				Vector zero(m);
				Matrix matrix(m, n);
				for(size_t i = 0; i < m; ++i) {
					zero(i) = 0;
				}
				for(size_t i = 0; i < m; ++i) {
					for(size_t j = 0; j < D; ++j) {
						matrix(i, j) = splitPoints[i][j];
					}
					matrix(i, D) = 1;
				}

				Eigen::JacobiSVD<Matrix> svd(matrix, Eigen::ComputeFullV);
				Vector vector = svd.matrixV().col(D);

				dVector<D, Precision> orthogonal;
				for(size_t j = 0; j < D; ++j) {
					orthogonal[j] = vector(j);
				}
				
				auto normal = 1/len(orthogonal) * orthogonal;
				auto offset = -vector(D)/scalarProduct(normal, orthogonal);

				// separate points
				Point_Ids leftPoints;
				Point_Ids rightPoints;
				for(auto id : pointIds) {
					auto& coords = points[id].coords;
					if(distanceToPlane(coords, normal, offset) > 0) {
						leftPoints.insert(id);
					} else {
						rightPoints.insert(id);
					}
				}

				// build partition tree
				auto leftIC =
					std::make_unique<HyperplaneIntersectionChecker<D, Precision>>(bounds,
																				  normal,	
																				  offset + 1e-6);
				auto leftTree = makePartitioning(bounds, points, std::move(leftPoints),
				                                 std::move(leftIC), numPartitions / 2);
				
				auto rightIC =
					std::make_unique<HyperplaneIntersectionChecker<D, Precision>>(bounds,
																				  -normal,
																				  -offset + 1e-6);
				auto rightTree = makePartitioning(bounds, points, std::move(rightPoints),
				                                  std::move(rightIC), numPartitions / 2);
				result.intersectionChecker = std::move(intersectionChecker);
				result.attachment = typename decltype(result)::ChildContainer {
					std::move(leftTree), std::move(rightTree)};
			} else {
				result.intersectionChecker = std::move(intersectionChecker);
				result.attachment = pointIds;
			}
			return result;
		}
         
        Sampling<D, Precision> mSampling;
        size_t mPartitionSize;
        Sampler<D, Precision> mSampler;
    };
}


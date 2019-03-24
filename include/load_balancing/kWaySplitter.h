#pragma once

#include "Partitioner.h"
#include "PartitionTree.h"

namespace LoadBalancing
{ 
    template <uint D, typename Precision, typename MonitorT>
    struct KWaySplitter : public Partitioner<D, Precision, MonitorT>
    {
	    KWaySplitter(size_t numPartitions)
		    : splits(generateSplits(numPartitions))
        {
        }

        PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds,
                                        MonitorT& /*monitor*/) override
        {
			PartitionTree<D, Precision> tree;
			tree.attachment = partitionRecursively(points, pointIds, 0);
			tree.intersectionChecker =
				std::make_unique<BoundsIntersectionChecker<D, Precision>>(bounds);
			return tree;
        }
        
        std::string info() const override
        {
            return "k-way cycle";
        }

	private:
        using ChildContainer = typename PartitionTree<D, Precision>::ChildContainer;

		static std::vector<size_t> generateSplits(size_t numPartitions) {
	        assert(isPow2(numPartitions));

			size_t log2P = log2(numPartitions);
			std::vector<size_t> result(D, 1);
			for(size_t i = 0; i < log2P; ++i) {
				result[i%D] *= 2;
			}
			return result;
		}

		std::vector<size_t> splits;

		ChildContainer partitionRecursively(const dPoints<D, Precision>& points,
		                                            Point_Ids pointIds,
		                                            size_t depth) {
			ChildContainer result;
			auto stats = getPointStatsSeq(pointIds.begin(), pointIds.end(), points);
			dBox<D, Precision> bounds(stats.min, stats.max);

			if(depth < splits.size()) {
				size_t dim = depth % D;
				float min = stats.min[dim];
				float max = stats.max[dim];
				size_t numPartitions = splits[depth];

				std::vector<Point_Ids> partitioning(numPartitions);

				for(auto id : pointIds) {
					if(!dPoint<D, Precision>::isFinite(id))
						continue; // skip infinite points, they will be added later

					float normalizedAxisPosition = (points[id].coords[dim] - min)/(max - min);
					int partition = (size_t)std::floor(normalizedAxisPosition * numPartitions);
					partition = std::clamp(partition, 0, (int)numPartitions - 1);

					partitioning[partition].insert(id);
				}

				for(auto&& partition : partitioning) {
					auto subresult = partitionRecursively(points, std::move(partition), depth + 1);
					result.reserve(result.size() + subresult.size());
					std::move(std::begin(subresult), std::end(subresult),
					          std::back_inserter(result));
				}
			} else {
				PartitionTree<D, Precision> leave;
				for (tIdType k = dPoint<D, Precision>::cINF; k != 0; ++k) {
					pointIds.insert(k);
				}
				leave.attachment = std::move(pointIds);
				leave.intersectionChecker =
					std::make_unique<BoundsIntersectionChecker<D, Precision>>(bounds);
				result.push_back(std::move(leave));
			}
			return result;
		}
	};
}


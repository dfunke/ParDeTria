#pragma once
#include "Geometry.h"
#include "load_balancing/IntersectionChecker.h"
#include "load_balancing/BoxUtils.h"

namespace LoadBalancing
{
	template <uint D, typename Precision>
	struct BoundsIntersectionChecker : IntersectionChecker<D, Precision>
	{
		BoundsIntersectionChecker(dBox<D, Precision> bounds)
			: IntersectionChecker<D, Precision>(std::move(bounds)) {
		}

		virtual bool intersects(const dSphere<D, Precision>& sphere) const override
		{
			return IntersectionChecker<D, Precision>::bounds().intersects(sphere);
		}

		std::unique_ptr<IntersectionChecker<D, Precision>> copy() const override
		{
			return std::make_unique<BoundsIntersectionChecker>(*this);
		}
	};

	template <uint D, typename Precision>
	struct BoundsIntersectionPartitionMaker : IntersectionPartitionMaker<D, Precision>
	{
		template <typename PartitionAssigner>
		std::vector<IntersectionPartition<D, Precision>>
		operator()(const dPoints<D, Precision>& points,
		           const Point_Ids& ids, size_t partitions, PartitionAssigner part)
		{
			std::vector<Concurrent_Growing_Point_Ids> idsets;
			idsets.reserve(partitions);
			for(size_t i = 0; i < partitions; ++i) {
				idsets.emplace_back(ids.size()/partitions);
			}
			using Handle = tbb::enumerable_thread_specific<hConcurrent_Growing_Point_Ids,
			      tbb::cache_aligned_allocator<hConcurrent_Growing_Point_Ids>,
			      tbb::ets_key_usage_type::ets_key_per_instance>;
	                
			tbb::parallel_for(ids.range(), [&](const auto &range) {
			    for(auto id : range) {
					if(dPoint<D, Precision>::isFinite(id)) {
						auto partition = part(id);
						
						Handle tsPointHandle = std::ref(idsets[partition]);
						auto& pointsHandle = tsPointHandle.local();
						pointsHandle.insert(id);
					}
				}
			});
			
			std::vector<Point_Ids> idss;
			for(auto& ids : idsets) {
				idss.emplace_back(std::move(ids.data()));
			}
			std::vector<dBox<D, Precision>> boundss(partitions);
			std::vector<bool> initialized(partitions);
			
			tbb::parallel_for(size_t(0), partitions, [&](auto partition) {
			    for(auto id : idss[partition]) {
					const auto& coords = points[id].coords;
					if(initialized[partition]) {
						enlargeBoxAroundVector<D, Precision>(boundss[partition], coords);
					} else {
						boundss[partition].low = coords;
						boundss[partition].high = boundss[partition].low;
						initialized[partition] = true;
					}
				}
			});
			
			for (tIdType k = dPoint<D, Precision>::cINF; k != 0; ++k) {
				for(auto& e : idss) {
					e.insert(k);
				}
			}

			std::vector<IntersectionPartition<D, Precision>> result(idss.size());
			for(size_t i = 0; i < result.size(); ++i) {
				result[i].pointIds = std::move(idss[i]);
				using Checker = BoundsIntersectionChecker<D, Precision>;
				result[i].intersectionChecker = std::make_unique<Checker>(std::move(boundss[i]));
			}
			
			return result;
		}
	};
}

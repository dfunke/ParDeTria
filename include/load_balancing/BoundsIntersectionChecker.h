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
	};

	template <uint D, typename Precision>
	struct BoundsIntersectionPartitionMaker : IntersectionPartitionMaker<D, Precision>
	{
		template <typename PartitionAssigner>
		std::vector<IntersectionPartition<D, Precision>> operator()(const dPoints<D, Precision>& points, const Point_Ids& ids, PartitionAssigner part)
		{
			std::vector<Point_Ids> idss;
			std::vector<dBox<D, Precision>> boundss;
			std::vector<bool> initialized;
			for(auto id : ids) {
				if(dPoint<D, Precision>::isFinite(id)) {
					auto partition = part(id);
					
					if(idss.size() <= partition) {
						idss.resize(partition + 1);
						boundss.resize(partition + 1);
						initialized.resize(partition + 1, false);
					}
					idss[partition].insert(id);
                	const auto& coords = points[id].coords;
					if(initialized[partition]) {
						enlargeBoxAroundVector<D, Precision>(boundss[partition], coords);
					} else {
						boundss[partition].low = coords;
						boundss[partition].high = boundss[partition].low;
						initialized[partition] = true;
					}
				}
			}
			
			for (tIdType k = dPoint<D, Precision>::cINF; k != 0; ++k) {
				for(auto& e : idss) {
					e.insert(k);
				}
			}

			std::vector<IntersectionPartition<D, Precision>> result(idss.size());
			for(size_t i = 0; i < result.size(); ++i) {
				result[i].pointIds = std::move(idss[i]);
				result[i].intersectionChecker = std::make_unique<BoundsIntersectionChecker<D, Precision>>(std::move(boundss[i]));
			}
			
			return result;
		}
	};
}

#pragma once
#include <algorithm>
#include <unordered_set>
#include <functional>
#include "load_balancing/IntersectionChecker.h"
#include "load_balancing/Grid.h"

namespace LoadBalancing
{
	template <uint D, typename Precision, typename IndexPrecision>
	struct GridIntersectionChecker : IntersectionChecker<D, Precision>
	{
		template <typename ForwardIt>
		GridIntersectionChecker(dBox<D, Precision> bounds, Grid<D, Precision, IndexPrecision> grid, ForwardIt cellsBegin, ForwardIt cellsEnd)
			: IntersectionChecker<D, Precision>(std::move(bounds)), grid(std::move(grid)), cells(cellsBegin, cellsEnd)
		{}
		
		virtual bool intersects(const dSphere<D, Precision>& sphere) const override
		{
			auto discreteSphere = grid.intersectingIndices(sphere);
			return std::any_of(discreteSphere.begin(), discreteSphere.end(), [this] (const auto& i) {
			                   return cells.count(i) == 1;
							   });
		}

	private:
		Grid<D, Precision, IndexPrecision> grid;
		std::unordered_set<dIndex<D, IndexPrecision>, dIndexHasher<D, IndexPrecision>> cells;
	};
	
	template <uint D, typename Precision, typename IndexPrecision = int64_t>
	struct GridIntersectionPartitionMaker : IntersectionPartitionMaker<D, Precision>
	{
		GridIntersectionPartitionMaker(Grid<D, Precision, IndexPrecision> grid)
			: grid(std::move(grid)) {
		}

		template <typename PartitionAssigner>
		std::vector<IntersectionPartition<D, Precision>> operator()(const dPoints<D, Precision>& points, const Point_Ids& ids, PartitionAssigner part)
		{
			std::vector<Point_Ids> idss;
			std::vector<dBox<D, Precision>> boundss;
			std::vector<bool> initialized;
			std::unordered_set<dIndex<D, IndexPrecision>, dIndexHasher<D, IndexPrecision>> indices;

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
					indices.insert(grid.indexAt(coords));
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
				result[i].intersectionChecker = std::make_unique<GridIntersectionChecker<D, Precision, IndexPrecision>>(std::move(boundss[i]),
																										grid, indices.begin(),
																										indices.end());
			}
			
			return result;
		}
		
	private:
		Grid<D, Precision, IndexPrecision> grid;
	};
}


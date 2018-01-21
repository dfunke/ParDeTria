#pragma once
#include <algorithm>
#include <unordered_set>
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
			auto discreteSphere = grid.intersects(sphere);
			return std::any_of(discreteSphere.begin(), discreteSphere.end(), [this] (const auto& i) {
			                   return cells.count(i) == 1;
							   });
		}

	private:
		Grid<D, Precision, IndexPrecision> grid;
		std::unordered_set<dIndex<D, IndexPrecision>> cells;
	};
}

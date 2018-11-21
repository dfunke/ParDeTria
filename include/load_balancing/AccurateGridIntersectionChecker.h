#pragma once
#include <algorithm>
#include <unordered_set>
#include <functional>
#include "load_balancing/IntersectionChecker.h"
#include "load_balancing/Grid.h"
#include "datastructures/AABBTree.h"
#include "load_balancing/GridIntersectionChecker.h"


namespace LoadBalancing
{
	template <uint D, typename Precision, typename IndexPrecision = int64_t>
	struct AccurateGridIntersectionChecker : GridIntersectionChecker<D, Precision>
	{
		struct CellData {
			std::vector<tIdType> pointIds;
		};
		
		using CellsWithPoints = std::unordered_map<dIndex<D>,
			CellData, dIndexHasher<D, int64_t>>;

		private:
		const dPoints<D, Precision>* points;
		std::shared_ptr<CellsWithPoints> mCells;


		public:
		template <typename ForwardIt>
		AccurateGridIntersectionChecker(dBox<D, Precision> bounds,
		                                Grid<D, Precision, IndexPrecision> grid,
		                                ForwardIt cellsBegin,
		                                ForwardIt cellsEnd)
			: GridIntersectionChecker<D, Precision>(std::move(bounds),
			                                        std::move(grid),
			                                        std::move(cellsBegin),
			                                        std::move(cellsEnd))
		{
		}

		AccurateGridIntersectionChecker(GridIntersectionChecker<D, Precision> gic,
		                                std::shared_ptr<CellsWithPoints> cells,
		                                const dPoints<D, Precision>& points)
			: GridIntersectionChecker<D, Precision>(std::move(gic)),
			  points(&points),
			  mCells(std::move(cells)) { }
		
		virtual bool intersects(const dSphere<D, Precision>& sphere) const override
		{
			auto cells = GridIntersectionChecker<D, Precision>::intersectingCells(sphere);
			for(const auto& cell : cells) {
				auto it = mCells->find(cell);
				assert(it != mCells->end());
				for(tIdType pointId : it->second.pointIds) {
					Precision distSquared = lenSquared(sphere.center - (*points)[pointId].coords);
					if(distSquared < sphere.radius * sphere.radius)
						return true;
				}
			}
			
			return false;
		}

		std::unique_ptr<IntersectionChecker<D, Precision>> copy() const override
		{
			return std::make_unique<AccurateGridIntersectionChecker>(*this);
		}
	};
	
	template <uint D, typename Precision, typename IndexPrecision = int64_t>
	struct AccurateGridIntersectionPartitionMaker : GridIntersectionPartitionMaker<D, Precision>
	{
		AccurateGridIntersectionPartitionMaker(Grid<D, Precision, IndexPrecision> mGrid)
			: GridIntersectionPartitionMaker<D, Precision>(std::move(mGrid)) {
		}

		template <typename PartitionAssigner>
		std::vector<IntersectionPartition<D, Precision>>
		operator()(const dPoints<D, Precision>& points,
		           const Point_Ids& ids, size_t partitions, PartitionAssigner part)
		{

		    VTUNE_TASK(AssignPointsAccurately);
		   
		   	using GIC = GridIntersectionChecker<D, Precision>;
			using AGIC = AccurateGridIntersectionChecker<D, Precision>;

			std::vector<IntersectionPartition<D, Precision>> gridIntersectionPartitioning =
				GridIntersectionPartitionMaker<D, Precision>::operator()(points, ids, partitions,
																		  std::move(part));
			using Cells = typename AGIC::CellsWithPoints;
			Cells cells;
			for(tIdType id : ids) {
				auto index = GridIntersectionPartitionMaker<D, Precision>::mGrid.indexAt(points[id].coords);
				auto it = cells.find(index);
				if(it != cells.end()) {
					it->second.pointIds.push_back(id);
				} else {
					cells[index] = typename AGIC::CellData{{id}};
				}
			}

			auto cellsPtr = std::make_shared<decltype(cells)>(std::move(cells));
			
			std::vector<IntersectionPartition<D, Precision>> result;
			for(auto& partition : gridIntersectionPartitioning) {
				auto intersectionChecker = std::move(partition.intersectionChecker);
				auto gridIntersectionChecker =
					*(dynamic_cast<GIC*>(intersectionChecker.release()));
				partition.intersectionChecker =
					std::make_unique<AGIC>(std::move(gridIntersectionChecker), cellsPtr, points);
			}

			return gridIntersectionPartitioning;
		}
	};
}


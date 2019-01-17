#pragma once
#include <cassert>
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
		CellsWithPoints mCells;


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
		                                CellsWithPoints cells,
		                                const dPoints<D, Precision>& points)
			: GridIntersectionChecker<D, Precision>(std::move(gic)),
			  points(&points),
			  mCells(std::move(cells)) { }
		
		virtual bool intersects(const dSphere<D, Precision>& sphere,
		                        const dSimplex<D, Precision>& simplex) const override
		{
			auto cells = GridIntersectionChecker<D, Precision>::intersectingCells(sphere);
			for(const auto& cell : cells) {
				auto it = mCells.find(cell);
				assert(it != mCells.end());
				for(tIdType pointId : it->second.pointIds) {
					if(checkCircumsphereContainsPoint((*points)[pointId], sphere, simplex))
						return true;
				}
			}
			return false;
		}

		std::unique_ptr<IntersectionChecker<D, Precision>> copy() const override
		{
			return std::make_unique<AccurateGridIntersectionChecker>(*this);
		}

	private:
		bool checkCircumsphereContainsPoint(const dPoint<2, Precision>& point,
		                                    const dSphere<2, Precision>&,
		                                    const dSimplex<2, Precision>& simplex) const {
			return GeometryCore<2, Precision>::inSphere(point,
			                                            (*points)[simplex.vertices[0]],
			                                            (*points)[simplex.vertices[1]],
			                                            (*points)[simplex.vertices[2]]);
		}
		
		bool checkCircumsphereContainsPoint(const dPoint<3, Precision>& point,
		                                    const dSphere<3, Precision>&,
		                                    const dSimplex<3, Precision>& simplex) const {
			return GeometryCore<3, Precision>::inSphere(point,
			                                            (*points)[simplex.vertices[0]],
			                                            (*points)[simplex.vertices[1]],
			                                            (*points)[simplex.vertices[2]],
			                                            (*points)[simplex.vertices[3]]);
		}
		template <uint D2>
		bool checkCircumsphereContainsPoint(const dPoint<D2, Precision>& point,
		                                    const dSphere<D2, Precision>& sphere,
		                                    const dSimplex<D2, Precision>&) const {
			Precision distSquared =	lenSquared(sphere.center - point.coords);
			return distSquared < sphere.radius * sphere.radius;
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
			std::vector<IntersectionPartition<D, Precision>> result;
			for(auto& intersectionPartition : gridIntersectionPartitioning) {
				using Cells = typename AGIC::CellsWithPoints;
				Cells cells;
				for(tIdType id : intersectionPartition.pointIds) {
					auto index = GridIntersectionPartitionMaker<D, Precision>::mGrid
						.indexAt(points[id].coords);
					auto it = cells.find(index);
					if(it != cells.end()) {
						it->second.pointIds.push_back(id);
					} else {
						cells[index] = typename AGIC::CellData{{id}};
					}
				}

				auto intersectionChecker = std::move(intersectionPartition.intersectionChecker);
				auto gridIntersectionChecker =
					*(dynamic_cast<GIC*>(intersectionChecker.release()));
				intersectionPartition.intersectionChecker =
					std::make_unique<AGIC>(std::move(gridIntersectionChecker), std::move(cells), points);
			}

			return gridIntersectionPartitioning;
		}
	};
}


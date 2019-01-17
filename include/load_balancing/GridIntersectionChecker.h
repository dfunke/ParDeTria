#pragma once
#include <algorithm>
#include <unordered_set>
#include <functional>
#include "load_balancing/IntersectionChecker.h"
#include "load_balancing/Grid.h"
#include "datastructures/AABBTree.h"

namespace LoadBalancing
{
	template <uint D, typename Precision, typename IndexPrecision = int64_t>
	struct GridIntersectionChecker : IntersectionChecker<D, Precision>
	{
		template <typename ForwardIt>
		GridIntersectionChecker(dBox<D, Precision> bounds, Grid<D, Precision, IndexPrecision> grid, ForwardIt cellsBegin, ForwardIt cellsEnd)
			: IntersectionChecker<D, Precision>(std::move(bounds)), mGrid(std::move(grid)), cells(cellsBegin, cellsEnd), aabb(bounds, grid.cellWidth())
		{

//			printf("\n creating %.0f-%.0f %.0f-%.0f %.0f-%.0f\n",bounds.low[0], bounds.high[0], bounds.low[1], bounds.high[1], bounds.low[2], bounds.high[2]);
            VTUNE_TASK(GridIntersectionChecker_CTOR);
            for(uint i = 0; i < cells.size(); ++i) {
				const auto & c = cells[i];
				auto b = grid.boundsOf(c);
				aabb.insert(b);

//				printf("%u: %lu %lu %lu - %.0f-%.0f %.0f-%.0f %.0f-%.0f\n", i, c[0], c[1], c[2], b.low[0], b.high[0], b.low[1], b.high[1], b.low[2], b.high[2]);
				//RASSERT(c[0] != 32 && c[1] != 31 && c[2] != 43);

//				if(!aabb.insert(b)) {
//					printf("\tfailed\n");
//					RASSERT(false);
//				}
			}

		}
		
		virtual bool intersects(const dSphere<D, Precision>& sphere,
		                        const dSimplex<D, Precision>&) const override
		{
			ASSERT(aabb.intersects(sphere) == std::any_of(cells.begin(), cells.end(), [this, &sphere](const auto& i) {
				return intersectsWith<D, Precision, IndexPrecision>(mGrid, sphere, i);
			}));


			return aabb.intersects(sphere);
		}

		using Intersection = typename AABBTree<D, Precision>::Intersection;
		Intersection intersectingCells(const dSphere<D, Precision>& sphere) const {
			return aabb.intersectingCells(sphere);
		}

		std::unique_ptr<IntersectionChecker<D, Precision>> copy() const override
		{
			return std::make_unique<GridIntersectionChecker>(*this);
		}
		
		Grid<D, Precision, IndexPrecision> grid() const {
			return mGrid;
		}
		
		auto cellsBegin() const {
			return cells.cbegin();
		}
		
		auto cellsEnd() const {
			return cells.cend();
		}

	private:
		Grid<D, Precision, IndexPrecision> mGrid;
		std::vector<dIndex<D, IndexPrecision>> cells;
		AABBTree<D, Precision> aabb;
	};
	
	template <uint D, typename Precision, typename IndexPrecision = int64_t>
	struct GridIntersectionPartitionMaker : IntersectionPartitionMaker<D, Precision>
	{
		GridIntersectionPartitionMaker(Grid<D, Precision, IndexPrecision> mGrid)
			: mGrid(std::move(mGrid)) {
		}

		template <typename PartitionAssigner>
		std::vector<IntersectionPartition<D, Precision>>
		operator()(const dPoints<D, Precision>& points,
		           const Point_Ids& ids, size_t partitions, PartitionAssigner part)
		{

		    VTUNE_TASK(AssignPoints);

			std::vector<Concurrent_Growing_Point_Ids> idsets;
			idsets.reserve(partitions);

			using Handle  =tbb::enumerable_thread_specific<hConcurrent_Growing_Point_Ids,
					tbb::cache_aligned_allocator<hConcurrent_Growing_Point_Ids>,
					tbb::ets_key_usage_type::ets_key_per_instance>;

			std::vector<Handle*> tsPointHandles;
			tsPointHandles.resize(partitions);

			for(size_t i = 0; i < partitions; ++i) {
				idsets.emplace_back(ids.size()/partitions);
				tsPointHandles[i] = new Handle(std::ref(idsets[i]));
			}

	                
			tbb::parallel_for(ids.range(), [&](const auto &range) {
			    for(auto id : range) {
					if(dPoint<D, Precision>::isFinite(id)) {
						auto partition = part(id);

						auto& pointsHandle = tsPointHandles[partition]->local();
						pointsHandle.insert(id);
					}
				}
			});

			for(size_t i = 0; i < partitions; ++i) {
				delete tsPointHandles[i];
			}

			std::vector<Point_Ids> idss;
			for(auto& ids : idsets) {
				idss.emplace_back(std::move(ids.data()));
			}
			std::vector<dBox<D, Precision>> boundss(partitions);
			std::vector<bool> initialized(partitions);
			std::vector<std::unordered_set<dIndex<D, IndexPrecision>,
				dIndexHasher<D, IndexPrecision>>> indices(partitions);
			
			tbb::parallel_for(size_t(0), partitions, [&](auto partition) {
			    for(auto id : idss[partition]) {
					const auto& coords = points[id].coords;
					indices[partition].insert(mGrid.indexAt(coords));
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
			tbb::parallel_for(size_t(0), partitions, [&](auto i) {
				result[i].pointIds = std::move(idss[i]);
				using Checker = GridIntersectionChecker<D, Precision, IndexPrecision>;
				result[i].intersectionChecker = std::make_unique<Checker>(std::move(boundss[i]), mGrid,
																   indices[i].begin(),
																   indices[i].end());
			});
			
			return result;
		}
		
	protected:
		Grid<D, Precision, IndexPrecision> mGrid;
	};
}


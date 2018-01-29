#pragma once
#include <cassert>
#include <cstdint>
#include <vector>
#include <array>
#include <type_traits>
#include <algorithm>
#include "Geometry.h"
#include "load_balancing/VectorOperations.h"
#include "load_balancing/dIndex.h"

namespace LoadBalancing
{

template <uint D, typename Precision, typename IndexPrecision = int64_t>
struct Grid
{
	Grid(Precision cellWidth);
	Precision cellWidth() const;
	dIndex<D, IndexPrecision> indexAt(const dVector<D, Precision>& v) const;
	dVector<D, Precision> centerOf(const dIndex<D, IndexPrecision>& i) const;
	dBox<D, Precision> boundsOf(const dIndex<D, IndexPrecision>& i) const;

	bool next(dIndex<D, IndexPrecision>& gridIter, IndexPrecision end) const;
	Precision diagonalCellLength() const;

private:
	Precision mCellWidth;

	IndexPrecision fitOneDimension(Precision x) const;
};

template <uint D, typename Precision, typename IndexPrecision>
std::vector<dIndex<D, IndexPrecision>> intersectingIndices(const dSphere<D, Precision>& sphere, const Grid<D, Precision, IndexPrecision>& grid);


template <uint D, typename Precision, typename IndexPrecision>
bool intersectsWith(const dSphere<D, Precision>& sphere, const dIndex<D, IndexPrecision>& i);

template <uint D, typename Precision, typename IndexPrecision>
bool intersectsWith(const Grid<D, Precision, IndexPrecision>& grid, const dSphere<D, Precision>& sphere, const dIndex<D, IndexPrecision>& i)
{
	auto cell = i;
	auto center = grid.indexAt(sphere.center);
	IndexPrecision discreteRadius = std::ceil(sphere.radius/grid.cellWidth());
	std::transform(cell.begin(), cell.end(), center.begin(), cell.begin(), std::minus<IndexPrecision>());
	bool result = true;
	for(const auto& c : cell)
		if(std::abs(c) > discreteRadius)
			result = false;
	
	return result;
}

template <uint D, typename Precision, typename IndexPrecision>
Grid<D, Precision, IndexPrecision>::Grid(Precision cellWidth)
	: mCellWidth(cellWidth)
{
	assert(mCellWidth > 0);
}

template <uint D, typename Precision, typename IndexPrecision>
Precision Grid<D, Precision, IndexPrecision>::cellWidth() const
{
	return mCellWidth;
}

template <uint D, typename Precision, typename IndexPrecision>
dIndex<D, IndexPrecision> Grid<D, Precision, IndexPrecision>::indexAt(const dVector<D, Precision>& v) const
{
	dIndex<D, IndexPrecision> result;
	std::transform(begin(v), end(v), begin(result), [this] (auto x) -> auto { return fitOneDimension(x); });
	return result;
}

template <uint D, typename Precision, typename IndexPrecision>
dVector<D, Precision> Grid<D, Precision, IndexPrecision>::centerOf(const dIndex<D, IndexPrecision>& i) const
{
	dVector<D, Precision> result;
	std::transform(i.begin(), i.end(), result.begin(), [this] (auto x) { return x * mCellWidth; });
	return result;
}
	
template <uint D, typename Precision, typename IndexPrecision>
dBox<D, Precision> Grid<D, Precision, IndexPrecision>::boundsOf(const dIndex<D, IndexPrecision>& i) const
{
	dBox<D, Precision> result;
	auto center = centerOf(i);
	std::transform(center.begin(), center.end(), result.low.begin(), [this] (auto x) { return x - mCellWidth/2; });
	std::transform(center.begin(), center.end(), result.high.begin(), [this] (auto x) { return x + mCellWidth/2; });
	return result;
}

template <uint D, typename Precision, typename IndexPrecision>
std::vector<dIndex<D, IndexPrecision>> intersectingIndices(const dSphere<D, Precision>& sphere, const Grid<D, Precision, IndexPrecision>& grid)
{
	std::vector<dIndex<D, IndexPrecision>> result;
	
	IndexPrecision discreteRadius = std::ceil(sphere.radius/grid.cellWidth());
	auto centerIndex = grid.indexAt(sphere.center);	
	auto maxDist = grid.diagonalCellLength() / 2 + sphere.radius;

	dIndex<D, IndexPrecision> gridIter;
	gridIter.fill(-discreteRadius);
	do {
		dIndex<D, IndexPrecision> cell;
		std::transform(centerIndex.begin(), centerIndex.end(), gridIter.begin(), cell.begin(), std::plus<>());
		auto cellCenter = grid.centerOf(cell);

		if(lenSquared(cellCenter - sphere.center) <= maxDist * maxDist) {
			result.push_back(cell);
		}
		
	} while(grid.next(gridIter, discreteRadius));

	return result;
}

template <uint D, typename Precision, typename IndexPrecision>
bool Grid<D, Precision, IndexPrecision>::next(dIndex<D, IndexPrecision>& gridIter, IndexPrecision width) const
{
		auto it = std::find_if(gridIter.begin(), gridIter.end(), [width] (auto x) { return x < width; });
		std::fill(gridIter.begin(), it, -width);
		if(it != gridIter.end())
			++(*it);
		
		return it != gridIter.end();
}


template <typename Precision, typename IndexPrecision>
std::vector<dIndex<2, IndexPrecision>> intersectingIndices(const dSphere<2, Precision>& sphere, const Grid<2, Precision, IndexPrecision>& grid)
{
	std::vector<dIndex<2, IndexPrecision>> result;
	
	IndexPrecision discreteRadius = std::ceil(sphere.radius/grid.cellWidth());
	auto centerIndex = grid.indexAt(sphere.center);

	for(IndexPrecision x = centerIndex[0]-discreteRadius; x <= centerIndex[0]+discreteRadius; ++x) {
		for(IndexPrecision y = centerIndex[1]-discreteRadius; y <= centerIndex[1]+discreteRadius; ++y) {
			dIndex<2, IndexPrecision> cell{x, y};
			result.push_back(std::move(cell));
		}
	}
	return result;
}

template <uint D, typename Precision, typename IndexPrecision>
Precision Grid<D, Precision, IndexPrecision>::diagonalCellLength() const
{
	return sqrt(D) * mCellWidth;
}

template <uint D, typename Precision, typename IndexPrecision>
IndexPrecision Grid<D, Precision, IndexPrecision>::fitOneDimension(Precision x) const
{
	return std::floor(x/mCellWidth + 0.5);
}

}

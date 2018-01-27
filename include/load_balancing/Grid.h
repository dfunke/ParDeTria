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
	dIndex<D, IndexPrecision> indexAt(const dVector<D, Precision>& v) const;
	dVector<D, Precision> centerOf(const dIndex<D, IndexPrecision>& i) const;
	dBox<D, Precision> boundsOf(const dIndex<D, IndexPrecision>& i) const;
	std::vector<dIndex<D, IndexPrecision>> intersectingIndices(const dSphere<D, Precision>& sphere) const;

private:
	Precision cellWidth;

	IndexPrecision fitOneDimension(Precision x) const;
	bool next(dIndex<D, IndexPrecision>& gridIter, IndexPrecision end) const;
	Precision diagonalCellLength() const;
};


template <uint D, typename Precision, typename IndexPrecision>
Grid<D, Precision, IndexPrecision>::Grid(Precision cellWidth)
	: cellWidth(cellWidth)
{
	assert(cellWidth > 0);
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
	std::transform(i.begin(), i.end(), result.begin(), [this] (auto x) { return x * cellWidth; });
	return result;
}
	
template <uint D, typename Precision, typename IndexPrecision>
dBox<D, Precision> Grid<D, Precision, IndexPrecision>::boundsOf(const dIndex<D, IndexPrecision>& i) const
{
	dBox<D, Precision> result;
	auto center = centerOf(i);
	std::transform(center.begin(), center.end(), result.low.begin(), [this] (auto x) { return x - cellWidth/2; });
	std::transform(center.begin(), center.end(), result.high.begin(), [this] (auto x) { return x + cellWidth/2; });
	return result;
}

template <uint D, typename Precision, typename IndexPrecision>
std::vector<dIndex<D, IndexPrecision>> Grid<D, Precision, IndexPrecision>::intersectingIndices(const dSphere<D, Precision>& sphere) const
{
	std::vector<dIndex<D, IndexPrecision>> result;
	
	IndexPrecision discreteRadius = std::ceil(sphere.radius/cellWidth + 1);
	auto centerIndex = indexAt(sphere.center);	
	auto maxDist = diagonalCellLength() / 2 + sphere.radius;

	dIndex<D, IndexPrecision> gridIter;
	gridIter.fill(-discreteRadius);
	do {
		dIndex<D, IndexPrecision> cell;
		std::transform(centerIndex.begin(), centerIndex.end(), gridIter.begin(), cell.begin(), std::plus<>());
		auto cellCenter = centerOf(cell);

		if(lenSquared(cellCenter - sphere.center) <= maxDist * maxDist) {
			result.push_back(cell);
		}
		
	} while(next(gridIter, discreteRadius));

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

template <uint D, typename Precision, typename IndexPrecision>
Precision Grid<D, Precision, IndexPrecision>::diagonalCellLength() const
{
	return sqrt(D) * cellWidth;
}

template <uint D, typename Precision, typename IndexPrecision>
IndexPrecision Grid<D, Precision, IndexPrecision>::fitOneDimension(Precision x) const
{
	return std::floor(x/cellWidth + 0.5);
}

}

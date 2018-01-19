#pragma once
#include <cassert>
#include <cstdint>
#include <vector>
#include <array>
#include <type_traits>
#include <algorithm>
#include "Geometry.h"
#include "load_balancing/dIndex.h"

namespace LoadBalancing
{

template <uint D, typename Precision, typename IndexPrecision = int64_t>
struct Grid
{
	Grid(Precision cellWidth);
	dIndex<D, IndexPrecision> indexAt(const dVector<D, Precision>& v);
	std::vector<dIndex<D, IndexPrecision>> intersectingIndices(const dSphere<D, Precision>& sphere);

private:
	Precision cellWidth;

	IndexPrecision fitOneDimension(Precision x) const;
};


template <uint D, typename Precision, typename IndexPrecision>
Grid<D, Precision, IndexPrecision>::Grid(Precision cellWidth)
	: cellWidth(cellWidth)
{
	assert(cellWidth > 0);
}

template <uint D, typename Precision, typename IndexPrecision>
dIndex<D, IndexPrecision> Grid<D, Precision, IndexPrecision>::indexAt(const dVector<D, Precision>& v)
{
	dIndex<D, IndexPrecision> result;
	std::transform(begin(v), end(v), begin(result), [this] (auto x) -> auto { return fitOneDimension(x); });
	return result;
}

template <uint D, typename Precision, typename IndexPrecision>
std::vector<dIndex<D, IndexPrecision>> Grid<D, Precision, IndexPrecision>::intersectingIndices(const dSphere<D, Precision>& sphere)
{
	assert(false);
	return std::vector<dIndex<D, IndexPrecision>>();
}

template <uint D, typename Precision, typename IndexPrecision>
IndexPrecision Grid<D, Precision, IndexPrecision>::fitOneDimension(Precision x) const
{
	return std::floor((x + 0.5)/cellWidth);
}

}

#pragma once
#include <array>
template <uint D, typename IndexPrecision = int64_t>
//static_assert(std::is_integral<IndexPrecision::value_type>::value, "IndexPrecision is required to be integral.");
using dIndex = std::array<IndexPrecision, D>;

template <uint D, typename IndexPrecision>

struct dIndexHasher
{
	std::size_t operator()(const dIndex<D, IndexPrecision>& key) const
	{
		size_t result = 0;
		for(auto k : key) {
			result = result * 31 + std::hash<size_t>{}(k);
		}
		return result;
	}
};


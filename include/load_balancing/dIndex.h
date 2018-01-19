#pragma once
#include <array>
template <uint D, typename IndexPrecision = int64_t>
//static_assert(std::is_integral<IndexPrecision::value_type>::value, "IndexPrecision is required to be integral.");
using dIndex = std::array<IndexPrecision, D>;


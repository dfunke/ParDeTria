#include <iostream>
#include "load_balancing/Grid.h"

namespace lb = LoadBalancing;

template <uint D, typename IndexPrecision>
std::ostream& out(std::ostream& stream, const dIndex<D, IndexPrecision>& i)
{
	stream << "(";
	copy(i.cbegin(), i.cend(), std::ostream_iterator<IndexPrecision>(stream, " "));
	stream << ")";
	return stream;
}

template <uint D, typename IndexPrecision>
bool test(const dIndex<D, IndexPrecision>& real, const dIndex<D, IndexPrecision>& expected)
{
	bool result = real == expected;
	if(!result) {
		std::cout << "Error: expected ";
		out<D, IndexPrecision>(std::cout, expected) << " but got ";
		out<D, IndexPrecision>(std::cout, real) << "\n";
	}
	return result;
}

int main() {
	bool result = true;

	constexpr auto D = 2;
	using Precision = double;
	
	lb::Grid<D, Precision> grid(1.0);
	std::cout << "Test:\n";
	result = result && test<D, int64_t>(grid.indexAt({0.0, 0.0}), {0, 0});
	result = result && test<D, int64_t>(grid.indexAt({0.0, 1.0}), {0, 1});
	result = result && test<D, int64_t>(grid.indexAt({-1.0, 0.0}), {-1, 0});
	result = result && test<D, int64_t>(grid.indexAt({0.49, -0.49}), {0, 0});
	result = result && test<D, int64_t>(grid.indexAt({0.5, -0.5}), {1, 0});
	result = result && test<D, int64_t>(grid.indexAt({0.5, -0.5001}), {1, -1});
	result = result && test<D, int64_t>(grid.indexAt({-1.49, 1.49}), {-1, 1});
	result = result && test<D, int64_t>(grid.indexAt({1.5, -1.49}), {2, -1});

	return result ? 0 : 1;
}

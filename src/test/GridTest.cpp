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
bool testEqual(const dIndex<D, IndexPrecision>& real, const dIndex<D, IndexPrecision>& expected)
{
	bool result = real == expected;
	if(!result) {
		std::cout << "Error: expected ";
		out<D, IndexPrecision>(std::cout, expected) << " but got ";
		out<D, IndexPrecision>(std::cout, real) << "\n";
	}
	return result;
}

template <uint D, typename IndexPrecision>
bool testEqualSets(const std::vector<dIndex<D, IndexPrecision>>& real, const std::vector<dIndex<D, IndexPrecision>>& expected) {
	std::vector<dIndex<D, IndexPrecision>> superfluous;
	std::copy_if(real.begin(), real.end(), std::back_inserter(superfluous), [&expected] (const auto& i) -> bool {
	                                 return std::find(expected.begin(), expected.end(), i) == expected.end();
									 });

	std::vector<dIndex<D, IndexPrecision>> missing;
	std::copy_if(expected.begin(), expected.end(), std::back_inserter(missing), [&real] (const auto& i) -> bool {
	                                 return std::find(real.begin(), real.end(), i) == real.end();
									 });

	bool result = missing.empty() && superfluous.empty();
	if(!result) {
		std::cout << "Error:\n";
		for(const auto& i : missing) {
			out<D, IndexPrecision>(std::cout << "\tmissing: ", i) << "\n";
		}
		
		for(const auto& i : superfluous) {
			out<D, IndexPrecision>(std::cout << "\tsuperfluous: ", i) << "\n";
		}
	}
	
	return result;
}

void info() {
	static size_t k;
	std::cout << "Test " << ++k << "..\n";
}

int main() {
	bool result = true;

	constexpr auto D = 2;
	using Precision = double;
	
	std::cout << "2D grid with r=1\n";
	lb::Grid<D, Precision> grid(1.0);
	info();
	result = result & testEqual<D, int64_t>(grid.indexAt({0.0, 0.0}), {0, 0});
	info();
	result = result & testEqual<D, int64_t>(grid.indexAt({0.0, 1.0}), {0, 1});
	info();
	result = result & testEqual<D, int64_t>(grid.indexAt({-1.0, 0.0}), {-1, 0});
	info();
	result = result & testEqual<D, int64_t>(grid.indexAt({0.49, -0.49}), {0, 0});
	info();
	result = result & testEqual<D, int64_t>(grid.indexAt({0.5, -0.5}), {1, 0});
	info();
	result = result & testEqual<D, int64_t>(grid.indexAt({0.5, -0.5001}), {1, -1});
	info();
	result = result & testEqual<D, int64_t>(grid.indexAt({-1.49, 1.49}), {-1, 1});
	info();
	result = result & testEqual<D, int64_t>(grid.indexAt({1.5, -1.49}), {2, -1});

	info();
	result = result & testEqualSets<D, int64_t>(grid.intersectingIndices(dSphere<D, Precision>{{0.0, 0.0}, 0.2}), {{0, 0}});
	info();
	result = result & testEqualSets<D, int64_t>(grid.intersectingIndices(dSphere<D, Precision>{{5.0, 2.0}, 0.2}), {{5, 2}});
	info();
	result = result & testEqualSets<D, int64_t>(grid.intersectingIndices(dSphere<D, Precision>{{0.0, 0.0}, 0.6}), {{0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}});
	info();
	result = result & testEqualSets<D, int64_t>(grid.intersectingIndices(dSphere<D, Precision>{{0.0, 0.0}, 0.8}), {{0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}, {1, 1}, {-1, 1}, {-1, -1}, {1, -1}});

	constexpr double pi = atan(1)*4;
	std::cout << "2D grid with r=" << pi << "\n";
	lb::Grid<D, Precision> piGrid(pi);
	info();
	result = result & testEqual<D, int64_t>(piGrid.indexAt({0.0, 0.0}), {0, 0});
	info();
	result = result & testEqual<D, int64_t>(piGrid.indexAt({0.0, pi}), {0, 1});
	info();
	result = result & testEqual<D, int64_t>(piGrid.indexAt({pi/2 - 0.1, -pi/2 + 0.1}), {0, 0});	
	info();
	result = result & testEqual<D, int64_t>(piGrid.indexAt({pi/2, -pi/2}), {1, 0});
	info();
	result = result & testEqual<D, int64_t>(piGrid.indexAt({pi/2, -pi/2 - 0.1}), {1, -1});
	info();
	result = result & testEqual<D, int64_t>(piGrid.indexAt({5, 2}), {2, 1});

	info();
	result = result & testEqualSets<D, int64_t>(piGrid.intersectingIndices(dSphere<D, Precision>{{0.0, 0.0}, 0.2*pi}), {{0, 0}});
	info();
	result = result & testEqualSets<D, int64_t>(piGrid.intersectingIndices(dSphere<D, Precision>{{5.0*pi, 2.0*pi}, 0.2}), {{5, 2}});
	info();
	result = result & testEqualSets<D, int64_t>(piGrid.intersectingIndices(dSphere<D, Precision>{{0.0, 0.0}, 0.6}), {{0, 0}});
	info();
	result = result & testEqualSets<D, int64_t>(piGrid.intersectingIndices(dSphere<D, Precision>{{0.0, 0.0}, 0.6*pi}), {{0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}});


	std::cout << "3D grid with r=1.0\n";
	lb::Grid<3, Precision> grid3(1.0);
	info();
	result = result & testEqual<3, int64_t>(grid3.indexAt({0, 0, 0}), {0, 0, 0});
	info();
	result = result & testEqual<3, int64_t>(grid3.indexAt({0.5, 0.5, 0.5}), {1, 1, 1});
	info();
	result = result & testEqual<3, int64_t>(grid3.indexAt({0.5, -0.5, 0.5}), {1, 0, 1});

	info();
	result = result & testEqualSets<3, int64_t>(grid3.intersectingIndices(dSphere<3, Precision>{{0.0, 0.0, 0.0}, 0.1}), {{0, 0, 0}});
	info();
	result = result & testEqualSets<3, int64_t>(grid3.intersectingIndices(dSphere<3, Precision>{{0.0, 0.0, 0.0}, 0.2}), {{0, 0, 0},
																		   {1, 0, 0}, {-1, 0, 0},
																		   {0, 1, 0}, {0, -1, 0},
																		   {0, 0, 1}, {0, 0, -1}});
	info();
	result = result & testEqualSets<3, int64_t>(grid3.intersectingIndices(dSphere<3, Precision>{{0.0, 0.0, 0.0}, 0.6}), {{0, 0, 0},
																		   {1, 0, 0}, {-1, 0, 0},
																		   {0, 1, 0}, {0, -1, 0},
																		   {0, 0, 1}, {0, 0, -1},
																		   {1, 1, 0}, {-1, 1, 0}, {1, -1, 0}, {-1, -1, 0},
																		   {1, 0, 1}, {-1, 0, 1}, {1, 0, -1}, {-1, 0, -1},
																		   {0, 1, 1}, {0, -1, 1}, {0, 1, -1}, {0, -1, -1}
																		   });
	info();
	result = result & testEqualSets<3, int64_t>(grid3.intersectingIndices(dSphere<3, Precision>{{0.0, 0.0, 0.0}, 0.9}), {{0, 0, 0},
																		   {1, 0, 0}, {-1, 0, 0},
																		   {0, 1, 0}, {0, -1, 0},
																		   {0, 0, 1}, {0, 0, -1},
																		   {1, 1, 0}, {-1, 1, 0}, {1, -1, 0}, {-1, -1, 0},
																		   {1, 0, 1}, {-1, 0, 1}, {1, 0, -1}, {-1, 0, -1},
																		   {0, 1, 1}, {0, -1, 1}, {0, 1, -1}, {0, -1, -1},
																		   {1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {-1, -1, 1}, {1, 1, -1}, {-1, 1, -1}, {1, -1, -1}, {-1, -1, -1}
																		   });


	std::cout << "3D grid with r=" << pi << "\n";
	lb::Grid<3, Precision> piGrid3(pi);
	info();
	result = result & testEqual<3, int64_t>(piGrid3.indexAt({0, 0, 0}), {0, 0, 0});
	info();
	result = result & testEqual<3, int64_t>(piGrid3.indexAt({pi*0.5, pi*0.5, pi*0.5}), {1, 1, 1});
	info();
	result = result & testEqual<3, int64_t>(piGrid3.indexAt({pi*0.5, pi*(-0.5), pi*0.5}), {1, 0, 1});
	return result ? 0 : 1;
}

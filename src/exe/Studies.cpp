// own
#include "Geometry.h"
#include "DCTriangulator.h"
#include "Partitioner.h"

#include "utils/Random.h"
#include "utils/Generator.h"
#include "utils/Logger.h"
#include "utils/CSV.h"
#include "utils/Timings.h"
#include "utils/System.h"
#include "utils/ASSERT.h"
#include "utils/DBConnection.h"
#include "utils/Serialization.hxx"

#include <boost/program_options.hpp>
#include <tbb/task_scheduler_init.h>

//**************************

// CGAL
#define CGAL_LINKED_WITH_TBB

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/hilbert_sort.h>

// boost
#include <boost/iterator/transform_iterator.hpp>
#include <iomanip>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

template<uint D, typename Precision, uint I>
struct LessStruct {
    bool operator()(const dPoint<D, Precision> &p, const dPoint<D, Precision> &q) const {
        return p.coords[I] < q.coords[I];
    }
};

template<uint D, typename Precision>
struct SpatialSortingTraits;

template<typename Precision>
struct SpatialSortingTraits<2, Precision> {

    typedef dPoint<2, Precision> Point_2;
    typedef LessStruct<2, Precision, 0> Less_x_2;
    typedef LessStruct<2, Precision, 1> Less_y_2;

    Less_x_2 less_x_2_object() const {
        return Less_x_2();
    }

    Less_y_2 less_y_2_object() const {
        return Less_y_2();
    }
};

template<typename Precision>
struct SpatialSortingTraits<3, Precision> {

    typedef dPoint<3, Precision> Point_3;
    typedef LessStruct<3, Precision, 0> Less_x_3;
    typedef LessStruct<3, Precision, 1> Less_y_3;
    typedef LessStruct<3, Precision, 2> Less_z_3;

    Less_x_3 less_x_3_object() const {
        return Less_x_3();
    }

    Less_y_3 less_y_3_object() const {
        return Less_y_3();
    }

    Less_z_3 less_z_3_object() const {
        return Less_z_3();
    }
};

//**************************

#define D 3
#define Precision double

bool operator==(const dPoint<3, Precision> &a, const K::Point_3 &b) {
    return a.coords[0] == b.x() && a.coords[1] == b.y() && a.coords[2] == b.z();
}

int main() {

    dBox<D, Precision> bounds(
            dVector<D, Precision>( {{ 0, 0, 0 }}),
    dVector < D, Precision > ( {{ 100, 100, 100 }}));
    const uint N = 1e2;

    //use same start seed for all experiment runs
    tGenerator gen(START_SEED);
    auto dice = RandomFactory<Precision>::make('u', gen);
    dPoints<D, Precision> original = genPoints(N, bounds, dice);

    SpatialSortingTraits<D, Precision> sst;

    dPoints<D, Precision> hilbertSort = original;
    CGAL::hilbert_sort(hilbertSort.begin() + 1, hilbertSort.end(), sst);

    dPoints<D, Precision> spatialSort = original;
    CGALSpatialSorter<D, Precision> sorter;
    sorter.sort(spatialSort);

    std::vector<K::Point_3> cgalHilbert;
    cgalHilbert.reserve(N);
    for (const auto &p : original)
        cgalHilbert.emplace_back(p.coords[0], p.coords[1], p.coords[2]);

    CGAL::hilbert_sort(cgalHilbert.begin() + 1, cgalHilbert.end());

    std::vector<K::Point_3> cgalSpatial;
    cgalSpatial.reserve(N);
    for (const auto &p : original)
        cgalSpatial.emplace_back(p.coords[0], p.coords[1], p.coords[2]);

    CGAL::spatial_sort(cgalSpatial.begin() + 1, cgalSpatial.end());

    std::cout << std::setw(4) << "idx"
    << std::setw(30) << "Unsorted"

    << std::setw(30) << "HilbertSort"
    << std::setw(30) << "cgalHilbert"

    << std::setw(30) << "SpatialSort"
    << std::setw(30) << "cgalSpatial"
    << std::endl;

    for (uint i = 1; i <= N; ++i) {
        if (!(hilbertSort[i] == cgalHilbert[i]))
            std::cerr << "[" << i << "] Invalid HilbertSort" << std::endl;
        if (!(spatialSort[i] == cgalSpatial[i]))
            std::cerr << "[" << i << "] Invalid SpatialSort" << std::endl;

        std::cout << std::setw(4) << i
        << std::setw(30) << original[i]

        << std::setw(30) << hilbertSort[i]
        << std::setw(30) << cgalHilbert[i]

        << std::setw(30) << spatialSort[i]
        << std::setw(30) << cgalSpatial[i]

        << std::endl;


    }

    //*************************************************************

    std::vector<tIdType> allPoints;
    allPoints.reserve(N);
    for (uint i = 1; i <= N; ++i)
        allPoints.emplace_back(i);

    auto stats = getPointStats(allPoints.begin(), allPoints.end(), spatialSort);
    const uint k = 0;

    std::vector<std::vector<tIdType>> spatialPart;
    spatialPart.resize(2);

    std::vector<std::vector<tIdType>> hilbertPart;
    hilbertPart.resize(2);

    for (uint i = 1; i <= N; ++i) {
        uint part = (spatialSort[i].coords[k] > stats.mid[k]);
        spatialPart[part].push_back(i);

        part = (hilbertSort[i].coords[k] > stats.mid[k]);
        hilbertPart[part].push_back(i);
    }

    std::cout
    << std::setw(10) << "spatialPart. 0"
    << std::setw(10) << "spatialPart. 1"
    << std::setw(1) << "|"
    << std::setw(10) << "hilbertPart. 0"
    << std::setw(10) << "hilbertPart. 1"
    << std::endl;

    auto spatialIT0 = spatialPart[0].begin();
    auto spatialIT1 = spatialPart[1].begin();
    auto hilbertIT0 = hilbertPart[0].begin();
    auto hilbertIT1 = hilbertPart[1].begin();

    while (spatialIT0 != spatialPart[0].end() && spatialIT1 != spatialPart[1].end()
           && hilbertIT0 != hilbertPart[0].end() && hilbertIT1 != hilbertPart[1].end()) {
        std::cout
        << std::setw(10) << (spatialIT0 != spatialPart[0].end() ? *(spatialIT0++) : 0)
        << std::setw(10) << (spatialIT1 != spatialPart[1].end() ? *(spatialIT1++) : 0)
        << std::setw(1) << "|"
        << std::setw(10) << (hilbertIT0 != hilbertPart[0].end() ? *(hilbertIT0++) : 0)
        << std::setw(10) << (hilbertIT1 != hilbertPart[1].end() ? *(hilbertIT1++) : 0)
        << std::endl;
    }

}
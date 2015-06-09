#include "SpatialSort.h"

#ifndef NDEBUG

#include <csignal>

#endif

#include <atomic>
#include <type_traits>

#include "utils/ASSERT.h"
#include "utils/StaticSort.h"
#include "utils/VTuneAdapter.h"

// CGAL
#define CGAL_LINKED_WITH_TBB

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/spatial_sort.h>

// boost
#include <boost/iterator/transform_iterator.hpp>

#include <CGAL/Unique_hash_map.h>


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

// partial specializations

template<uint D, typename Precision>
void CGALSpatialSorter<D, Precision>::sort(dPoints<D, Precision> &points
        /*, bool filterInfinite = false */) {

    SpatialSortingTraits<D, Precision> sst;

    VTUNE_TASK(SpatialSort);
    CGAL::spatial_sort(points.begin(), points.end(), sst);
}

// specializations
template
class CGALSpatialSorter<2, float>;

template
class CGALSpatialSorter<3, float>;

template
class CGALSpatialSorter<2, double>;

template
class CGALSpatialSorter<3, double>;
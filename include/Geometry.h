/*
 * Geometry.h
 *
 *  Created on: Oct 28, 2014
 *      Author: dfunke
 */

#pragma once

#include <iostream>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <array>
#include <cmath>
#include <functional>
#include <set>

#include <datastructures/Growing_LP.hxx>
//#include "datastructures/IndexedVector.hxx"
//#include "datastructures/VectorAdapter.hxx"
#include "datastructures/VectorAdapter2.hxx"
#include "datastructures/LP_MultiMap.hxx"
//#include "datastructures/BlockedArray.hxx"
#include "datastructures/BlockedArray2.hxx"
//#include "datastructures/Bit_Set.hxx"

#include "utils/Logger.h"
#include "utils/ASSERT.h"
#include "utils/Timings.h"
#include "utils/Misc.h"
#include "utils/StaticSort.h"

typedef uint64_t tHashType;
typedef uint64_t tIdType;

typedef LP_Set<tIdType, true> Simplex_Ids;
typedef Concurrent_LP_Set<tIdType, false> Concurrent_Fixed_Simplex_Ids;
typedef GrowingHashTable<Concurrent_LP_Set<tIdType, true>> Concurrent_Growing_Simplex_Ids;
typedef GrowingHashTableHandle<Concurrent_LP_Set<tIdType, true>> hConcurrent_Growing_Simplex_Ids;

typedef LP_Set<tIdType, true> Point_Ids;
typedef Concurrent_LP_Set<tIdType, false> Concurrent_Fixed_Point_Ids;
typedef GrowingHashTable<Concurrent_LP_Set<tIdType, true>> Concurrent_Growing_Point_Ids;
typedef GrowingHashTableHandle<Concurrent_LP_Set<tIdType, true>> hConcurrent_Growing_Point_Ids;

// basic data structure for d-dimensional data
template<std::size_t D, typename Precision>
using dVector = std::array<Precision, D>;

//#############################################################################

// sphere
template<uint D, typename Precision>
struct dSphere {
    dVector<D, Precision> center;
    Precision radius;
};

//#############################################################################

// (hyper)-rectangle
template<uint D, typename Precision>
struct dBox {
    dVector<D, Precision> low;
    dVector<D, Precision> high;

    dBox() { }

    dBox(const dVector<D, Precision> &_low, const dVector<D, Precision> &_high)
            : low(_low), high(_high) { }

    bool contains(const dVector<D, Precision> &p) const {
        PROFILER_INC("dBox_containsPoint");

        for (uint d = 0; d < D; ++d) {
            if (!(low[d] <= p[d] && p[d] <= high[d]))
                return false;
        }

        return true;
    }

    /* tests whether sphere is FULLY contained in box */
    bool contains(const dSphere<D, Precision> &sphere) const {
        PROFILER_INC("dBox_containsSphere");

        Precision r2 = sphere.radius * sphere.radius;

        for (uint d = 0; d < D; ++d) {
            if (!(low[d] <= sphere.center[d] && sphere.center[d] <= high[d])) {
                return false;
            }
        }

        // the center of the sphere is within the box
        for (uint d = 0; d < D; ++d) {
            auto p = sphere.center;
            // project p to the boundary of box in dimension d closest to center of
            // the
            // sphere
            p[d] = sphere.center[d] < (high[d] + low[d]) / 2 ? low[d] : high[d];

            Precision dist = 0;
            for (uint i = 0; i < D; ++i)
                dist += (sphere.center[i] - p[i]) * (sphere.center[i] - p[i]);

            if (dist < r2)
                return false;
        }

        return true;
    }

    /* Tests whehter a sphere intersects with the box */
    bool intersects(const dSphere<D, Precision> &sphere) const {
        PROFILER_INC("dBox_intersectsSphere");

        Precision r2 = sphere.radius * sphere.radius;

        Precision dist = 0;
        for (uint d = 0; d < D; ++d) {
            Precision e = std::max(low[d] - sphere.center[d], (Precision) 0) +
                          std::max(sphere.center[d] - high[d], (Precision) 0);

            if (sphere.radius < e)
                return false;
            dist += e * e;
        }

        if (dist <= r2)
            return true;

        return false;
    }

    dBox<D, Precision> operator-(const dBox<D, Precision> & other) const{
        dBox<D, Precision> result(low, high);

        // we can only substract a rectangle other that is equal to ours, except in one dimension

#ifndef NDEBUG
        uint mods = 0;
#endif
        for (uint d = 0; d < D; ++d) {
            if (low[d] != other.low[d]) {
                result.high[d] = other.low[d];
#ifndef NDEBUG
                ++mods;
#endif
            }

            if (high[d] != other.high[d]) {
                result.low[d] = other.high[d];
#ifndef NDEBUG
                ++mods;
#endif
            }
        }

        ASSERT(mods == 1);

        return result;
    };

    bool operator==(const dBox<D, Precision> &other) const {
        return low == other.low && high == other.high;
    }
};

//#############################################################################

// point with ID and simplices list

template<uint D, typename Precision>
class dPoint {

public:
    dPoint() { }

    dPoint(const dVector<D, Precision> &_coords)
            : coords(_coords) { }

public:
    dVector<D, Precision> coords;

public:
    static inline bool isFinite(const tIdType &i) {
        PROFILER_INC("dPoint_staticIsFinite");
        return !((i & cINF) == cINF);
    }

    static inline tIdType infIndex(const tIdType &i) {
        ASSERT(!isFinite(i));
        return ~cINF & i;
    }

    static constexpr tIdType cINF = ~tIdType(0) ^((1 << D) - 1);
    static constexpr tIdType nINF = ~tIdType(0) - cINF + 1;
};

//#############################################################################

// IndexedVector of dPoints

template<uint D, typename Precision>
class dPoints : public VectorAdapter2<dPoint<D, Precision>> {

    typedef VectorAdapter2<dPoint<D, Precision>> super;

public:
    dPoints() : VectorAdapter2<dPoint<D, Precision>>() { }

    dPoints(const VectorAdapter2<dPoint<D, Precision>> &other)
            : VectorAdapter2<dPoint<D, Precision>>(other) { }

    bool operator==(const Concurrent_Growing_Point_Ids &other) const {
        return operator==(other.handle());
    }

    bool contains(const tIdType &i) const {
        PROFILER_INC("dPoints_contains");

        //we have to check whether i is GREATER than offset, because zero is reserved

        return (__builtin_expect(dPoint<D, Precision>::isFinite(i), true))
               ? (i > super::offset() && super::_finIdx(i) < super::finite_size())
               : (super::_infIdx(i) < dPoint<D, Precision>::nINF);
    }

    auto range() const {
        // we have to adapt the range to exclude 0
        return tbb::blocked_range<std::size_t>(super::offset() + 1, super::offset() + super::finite_size());
    }

    template<class Container>
    bool operator==(const Container &other) const {
        PROFILER_INC("dPoints_compare");

        //account for 0 not being used
        if (this->size() - 1 != other.size())
            return false;

        for (const auto &p : other) {
            if (!this->contains(p))
                return false;
        }

        return true;
    }

    template<class Container>
    bool operator!=(const Container &other) const { return !operator==(other); }
};

//#############################################################################

template<uint D, typename Precision>
class dSimplex;

template<uint D, typename Precision>
class GeometryCore {
public:
    static Precision orientation(const dPoint<D, Precision> &s0, const dPoint<D, Precision> &s1,
                                 const dPoint<D, Precision> &s2);

    static Precision orientation(const dPoint<D, Precision> &s0, const dPoint<D, Precision> &s1,
                                 const dPoint<D, Precision> &s2, const dPoint<D, Precision> &s3);

    static bool inSphere(const dPoint<D, Precision> &p,
                         const dPoint<D, Precision> &s0, const dPoint<D, Precision> &s1,
                         const dPoint<D, Precision> &s2);

    static bool inSphere(const dPoint<D, Precision> &p,
                         const dPoint<D, Precision> &s0, const dPoint<D, Precision> &s1, const dPoint<D, Precision> &s2,
                         const dPoint<D, Precision> &s3);

    static dSphere<D, Precision>
            circumsphere(const dPoint<D, Precision> &s0, const dPoint<D, Precision> &s1,
                         const dPoint<D, Precision> &s2);

    static dSphere<D, Precision>
            circumsphere(const dPoint<D, Precision> &s0, const dPoint<D, Precision> &s1, const dPoint<D, Precision> &s2,
                         const dPoint<D, Precision> &s3);
};

template<uint D, typename Precision>
class GeometryHelper {
public:

    template<class Collection>
    static Precision orientation(const dSimplex<D, Precision> &simplex,
                                 const Collection &points);

    template<class Collection>
    static bool inSphere(const dSimplex<D, Precision> &simplex,
                         const dPoint<D, Precision> &p,
                         const Collection &points);

    template<class Collection>
    static dSphere<D, Precision>
            circumsphere(const dSimplex<D, Precision> &simplex,
                         const Collection &points);
};

// d-Simplex

template<uint D, typename Precision>
class dSimplex {

public:
    dSimplex() : id(dSimplex<D, Precision>::cINF) { }

    dSimplex(const std::array<tIdType, D + 1> &_vertices)
            : id(dSimplex<D, Precision>::cINF), vertices(_vertices) { }

    bool operator==(const dSimplex<D, Precision> &a) const {
        PROFILER_INC("dSimplex_compare");

        if (!((id == cINF) ^ (a.id == cINF))) {
            // either none or both simplices are infinity
            // compare ids
            PROFILER_INC("dSimplex_compare-Id");

            return id == a.id;
        } else {
            // either one is infinity
            // compare vertices

            PROFILER_INC("dSimplex_compare-Points");

            return equalVertices(a);
        }
    }

    // for use as map key
    bool operator<(const dSimplex<D, Precision> &a) const { return id < a.id; }

    bool equalVertices(const dSimplex<D, Precision> &a) const {
        PROFILER_INC("dSimplex_equalVertices");

        // first compare fingerprints
        if (fingerprint() != a.fingerprint()) {
            PROFILER_INC("dSimplex_equalVertices-Fingerprint");

            return false;
        }

        // fingerprints are equal => compare vertices
        // both vertex arrays are sorted by vertex id
        PROFILER_INC("dSimplex_equalVertices-Points");

        for (uint i = 0; i < D + 1; ++i) {
            if (vertices[i] != a.vertices[i]) {
                return false;
            }
        }

        return true;
    }

    uint countSharedVertices(const dSimplex<D, Precision> &a) const {
        PROFILER_INC("dSimplex_sharedVertices");

        uint sharedVertices = 0;

        uint i = 0;
        uint j = 0;
        uint tmp = 0;
        while (i < D + 1 && j < D + 1) {
            sharedVertices += vertices[i] == a.vertices[j];

            tmp = i;
            i += vertices[i] <= a.vertices[j];
            j += vertices[tmp] >= a.vertices[j];
        }

        ASSERT(0 <= sharedVertices && sharedVertices <= D + 1);

        return sharedVertices;
    }

    bool equalNeighbors(const dSimplex<D, Precision> &a) const {
        PROFILER_INC("dSimplex_equalNeighbors");

        return neighbors != a.neighbors;
    }

    bool operator==(const tIdType a) const { return id == a; }

    bool contains(tIdType p) const {
        PROFILER_INC("dSimplex_containsId");

        // vertices are sorted by point id
        for (uint d = 0; d < D + 1 && vertices[d] <= p; ++d) {
            if (p == vertices[d])
                return true;
        }

        return false;
    }

    /*bool contains(const dPoint<D, Precision> &p) const {
        PROFILER_INC("dSimplex_containsPoint");

        for (uint d = 0; d < D + 1; ++d) {
            if (p == vertices[d])
                return true;
        }

        return false;
    }*/

    template<typename Container>
    bool containsAll(const Container &p) const {
        PROFILER_INC("dSimplex_containsAll");

        for (const auto &x : p) {
            if (!contains(x))
                return false;
        }

        return true;
    }

    template<typename Container>
    bool containsAny(const Container &p) const {
        PROFILER_INC("dSimplex_containsAny");

        for (const auto &x : p) {
            if (contains(x))
                return true;
        }

        return false;
    }

    bool isFinite() const {
        PROFILER_INC("dSimplex_isFinite");

        //the vertices are sorted by ID -> infinite points are at the end
        for (int i = D; i >= 0; --i) {
            if (!dPoint<D, Precision>::isFinite(vertices[i]))
                return false;
        }

        return true;
    }

    template<class Collection>
    Precision orientation(const Collection &points) const {
        PROFILER_INC("dSimplex_orientation");

        return GeometryHelper<D, Precision>::orientation(*this, points);
    }

    template<class Collection>
    bool inSphere(const dPoint<D, Precision> &p, const Collection &points) const {
        PROFILER_INC("dSimplex_inSphere");

        return GeometryHelper<D, Precision>::inSphere(*this, p, points);
    }

    template<class Collection>
    dSphere<D, Precision> circumsphere(
            const Collection &points) const {
        PROFILER_INC("dSimplex_circumsphere");

        return GeometryHelper<D, Precision>::circumsphere(*this, points);
    }

    bool isNeighbor(const dSimplex<D, Precision> &other) const {
        PROFILER_INC("dSimplex_isNeighbor");

        uint sharedVertices = countSharedVertices(other);

        ASSERT(sharedVertices <= D || id == other.id);

        return sharedVertices == D;
    }

    void sortVertices() {
        static_insertion_sort_with_data(vertices, neighbors);
    }

    tHashType genFingerprint() {
        PROFILER_INC("dSimplex_genFingerprint");

        vFingerprint = 0;
        for (uint i = 0; i < D + 1; ++i) {
            vFingerprint ^= _rol(vertices[i]);
        }

        return vFingerprint;
    }

    inline tHashType fingerprint() const {
        PROFILER_INC("dSimplex_fingerprint");
        return vFingerprint;
    }

    inline tHashType faceFingerprint(const uint i) const {
        PROFILER_INC("dSimplex_faceFingerprint");
        return (vFingerprint ^ _rol(vertices[i])) | 1; //guarantee that hash does _not_ become zero
    }

private:

    inline tHashType _rol(const tHashType x) const {
        return _rol(x, x & (sizeof(tHashType) * 8 - 1));
    }

    inline tHashType _rol(const tHashType x, const uint r) const {
        return (x << r) | (x >> (sizeof(tHashType) * 8 - r));
    }

public:
    tIdType id;
    std::array<tIdType, D + 1> vertices;
    tHashType vFingerprint;
    std::array<tIdType, D + 1> neighbors;

public:
    static bool isFinite(const tIdType &i) {
        PROFILER_INC("dSimplex_staticIsFinite");

        return i != cINF;
    }

    static constexpr tIdType cINF = ~tIdType(0);
};

template<std::size_t D, typename Precision>
std::string to_string(const dVector<D, Precision> &p);

template<uint D, typename Precision>
std::string to_string(const dPoint<D, Precision> &p);

template<uint D, typename Precision>
std::string to_string(const dSimplex<D, Precision> &p);

template<uint D, typename Precision>
std::string to_string(const dSphere<D, Precision> &p);

template<uint D, typename Precision>
std::string to_string(const dBox<D, Precision> &b);

template<std::size_t D, typename Precision>
std::ostream &operator<<(std::ostream &o, const dVector<D, Precision> &p);

template<uint D, typename Precision>
std::ostream &operator<<(std::ostream &o, const dPoint<D, Precision> &p);

template<uint D, typename Precision>
std::ostream &operator<<(std::ostream &o, const dSphere<D, Precision> &p);

template<uint D, typename Precision>
std::ostream &operator<<(std::ostream &o, const dSimplex<D, Precision> &p);

template<uint D, typename Precision>
std::ostream &operator<<(std::ostream &o, const dBox<D, Precision> &b);

template<uint D, typename Precision>
struct CrossCheckReport;

template<uint D, typename Precision>
struct VerificationReport;

//class cConvexHull : private tbb::concurrent_unordered_set<tIdType> {
//
//private:
//    typedef tbb::concurrent_unordered_set<tIdType> base;
//
//public:
//
//    auto insert(const tIdType &id) {
//        PROFILER_INC("dSimplices_convexHull_insert");
//
//        return base::insert(id);
//    }
//
//    template<typename IT>
//    void insert(IT first, IT last) {
//        PROFILER_ADD("dSimplices_convexHull_insert", std::distance(first, last));
//
//        base::insert(first, last);
//    }
//
//    auto begin() {
//        PROFILER_INC("dSimplices_convexHull_begin");
//
//        return base::begin();
//    }
//
//    auto end() {
//        return base::end();
//    }
//
//    auto begin() const {
//        PROFILER_INC("dSimplices_convexHull_begin");
//
//        return base::begin();
//    }
//
//    auto end() const {
//        return base::end();
//    }
//
//    auto begin(std::size_t i) {
//        PROFILER_INC("dSimplices_convexHull_localBegin");
//
//        return base::unsafe_begin(i);
//    }
//
//    auto end(std::size_t i) {
//
//        return base::unsafe_end(i);
//    }
//
//    auto begin(std::size_t i) const {
//        PROFILER_INC("dSimplices_convexHull_localBegin");
//
//        return base::unsafe_begin(i);
//    }
//
//    auto end(std::size_t i) const {
//
//        return base::unsafe_end(i);
//    }
//
//    void reserve(std::size_t s) {
//        PROFILER_INC("dSimplices_convexHull_reserve");
//
//        base::rehash(nextPow2(s));
//    }
//
//    auto size() const {
//        return base::size();
//    }
//
//    auto bucket_count() const {
//        return base::unsafe_bucket_count();
//    }
//
//    auto count(tIdType k) const {
//        PROFILER_INC("dSimplices_convexHull_count");
//
//        return base::count(k);
//    }
//
//    auto clear() {
//        return base::clear();
//    }
//
//};

typedef GrowingHashTable<Concurrent_LP_MultiMap<tIdType, tIdType>> cWuFaces;
typedef GrowingHashTableHandle<Concurrent_LP_MultiMap<tIdType, tIdType>> hcWuFaces;
//typedef tbb::concurrent_unordered_multimap<tHashType, tIdType> cWuFaces;

//class cWuFaces : private tbb::concurrent_unordered_multimap<tHashType, tIdType> {
//
//private:
//    typedef tbb::concurrent_unordered_multimap<tHashType, tIdType> base;
//
//public:
//
//    void reserve(std::size_t s) {
//        PROFILER_INC("dSimplices_wuFaces_reserve");
//
//        base::rehash(nextPow2(s));
//    }
//
//    auto begin() {
//        PROFILER_INC("dSimplices_wuFaces_begin");
//
//        return base::begin();
//    }
//
//    auto end() {
//        return base::end();
//    }
//
//    auto begin() const {
//        PROFILER_INC("dSimplices_wuFaces_begin");
//
//        return base::begin();
//    }
//
//    auto end() const {
//        return base::end();
//    }
//
//    auto begin(std::size_t i) {
//        PROFILER_INC("dSimplices_wuFaces_localBegin");
//
//        return base::unsafe_begin(i);
//    }
//
//    auto end(std::size_t i) {
//
//        return base::unsafe_end(i);
//    }
//
//    auto begin(std::size_t i) const {
//        PROFILER_INC("dSimplices_wuFaces_localBegin");
//
//        return base::unsafe_begin(i);
//    }
//
//    auto end(std::size_t i) const {
//
//        return base::unsafe_end(i);
//    }
//
//    template<typename Pair>
//    auto insert(Pair &&__x) {
//        PROFILER_INC("dSimplices_wuFaces_insert");
//
//        return base::insert(std::forward<Pair>(__x));
//    }
//
//    template<typename... Args>
//    auto emplace(Args &&... __args) {
//        PROFILER_INC("dSimplices_wuFaces_emplace");
//
//        return base::emplace(std::forward<Args>(__args)...);
//    }
//
//    template<typename IT>
//    void insert(IT first, IT last) {
//        PROFILER_ADD("dSimplices_wuFaces_insert", std::distance(first, last));
//
//        base::insert(first, last);
//    }
//
//    auto equal_range(const tHashType &key) const {
//        PROFILER_INC("dSimplices_wuFaces_equal_range");
//
//        return base::equal_range(key);
//    }
//
//    auto equal_range(const tHashType &key) {
//        PROFILER_INC("dSimplices_wuFaces_equal_range");
//
//        return base::equal_range(key);
//    }
//
//    auto size() const {
//        return base::size();
//    }
//
//    auto bucket_count() const {
//        return base::unsafe_bucket_count();
//    }
//
//};

template<uint D, typename Precision>
class dSimplices : public BlockedArray2<dSimplex<D, Precision>, tIdType> {

private:
    typedef BlockedArray2<dSimplex<D, Precision>, tIdType> base;

public:
    typedef std::array<uint, 256> tHash;

public:
    dSimplices() : base(0, 1), convexHull(1) { }

    dSimplices(const tIdType min, const tIdType max) : base(min, max),
                                                       convexHull((max-min) / 4) { }

    dSimplices(dSimplices &&other) : base(std::move(other)),
                                     convexHull(std::move(other.convexHull)) { }

    dSimplices &operator=(dSimplices &&other) {
        base::operator=(std::move(other));
        convexHull = std::move(other.convexHull);

        return *this;
    }

    VerificationReport<D, Precision> verify(const dPoints<D, Precision> &points, const Point_Ids *pointIds = nullptr) const;

    CrossCheckReport<D, Precision> crossCheck(const dSimplices<D, Precision> &realSimplices) const;

    tHash genFingerprint() const;

//    uint countDuplicates(const Simplex_Ids &simplices) const;

public:

    auto &addBlock(const tIdType min, const tIdType max) {
        return base::addBlock(min, max);
    }

    bool unsafe_contains(const dSimplex<D, Precision> &s) const {
        return unsafe_contains(s.id);
    }

    bool unsafe_contains(const tIdType id) const {
        if (base::unsafe_contains(id))
            return _valid(this->unsafe_at(id));

        return false;
    }

    bool contains(const dSimplex<D, Precision> &s, uint & hint) const {
        return contains(s.id, hint);
    }

    bool contains(const tIdType id, uint & hint) const {
        if (base::contains(id, hint))
            return _valid(this->at(id, hint));

        return false;
    }

    bool unsafe_inRange(const dSimplex<D, Precision> &s) const {
        return unsafe_inRange(s.id);
    }

    bool unsafe_inRange(const tIdType id) const {
        return base::unsafe_contains(id);
    }

    bool inRange(const dSimplex<D, Precision> &s, uint & hint) const {
        return inRange(s.id, hint);
    }

    bool inRange(const tIdType id, uint & hint) const {
        return base::contains(id, hint);
    }

    void merge(dSimplices<D, Precision> && other){
        base::merge(std::move(other));
        convexHull.merge(std::move(other.convexHull));
    }

    template <class Filter>
    void merge(dSimplices<D, Precision> && other, const Filter & filter, const bool cmp = false){
        base::merge(std::move(other));
        convexHull.merge(std::move(other.convexHull), filter, cmp);
    }

//    void offset(const tIdType offset){
//        base::offset(offset);
//        simplexID += offset;
//    }

public:
    typedef _detail::filtered_block_iterator<dSimplices<D, Precision>, dSimplex<D, Precision>> iterator;
    friend iterator;

    typedef _detail::filtered_block_iterator<const dSimplices<D, Precision>, dSimplex<D, Precision>> const_iterator;
    friend const_iterator;

    typedef _detail::range_type<dSimplices<D, Precision>, iterator> range_type;
    typedef _detail::range_type<const dSimplices<D, Precision>, const_iterator> const_range_type;

    const_iterator begin() const {
        return const_iterator(*this, 0, this->lowerBound(), true);
    }

    const_iterator end() const {
        return const_iterator(*this, this->m_blocks.size() - 1, this->upperBound(), false);
    }

    iterator begin() {
        return iterator(*this, 0, this->lowerBound(), true);
    }

    iterator end() {
        return iterator(*this, this->m_blocks.size() - 1, this->upperBound(), false);
    }

    const_range_type range() const {
        return const_range_type(*this);
    }

    range_type range() {
        return range_type(*this);
    }

private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, __attribute((unused)) const unsigned int version) {
        ar & boost::serialization::base_object<base>(*this);
        ar & convexHull;
    }


private:
    bool _valid(const dSimplex<D, Precision> &a) const {
        return dSimplex<D, Precision>::isFinite(a.id);
    }

public:
    Simplex_Ids convexHull;

public:
    static std::atomic<tIdType> simplexID;

};

template<uint D, typename Precision>
using dSimplicesHandle = BlockedArray2Handle<dSimplices<D, Precision>, tIdType>;

template<uint D, typename Precision>
using dSimplicesConstHandle = BlockedArray2Handle<const dSimplices<D, Precision>, tIdType>;

template<uint D, typename Precision>
struct CrossCheckReport {
    bool valid;
    std::vector<dSimplex<D, Precision>> missing;
    std::vector<dSimplex<D, Precision>> invalid;
};

template<uint D, typename Precision>
struct VerificationReport {
    bool valid;
    std::map<dSimplex<D, Precision>, std::set<tIdType>> inCircle;
    std::vector<std::pair<dSimplex<D, Precision>, dSimplex<D, Precision>>>
            wrongNeighbors;
};

//#############################################################################

template<typename Precision>
class GeometryHelper<2, Precision> {
public:

    template<class Collection>
    static Precision orientation(const dSimplex<2, Precision> &simplex,
                                 const Collection &points) {
        return GeometryCore<2, Precision>::orientation(points[simplex.vertices[0]], points[simplex.vertices[1]],
                                                       points[simplex.vertices[2]]);
    }

    template<class Collection>
    static bool inSphere(const dSimplex<2, Precision> &simplex,
                         const dPoint<2, Precision> &p,
                         const Collection &points) {
        return GeometryCore<2, Precision>::inSphere(p, points[simplex.vertices[0]], points[simplex.vertices[1]],
                                                    points[simplex.vertices[2]]);
    }

    template<class Collection>
    static dSphere<2, Precision>
    circumsphere(const dSimplex<2, Precision> &simplex,
                 const Collection &points) {
        return GeometryCore<2, Precision>::circumsphere(points[simplex.vertices[0]], points[simplex.vertices[1]],
                                                        points[simplex.vertices[2]]);
    }
};

template<typename Precision>
class GeometryHelper<3, Precision> {
public:

    template<class Collection>
    static Precision orientation(const dSimplex<3, Precision> &simplex,
                                 const Collection &points) {
        return GeometryCore<3, Precision>::orientation(points[simplex.vertices[0]], points[simplex.vertices[1]],
                                                       points[simplex.vertices[2]], points[simplex.vertices[3]]);
    }

    template<class Collection>
    static bool inSphere(const dSimplex<3, Precision> &simplex,
                         const dPoint<3, Precision> &p,
                         const Collection &points) {

        return GeometryCore<3, Precision>::inSphere(p, points[simplex.vertices[0]], points[simplex.vertices[1]],
                                                    points[simplex.vertices[2]], points[simplex.vertices[3]]);
    }

    template<class Collection>
    static dSphere<3, Precision>
    circumsphere(const dSimplex<3, Precision> &simplex,
                 const Collection &points) {
        return GeometryCore<3, Precision>::circumsphere(points[simplex.vertices[0]], points[simplex.vertices[1]],
                                                        points[simplex.vertices[2]], points[simplex.vertices[3]]);
    }
};
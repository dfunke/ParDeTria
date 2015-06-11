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

#include "datastructures/IndexedVector.hxx"
#include "datastructures/VectorAdapter.hxx"
#include "datastructures/LP_MultiMap.hxx"
#include "datastructures/BlockedArray.hxx"


#include "utils/Logger.h"
#include "utils/ASSERT.h"
#include "utils/Timings.h"
#include "utils/Misc.h"

typedef uint tHashType;
typedef uint tIdType;

class Ids : private tbb::concurrent_unordered_set<tIdType> {

private:
    typedef tbb::concurrent_unordered_set<tIdType> base;

public:

    Ids() : base() { }

    Ids(const Ids &other) : base(other) {
        PROFILER_INC("Ids_copy");
    }

    Ids(const Ids &&other) : base(std::move(other)) {
        PROFILER_INC("Ids_move");
    }

    Ids(const base &other) : base(other) {
        PROFILER_INC("Ids_copy");
    }

    Ids(const base &&other) : base(std::move(other)) {
        PROFILER_INC("Ids_move");
    }

    auto operator=(const Ids &other) {
        PROFILER_INC("Ids_copy");

        return base::operator=(other);
    }

    auto operator=(const base &other) {
        PROFILER_INC("Ids_copy");

        return base::operator=(other);
    }

    auto operator=(const Ids &&other) {
        PROFILER_INC("Ids_move");

        return base::operator=(std::move(other));
    }

    auto operator=(const base &&other) {
        PROFILER_INC("Ids_move");

        return base::operator=(std::move(other));
    }

    auto insert(const tIdType &id) {
        PROFILER_INC("Ids_insert");

        return base::insert(id);
    }

    template<typename IT>
    void insert(IT first, IT last) {
        PROFILER_ADD("Ids_insert", std::distance(first, last));

        base::insert(first, last);
    }

    auto begin() {
        PROFILER_INC("Ids_begin");

        return base::begin();
    }

    auto end() {
        return base::end();
    }

    auto begin() const {
        PROFILER_INC("Ids_begin");

        return base::begin();
    }

    auto end() const {
        return base::end();
    }

    auto range() {
        PROFILER_INC("Ids_begin");

        return base::range();
    }

    auto range() const {
        PROFILER_INC("Ids_begin");

        return base::range();
    }

    auto begin(std::size_t i) {
        PROFILER_INC("Ids_localBegin");

        return base::unsafe_begin(i);
    }

    auto end(std::size_t i) {

        return base::unsafe_end(i);
    }

    auto begin(std::size_t i) const {
        PROFILER_INC("Ids_localBegin");

        return base::unsafe_begin(i);
    }

    auto end(std::size_t i) const {

        return base::unsafe_end(i);
    }

    void reserve(std::size_t s) {
        PROFILER_INC("Ids_reserve");

        base::rehash(nextPow2(s));
    }

    auto size() const {
        return base::size();
    }

    auto bucket_count() const {
        return base::unsafe_bucket_count();
    }

    auto count(tIdType k) const {
        PROFILER_INC("Ids_count");

        return base::count(k);
    }

    auto clear() {
        return base::clear();
    }

    void erase(tIdType k) {
        base::unsafe_erase(k);
    }

};

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
};

//#############################################################################

// point with ID and simplices list

template<uint D, typename Precision>
class dPoint {

public:
    dPoint() { }

    dPoint(const dVector<D, Precision> &_coords)
            : coords(_coords) { }

    /*dPoint(const dPoint<D, Precision> &a) : id(a.id), coords(a.coords) {}

    dPoint &operator=(const dPoint<D, Precision> &a) {
      id = a.id;
      coords = a.coords;

      return *this;
    }*/

    /*bool operator==(const dPoint<D, Precision> &a) const {
        PROFILER_INC("dPoint_compare");

        for (uint i = 0; i < D; ++i) {
            if (coords[i] != a.coords[i]) {
                return false;
            }
        }
        return true;
    }*/

    //bool operator==(const tIdType &a) const { return id == a; }

    /*inline bool isFinite() const {
        PROFILER_INC("dPoint_isFinite");

        return isFinite(id);
    }*/

public:
    //tIdType id;
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

    static constexpr tIdType cINF = ~(0) ^((1 << D) - 1);
    static constexpr tIdType nINF = ~(0) - cINF + 1;
};

//#############################################################################

// IndexedVector of dPoints

template<uint D, typename Precision>
class dPoints : public VectorAdapter<dPoint<D, Precision>> {

public:
    dPoints() : VectorAdapter<dPoint<D, Precision>>() { }

    dPoints(const VectorAdapter<dPoint<D, Precision>> &other)
            : VectorAdapter<dPoint<D, Precision>>(other) { }

    bool operator==(const Ids &other) const {
        PROFILER_INC("dPoints_compare");

        if (this->size() != other.size())
            return false;

        for (const auto &p : other) {
            if (!this->contains(p))
                return false;
        }

        return true;
    }

    bool operator!=(const Ids &other) const { return !operator==(other); }
};

//#############################################################################

// d-Simplex

template<uint D, typename Precision>
class dSimplex {

public:
    dSimplex() : id(dSimplex<D, Precision>::cINF) { }

    dSimplex(const std::array<uint, D + 1> &_vertices)
            : id(dSimplex<D, Precision>::cINF), vertices(_vertices) { }

    /*dSimplex(const dSimplex<D, Precision> &a)
        : id(a.id), vertices(a.vertices), vertexFingerprint(a.vertexFingerprint),
          neighbors(a.neighbors) {}

    dSimplex &operator=(const dSimplex<D, Precision> &a) {
      id = a.id;
      vertices = a.vertices;
      vertexFingerprint = a.vertexFingerprint;
      neighbors = a.neighbors;

      return *this;
    }*/

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

    bool operator==(const uint &a) const { return id == a; }

    bool contains(uint p) const {
        PROFILER_INC("dSimplex_containsId");

        for (uint d = 0; d < D + 1; ++d) {
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

    // dimension specific implementations in cpp file
    Precision orientation(const dPoints<D, Precision> &points) const;

    bool inSphere(const dPoint<D, Precision> &p,
                  const dPoints<D, Precision> &points) const;

    dSphere<D, Precision> circumsphere(const dPoints<D, Precision> &points) const;

    bool isNeighbor(const dSimplex<D, Precision> &other) const {
        PROFILER_INC("dSimplex_isNeighbor");

        uint sharedVertices = countSharedVertices(other);

        ASSERT(sharedVertices <= D || id == other.id);

        return sharedVertices == D;
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
    uint id;
    std::array<tIdType, D + 1> vertices;
    tHashType vFingerprint;
    std::array<tIdType, D + 1> neighbors;

public:
    static bool isFinite(const uint &i) {
        PROFILER_INC("dSimplex_staticIsFinite");

        return i != cINF;
    }

    static constexpr uint cINF = ~(0);
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

typedef GrowingHashTable<Concurrent_LP_MultiMap> cWuFaces;
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

//forward declare
struct PartialTriangulation;

template<uint D, typename Precision>
class dSimplices : public Concurrent_BlockedArray<dSimplex<D, Precision>> {

public:
    dSimplices() : Concurrent_BlockedArray<dSimplex<D, Precision>>(),
                   tetrahedronID(1) { }

    dSimplices(dSimplices &&other) : Concurrent_BlockedArray<dSimplex<D, Precision>>(std::move(other)),
                                     tetrahedronID(other.tetrahedronID.load()) { }

    VerificationReport<D, Precision>
            verify(const PartialTriangulation &pt,
                   const dPoints<D, Precision> &points) const;

    CrossCheckReport<D, Precision>
            crossCheck(const PartialTriangulation &pt,
                       const dSimplices<D, Precision> &realSimplices,
                       const PartialTriangulation &realPT) const;

public:
    std::atomic<uint> tetrahedronID;

};

template<uint D, typename Precision>
struct CrossCheckReport {
    bool valid;
    std::vector<dSimplex<D, Precision>> missing;
    std::vector<dSimplex<D, Precision>> invalid;
};

template<uint D, typename Precision>
struct VerificationReport {
    bool valid;
    std::map<dSimplex<D, Precision>, Ids> inCircle;
    std::vector<std::pair<dSimplex<D, Precision>, dSimplex<D, Precision>>>
            wrongNeighbors;
};

//#############################################################################
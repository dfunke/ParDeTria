/*
 * Geometry.h
 *
 *  Created on: Oct 28, 2014
 *      Author: dfunke
 */

#pragma once

#include <iostream>
#include <unordered_set>
#include <map>
#include <array>
#include <cmath>

#include <tbb/spin_mutex.h>
#include <tbb/spin_rw_mutex.h>

#include "utils/IndexedVector.hxx"
#include "utils/VectorAdapter.hxx"
#include "utils/Logger.h"
#include "utils/ASSERT.h"

typedef std::unordered_set<uint> Ids;

// dimensionality of our problem
// const uint D = 2;

// basic data structure for d-dimensional data
template <std::size_t D, typename Precision>
using dVector = std::array<Precision, D>;

//#############################################################################

// sphere
template <uint D, typename Precision> struct dSphere {
  dVector<D, Precision> center;
  Precision radius;
};

//#############################################################################

// (hyper)-rectangle
template <uint D, typename Precision> struct dBox {
  dVector<D, Precision> low;
  dVector<D, Precision> high;

  bool contains(const dVector<D, Precision> &p) const {
    for (uint d = 0; d < D; ++d) {
      if (!(low[d] <= p[d] && p[d] <= high[d]))
        return false;
    }

    return true;
  }

  /* tests whether sphere is FULLY contained in box */
  bool contains(const dSphere<D, Precision> &sphere) const {

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

    Precision r2 = sphere.radius * sphere.radius;

    Precision dist = 0;
    for (uint d = 0; d < D; ++d) {
      Precision e = std::max(low[d] - sphere.center[d], (Precision)0) +
                    std::max(sphere.center[d] - high[d], (Precision)0);

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

template <uint D, typename Precision> class dPoint {

public:
  dPoint() : id(dPoint<D, Precision>::cINF) {}

  dPoint(const dVector<D, Precision> &_coords)
      : id(dPoint<D, Precision>::cINF), coords(_coords) {}

  dPoint(const dPoint<D, Precision> &a) : id(a.id), coords(a.coords) {}

  dPoint &operator=(const dPoint<D, Precision> &a) {
    id = a.id;
    coords = a.coords;

    return *this;
  }

  bool operator==(const dPoint<D, Precision> &a) const {

    // COUT << "Comparing POINTS THIS " << *this << " and OTHER " << a << ": ";

    if (!(isFinite() ^ a.isFinite())) {
      // either none or both points are infinity
      // compare ids
      // std::cout << (this->id == a.id) << std::endl);
      return id == a.id;
    } else {
      // either one is infinity
      // compare coordinates
      for (uint i = 0; i < D; ++i) {
        if (coords[i] != a.coords[i]) {
          // std::cout << false << std::endl);
          return false;
        }
      }
      // std::cout << true << std::endl);
      return true;
    }
  }

  bool operator==(const uint &a) const { return id == a; }

  inline bool isFinite() const { return isFinite(id); }

public:
  uint id;
  dVector<D, Precision> coords;

public:
  static inline bool isFinite(const uint &i) {
    // INDENT
    // PLOG(i << " -> " << (i & cINF) << " != " << cINF << ": " << !((i &
    // cINF) == cINF) << std::endl);
    // DEDENT
    return !((i & cINF) == cINF);
  }

  static inline uint infIndex(const uint &i) {
    ASSERT(!isFinite(i));
    return ~cINF & i;
  }

  static constexpr uint cINF = ~(0) ^ ((1 << D) - 1);
  static constexpr uint nINF = ~(0) - cINF + 1;
};

//#############################################################################

// IndexedVector of dPoints

template <uint D, typename Precision>
class dPoints : public VectorAdapter<dPoint<D, Precision>> {

public:
  dPoints() : VectorAdapter<dPoint<D, Precision>>() {}

  dPoints(const VectorAdapter<dPoint<D, Precision>> &other)
      : VectorAdapter<dPoint<D, Precision>>(other) {}

  bool operator==(const Ids &other) const {
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

template <uint D, typename Precision> class dSimplex {

public:
  dSimplex() : id(dSimplex<D, Precision>::cINF) {}

  dSimplex(const std::array<uint, D + 1> &_vertices)
      : id(dSimplex<D, Precision>::cINF), vertices(_vertices) {}

  dSimplex(const dSimplex<D, Precision> &a)
      : id(a.id), vertices(a.vertices), vertexFingerprint(a.vertexFingerprint),
        neighbors(a.neighbors) {}

  dSimplex &operator=(const dSimplex<D, Precision> &a) {
    // acquire lock, as we are modifying the point
    tbb::spin_mutex::scoped_lock lock(mtx);

    id = a.id;
    vertices = a.vertices;
    vertexFingerprint = a.vertexFingerprint;
    neighbors = a.neighbors;

    return *this;
  }

  bool operator==(const dSimplex<D, Precision> &a) const {

    // COUT << "Comparing SIMPLICES THIS " << *this << " and OTHER " << a << ":
    // ";

    if (!((id == cINF) ^ (a.id == cINF))) {
      // either none or both simplices are infinity
      // compare ids
      // std::cout << (this->id == a.id) << std::endl);
      return id == a.id;
    } else {
      // either one is infinity
      // compare vertices
      return equalVertices(a);
    }
  }

  // for use as map key
  bool operator<(const dSimplex<D, Precision> &a) const { return id < a.id; }

  bool equalVertices(const dSimplex<D, Precision> &a) const {
    // compare vertices
    // both vertex arrays are sorted by vertex id

    for (uint i = 0; i < D + 1; ++i) {
      if (vertices[i] != a.vertices[i]) {
        return false;
      }
    }

    return true;
  }

  uint countSharedVertices(const dSimplex<D, Precision> &a) const {
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
    return neighbors != a.neighbors;
  }

  bool operator==(const uint &a) const { return id == a; }

  bool contains(uint p) const {
    for (uint d = 0; d < D + 1; ++d) {
      if (p == vertices[d])
        return true;
    }

    return false;
  }

  bool contains(const dPoint<D, Precision> &p) const {
    for (uint d = 0; d < D + 1; ++d) {
      if (p == vertices[d])
        return true;
    }

    return false;
  }

  template <typename Container> bool containsAll(const Container &p) const {
    for (const auto &x : p) {
      if (!contains(x))
        return false;
    }

    return true;
  }

  template <typename Container> bool containsAny(const Container &p) const {
    for (const auto &x : p) {
      if (contains(x))
        return true;
    }

    return false;
  }

  bool isFinite() const {
    bool finite = true;

    for (uint i = 0; i < D + 1; ++i) {
      if (!dPoint<D, Precision>::isFinite(vertices[i]))
        finite = false;
    }

    return finite;
  }

  // dimension specific implementations in cpp file
  Precision orientation(const dPoints<D, Precision> &points) const;
  bool inSphere(const dPoint<D, Precision> &p,
                const dPoints<D, Precision> &points) const;
  dSphere<D, Precision> circumsphere(const dPoints<D, Precision> &points) const;

  bool isNeighbor(const dSimplex<D, Precision> &other) const {
    uint sharedVertices = countSharedVertices(other);

    ASSERT(sharedVertices <= D || id == other.id);

    // if(sharedVertices > 2 && id != other.id)
    //	PLOG("Distinct vertices " << *this " and " << other << " share more
    // than 2 hits" << std::endl);

    return sharedVertices == D;
  }

  uint fingerprint() {
    vertexFingerprint = 0;
    for (uint i = 0; i < D + 1; ++i) {
      vertexFingerprint ^= vertices[i];
    }

    return vertexFingerprint;
  }

public:
  uint id;
  std::array<uint, D + 1> vertices;
  uint vertexFingerprint;
  Ids neighbors;
  tbb::spin_mutex mtx;

public:
  static bool isFinite(const uint &i) { return i != cINF; }

  static constexpr uint cINF = ~(0);
};

template <std::size_t D, typename Precision>
std::string to_string(const dVector<D, Precision> &p);

template <uint D, typename Precision>
std::string to_string(const dPoint<D, Precision> &p);

template <uint D, typename Precision>
std::string to_string(const dSimplex<D, Precision> &p);

template <uint D, typename Precision>
std::string to_string(const dSphere<D, Precision> &p);

template <uint D, typename Precision>
std::string to_string(const dBox<D, Precision> &b);

template <std::size_t D, typename Precision>
std::ostream &operator<<(std::ostream &o, const dVector<D, Precision> &p);

template <uint D, typename Precision>
std::ostream &operator<<(std::ostream &o, const dPoint<D, Precision> &p);

template <uint D, typename Precision>
std::ostream &operator<<(std::ostream &o, const dSphere<D, Precision> &p);

template <uint D, typename Precision>
std::ostream &operator<<(std::ostream &o, const dSimplex<D, Precision> &p);

template <uint D, typename Precision>
std::ostream &operator<<(std::ostream &o, const dBox<D, Precision> &b);

template <uint D, typename Precision> struct CrossCheckReport;

template <uint D, typename Precision> struct VerificationReport;

template <uint D, typename Precision>
class dSimplices : public IndexedVector<dSimplex<D, Precision>> {

public:
  dSimplices() : IndexedVector<dSimplex<D, Precision>>() {}

  dSimplices(const IndexedVector<dSimplex<D, Precision>> &other)
      : IndexedVector<dSimplex<D, Precision>>(other) {}

  bool operator==(const dSimplices<D, Precision> &other) const {
    if (this->size() != other.size()) {
      PLOG("my size: " << this->size() << " other size: " << other.size()
                       << std::endl);
      return false;
    }

    for (const auto &otherSimplex : other) {
      // find my simplex, compares simplex id or vertices ids
      auto mySimplex = std::find(this->begin(), this->end(), otherSimplex);

      if (mySimplex == this->end()) {
        PLOG("did not find simplex" << otherSimplex << std::endl);
        return false;
      }

      // check neighbors
      for (const auto &n : mySimplex->neighbors) {
        // TODO handle infinite neighbors better
        if (dSimplex<D, Precision>::isFinite(n) &&
            otherSimplex.neighbors.count(n) != 1) {
          // the other triangulation does not contain n as neighbor
          PLOG("wrong neighbors: " << *mySimplex << " -- " << otherSimplex
                                   << std::endl);
          return false;
        }
      }
    }

    return true;
  }

  bool operator!=(const dSimplices<D, Precision> &other) const {
    return !operator==(other);
  }

  VerificationReport<D, Precision>
  verify(const Ids &partitionPoints, const dPoints<D, Precision> &points) const;

  CrossCheckReport<D, Precision>
  crossCheck(const dSimplices<D, Precision> &realDT) const;

  uint countDuplicates() const;

  template <typename Container>
  dSimplices<D, Precision> findSimplices(const Container &points,
                                         const bool all = false) const {

    dSimplices<D, Precision> result;

    if (all) {
      for (const auto &s : *this) {
        if (s.containsAll(points))
          result.insert(s);
      }
    } else {
      for (const auto &s : *this) {
        if (s.containsAny(points))
          result.insert(s);
      }
    }

    return result;
  }

public:
  tbb::concurrent_unordered_map<uint, tbb::concurrent_vector<uint>> whereUsed;
};

template <uint D, typename Precision> struct CrossCheckReport {
  bool valid;
  dSimplices<D, Precision> missing;
  dSimplices<D, Precision> invalid;
};

template <uint D, typename Precision> struct VerificationReport {
  bool valid;
  std::map<dSimplex<D, Precision>, Ids> inCircle;
  std::vector<std::pair<dSimplex<D, Precision>, dSimplex<D, Precision>>>
      wrongNeighbors;
};

//#############################################################################

/*
 * Geometry.h
 *
 *  Created on: Oct 28, 2014
 *      Author: dfunke
 */

#pragma once

#include <iostream>
#include <set>
#include <map>
#include <array>
#include <cmath>

#include "utils/IndexedVector.hxx"
#include "utils/Logger.h"
#include "utils/ASSERT.h"

typedef float tCoordinate;
typedef std::set<uint> Ids;

// dimensionality of our problem
// const uint D = 2;

// basic data structure for d-dimensional data
template <uint D> using dVector = std::array<tCoordinate, D>;

//#############################################################################

// sphere
template <uint D> struct dSphere {
  dVector<D> center;
  tCoordinate radius;
};

//#############################################################################

// (hyper)-rectangle
template <uint D> struct dBox {
  dVector<D> coords;
  dVector<D> dim;

  /* tests whether sphere is FULLY contained in box */
  bool contains(const dSphere<D> &sphere) const {

    tCoordinate r2 = sphere.radius * sphere.radius;

    for (uint d = 0; d < D; ++d) {
      if (!(coords[d] <= sphere.center[d] &&
            sphere.center[d] <= coords[d] + dim[d])) {
        return false;
      }
    }

    // the center of the sphere is within the box
    for (uint d = 0; d < D; ++d) {
      auto p = sphere.center;
      // project p to the boundary of box in dimension d closest to center of
      // the
      // sphere
      p[d] = sphere.center[d] < coords[d] + (dim[d] / 2) ? coords[d]
                                                         : coords[d] + dim[d];

      tCoordinate dist = 0;
      for (uint i = 0; i < D; ++i)
        dist += (sphere.center[i] - p[i]) * (sphere.center[i] - p[i]);

      if (dist < r2)
        return false;
    }

    return true;
  }
};

//#############################################################################

// point with ID and simplices list

template <uint D> class dPoint {

public:
  bool operator==(const dPoint<D> &a) const {

    // COUT << "Comparing POINTS THIS " << *this << " and OTHER " << a << ": ";

    if (!(isFinite() ^ a.isFinite())) {
      // either none or both points are infinity
      // compare ids
      // std::cout << (this->id == a.id) << std::endl;
      return id == a.id;
    } else {
      // either one is infinity
      // compare coordinates
      for (uint i = 0; i < D; ++i) {
        if (coords[i] != a.coords[i]) {
          // std::cout << false << std::endl;
          return false;
        }
      }
      // std::cout << true << std::endl;
      return true;
    }
  }

  bool operator==(const uint &a) const { return id == a; }

  inline bool isFinite() const { return isFinite(id); }

public:
  uint id;
  dVector<D> coords;
  Ids simplices;

public:
  static inline bool isFinite(const uint &i) {
    // INDENT
    // PLOG << i << " -> " << (i & cINF) << " != " << cINF << ": " << !((i &
    // cINF) == cINF) << std::endl;
    // DEDENT
    return !((i & cINF) == cINF);
  }

  static constexpr uint cINF = ~(0) ^ ((1 << D) - 1);
};

//#############################################################################

// IndexedVector of dPoints

template <uint D> class dPoints : public IndexedVector<dPoint<D>> {

public:
  dPoints() : IndexedVector<dPoint<D>>() {}

  dPoints(const IndexedVector<dPoint<D>> &other)
      : IndexedVector<dPoint<D>>(other) {}

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

template <uint D> class dSimplex {

public:
  bool operator==(const dSimplex<D> &a) const {

    // COUT << "Comparing SIMPLICES THIS " << *this << " and OTHER " << a << ":
    // ";

    if (!((id == cINF) ^ (a.id == cINF))) {
      // either none or both simplices are infinity
      // compare ids
      // std::cout << (this->id == a.id) << std::endl;
      return id == a.id;
    } else {
      // either one is infinity
      // compare vertices
      for (uint i = 0; i < D + 1; ++i) {
        bool found = false;
        for (uint j = 0; j < D + 1; ++j) {
          if (vertices[i] == a.vertices[j]) {
            // std::cout << this->vertices[i] << " == " << a.vertices[j];
            found = true;
            break;
          }
        }
        if (!found) {
          // std::cout << false << std::endl;
          return false;
        }
      }
      // std::cout << true << std::endl;
      return true;
    }
  }

  // for use as map key
  bool operator<(const dSimplex<D> &a) const { return id < a.id; }

  bool equalVertices(const dSimplex<D> &a) const {
    // compare vertices
    for (uint i = 0; i < D + 1; ++i) {
      bool found = false;
      for (uint j = 0; j < D + 1; ++j) {
        if (vertices[i] == a.vertices[j]) {
          // std::cout << this->vertices[i] << " == " << a.vertices[j];
          found = true;
          break;
        }
      }
      if (!found) {
        // std::cout << false << std::endl;
        return false;
      }
    }

    return true;
  }

  uint countSharedVertices(const dSimplex<D> &a) const {
    uint count = 0;
    for (uint i = 0; i < D + 1; ++i) {
      bool found = false;
      for (uint j = 0; j < D + 1; ++j) {
        if (vertices[i] == a.vertices[j]) {
          // std::cout << this->vertices[i] << " == " << a.vertices[j];
          found = true;
          break;
        }
      }
      count += found;
    }

    return count;
  }

  bool equalNeighbors(const dSimplex<D> &a) const {
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

  bool contains(const dPoint<D> &p) const {
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
      if (!dPoint<D>::isFinite(vertices[i]))
        finite = false;
    }

    return finite;
  }

  // dimension specific implementations in cpp file
  tCoordinate orientation(const dPoints<D> &points) const;
  bool inSphere(const dPoint<D> &p, const dPoints<D> &points) const;
  dSphere<D> circumsphere(const dPoints<D> &points) const;

  bool isNeighbor(const dSimplex<D> &other) const {
    uint sharedVertices = 0;

    for (uint i = 0; i < D + 1; ++i) {
      for (uint j = 0; j < D + 1; ++j) {
        sharedVertices += (vertices[i] == other.vertices[j]);
      }
    }

    ASSERT(sharedVertices <= D || id == other.id);

    // if(sharedVertices > 2 && id != other.id)
    //	PLOG << "Distinct vertices " << *this " and " << other << " share more
    // than 2 hits";

    return sharedVertices == D;
  }

public:
  uint id;
  std::array<uint, D + 1> vertices;
  Ids neighbors;

public:
  static bool isFinite(const uint &i) { return i != cINF; }

  static constexpr uint cINF = ~(0);
};

template <uint D> std::string to_string(const dPoint<D> &p);

template <uint D> std::string to_string(const dSimplex<D> &p);

template <uint D> std::string to_string(const dSphere<D> &p);

template <uint D> std::string to_string(const dBox<D> &b);

template <uint D> std::ostream &operator<<(std::ostream &o, const dPoint<D> &p);

template <uint D>
std::ostream &operator<<(std::ostream &o, const dSphere<D> &p);

template <uint D>
std::ostream &operator<<(std::ostream &o, const dSimplex<D> &p);

template <uint D> std::ostream &operator<<(std::ostream &o, const dBox<D> &b);

template <uint D> struct CrossCheckReport;

template <uint D> struct VerificationReport;

template <uint D> class dSimplices : public IndexedVector<dSimplex<D>> {

public:
  dSimplices() : IndexedVector<dSimplex<D>>() {}

  dSimplices(const IndexedVector<dSimplex<D>> &other)
      : IndexedVector<dSimplex<D>>(other) {}

  bool operator==(const dSimplices<D> &other) const {
    if (this->size() != other.size()) {
      PLOG << "my size: " << this->size() << " other size: " << other.size()
           << std::endl;
      return false;
    }

    for (const auto &otherSimplex : other) {
      // find my simplex, compares simplex id or vertices ids
      auto mySimplex = std::find(this->begin(), this->end(), otherSimplex);

      if (mySimplex == this->end()) {
        PLOG << "did not find simplex" << otherSimplex << std::endl;
        return false;
      }

      // check neighbors
      for (const auto &n : mySimplex->neighbors) {
        // TODO handle infinite neighbors better
        if (dSimplex<D>::isFinite(n) && otherSimplex.neighbors.count(n) != 1) {
          // the other triangulation does not contain n as neighbor
          PLOG << "wrong neighbors: " << *mySimplex << " -- " << otherSimplex
               << std::endl;
          return false;
        }
      }
    }

    return true;
  }

  bool operator!=(const dSimplices<D> &other) const {
    return !operator==(other);
  }

  VerificationReport<D> verify(const dPoints<D> &points) const;

  CrossCheckReport<D> crossCheck(const dSimplices<D> &realDT) const;

  uint countDuplicates() const {
    uint duplicates = 0;

    for (const auto &s : *this) {
      auto simplices = findSimplices(s.vertices, true);
      duplicates += simplices.size() - 1;
    }

    PLOG << "Found " << duplicates << " duplicates" << std::endl;

    return duplicates;
  }

  template <typename Container>
  dSimplices<D> findSimplices(const Container &points,
                              const bool all = false) const {

    dSimplices<D> result;

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
};

template <uint D> struct CrossCheckReport {
  bool valid;
  dSimplices<D> missing;
  dSimplices<D> invalid;
};

template <uint D> struct VerificationReport {
  bool valid;
  std::map<dSimplex<D>, Ids> inCircle;
};

//#############################################################################

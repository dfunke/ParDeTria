#pragma once

#include "Geometry.h"

#include "utils/ASSERT.h"

template <uint D> class Partition {

public:
  bool contains(const uint &p) const { return points.count(p) == 1; }

  bool contains(const dPoint<D> &p) const { return points.count(p.id) == 1; }

  bool contains(const dSimplex<D> &s, bool partially = false) const {
    for (const auto &p : s.vertices) {
      // store contains result
      bool t = contains(p);

      if (!partially && !t)
        // we are looking for complete containment, but one point wasn't found
        return false;
      if (partially && t)
        // we are looking for partial containment and found a point
        return true;
    }

    return !partially;
  }

public:
  uint id;
  Ids points;
  dBox<D> bounds;
};

template <uint D> std::string to_string(const Partition<D> &p);

template <uint D>
std::ostream &operator<<(std::ostream &o, const Partition<D> &p);

template <uint D> class Partitioning : public IndexedVector<Partition<D>> {

public:
  uint partition(const uint p) const {
    for (const auto &part : *this) {
      if (part.points.count(p) > 0)
        return part.id;
    }

    throw std::out_of_range("Partition of point " + std::to_string(p) +
                            "not found");
  }

  uint partition(const dPoint<D> &p) const { return partition(p.id); }
};

template <uint D> struct dPointStats {
  dPoint<D> min;
  dPoint<D> mid;
  dPoint<D> max;
};

template <uint D, class InputIt>
dPointStats<D> getPointStats(const InputIt &first, const InputIt &last,
                             const dPoints<D> &points,
                             const bool ignoreInfinite = true) {

  dPointStats<D> stats;

  for (uint dim = 0; dim < D; ++dim) {
    stats.min.coords[dim] = std::numeric_limits<tCoordinate>::max();
    stats.max.coords[dim] = std::numeric_limits<tCoordinate>::min();

    for (auto it = first; it != last; ++it) {

      uint id = *it;
      ASSERT(points.contains(id));

      if (points[id].isFinite() || !ignoreInfinite) {
        stats.min.coords[dim] =
            std::min(stats.min.coords[dim], points[id].coords[dim]);
        stats.max.coords[dim] =
            std::max(stats.max.coords[dim], points[id].coords[dim]);
      }
    }

    stats.mid.coords[dim] = (stats.max.coords[dim] + stats.min.coords[dim]) / 2;
  }

  return stats;
}

template <uint D> class Partitioner {

public:
  virtual ~Partitioner() {}

  virtual Partitioning<D> partition(const Ids &ids, const dPoints<D> &points,
                                    const std::string &provenance) const = 0;
};

template <uint D> class dPartitioner : public Partitioner<D> {

public:
  Partitioning<D> partition(const Ids &ids, const dPoints<D> &points,
                            const std::string &provenance) const;
};

template <uint D> class kPartitioner : public Partitioner<D> {

public:
  kPartitioner(uint _k) : k(_k) {}

  Partitioning<D> partition(const Ids &ids, const dPoints<D> &points,
                            const std::string &provenance) const;

private:
  uint k; // dimension to partition
};

template <uint D> class CyclePartitioner : public Partitioner<D> {

public:
  Partitioning<D> partition(const Ids &ids, const dPoints<D> &points,
                            const std::string &provenance) const;
};

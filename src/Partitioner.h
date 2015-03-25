#pragma once

#include "Geometry.h"

#include "utils/ASSERT.h"

template <uint D, typename Precision> class Partition {

public:
  bool contains(const uint &p) const { return points.count(p) == 1; }

  bool contains(const dPoint<D, Precision> &p) const {
    assert(points.count(p.id) == 0 || bounds.contains(p.coords));
    return points.count(p.id) == 1 && bounds.contains(p.coords);
  }

  bool contains(const dSimplex<D, Precision> &s, bool partially = false) const {
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
  dBox<D, Precision> bounds;
};

template <uint D, typename Precision>
std::string to_string(const Partition<D, Precision> &p);

template <uint D, typename Precision>
std::ostream &operator<<(std::ostream &o, const Partition<D, Precision> &p);

template <uint D, typename Precision>
class Partitioning : public IndexedVector<Partition<D, Precision>> {

public:
  uint partition(const uint p) const {
    for (const auto &part : *this) {
      if (part.points.count(p) > 0)
        return part.id;
    }

    throw std::out_of_range("Partition of point " + std::to_string(p) +
                            "not found");
  }

  uint partition(const dPoint<D, Precision> &p) const {
    return partition(p.id);
  }
};

template <uint D, typename Precision> struct dPointStats {
  dVector<D, Precision> min;
  dVector<D, Precision> mid;
  dVector<D, Precision> max;
};

template <uint D, typename Precision>
dPointStats<D, Precision> getPointStats(const uint &first, const uint &last,
                                        const dPoints<D, Precision> &points,
                                        const bool ignoreInfinite = true) {

  dPointStats<D, Precision> stats;

  for (uint dim = 0; dim < D; ++dim) {
    stats.min[dim] = std::numeric_limits<Precision>::max();
    stats.max[dim] = std::numeric_limits<Precision>::min();

    for (uint id = first; id < last; ++id) {

      ASSERT(points.contains(id));

      if (dPoint<D, Precision>::isFinite(id) || !ignoreInfinite) {
        stats.min[dim] = std::min(stats.min[dim], points[id].coords[dim]);
        stats.max[dim] = std::max(stats.max[dim], points[id].coords[dim]);
      }
    }

    stats.mid[dim] = (stats.max[dim] + stats.min[dim]) / 2;
  }

  return stats;
}

template <uint D, typename Precision, class InputIt>
dPointStats<D, Precision> getPointStats(const InputIt &first,
                                        const InputIt &last,
                                        const dPoints<D, Precision> &points,
                                        const bool ignoreInfinite = true) {

  dPointStats<D, Precision> stats;

  for (uint dim = 0; dim < D; ++dim) {
    stats.min[dim] = std::numeric_limits<Precision>::max();
    stats.max[dim] = std::numeric_limits<Precision>::min();

    for (auto it = first; it != last; ++it) {

      uint id = *it;
      ASSERT(points.contains(id));

      if (dPoint<D, Precision>::isFinite(id) || !ignoreInfinite) {
        stats.min[dim] = std::min(stats.min[dim], points[id].coords[dim]);
        stats.max[dim] = std::max(stats.max[dim], points[id].coords[dim]);
      }
    }

    stats.mid[dim] = (stats.max[dim] + stats.min[dim]) / 2;
  }

  return stats;
}

template <uint D, typename Precision> class Partitioner {

public:
  virtual ~Partitioner() {}

  virtual Partitioning<D, Precision>
  partition(const Ids &ids, const dPoints<D, Precision> &points,
            const std::string &provenance) const = 0;
};

template <uint D, typename Precision>
class dPartitioner : public Partitioner<D, Precision> {

public:
  Partitioning<D, Precision> partition(const Ids &ids,
                                       const dPoints<D, Precision> &points,
                                       const std::string &provenance) const;
};

template <uint D, typename Precision>
class kPartitioner : public Partitioner<D, Precision> {

public:
  kPartitioner(uint _k) : k(_k) {}

  Partitioning<D, Precision> partition(const Ids &ids,
                                       const dPoints<D, Precision> &points,
                                       const std::string &provenance) const;

private:
  uint k; // dimension to partition
};

template <uint D, typename Precision>
class CyclePartitioner : public Partitioner<D, Precision> {

public:
  Partitioning<D, Precision> partition(const Ids &ids,
                                       const dPoints<D, Precision> &points,
                                       const std::string &provenance) const;
};

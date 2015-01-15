#pragma once

#include "Geometry.h"

template <uint D> class Partition {

public:
  uint id;
  Ids points;
  dBox<D> bounds;
};

template <uint D> std::string to_string(const Partition<D> &p);

template <uint D>
std::ostream &operator<<(std::ostream &o, const Partition<D> &p);

template <uint D> using Partitioning = IndexedVector<Partition<D>>;

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
      assert(points.contains(id));

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

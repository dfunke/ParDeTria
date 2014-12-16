#pragma once

#include "Geometry.h"

class Partition {

public:
  uint id;
  Ids points;
  dBox bounds;
};

std::string to_string(const Partition &p);
std::ostream &operator<<(std::ostream &o, const Partition &p);

typedef IndexedVector<Partition> Partitioning;

struct dPointStats {
  dPoint min;
  dPoint mid;
  dPoint max;
};

template <class InputIt>
dPointStats getPointStats(const InputIt &first, const InputIt &last,
                          const dPoints &points,
                          const bool ignoreInfinite = true) {

  dPointStats stats;

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

class Partitioner {

public:
  virtual ~Partitioner() {}

  virtual Partitioning partition(const Ids &ids, const dPoints &points,
                                 const std::string &provenance) const = 0;
};

class dPartitioner : public Partitioner {

public:
  Partitioning partition(const Ids &ids, const dPoints &points,
                         const std::string &provenance) const;
};

class kPartitioner : public Partitioner {

public:
  kPartitioner(uint _k) : k(_k) {}

  Partitioning partition(const Ids &ids, const dPoints &points,
                         const std::string &provenance) const;

private:
  uint k; // dimension to partition
};

class CyclePartitioner : public Partitioner {

public:
  Partitioning partition(const Ids &ids, const dPoints &points,
                         const std::string &provenance) const;
};

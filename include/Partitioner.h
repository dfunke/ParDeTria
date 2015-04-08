#pragma once

#include "Geometry.h"

#include <type_traits>

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

namespace __detail {
template <class T,
          typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
uint _id(const T &t) {
  return t;
}

template <class T,
          typename std::enable_if<!std::is_integral<T>::value, int>::type = 0>
uint _id(const T &t) {
  return *t;
}
}

template <uint D, typename Precision, typename InputIt>
dPointStats<D, Precision> getPointStats(const InputIt &first,
                                        const InputIt &last,
                                        const dPoints<D, Precision> &points,
                                        const bool ignoreInfinite = true) {

  dPointStats<D, Precision> stats;

  for (uint dim = 0; dim < D; ++dim) {
    stats.min[dim] = std::numeric_limits<Precision>::max();
    stats.max[dim] = std::numeric_limits<Precision>::min();

    for (auto it = first; it != last; ++it) {

      const uint id = __detail::_id(it);
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

public:
    static std::unique_ptr<Partitioner<D, Precision>> && make(const unsigned char type);
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


template<uint D, typename Precision>
std::unique_ptr<Partitioner<D, Precision>> && Partitioner<D, Precision>::make(const unsigned char type) {

  std::unique_ptr<Partitioner<D, Precision>> partitioner_ptr;
  switch (type) {
    case 'd':
      partitioner_ptr = std::make_unique<dPartitioner<D, Precision>>();
          break;
    case 'c':
      partitioner_ptr = std::make_unique<CyclePartitioner<D, Precision>>();
          break;
    default:
      // p must be a dimension - subtract '0' to get integer value
      int d = type - '0';
          ASSERT(0 <= d && d < D);
          partitioner_ptr = std::make_unique<kPartitioner<D, Precision>>(d);
          break;
  }

  return std::move(partitioner_ptr);
}

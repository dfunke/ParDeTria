#pragma once

#include <memory>

#include "Geometry.h"

#include <type_traits>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>


#include "utils/ASSERT.h"

template<uint D, typename Precision>
class Partition {

public:
    Partition(const std::size_t &nPoints) : points(nPoints) { }
    Partition(Point_Ids &&_points) : points(std::move(_points)) { }

    Partition(Partition &&other)
            : id(other.id),
              points(std::move(other.points)),
              bounds(other.bounds) { }

public:
    bool contains(const tIdType &p) const { return points.count(p) == 1; }

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
    Point_Ids points;
    dBox<D, Precision> bounds;
};

template<uint D, typename Precision>
std::string to_string(const Partition<D, Precision> &p);

template<uint D, typename Precision>
std::ostream &operator<<(std::ostream &o, const Partition<D, Precision> &p);

template<uint D, typename Precision>
class Partitioning : public std::vector<Partition<D, Precision>> {

public:
    uint partition(const tIdType p) const {
        for (const auto &part : *this) {
            if (part.points.count(p) > 0)
                return part.id;
        }

        throw std::out_of_range("Partition of point " + std::to_string(p) + " not found");
    }
};

template<uint D, typename Precision>
struct dPointStats {

    dPointStats() {

        for (uint dim = 0; dim < D; ++dim) {
            min[dim] = std::numeric_limits<Precision>::max();
            mid[dim] = 0;
            max[dim] = std::numeric_limits<Precision>::min();
        }
    }

    dVector<D, Precision> min;
    dVector<D, Precision> mid;
    dVector<D, Precision> max;

    bool operator==(const dPointStats<D, Precision> & other) const {
        return min == other.min && mid == other.mid && max == other.max;
    }
};

template<uint D, typename Precision>
std::string to_string(const dPointStats<D, Precision> &p);

template<uint D, typename Precision>
std::ostream &operator<<(std::ostream &o, const dPointStats<D, Precision> &p);

namespace __detail {
    template<class T,
            typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    tIdType _id(const T &t) {
        return t;
    }

    template<class T,
            typename std::enable_if<!std::is_integral<T>::value, int>::type = 0>
    tIdType _id(const T &t) {
        return *t;
    }
}

template<uint D, typename Precision, typename InputIt>
dPointStats<D, Precision> getPointStatsSeq(const InputIt &first,
                                           const InputIt &last,
                                           const dPoints<D, Precision> &points,
                                           const bool ignoreInfinite = true) {

    dPointStats<D, Precision> stats;

    for (uint dim = 0; dim < D; ++dim) {
        stats.min[dim] = std::numeric_limits<Precision>::max();
        stats.max[dim] = std::numeric_limits<Precision>::min();

        for (auto it = first; it != last; ++it) {

            const tIdType id = __detail::_id(it);
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

template<uint D, typename Precision, typename Container>
dPointStats<D, Precision> getPointStats(const Container &ids,
                                        const dPoints<D, Precision> &points,
                                        const bool ignoreInfinite = true) {

    tbb::enumerable_thread_specific<dPointStats<D, Precision>,
            tbb::cache_aligned_allocator<dPointStats<D, Precision>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsPointStats;

    tbb::parallel_for(ids.range(), [&tsPointStats, &points, ignoreInfinite](const auto &r) {

        auto &stats = tsPointStats.local();

        for (auto it = r.begin(); it != r.end(); ++it) {

            const tIdType id = __detail::_id(it);

            if (dPoint<D, Precision>::isFinite(id) || !ignoreInfinite) {

                ASSERT(id != 0);
                ASSERT(points.contains(id));
                const auto &coords = points[id].coords;

                for (uint dim = 0; dim < D; ++dim) {
                    stats.min[dim] = std::min(stats.min[dim], coords[dim]);
                    stats.max[dim] = std::max(stats.max[dim], coords[dim]);
                }
            }
        }
    });

    auto stats = tsPointStats.combine([](const dPointStats<D, Precision> &a, const dPointStats<D, Precision> &b) {
        dPointStats<D, Precision> s;

        for (uint dim = 0; dim < D; ++dim) {
            s.min[dim] = std::min(a.min[dim], b.min[dim]);
            s.max[dim] = std::max(a.max[dim], b.max[dim]);
        }

        return s;
    });

    for (uint dim = 0; dim < D; ++dim) {
        stats.mid[dim] = (stats.max[dim] + stats.min[dim]) / 2;
    }

    return stats;
}

template<uint D, typename Precision>
class Partitioner {

public:
    virtual ~Partitioner() = default;

    virtual Partitioning<D, Precision>
            partition(const dPointStats<D, Precision> &stats,
                      const Point_Ids &ids,
                      const dPoints<D, Precision> &points,
                      const std::string &provenance,
                      bool ignoreInfinte = false) const = 0;

    Partitioning<D, Precision>
    partition(const Point_Ids &ids,
              const dPoints<D, Precision> &points,
              const std::string &provenance,
              bool ignoreInfinte = false) const {

        return partition(getPointStats(ids, points),
                         ids, points, provenance, ignoreInfinte);
    }

public:
    static std::unique_ptr<Partitioner<D, Precision>> make(const unsigned char type);
};

template<uint D, typename Precision>
class dWayPartitioner : public Partitioner<D, Precision> {

public:
    Partitioning<D, Precision> partition(const dPointStats<D, Precision> &stats,
                                         const Point_Ids &ids,
                                         const dPoints<D, Precision> &points,
                                         const std::string &provenance,
                                         bool ignoreInfinte = false) const;
};

template<uint D, typename Precision>
class OneDimPartitioner : public Partitioner<D, Precision> {

public:
    OneDimPartitioner(uint _k) : k(_k) { }

    Partitioning<D, Precision> partition(const dPointStats<D, Precision> &stats,
                                         const Point_Ids &ids,
                                         const dPoints<D, Precision> &points,
                                         const std::string &provenance,
                                         bool ignoreInfinte = false) const;

private:
    uint k; // dimension to partition
};

template<uint D, typename Precision>
class CyclePartitioner : public Partitioner<D, Precision> {

public:
    Partitioning<D, Precision> partition(const dPointStats<D, Precision> &stats,
                                         const Point_Ids &ids,
                                         const dPoints<D, Precision> &points,
                                         const std::string &provenance,
                                         bool ignoreInfinte = false) const;
};


template<uint D, typename Precision>
std::unique_ptr<Partitioner<D, Precision>> Partitioner<D, Precision>::make(const unsigned char type) {

    std::unique_ptr<Partitioner<D, Precision>> partitioner_ptr;
    switch (type) {
        case 'd':
            partitioner_ptr = std::make_unique<dWayPartitioner<D, Precision>>();
            break;
        case 'c':
            partitioner_ptr = std::make_unique<CyclePartitioner<D, Precision>>();
            break;
        default:
            // p must be a dimension - subtract '0' to get integer value
            int d = type - '0';
            ASSERT(0 <= d && (uint) d < D);
            partitioner_ptr = std::make_unique<OneDimPartitioner<D, Precision>>(d);
            break;
    }

    ASSERT(partitioner_ptr);
    return partitioner_ptr;
}

template<uint D, typename Precision>
std::string to_string(const dPointStats<D, Precision> &p);

template<uint D, typename Precision>
std::ostream &operator<<(std::ostream &o, const dPointStats<D, Precision> &p);

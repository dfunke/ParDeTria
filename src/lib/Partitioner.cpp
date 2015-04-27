#include "Partitioner.h"

template<uint D, typename Precision>
Partitioning<D, Precision> dPartitioner<D, Precision>::partition(
        const Ids &ids, const dPoints<D, Precision> &points,
        __attribute((unused)) const std::string &provenance) const {
    // do mid-point based partitioning for now
    auto stats = getPointStats(ids.begin(), ids.end(), points);

    PLOG("Midpoint is " << stats.mid << std::endl);

    Partitioning<D, Precision> partitioning;
    partitioning.resize(pow(2, D));

    for (uint i = 0; i < pow(2, D); ++i) {
        partitioning[i].id = i;

        for (uint d = 0; d < D; ++d) {
            partitioning[i].bounds.low[d] =
                    i & (1 << d) ? stats.mid[d] : stats.min[d];
            partitioning[i].bounds.high[d] =
                    i & (1 << d) ? stats.max[d] : stats.mid[d];
        }
    }

    for (auto &id : ids) {

        if (!dPoint<D, Precision>::isFinite(id))
            continue; // skip infinite points, they will be added later

        ASSERT(points.contains(id));
        const auto &p = points[id];

        uint part = 0;
        for (uint dim = 0; dim < D; ++dim) {
            part |= (p.coords[dim] > stats.mid[dim]) << dim;
        }

        ASSERT(partitioning[part].bounds.contains(p.coords));

        // LOG("Adding " << p << " to " << part << std::endl);
        partitioning[part].points.insert(p.id);
    }

    // add infinite points
    for (uint i = 0; i < partitioning.size(); ++i) {
        for (uint k = dPoint<D, Precision>::cINF; k != 0; ++k) {
            partitioning[i].points.insert(k);
        }
    }

    return partitioning;
}

template<uint D, typename Precision>
Partitioning<D, Precision> kPartitioner<D, Precision>::partition(
        const Ids &ids, const dPoints<D, Precision> &points,
        __attribute((unused)) const std::string &provenance) const {
    // do mid-point based partitioning for now
    auto stats = getPointStats(ids.begin(), ids.end(), points);

    PLOG("Midpoint is " << stats.mid << std::endl);

    Partitioning<D, Precision> partitioning;
    partitioning.resize(2);

    partitioning[0].id = 0;
    for (uint d = 0; d < D; ++d) {
        partitioning[0].bounds.low[d] = stats.min[d];
        partitioning[0].bounds.high[d] = d == k ? stats.mid[d] : stats.max[d];
    }

    partitioning[1].id = 1;
    for (uint d = 0; d < D; ++d) {
        partitioning[1].bounds.low[d] = d == k ? stats.mid[d] : stats.min[d];
        partitioning[1].bounds.high[d] = stats.max[d];
    }

    for (auto &id : ids) {

        if (!dPoint<D, Precision>::isFinite(id))
            continue; // skip infinite points, they will be added later

        ASSERT(points.contains(id));
        const auto &p = points[id];

        uint part = (p.coords[k] > stats.mid[k]);

        ASSERT(partitioning[part].bounds.contains(p.coords));

        // LOG("Adding " << p << " to " << part << std::endl);
        partitioning[part].points.insert(p.id);
    }

    // add infinite points
    for (uint i = 0; i < partitioning.size(); ++i) {
        for (uint k = dPoint<D, Precision>::cINF; k != 0; ++k) {
            partitioning[i].points.insert(k);
        }
    }

    return partitioning;
}

template<uint D, typename Precision>
Partitioning<D, Precision>
CyclePartitioner<D, Precision>::partition(const Ids &ids,
                                          const dPoints<D, Precision> &points,
                                          const std::string &provenance) const {
    // do mid-point based partitioning for now
    auto stats = getPointStats(ids.begin(), ids.end(), points);

    // cycle is lenght of provenance - 1 modulo D
    uint k = (provenance.size() - 1) % D;

    PLOG("Midpoint is " << stats.mid << std::endl);
    PLOG("Splitting dimension is " << k << std::endl);

    Partitioning<D, Precision> partitioning;
    partitioning.resize(2);

    partitioning[0].id = 0;
    for (uint d = 0; d < D; ++d) {
        partitioning[0].bounds.low[d] = stats.min[d];
        partitioning[0].bounds.high[d] = d == k ? stats.mid[d] : stats.max[d];
    }

    partitioning[1].id = 1;
    for (uint d = 0; d < D; ++d) {
        partitioning[1].bounds.low[d] = d == k ? stats.mid[d] : stats.min[d];
        partitioning[1].bounds.high[d] = stats.max[d];
    }

    for (auto &id : ids) {
        if (!dPoint<D, Precision>::isFinite(id))
            continue; // skip infinite points, they will be added later

        ASSERT(points.contains(id));
        const auto &p = points[id];

        uint part = (p.coords[k] > stats.mid[k]);

        ASSERT(partitioning[part].bounds.contains(p.coords));

        // LOG("Adding " << p << " to " << part << std::endl);
        partitioning[part].points.insert(p.id);
    }

    // add infinite points
    for (uint i = 0; i < partitioning.size(); ++i) {
        for (uint k = dPoint<D, Precision>::cINF; k != 0; ++k) {
            partitioning[i].points.insert(k);
        }
    }

    return partitioning;
}

// specializations
template
class dPartitioner<2, float>;

template
class kPartitioner<2, float>;

template
class CyclePartitioner<2, float>;

template
class dPartitioner<3, float>;

template
class kPartitioner<3, float>;

template
class CyclePartitioner<3, float>;

template
class dPartitioner<2, double>;

template
class kPartitioner<2, double>;

template
class CyclePartitioner<2, double>;

template
class dPartitioner<3, double>;

template
class kPartitioner<3, double>;

template
class CyclePartitioner<3, double>;

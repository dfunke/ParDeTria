#include "Partitioner.h"

template<uint D, typename Precision>
Partitioning<D, Precision> dPartitioner<D, Precision>::partition(
        const Ids &ids, const dPoints<D, Precision> &points,
        __attribute((unused)) const std::string &provenance) const {
    // do mid-point based partitioning for now
    auto idsHandle = ids.handle();
    auto stats = getPointStats(idsHandle.begin(), idsHandle.end(), points);

    PLOG("Midpoint is " << stats.mid << std::endl);

    Partitioning<D, Precision> partitioning;
    partitioning.reserve(pow(2, D));

    std::vector<GrowingHashTableHandle<Concurrent_LP_Set>> pointHandles;
    pointHandles.reserve(partitioning.size());

    for (uint i = 0; i < pow(2, D); ++i) {
        partitioning.emplace_back(idsHandle.size() / (pow(2, D)));
        partitioning[i].id = i;

        for (uint d = 0; d < D; ++d) {
            partitioning[i].bounds.low[d] =
                    i & (1 << d) ? stats.mid[d] : stats.min[d];
            partitioning[i].bounds.high[d] =
                    i & (1 << d) ? stats.max[d] : stats.mid[d];
        }

        pointHandles.emplace_back(partitioning[i].points.handle());
    }

    for (auto id : idsHandle) {

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
        pointHandles[part].insert(id);
    }

    // add infinite points
    for (uint i = 0; i < partitioning.size(); ++i) {
        for (uint k = dPoint<D, Precision>::cINF; k != 0; ++k) {
            pointHandles[i].insert(k);
        }
    }

    return partitioning;
}

template<uint D, typename Precision>
Partitioning<D, Precision> kPartitioner<D, Precision>::partition(
        const Ids &ids, const dPoints<D, Precision> &points,
        __attribute((unused)) const std::string &provenance) const {
    // do mid-point based partitioning for now
    auto idsHandle = ids.handle();
    auto stats = getPointStats(idsHandle.begin(), idsHandle.end(), points);

    PLOG("Midpoint is " << stats.mid << std::endl);

    Partitioning<D, Precision> partitioning;
    partitioning.reserve(2);

    std::vector<GrowingHashTableHandle<Concurrent_LP_Set>> pointHandles;
    pointHandles.reserve(partitioning.size());

    partitioning.emplace_back(idsHandle.size() / 2);
    partitioning[0].id = 0;
    for (uint d = 0; d < D; ++d) {
        partitioning[0].bounds.low[d] = stats.min[d];
        partitioning[0].bounds.high[d] = d == k ? stats.mid[d] : stats.max[d];
    }
    pointHandles.emplace_back(partitioning[0].points.handle());

    partitioning.emplace_back(idsHandle.size() / 2);
    partitioning[1].id = 1;
    for (uint d = 0; d < D; ++d) {
        partitioning[1].bounds.low[d] = d == k ? stats.mid[d] : stats.min[d];
        partitioning[1].bounds.high[d] = stats.max[d];
    }
    pointHandles.emplace_back(partitioning[1].points.handle());

    for (auto id : idsHandle) {

        if (!dPoint<D, Precision>::isFinite(id))
            continue; // skip infinite points, they will be added later

        ASSERT(points.contains(id));
        const auto &p = points[id];

        uint part = (p.coords[k] > stats.mid[k]);

        ASSERT(partitioning[part].bounds.contains(p.coords));

        // LOG("Adding " << p << " to " << part << std::endl);
        pointHandles[part].insert(id);
    }

    // add infinite points
    for (uint i = 0; i < partitioning.size(); ++i) {
        for (uint k = dPoint<D, Precision>::cINF; k != 0; ++k) {
            pointHandles[i].insert(k);
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
    auto idsHandle = ids.handle();
    auto stats = getPointStats(idsHandle.begin(), idsHandle.end(), points);

    // cycle is lenght of provenance - 1 modulo D
    uint k = (provenance.size() - 1) % D;

    PLOG("Midpoint is " << stats.mid << std::endl);
    PLOG("Splitting dimension is " << k << std::endl);

    Partitioning<D, Precision> partitioning;
    partitioning.reserve(2);

    std::vector<GrowingHashTableHandle<Concurrent_LP_Set>> pointHandles;
    pointHandles.reserve(partitioning.size());

    partitioning.emplace_back(idsHandle.size() / 2);
    partitioning[0].id = 0;
    for (uint d = 0; d < D; ++d) {
        partitioning[0].bounds.low[d] = stats.min[d];
        partitioning[0].bounds.high[d] = d == k ? stats.mid[d] : stats.max[d];
    }
    pointHandles.emplace_back(partitioning[0].points.handle());

    partitioning.emplace_back(idsHandle.size() / 2);
    partitioning[1].id = 1;
    for (uint d = 0; d < D; ++d) {
        partitioning[1].bounds.low[d] = d == k ? stats.mid[d] : stats.min[d];
        partitioning[1].bounds.high[d] = stats.max[d];
    }
    pointHandles.emplace_back(partitioning[1].points.handle());

    for (auto id : idsHandle) {
        if (!dPoint<D, Precision>::isFinite(id))
            continue; // skip infinite points, they will be added later

        ASSERT(points.contains(id));
        const auto &p = points[id];

        uint part = (p.coords[k] > stats.mid[k]);

        ASSERT(partitioning[part].bounds.contains(p.coords));

        // LOG("Adding " << p << " to " << part << std::endl);
        pointHandles[part].insert(id);
    }

    // add infinite points
    for (uint i = 0; i < partitioning.size(); ++i) {
        for (uint k = dPoint<D, Precision>::cINF; k != 0; ++k) {
            pointHandles[i].insert(k);
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

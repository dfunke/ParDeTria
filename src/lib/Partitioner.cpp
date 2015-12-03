#include "Partitioner.h"

template<uint D, typename Precision>
Partitioning<D, Precision> dWayPartitioner<D, Precision>::partition(
        const dPointStats<D, Precision> &stats,
        const Point_Ids &ids,
        const dPoints<D, Precision> &points,
        __attribute((unused)) const std::string &provenance,
	bool ignoreInfinte) const {
    // do mid-point based partitioning for now

    PLOG("Midpoint is " << stats.mid << std::endl);

    Partitioning<D, Precision> partitioning;
    partitioning.reserve(pow(2, D));

    for (uint i = 0; i < pow(2, D); ++i) {
        partitioning.emplace_back(ids.size() / (pow(2, D)));
        partitioning[i].id = i;

        for (uint d = 0; d < D; ++d) {
            partitioning[i].bounds.low[d] =
                    i & (1 << d) ? stats.mid[d] : stats.min[d];
            partitioning[i].bounds.high[d] =
                    i & (1 << d) ? stats.max[d] : stats.mid[d];
        }
    }

    VTUNE_TASK(PartitioningDistribute);
    for (auto id : ids) {

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
        partitioning[part].points.insert(id);
    }
    VTUNE_END_TASK(PartitioningDistribute);

    // add infinite points
    if(!ignoreInfinte) {
        for (uint i = 0; i < partitioning.size(); ++i) {
            for (tIdType k = dPoint<D, Precision>::cINF; k != 0; ++k) {
                partitioning[i].points.insert(k);
            }
        }
    }

    return partitioning;
}

template<uint D, typename Precision>
Partitioning<D, Precision> OneDimPartitioner<D, Precision>::partition(
        const dPointStats<D, Precision> &stats,
        const Point_Ids &ids,
        const dPoints<D, Precision> &points,
        __attribute((unused)) const std::string &provenance,
	bool ignoreInfinte) const {
    // do mid-point based partitioning for now

    PLOG("Midpoint is " << stats.mid << std::endl);

    Partitioning<D, Precision> partitioning;
    partitioning.reserve(2);

    partitioning.emplace_back(ids.size() / 2);
    partitioning[0].id = 0;
    for (uint d = 0; d < D; ++d) {
        partitioning[0].bounds.low[d] = stats.min[d];
        partitioning[0].bounds.high[d] = d == k ? stats.mid[d] : stats.max[d];
    }

    partitioning.emplace_back(ids.size() / 2);
    partitioning[1].id = 1;
    for (uint d = 0; d < D; ++d) {
        partitioning[1].bounds.low[d] = d == k ? stats.mid[d] : stats.min[d];
        partitioning[1].bounds.high[d] = stats.max[d];
    }

    VTUNE_TASK(PartitioningDistribute);
    for (auto id : ids) {

        if (!dPoint<D, Precision>::isFinite(id))
            continue; // skip infinite points, they will be added later

        ASSERT(points.contains(id));
        const auto &p = points[id];

        uint part = (p.coords[k] > stats.mid[k]);

        ASSERT(partitioning[part].bounds.contains(p.coords));

        // LOG("Adding " << p << " to " << part << std::endl);
        partitioning[part].points.insert(id);
    }
    VTUNE_END_TASK(PartitioningDistribute);

    // add infinite points
    if(!ignoreInfinte) {
        for (uint i = 0; i < partitioning.size(); ++i) {
            for (tIdType k = dPoint<D, Precision>::cINF; k != 0; ++k) {
                partitioning[i].points.insert(k);
            }
        }
    }

    return partitioning;
}

template<uint D, typename Precision>
Partitioning<D, Precision>
CyclePartitioner<D, Precision>::partition(const dPointStats<D, Precision> &stats,
                                          const Point_Ids &ids,
                                          const dPoints<D, Precision> &points,
                                          const std::string &provenance,
					  bool ignoreInfinte) const {

    // cycle is lenght of provenance - 1 modulo D
    uint k = (provenance.size() - 1) % D;

    PLOG("Midpoint is " << stats.mid << std::endl);
    PLOG("Splitting dimension is " << k << std::endl);

    std::vector<Concurrent_Fixed_Point_Ids> cPartPoints;
    cPartPoints.reserve(2);
    cPartPoints.emplace_back(ids.size());
    cPartPoints.emplace_back(ids.size());


    VTUNE_TASK(PartitioningDistribute);
    tbb::parallel_for(ids.range(), [&stats, &points, &cPartPoints, k] (const auto &r) {

        for(auto id : r) {
            if (!dPoint<D, Precision>::isFinite(id))
                continue; // skip infinite points, they will be added later

            ASSERT(points.contains(id));
            const auto &p = points[id];

            uint part = (p.coords[k] > stats.mid[k]);

            //ASSERT(partitioning[part].bounds.contains(p.coords));

            // LOG("Adding " << p << " to " << part << std::endl);
#ifndef NDEBUG
            auto inserted =
#endif
                    cPartPoints[part].insert(id);
            ASSERT(inserted);
        }
    });
    VTUNE_END_TASK(PartitioningDistribute);

    Partitioning<D, Precision> partitioning;
    partitioning.reserve(2);

    partitioning.emplace_back(std::move(cPartPoints[0]));
    partitioning[0].id = 0;
    for (uint d = 0; d < D; ++d) {
        partitioning[0].bounds.low[d] = stats.min[d];
        partitioning[0].bounds.high[d] = d == k ? stats.mid[d] : stats.max[d];
    }

    partitioning.emplace_back(std::move(cPartPoints[1]));
    partitioning[1].id = 1;
    for (uint d = 0; d < D; ++d) {
        partitioning[1].bounds.low[d] = d == k ? stats.mid[d] : stats.min[d];
        partitioning[1].bounds.high[d] = stats.max[d];
    }

    // add infinite points
    if(!ignoreInfinte) {
        for (uint i = 0; i < partitioning.size(); ++i) {
            for (tIdType k = dPoint<D, Precision>::cINF; k != 0; ++k) {
                partitioning[i].points.insert(k);
            }
        }
    }

    return partitioning;
}

template<uint D, typename Precision>
Partitioning<D, Precision>
ExtendPartitioner<D, Precision>::partition(const dPointStats<D, Precision> &stats,
                                          const Point_Ids &ids,
                                          const dPoints<D, Precision> &points,
                                          __attribute__((unused)) const std::string &provenance,
                                          bool ignoreInfinte) const {

    // determine most vast dimension
    Precision maxExtend = 0;
    int maxDim = -1;

    for (uint d = 0; d < D; ++d) {
        if (stats.max[d] - stats.min[d] > maxExtend) {
            maxExtend = stats.max[d] - stats.min[d];
            maxDim = d;
        }
    }

    ASSERT(0 <= maxDim && maxDim < (signed) D);
    uint k = static_cast<uint>(maxDim);

    PLOG("Midpoint is " << stats.mid << std::endl);
    PLOG("Splitting dimension is " << k << std::endl);

    std::vector<Concurrent_Fixed_Point_Ids> cPartPoints;
    cPartPoints.reserve(2);
    cPartPoints.emplace_back(ids.size());
    cPartPoints.emplace_back(ids.size());


    VTUNE_TASK(PartitioningDistribute);
    tbb::parallel_for(ids.range(), [&stats, &points, &cPartPoints, k] (const auto &r) {

        for(auto id : r) {
            if (!dPoint<D, Precision>::isFinite(id))
                continue; // skip infinite points, they will be added later

            ASSERT(points.contains(id));
            const auto &p = points[id];

            uint part = (p.coords[k] > stats.mid[k]);

            //ASSERT(partitioning[part].bounds.contains(p.coords));

            // LOG("Adding " << p << " to " << part << std::endl);
#ifndef NDEBUG
            auto inserted =
#endif
                    cPartPoints[part].insert(id);
            ASSERT(inserted);
        }
    });
    VTUNE_END_TASK(PartitioningDistribute);

    Partitioning<D, Precision> partitioning;
    partitioning.reserve(2);

    partitioning.emplace_back(std::move(cPartPoints[0]));
    partitioning[0].id = 0;
    for (uint d = 0; d < D; ++d) {
        partitioning[0].bounds.low[d] = stats.min[d];
        partitioning[0].bounds.high[d] = d == k ? stats.mid[d] : stats.max[d];
    }

    partitioning.emplace_back(std::move(cPartPoints[1]));
    partitioning[1].id = 1;
    for (uint d = 0; d < D; ++d) {
        partitioning[1].bounds.low[d] = d == k ? stats.mid[d] : stats.min[d];
        partitioning[1].bounds.high[d] = stats.max[d];
    }

    // add infinite points
    if(!ignoreInfinte) {
        for (uint i = 0; i < partitioning.size(); ++i) {
            for (tIdType k = dPoint<D, Precision>::cINF; k != 0; ++k) {
                partitioning[i].points.insert(k);
            }
        }
    }

    return partitioning;
}

// specializations
template
class dWayPartitioner<2, float>;

template
class OneDimPartitioner<2, float>;

template
class CyclePartitioner<2, float>;

template
class ExtendPartitioner<2, float>;

template
class dWayPartitioner<3, float>;

template
class OneDimPartitioner<3, float>;

template
class CyclePartitioner<3, float>;

template
class ExtendPartitioner<3, float>;

template
class dWayPartitioner<2, double>;

template
class OneDimPartitioner<2, double>;

template
class CyclePartitioner<2, double>;

template
class ExtendPartitioner<2, double>;

template
class dWayPartitioner<3, double>;

template
class OneDimPartitioner<3, double>;

template
class CyclePartitioner<3, double>;

template
class ExtendPartitioner<3, double>;
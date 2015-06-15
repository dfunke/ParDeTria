#include "DCTriangulator.h"

// debug
#ifndef NDEBUG

#include <csignal>

#endif

// std library
#include <fstream>
#include <iostream>
#include <vector>
#include <functional>
#include <set>
#include <stdexcept>
#include <queue>
#include <atomic>
#include <thread>

#include <tbb/parallel_for.h>
// use own parallel_do
//#include <tbb/parallel_do.h>
#include "mods/parallel_do.h"
#include <tbb/parallel_sort.h>
#include <tbb/enumerable_thread_specific.h>

// tbb
#include <tbb/spin_mutex.h>

// own
#include "Geometry.h"
#include "CGALTriangulator.h"

#include "utils/Timings.h"
#include "utils/Logger.h"
#include "utils/ASSERT.h"
#include "utils/VTuneAdapter.h"

//**************************

template<uint D, typename Precision>
constexpr Precision DCTriangulator<D, Precision>::SAFETY;

template<uint D, typename Precision>
constexpr uint DCTriangulator<D, Precision>::BASE_CUTOFF;

//**************************

template<uint D, typename Precision>
DCTriangulator<D, Precision>::DCTriangulator(
        const dBox<D, Precision> &_bounds,
        dPoints<D, Precision> &_points,
        const uint _recursionDepth,
        const unsigned char splitter,
        const uint gridOccupancy,
        const bool parallelBaseSolver)
        : Triangulator<D, Precision>(_bounds, _points), recursionDepth(_recursionDepth) {


    // add infinite points to data set
    auto stats = getPointStats(std::size_t(0), this->points.size(), this->points);
    for (uint i = 0; i < pow(2, D); ++i) {
        VLOG("Point stats: " << stats.min << " - " << stats.mid << " - "
             << stats.max << std::endl);

        dPoint<D, Precision> p;
        //p.id = dPoint<D, Precision>::cINF | i;
        p.coords = stats.mid;

        for (uint d = 0; d < D; ++d)
            p.coords[d] +=
                    (i & (1 << d) ? 1 : -1) * 2 * SAFETY * (stats.max[d] - stats.min[d]);

        this->points.emplace_back(p);
    }

    if (parallelBaseSolver)
        baseTriangulator = std::make_unique<CGALTriangulator<D, Precision, true>>(this->baseBounds, this->points,
                                                                                  gridOccupancy);
    else
        baseTriangulator = std::make_unique<CGALTriangulator<D, Precision, false>>(this->baseBounds, this->points,
                                                                                   gridOccupancy);

    if (splitter != 0)
        partitioner = Partitioner<D, Precision>::make(splitter);

}

template<uint D, typename Precision>
void DCTriangulator<D, Precision>::getEdge(const PartialTriangulation &pt,
                                           const dSimplices<D, Precision> &simplices,
                                           const Partitioning<D, Precision> &partitioning, const uint &partition,
                                           Ids &edgePoints, Ids &edgeSimplices) {

    // infinite points to edge
    auto edgePointsHandle = edgePoints.handle();
    for (uint k = dPoint<D, Precision>::cINF; k != 0; ++k) {
        edgePointsHandle.insert(k);
    }

    auto convexHullHandle = pt.convexHull.handle();

    Ids wqa(convexHullHandle.capacity());
    auto wqaHandle = wqa.handle();

    wqaHandle.unsafe_copy(convexHullHandle); // set of already checked simplices
    wqaHandle.insert(dSimplex<D, Precision>::cINF); //we don't want to check the infinte vertex

    tbb::enumerable_thread_specific<GrowingHashTableHandle<Concurrent_LP_Set>,
            tbb::cache_aligned_allocator<GrowingHashTableHandle<Concurrent_LP_Set>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsWqaHandle(std::ref(wqa));

    tbb::enumerable_thread_specific<GrowingHashTableHandle<Concurrent_LP_Set>,
            tbb::cache_aligned_allocator<GrowingHashTableHandle<Concurrent_LP_Set>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsEdgePointsHandle(std::ref(edgePoints));

    tbb::enumerable_thread_specific<GrowingHashTableHandle<Concurrent_LP_Set>,
            tbb::cache_aligned_allocator<GrowingHashTableHandle<Concurrent_LP_Set>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsEdgeSimplicesHandle(std::ref(edgeSimplices));

    /* Walk along the neighbors,
     * testing for each neighbor whether its circumsphere is within the
     * partition or not
     */

    INDENT

    VTUNE_TASK(IdentifyEdge);
    tbb::parallel_do(convexHullHandle, [&](const uint id,
                                  tbb::parallel_do_feeder<uint> &feeder) {

        if(!dSimplex<D,Precision>::isFinite(id))
            return;

        auto wqaHandle = tsWqaHandle.local();
        auto edgePointsHandle = tsEdgePointsHandle.local();
        auto edgeSimplicesHandle = tsEdgeSimplicesHandle.local();

        const auto &simplex = simplices[id];
        const auto cs = simplex.circumsphere(this->points);
        bool intersectsBounds = false;
        for (uint i = 0; i < partitioning.size(); ++i) {
            if (i != partition && partitioning[i].bounds.intersects(cs))
                intersectsBounds = true;
        }
        if (intersectsBounds) {

            PLOG("Adding " << simplex
                 << " to edge -> circumcircle criterion"
                 << std::endl);
            edgeSimplicesHandle.insert(simplex.id);

            for (uint i = 0; i < D + 1; ++i) {
                edgePointsHandle.insert(simplex.vertices[i]);
            }

            for (const auto &n : simplex.neighbors) {
                if (wqaHandle.insert(n)) {
                    // n was not yet inspected
                    feeder.add(n);
                }
            }
        }
    });

    DEDENT
}

template<uint D, typename Precision>
cWuFaces DCTriangulator<D, Precision>::buildWhereUsed(const dSimplices<D, Precision> &DT,
                                                      const Ids &edgeSimplices) {

    auto edgeSimplicesHandle = edgeSimplices.handle();

    Ids wqa(edgeSimplicesHandle.size()); // set of already checked simplices
    tbb::enumerable_thread_specific<GrowingHashTableHandle<Concurrent_LP_Set>,
            tbb::cache_aligned_allocator<GrowingHashTableHandle<Concurrent_LP_Set>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsWqaHandle(std::ref(wqa));

    /* We need to build the wuFaces DS for the simplices of the first layer "inward" of the edge
     * plus one more layer, to re-find their neighbors
     */

    VTUNE_TASK(BuildWU);
    cWuFaces wuFaces((D + 1) * (D + 1) * (D + 1) * edgeSimplicesHandle.size());
    tbb::enumerable_thread_specific<GrowingHashTableHandle<Concurrent_LP_MultiMap>,
            tbb::cache_aligned_allocator<GrowingHashTableHandle<Concurrent_LP_MultiMap>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsWuFacesHandle(std::ref(wuFaces));

    tbb::parallel_for(edgeSimplicesHandle.range(), [&](const auto &range) {

        auto wqaHandle = tsWqaHandle.local();
        auto wuFacesHandle = tsWuFacesHandle.local();

        for (const auto &edgeSimplexID : range) {
            const auto &edgeSimplex = DT[edgeSimplexID];
            for (const auto &firstLayer : edgeSimplex.neighbors) {
                if (dSimplex<D, Precision>::isFinite(firstLayer) && !edgeSimplicesHandle.count(firstLayer)) {
                    // we have an "inward" neighbor, add it to the wuFaces DS
                    if (wqaHandle.insert(firstLayer)) {
                        for (uint i = 0; i < D + 1; ++i) {
                            auto facetteHash = DT[firstLayer].faceFingerprint(i);
                            wuFacesHandle.insert(facetteHash, firstLayer);
                        }
                    }
                    // now loop over its neighbors, adding the ones not belonging to the edge
                    // TODO maybe we can avoid this
                    for (const auto &secondLayer : DT[firstLayer].neighbors) {
                        if (dSimplex<D, Precision>::isFinite(secondLayer) && !edgeSimplicesHandle.count(secondLayer)
                            && wqaHandle.insert(secondLayer)) {

                            for (uint i = 0; i < D + 1; ++i) {
                                auto facetteHash = DT[secondLayer].faceFingerprint(i);
                                wuFacesHandle.insert(facetteHash, secondLayer);
                            }
                        }
                    }
                }
            }
        }
    });

    return wuFaces;
}

template<uint D, typename Precision>
void DCTriangulator<D, Precision>::updateNeighbors(
        dSimplices<D, Precision> &simplices,
        const PartialTriangulation &pt,
        const Ids &toCheck,
        const cWuFaces & wuFaces,
        __attribute__((unused)) const std::string &provenance) {

    INDENT
    const uint saveIndent = LOGGER.getIndent();

    auto toCheckHandle = toCheck.handle();
    Ids wqa(toCheckHandle.size());

    tbb::enumerable_thread_specific<ConstGrowingHashTableHandle<Concurrent_LP_Set>,
            tbb::cache_aligned_allocator<ConstGrowingHashTableHandle<Concurrent_LP_Set>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsSimplicesHandle(std::ref(pt.simplices));

    tbb::enumerable_thread_specific<GrowingHashTableHandle<Concurrent_LP_Set>,
            tbb::cache_aligned_allocator<GrowingHashTableHandle<Concurrent_LP_Set>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsWqaHandle(std::ref(wqa));

    tbb::enumerable_thread_specific<std::set<tIdType>,
            tbb::cache_aligned_allocator<std::set<tIdType>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsNeighborSet;

    auto wuFacesHandle = wuFaces.handle();

#ifndef NDEBUG
    std::atomic<uint> checked(0);
    std::atomic<uint> updated(0);
#endif

    VTUNE_TASK(UpdateNeighbors);
    tbb::parallel_do(toCheckHandle.range(), [&](const uint &id,
                                  tbb::parallel_do_feeder<uint> &feeder) {

        LOGGER.setIndent(saveIndent);
        auto wqa = tsWqaHandle.local();

        if (!wqa.insert(id)) {
            return;
        }

        dSimplex<D, Precision> &simplex = simplices[id];
        auto simplicesHandle = tsSimplicesHandle.local();

        PLOG("Checking neighbors of " << simplex << std::endl);
#ifndef NDEBUG
        ++checked;
#endif

        bool requiresUpdate = false;
        for (const auto &n : simplex.neighbors) {
            if (dSimplex<D, Precision>::isFinite(n) &&
                (!simplicesHandle.count(n) || !simplex.isNeighbor(simplices[n]))) {
                requiresUpdate = true;
                break;
            }
        }

        if (!requiresUpdate) {
            INDENT
            PLOG("No update necessary" << std::endl);
            DEDENT

            return;
        }

        PLOG("Updating neighbors of " << simplex << std::endl);
#ifndef NDEBUG
        ++updated;
#endif

        INDENT
        uint neighborIdx = 0;
        auto neighborSet = tsNeighborSet.local();
        neighborSet.clear();
        for (uint i = 0; i < D + 1; ++i) {

            auto facetteHash = simplex.faceFingerprint(i);
            auto range = wuFacesHandle.get(facetteHash);
            PLOG("Key: " << facetteHash << " #Results: " << std::distance(range.first, range.second) << std::endl;);
            for (auto it = range.first; it != range.second; ++it) {
                auto i = *it;
                if (i.second != simplex.id && dSimplex<D, Precision>::isFinite(i.second)
                    && simplicesHandle.count(i.second) && simplex.isNeighbor(simplices[i.second])) {
                    PLOG("Neighbor with " << simplices[i.second] << std::endl);

                    if (neighborSet.insert(i.second).second) {
                        simplex.neighbors[neighborIdx++] = i.second;
                        feeder.add(i.second);
                    }
                }
            }
        }
        DEDENT

        ASSERT(0 < neighborIdx && neighborIdx <= D + 1);

#ifndef NDEBUG
        if (!(neighborIdx > 0 && neighborIdx <= D + 1)) {
            LOG("Error: wrong number of neighbors for simplex " << simplex << std::endl);
        }
#endif

        // ASSERT(simplex.neighbors.size() > 0 && simplex.neighbors.size() <=
        // D+1);

        while (neighborIdx < D + 1) { // if it is a simplex at the border,
            // it might have one infinite
            // simplex as neighbor
            simplex.neighbors[neighborIdx++] = uint(dSimplex<D, Precision>::cINF);
        }
    });

#ifndef NDEBUG
    LOG("Updating neighbors - checked: " << checked << " updated: " << updated << std::endl);
#endif

    DEDENT
}

template<uint D, typename Precision>
PartialTriangulation DCTriangulator<D, Precision>::mergeTriangulation(std::vector<PartialTriangulation> &partialDTs,
                                                                      dSimplices<D, Precision> &DT,
                                                                      const Ids &edgeSimplices,
                                                                      const PartialTriangulation &edgeDT,
                                                                      const Partitioning<D, Precision> &partitioning,
                                                                      const std::string &provenance) {

    LOG("Merging partial DTs into one triangulation" << std::endl);

    VTUNE_TASK(MergeTriangulations);

    VTUNE_TASK(CombineTriangulations);
    PartialTriangulation pt = std::move(partialDTs[0]);
    auto simplicesHandle = pt.simplices.handle();
    auto convexHullHandle = pt.convexHull.handle();
    auto edgeSimplicesHandle = edgeSimplices.handle();

    for (uint i = 1; i < partialDTs.size(); ++i) {
        simplicesHandle.unsafe_merge(partialDTs[i].simplices.handle(), edgeSimplicesHandle);
        convexHullHandle.unsafe_merge(partialDTs[i].convexHull.handle(), edgeSimplicesHandle);
    }
    VTUNE_END_TASK(CombineTriangulations);

    //build the where used datastructure for the faces
    cWuFaces wuFaces = buildWhereUsed(DT, edgeSimplices);

    auto cmpFingerprint =
            [&DT](const tIdType &a, const tIdType &b) {
                return DT[a].fingerprint() < DT[b].fingerprint();
            };

    VTUNE_TASK(SortDeletedVertices);
    std::vector<tIdType> deletedSimplices(edgeSimplicesHandle.begin(), edgeSimplicesHandle.end());

    tbb::parallel_sort(deletedSimplices.begin(), deletedSimplices.end(), cmpFingerprint);
    VTUNE_END_TASK(SortDeletedVertices);


    // merge partial DTs and edge DT
    LOG("Merging triangulations" << std::endl);
    Ids insertedSimplices(edgeSimplicesHandle.size() / 2);

    VTUNE_TASK(AddBorderSimplices);

    tbb::enumerable_thread_specific<GrowingHashTableHandle<Concurrent_LP_Set>,
            tbb::cache_aligned_allocator<GrowingHashTableHandle<Concurrent_LP_Set>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsSimplicesHandle(std::ref(pt.simplices));

    tbb::enumerable_thread_specific<GrowingHashTableHandle<Concurrent_LP_Set>,
            tbb::cache_aligned_allocator<GrowingHashTableHandle<Concurrent_LP_Set>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsConvexHullHandle(std::ref(pt.convexHull));

    tbb::enumerable_thread_specific<GrowingHashTableHandle<Concurrent_LP_MultiMap>,
            tbb::cache_aligned_allocator<GrowingHashTableHandle<Concurrent_LP_MultiMap>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tswuFacesHandle(std::ref(wuFaces));

    tbb::enumerable_thread_specific<GrowingHashTableHandle<Concurrent_LP_Set>,
            tbb::cache_aligned_allocator<GrowingHashTableHandle<Concurrent_LP_Set>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsInsertedHandle(std::ref(insertedSimplices));

    tbb::parallel_for(edgeDT.simplices.handle().range(), [&](const auto &r) {

        auto simplicesHandle = tsSimplicesHandle.local();
        auto convexHullHandle = tsConvexHullHandle.local();
        auto wuFacesHandle = tswuFacesHandle.local();
        auto insertedHandle = tsInsertedHandle.local();

        for (const auto &i : r) {

            const dSimplex<D, Precision> &edgeSimplex = DT[i];

            // check whether edgeSimplex is completely contained in one partition
            uint p0 = partitioning.partition(edgeSimplex.vertices[0]);
            bool inOnePartition = partitioning[p0].contains(edgeSimplex);

            if (!inOnePartition) {
                simplicesHandle.insert(edgeSimplex.id);

                //convex hull treatment
                if (!edgeSimplex.isFinite())
                    convexHullHandle.insert(edgeSimplex.id);

                for (uint d = 0; d < D + 1; ++d) {
                    wuFacesHandle.insert((edgeSimplex.faceFingerprint(d)), edgeSimplex.id);
                }

                insertedHandle.insert(edgeSimplex.id);
            } else {
                // the simplex is completely in one partition -> it must have been found
                // before
                auto range =
                        std::equal_range(deletedSimplices.begin(), deletedSimplices.end(),
                                         edgeSimplex.id, cmpFingerprint);

                for (auto it = range.first; it != range.second; ++it) {
                    if (edgeSimplex.equalVertices(DT[*it])) {
                        simplicesHandle.insert(edgeSimplex.id);

                        //convex hull treatment
                        if (!edgeSimplex.isFinite())
                            convexHullHandle.insert(edgeSimplex.id);

                        for (uint d = 0; d < D + 1; ++d) {
                            wuFacesHandle.insert((edgeSimplex.faceFingerprint(d)),
                                               edgeSimplex.id);
                        }

                        insertedHandle.insert(edgeSimplex.id);
                    }
                }
            }
        }
    });
    VTUNE_END_TASK(AddBorderSimplices);

    //ASSERT(DT.countDuplicates() == 0);

    LOG("Updating neighbors" << std::endl);
    updateNeighbors(DT, pt, insertedSimplices, wuFaces, provenance);

    return pt;
}

template<uint D, typename Precision>
PartialTriangulation DCTriangulator<D, Precision>::_triangulateBase(dSimplices<D, Precision> &DT,
                                                                    const Ids &partitionPoints,
                                                                    const dBox<D, Precision> &bounds,
                                                                    const std::string provenance) {

    VTUNE_TASK(TriangulateBase);
    PROFILER_MEAS(provenance.find('e') == std::string::npos? "basecaseSize" : "edgeBasecaseSize", partitionPoints.handle().size());

    LOGGER.setIndent(provenance.length());
    LOG("triangulateBASE called on level " << provenance << " with "
        << partitionPoints.handle().size() << " points" << std::endl);

    INDENT
    auto pt = baseTriangulator->_triangulate(DT, partitionPoints, bounds, provenance);
    DEDENT

    return pt;
}

template<uint D, typename Precision>
PartialTriangulation DCTriangulator<D, Precision>::_triangulate(dSimplices<D, Precision> &DT,
                                                                const Ids &partitionPoints,
                                                                const dBox<D, Precision> &bounds,
                                                                const std::string provenance) {

    LOGGER.setIndent(provenance.length());

    PROFILER_MEAS("nPoints", partitionPoints.handle().size());

    LOG("triangulateDAC called on level " << provenance << " with "
        << partitionPoints.handle().size() << " points" << std::endl);


    auto partitionPointsHandle = partitionPoints.handle();
    if (provenance.length() - 1 < recursionDepth && partitionPointsHandle.size() > BASE_CUTOFF) {
        VTUNE_TASK(TriangulateRecursive);

        LOG("Recursive case" << std::endl);

        // partition input
        LOG("Partioning" << std::endl);
        INDENT

        VTUNE_TASK(Partitioning);
        const auto partioning =
                partitioner->partition(partitionPoints, this->points, provenance);
        VTUNE_END_TASK(Partitioning);

        DEDENT

        Ids edgePointIds(partitionPointsHandle.size() / 4);
        Ids edgeSimplexIds(partitionPointsHandle.size() / 4);

        std::vector<PartialTriangulation> partialDTs;
        partialDTs.resize(partioning.size());

        VTUNE_TASK(TriangulatePartitions);
        tbb::parallel_for(
                std::size_t(0), partioning.size(),
                [&](const uint i) {
                    LOGGER.setIndent(
                            provenance.length()); // new thread, initialize Logger indent
                    INDENT
                    LOG("Partition " << i << " on thread " << std::this_thread::get_id()
                        << std::endl);
                    partialDTs[i] = _triangulate(DT, partioning[i].points, partioning[i].bounds,
                                                 provenance + std::to_string(i));

                    VLOG("Partition " << i << ": extracting edge" << std::endl);
                    //TODO how to handle edge detection
                    getEdge(partialDTs[i], DT, partioning, i, edgePointIds, edgeSimplexIds);
                    DEDENT
                }
        );
        VTUNE_END_TASK(TriangulatePartitions);

        PROFILER_MEAS("edgePoints", edgePointIds.handle().size());
        PROFILER_MEAS("edgePointsPerc", edgePointIds.handle().size() / partitionPointsHandle.size());


        PROFILER_MEAS("edgeSimplices", edgeSimplexIds.handle().size());
        PROFILER_MEAS("edgeSimplicesPerc", edgeSimplexIds.handle().size() / ([&] () -> std::size_t {
            std::size_t size = 0;
            for(uint i = 0; i < partialDTs.size(); ++i)
                size += partialDTs[i].simplices.handle().size();
            return size;
        })());

        LOG("Edge has " << edgeSimplexIds.handle().size() << " simplices with "
            << edgePointIds.handle().size() << " points" << std::endl);

        INDENT

        VTUNE_TASK(TriangulateEdge);
        //TODO how to identify the edge triangulation
        PartialTriangulation edgeDT = _triangulate(DT, edgePointIds, bounds, provenance + "e");
        VTUNE_END_TASK(TriangulateEdge);

        PROFILER_MEAS("edgeDT", edgeDT.simplices.handle().size());
        PROFILER_MEAS("edgeDTPerc", edgeDT.simplices.handle().size()/ ([&] () -> std::size_t {
            std::size_t size = 0;
            for(uint i = 0; i < partialDTs.size(); ++i)
                size += partialDTs[i].simplices.handle().size();
            return size;
        })());

        DEDENT

        //TODO this is all wrong
        return mergeTriangulation(partialDTs, DT, edgeSimplexIds, edgeDT, partioning, provenance);

    } else { // base case
        return _triangulateBase(DT, partitionPoints, bounds, provenance);
    }
}

// specializations

template
class DCTriangulator<2, float>;

template
class DCTriangulator<3, float>;

template
class DCTriangulator<2, double>;

template
class DCTriangulator<3, double>;

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
#include <tbb/parallel_do.h>
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
        p.id = dPoint<D, Precision>::cINF | i;
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
void DCTriangulator<D, Precision>::getEdge(const Ids &convexHull,
                                           const dSimplices<D, Precision> &simplices,
                                           const Partitioning<D, Precision> &partitioning, const uint &partition,
                                           Ids &edgePoints, Ids &edgeSimplices) {

    // infinite points to edge
    for (uint k = dPoint<D, Precision>::cINF; k != 0; ++k) {
        edgePoints.insert(k);
    }

    Ids wqa;
    wqa.insert(convexHull.begin(), convexHull.end()); // set of already checked simplices
    wqa.insert(dSimplex<D, Precision>::cINF); //we don't want to check the infinte vertex
    std::deque<uint> wq(convexHull.begin(), convexHull.end());

    /* Walk along the neighbors,
     * testing for each neighbor whether its circumsphere is within the
     * partition or not
     */

    INDENT

    VTUNE_TASK(IdentifyEdge);
    while (!wq.empty()) {
        uint x = wq.front();
        wq.pop_front();

        const auto &simplex = simplices[x];
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
            edgeSimplices.insert(simplex.id);
            if(simplex.id == dSimplex<D, Precision>::cINF) {
                std::cerr << x << " " << simplices.size() << std::endl;
                raise(SIGINT);
            }

            for (uint i = 0; i < D + 1; ++i) {
                edgePoints.insert(simplex.vertices[i]);
            }

            for (const auto &n : simplex.neighbors) {
                if (wqa.insert(n).second) {
                    // n was not yet inspected
                    wq.push_back(n);
                }
            }
        }
    }

    DEDENT
}

template<uint D, typename Precision>
void DCTriangulator<D, Precision>::buildWhereUsed(const dSimplices<D, Precision> &DT,
                                                  const Ids &edgeSimplices,
                                                  cWuFaces &wuFaces) {

    Ids wqa; // set of already checked simplices

    /* We need to build the wuFaces DS for the simplices of the first layer "inward" of the edge
     * plus one more layer, to re-find their neighbors
     */


    VTUNE_TASK(BuildWU);
    wuFaces.reserve((D + 1) * (D + 1) * edgeSimplices.size());

    tbb::parallel_for(edgeSimplices.range(), [&](const auto &range) {

        for (const auto &edgeSimplexID : range) {
            const auto &edgeSimplex = DT[edgeSimplexID];
            for (const auto &firstLayer : edgeSimplex.neighbors) {
                if (dSimplex<D, Precision>::isFinite(firstLayer) && !edgeSimplices.count(firstLayer)) {
                    // we have an "inward" neighbor, add it to the wuFaces DS
                    if (wqa.insert(firstLayer).second) {
                        for (uint i = 0; i < D + 1; ++i) {
                            auto facetteHash = DT[firstLayer].faceFingerprint(i);
                            wuFaces.emplace(facetteHash, firstLayer);
                        }
                    }
                    // now loop over its neighbors, adding the ones not belonging to the edge
                    // TODO maybe we can avoid this
                    for (const auto &secondLayer : DT[firstLayer].neighbors) {
                        if (dSimplex<D, Precision>::isFinite(secondLayer) && !edgeSimplices.count(secondLayer)
                            && wqa.insert(secondLayer).second) {

                            for (uint i = 0; i < D + 1; ++i) {
                                auto facetteHash = DT[secondLayer].faceFingerprint(i);
                                wuFaces.emplace(facetteHash, secondLayer);
                            }
                        }
                    }
                }
            }
        }
    });
}

template<uint D, typename Precision>
void DCTriangulator<D, Precision>::updateNeighbors(
        dSimplices<D, Precision> &simplices,
        const PartialTriangulation &pt,
        const Ids &toCheck,
        __attribute__((unused)) const std::string &provenance) {

    INDENT
    const uint saveIndent = LOGGER.getIndent();

    tbb::concurrent_unordered_set<uint> wqa(toCheck.size());
    tbb::enumerable_thread_specific<std::set<uint>,
            tbb::cache_aligned_allocator<std::set<uint>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsNeighborSet;

#ifndef NDEBUG
    std::atomic<uint> checked(0);
    std::atomic<uint> updated(0);
#endif

    VTUNE_TASK(UpdateNeighbors);
    tbb::parallel_do(toCheck, [&](const uint &id,
                                  tbb::parallel_do_feeder<uint> &feeder) {

        LOGGER.setIndent(saveIndent);

        if (!wqa.insert(id).second) {
            return;
        }

        dSimplex<D, Precision> &simplex = simplices[id];

        PLOG("Checking neighbors of " << simplex << std::endl);
#ifndef NDEBUG
        ++checked;
#endif

        bool requiresUpdate = false;
        for (const auto &n : simplex.neighbors) {
            if (dSimplex<D, Precision>::isFinite(n) &&
                (!pt.simplices.count(n) || !simplex.isNeighbor(simplices[n]))) {
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
            auto range = simplices.wuFaces.equal_range(facetteHash);
            PLOG("Key: " << facetteHash << " #Results: " << std::distance(range.first, range.second) << std::endl;);
            for (auto it = range.first; it != range.second; ++it) {
                if (it->second != simplex.id && dSimplex<D, Precision>::isFinite(it->second)
                    && pt.simplices.count(it->second) && simplex.isNeighbor(simplices[it->second])) {
                    PLOG("Neighbor with " << simplices[it->second] << std::endl);

                    if (neighborSet.insert(it->second).second) {
                        simplex.neighbors[neighborIdx++] = it->second;
                        feeder.add(it->second);
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
    LOG("Updating neighbors - checked: " << checked << " updated: " << updated
        << " simplices: " << simplices.size()
        << std::endl);
#endif

    DEDENT
}

template<uint D, typename Precision>
PartialTriangulation DCTriangulator<D, Precision>::mergeTriangulation(std::vector<PartialTriangulation> &partialDTs,
                                                                      dSimplices<D, Precision> &DT,
                                                                      const Ids &edgeSimplices,
                                                                      const Ids &edgeDT,
                                                                      const Partitioning<D, Precision> &partitioning,
                                                                      const std::string &provenance) {

    LOG("Merging partial DTs into one triangulation" << std::endl);

    VTUNE_TASK(MergeTriangulations);

    VTUNE_TASK(CombineTriangulations);
    PartialTriangulation pt = std::move(partialDTs[0]);

    for (uint i = 1; i < partialDTs.size(); ++i) {
        pt.simplices.insert(partialDTs[i].simplices.begin(), partialDTs[i].simplices.end());
        pt.convexHull.insert(partialDTs[i].convexHull.begin(), partialDTs[i].convexHull.end());
    }

    // copy the simplices not belonging to the edge
    for (auto &edgeSimplex : edgeSimplices) {
        pt.simplices.erase(edgeSimplex);
        pt.convexHull.erase(edgeSimplex);
    }
    VTUNE_END_TASK(CombineTriangulations);

    //build the where used datastructure for the faces
    buildWhereUsed(DT, edgeSimplices, DT.wuFaces);

    auto cmpFingerprint =
            [&DT](const tIdType &a, const tIdType &b) {
                return DT[a].fingerprint() < DT[b].fingerprint();
            };

    VTUNE_TASK(SortDeletedVertices);
    std::vector<tIdType> deletedSimplices(edgeSimplices.begin(), edgeSimplices.end());

    tbb::parallel_sort(deletedSimplices.begin(), deletedSimplices.end(), cmpFingerprint);
    VTUNE_END_TASK(SortDeletedVertices);


    // merge partial DTs and edge DT
    LOG("Merging triangulations" << std::endl);
    Ids insertedSimplices;

    VTUNE_TASK(AddBorderSimplices);
    tbb::parallel_for(std::size_t(0), edgeDT.bucket_count(), [&](const uint i) {

        for (auto it = edgeDT.begin(i); it != edgeDT.end(i); ++it) {

            const dSimplex<D, Precision> &edgeSimplex = DT[*it];

            // check whether edgeSimplex is completely contained in one partition
            uint p0 = partitioning.partition(edgeSimplex.vertices[0]);
            bool inOnePartition = partitioning[p0].contains(edgeSimplex);

            if (!inOnePartition) {
                pt.simplices.insert(edgeSimplex.id);

                //convex hull treatment
                if (!edgeSimplex.isFinite())
                    pt.convexHull.insert(edgeSimplex.id);

                for (uint d = 0; d < D + 1; ++d) {
                    DT.wuFaces.emplace((edgeSimplex.faceFingerprint(d)), edgeSimplex.id);
                }

                insertedSimplices.insert(edgeSimplex.id);
            } else {
                // the simplex is completely in one partition -> it must have been found
                // before
                auto range =
                        std::equal_range(deletedSimplices.begin(), deletedSimplices.end(),
                                         edgeSimplex.id, cmpFingerprint);

                for (auto it = range.first; it != range.second; ++it) {
                    if (edgeSimplex.equalVertices(DT[*it])) {
                        pt.simplices.insert(edgeSimplex.id);

                        //convex hull treatment
                        if (!edgeSimplex.isFinite())
                            pt.convexHull.insert(edgeSimplex.id);

                        for (uint d = 0; d < D + 1; ++d) {
                            DT.wuFaces.emplace((edgeSimplex.faceFingerprint(d)),
                                               edgeSimplex.id);
                        }

                        insertedSimplices.insert(edgeSimplex.id);
                    }
                }
            }
        }
    });
    VTUNE_END_TASK(AddBorderSimplices);

    //ASSERT(DT.countDuplicates() == 0);

    LOG("Updating neighbors" << std::endl);
    updateNeighbors(DT, pt, insertedSimplices, provenance);

    return pt;
}

template<uint D, typename Precision>
PartialTriangulation DCTriangulator<D, Precision>::_triangulateBase(dSimplices<D, Precision> &DT,
                                                                    const Ids partitionPoints,
                                                                    const dBox<D, Precision> &bounds,
                                                                    const std::string provenance) {

    VTUNE_TASK(TriangulateBase);

    LOGGER.setIndent(provenance.length());

    LOG("triangulateBASE called on level " << provenance << " with "
        << partitionPoints.size() << " points" << std::endl);

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

    LOG("triangulateDAC called on level " << provenance << " with "
        << partitionPoints.size() << " points" << std::endl);


    if (provenance.length() - 1 < recursionDepth && partitionPoints.size() > BASE_CUTOFF) {
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

        Ids edgePointIds;
        Ids edgeSimplexIds;

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
                    getEdge(partialDTs[i].convexHull, DT, partioning, i, edgePointIds, edgeSimplexIds);
                    DEDENT
                }
        );
        VTUNE_END_TASK(TriangulatePartitions);

        LOG("Edge has " << edgeSimplexIds.size() << " simplices with "
            << edgePointIds.size() << " points" << std::endl);

        INDENT

        VTUNE_TASK(TriangulateEdge);
        //TODO how to identify the edge triangulation
        PartialTriangulation edgeDT = _triangulate(DT, edgePointIds, bounds, provenance + "e");
        VTUNE_END_TASK(TriangulateEdge);

        DEDENT

        //TODO this is all wrong
        return mergeTriangulation(partialDTs, DT, edgeSimplexIds, edgeDT.simplices, partioning, provenance);

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

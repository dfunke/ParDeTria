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
        const uint _threads,
        std::unique_ptr<Partitioner<D, Precision>> && _splitter,
        const uint gridOccupancy,
        const bool parallelBaseSolver,
        const bool _parallelEdgeTria,
        const bool addInfiniteVertices)
        : Triangulator<D, Precision>(_bounds, _points),
          recursionDepth(_splitter->getRecursionDepth(_threads)),
          parallelEdgeTria(_parallelEdgeTria),
          splitter(std::move(_splitter))
{


    LOG("Initializing DCTriangulator with " << _threads << " threads, recursion depth " << recursionDepth << ", splitter " << splitter->to_string() << std::endl);

    if (addInfiniteVertices) {
        // add infinite points to data set
        auto stats = getPointStats(this->points, this->points);
        for (uint i = 0; i < pow(2, D); ++i) {
            VLOG("Point stats: " << stats.min << " - " << stats.mid << " - "
                 << stats.max << std::endl);

            dPoint<D, Precision> p;
            //p.id = dPoint<D, Precision>::cINF | i;
            p.coords = stats.mid;

            for (uint d = 0; d < D; ++d)
                p.coords[d] =
                        (i & (1 << d) ? 1 : -1) * 2 * SAFETY * (stats.max[d] - stats.min[d]);

            this->points[dPoint<D, Precision>::cINF | i] = std::move(p);
        }
    }

    if (parallelBaseSolver)
        baseTriangulator = std::make_unique<CGALTriangulator<D, Precision, true>>(this->baseBounds, this->points,
                                                                                  gridOccupancy);
    else
        baseTriangulator = std::make_unique<CGALTriangulator<D, Precision, false>>(this->baseBounds, this->points,
                                                                                   gridOccupancy);
}

template<uint D, typename Precision>
void DCTriangulator<D, Precision>::getEdge(const dSimplices<D, Precision> &simplices,
                                           const Partitioning<D, Precision> &partitioning, const uint &partition,
                                           Concurrent_Growing_Point_Ids &edgePoints,
                                           Concurrent_Growing_Simplex_Ids &edgeSimplices) {

    // infinite points to edge
    auto edgePointsHandle = edgePoints.handle();
    for (tIdType k = dPoint<D, Precision>::cINF; k != 0; ++k) {
        edgePointsHandle.insert(k);
    }

    tbb::enumerable_thread_specific<hConcurrent_Growing_Point_Ids,
            tbb::cache_aligned_allocator<hConcurrent_Growing_Point_Ids>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsEdgePointsHandle(std::ref(edgePoints));

    tbb::enumerable_thread_specific<hConcurrent_Growing_Simplex_Ids,
            tbb::cache_aligned_allocator<hConcurrent_Growing_Simplex_Ids>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsEdgeSimplicesHandle(std::ref(edgeSimplices));

    tbb::enumerable_thread_specific<dSimplicesConstHandle<D, Precision>,
            tbb::cache_aligned_allocator<dSimplicesConstHandle<D, Precision>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsSimplicesHandle(std::ref(simplices));


    tMarkType queuedMark = ++simplices.mark;
    tMarkType doneMark = ++simplices.mark;
    ASSERT(queuedMark < doneMark);

    /* Walk along the neighbors,
     * testing for each neighbor whether its circumsphere is within the
     * partition or not
     */

    INDENT

    VTUNE_TASK(IdentifyEdge);
    tbb::parallel_do(simplices.convexHull, [&](const tIdType id,
                                               tbb::parallel_do_feeder<tIdType> &feeder) {

        if (!dSimplex<D, Precision>::isFinite(id))
            return;

        auto &edgePointsHandle = tsEdgePointsHandle.local();
        auto &edgeSimplicesHandle = tsEdgeSimplicesHandle.local();

        auto &simplicesHandle = tsSimplicesHandle.local();

        const auto &simplex = simplicesHandle[id];


        if (simplex.mark == doneMark) {
            return; //already checked
        }
        simplex.mark = doneMark;

        const auto cs = simplex.circumsphere(this->points);
        bool intersectsBounds = false;
        for (uint i = 0; i < partitioning.size() && !intersectsBounds; ++i) {
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
                if (dSimplex<D, Precision>::isFinite(n) && simplicesHandle[n].mark < queuedMark) {
                    // n was not yet inspected
                    simplicesHandle[n].mark = queuedMark;
                    feeder.add(n);
                }
            }
        }
    });

    DEDENT
}

template<uint D, typename Precision>
cWuFaces DCTriangulator<D, Precision>::buildWhereUsed(const dSimplices<D, Precision> &simplices,
                                                      const Simplex_Ids &edgeSimplices) {


    /* We need to build the wuFaces DS for the simplices of the first layer "inward" of the edge
     * plus one more layer, to re-find their neighbors
     */

    VTUNE_TASK(BuildWU);
    cWuFaces wuFaces((D + 1) * (D + 1) * (D + 1) * edgeSimplices.size());
    tbb::enumerable_thread_specific<hcWuFaces,
            tbb::cache_aligned_allocator<hcWuFaces>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsWuFacesHandle(std::ref(wuFaces));

    tbb::enumerable_thread_specific<dSimplicesConstHandle<D, Precision>,
            tbb::cache_aligned_allocator<dSimplicesConstHandle<D, Precision>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsSimplicesHandle(std::ref(simplices));

    tMarkType queuedMark = ++simplices.mark;
    tMarkType doneMark = ++simplices.mark;
    ASSERT(queuedMark < doneMark);

    tbb::parallel_for(edgeSimplices.range(), [&](const auto &range) {

        auto &wuFacesHandle = tsWuFacesHandle.local();

        auto &simplicesHandle = tsSimplicesHandle.local();

        for (const auto &edgeSimplexID : range) {
            const auto &edgeSimplex = simplicesHandle[edgeSimplexID];

            if (edgeSimplex.mark == doneMark) {
                continue; //already checked
            }
            edgeSimplex.mark = doneMark;

            for (const auto &firstLayer : edgeSimplex.neighbors) {
                if (dSimplex<D, Precision>::isFinite(firstLayer) && !edgeSimplices.count(firstLayer)) {
                    // we have an "inward" neighbor, add it to the wuFaces DS
                    if (simplicesHandle[firstLayer].mark < queuedMark) {
                        simplicesHandle[firstLayer].mark = queuedMark;

                        for (uint i = 0; i < D + 1; ++i) {
                            auto facetteHash = simplicesHandle[firstLayer].faceFingerprint(i);
                            wuFacesHandle.insert(facetteHash, firstLayer);
                        }
                    }
                    // now loop over its neighbors, adding the ones not belonging to the edge
                    // TODO maybe we can avoid this
//                    for (const auto &secondLayer : DT[firstLayer].neighbors) {
//                        if (dSimplex<D, Precision>::isFinite(secondLayer) && !edgeSimplices.count(secondLayer)
//                            && simplicesHandle[secondLayer].mark < gMark) {
//                            simplicesHandle[secondLayer].mark = gMark
//
//                            for (uint i = 0; i < D + 1; ++i) {
//                                auto facetteHash = DT[secondLayer].faceFingerprint(i);
//                                wuFacesHandle.insert(facetteHash, secondLayer);
//                            }
//                        }
//                    }
                }
            }
        }
    });

    return wuFaces;
}

template<uint D, typename Precision>
void DCTriangulator<D, Precision>::updateNeighbors(
        dSimplices<D, Precision> &simplices,
        const Concurrent_Id_Vector &toCheck,
        const cWuFaces &wuFaces,
        __attribute__((unused)) const std::string &provenance) {

    INDENT
    const uint saveIndent = LOGGER.getIndent();

    tbb::enumerable_thread_specific<dSimplicesHandle<D, Precision>,
            tbb::cache_aligned_allocator<dSimplicesHandle<D, Precision>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsSimplicesHandle(std::ref(simplices));

    tMarkType queuedMark = ++simplices.mark;
    tMarkType doneMark = ++simplices.mark;
    ASSERT(queuedMark < doneMark);


    auto wuFacesHandle = wuFaces.handle();

#ifndef NDEBUG
    std::atomic<tIdType> checked(0);
    std::atomic<tIdType> updated(0);
#endif

    VTUNE_TASK(UpdateNeighbors);
    tbb::parallel_do(toCheck.range(), [&](const tIdType &id,
                                          tbb::parallel_do_feeder<tIdType> &feeder) {

        LOGGER.setIndent(saveIndent);

        auto &simplicesHandle = tsSimplicesHandle.local();

        dSimplex<D, Precision> &simplex = simplicesHandle[id];

        if (simplex.mark == doneMark) {
            return; //already checked
        }
        simplex.mark = doneMark;

        PLOG("Checking neighbors of " << simplex << std::endl);
#ifndef NDEBUG
        ++checked;
#endif

        INDENT
        for (uint d = 0; d < D + 1; ++d) {

            /* We need to update the the neighbor if
             *  a) it is finite, but doesn't exist in the triangulation anymore or is not a neighbor anymore
             *  b) it is infinite
             */

            if (!dSimplex<D, Precision>::isFinite(simplex.neighbors[d]) // is infinite
                ||
                // not in triangulation anymore or not a neighbor anymore
                (!simplicesHandle.contains(simplex.neighbors[d]) || !simplex.isNeighbor(simplicesHandle[simplex.neighbors[d]]))
                    ) {

                // update needed
                simplex.neighbors[d] = dSimplex<D, Precision>::cINF; // reset to infinite
#ifndef NDEBUG
                bool set = false;
#endif

                auto facetteHash = simplex.faceFingerprint(d);
                auto range = wuFacesHandle.get(facetteHash);
                PLOG("Key: " << facetteHash << " #Results: " << std::distance(range.first, range.second) << std::endl;);
                for (auto it = range.first; it != range.second; ++it) {
                    auto cand = (*it).second;
                    if (cand != simplex.id //its not me
                        && dSimplex<D, Precision>::isFinite(cand) // its finite
                        && simplicesHandle.contains(cand) // its in the triangulation
                        && simplex.isNeighbor(simplicesHandle[cand]) // it is a neighbor
                        && !simplicesHandle[cand].contains(simplex.vertices[d]) // its neighboring at the correct edge
                            ) {
                        PLOG("Neighbor with " << simplicesHandle[cand] << std::endl);
                        ASSERT(!set || cand == simplex.neighbors[d]); //second term is to guard against double entries

                        simplex.neighbors[d] = cand;
                        if(simplicesHandle[cand].mark < queuedMark) {
                            simplicesHandle[cand].mark = queuedMark;
                            feeder.add(cand);
                        }
#ifndef NDEBUG
                        ++updated;
                        set = true;
#endif
                    }
                }
            }
        }
        DEDENT
    });

#ifndef NDEBUG
    VLOG("Updating neighbors - inspected " << checked << " simplices, updated " << updated << " neighbors - " <<
         ((double) updated) / ((double) checked) << " neighbors/simplex" << std::endl);
#endif

    DEDENT
}

template<uint D, typename Precision>
dSimplices<D, Precision> DCTriangulator<D, Precision>::mergeTriangulation(
        std::vector<dSimplices<D, Precision>> &&partialDTs,
        const Simplex_Ids &edgeSimplices,
        dSimplices<D, Precision> &&edgeDT,
        const Partitioning<D, Precision> &partitioning,
        const std::string &provenance) {

    LOG("Merging partial DTs into one triangulation" << std::endl);

    VTUNE_TASK(MergeTriangulations);

    VTUNE_TASK(CombineTriangulations);
    dSimplices<D, Precision> mergedDT = std::move(partialDTs[0]);

    for (uint i = 1; i < partialDTs.size(); ++i) {
        mergedDT.merge(std::move(partialDTs[i]));
    }
    mergedDT.filter(edgeSimplices);

    tbb::enumerable_thread_specific<dSimplicesHandle<D, Precision>,
            tbb::cache_aligned_allocator<dSimplicesHandle<D, Precision>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsMergedDTHandle(std::ref(mergedDT));

    tbb::enumerable_thread_specific<dSimplicesConstHandle<D, Precision>,
            tbb::cache_aligned_allocator<dSimplicesConstHandle<D, Precision>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsEdgeDTHandle(std::ref(edgeDT));

    VTUNE_END_TASK(CombineTriangulations);

    //build the where used datastructure for the faces
    cWuFaces wuFaces = buildWhereUsed(mergedDT, edgeSimplices);

    auto cmpSort =
            [&tsMergedDTHandle](const tIdType &a, const tIdType &b) {
                // contains() doesn't work because we deleted the simplices in mergedDT
                return tsMergedDTHandle.local()[a].fingerprint() < tsMergedDTHandle.local()[b].fingerprint();
            };

    VTUNE_TASK(SortDeletedVertices);
    Id_Vector deletedSimplices(edgeSimplices.begin(), edgeSimplices.end());

    tbb::parallel_sort(deletedSimplices.begin(), deletedSimplices.end(), cmpSort);
    VTUNE_END_TASK(SortDeletedVertices);


    // merge partial DTs and edge DT
    LOG("Merging triangulations" << std::endl);

    VTUNE_TASK(AddBorderSimplices);

    tbb::enumerable_thread_specific<hcWuFaces,
            tbb::cache_aligned_allocator<hcWuFaces>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsWuFacesHandle(std::ref(wuFaces));

    Concurrent_Id_Vector cConvexHull;
    Concurrent_Id_Vector insertedSimplices;

    auto cmpSearch =
            [&tsMergedDTHandle](const tIdType &a, const tHashType &hash) {
                // contains() doesn't work because we deleted the simplices in mergedDT
                return tsMergedDTHandle.local()[a].fingerprint() < hash;
            };

    tbb::parallel_for(edgeDT.range(), [&](auto &r) {

        auto &mergedDTHandle = tsMergedDTHandle.local();
        auto &wuFacesHandle = tsWuFacesHandle.local();

        // we need an explicit iterator loop here for the it < r.end() comparision
        // range-based for loop uses it != r.end() which doesn't work
        for (auto it = r.begin(); it < r.end(); ++it) {
            auto &edgeSimplex = *it;

            // check whether edgeSimplex is completely contained in one partition
            uint p0 = partitioning.partition(edgeSimplex.vertices[0]);
            bool insert = !partitioning[p0].contains(edgeSimplex);

            if (!insert) {
                // the simplex is completely in one partition -> it must have been found
                // before
                auto lb =
                        std::lower_bound(deletedSimplices.begin(), deletedSimplices.end(),
                                         edgeSimplex.fingerprint(), cmpSearch);
                RASSERT(lb == deletedSimplices.end() || *lb != 0);

                for (; lb != deletedSimplices.end() && mergedDTHandle[*lb].fingerprint() == edgeSimplex.fingerprint(); ++lb) {
                    if (edgeSimplex.equalVertices(mergedDTHandle[*lb])) {
                        insert = true;
                    }
                }
            }

            if(insert){
                insertedSimplices.push_back(edgeSimplex.id);

                if(!edgeSimplex.isFinite())
                    cConvexHull.push_back(edgeSimplex.id);

                for (uint d = 0; d < D + 1; ++d) {
                    wuFacesHandle.insert((edgeSimplex.faceFingerprint(d)), edgeSimplex.id);
                }
            } else {
                edgeSimplex.id = dSimplex<D, Precision>::cINF;
            }
        }
    });

    mergedDT.merge(std::move(edgeDT), false);
    mergedDT.convexHull.reserve(mergedDT.convexHull.size() + cConvexHull.size());
    std::move(cConvexHull.begin(), cConvexHull.end(), std::back_inserter(mergedDT.convexHull));

    VTUNE_END_TASK(AddBorderSimplices);

    //ASSERT(DT.countDuplicates() == 0);

    LOG("Updating neighbors" << std::endl);
    updateNeighbors(mergedDT, insertedSimplices, wuFaces, provenance);

    return mergedDT;
}

template<uint D, typename Precision>
dSimplices<D, Precision> DCTriangulator<D, Precision>::_triangulateBase(const Point_Ids &partitionPoints,
                                                                        const dBox<D, Precision> &bounds,
                                                                        const std::string provenance) {

    VTUNE_TASK(TriangulateBase);
    PROFILER_MEAS(provenance.find('e') == std::string::npos ? "basecaseSize" : "edgeBasecaseSize",
                  partitionPoints.size());

    LOGGER.setIndent(provenance.length());

    LOG("triangulateBASE called on level " << provenance << " with "
        << partitionPoints.size() << " points" << std::endl);

    INDENT
    auto dt = baseTriangulator->_triangulate(partitionPoints, bounds, provenance);
    DEDENT

    return dt;
}

template<uint D, typename Precision>
dSimplices<D, Precision> DCTriangulator<D, Precision>::_triangulate(const Point_Ids &partitionPoints,
                                                                    const dBox<D, Precision> &bounds,
                                                                    const std::string provenance) {

    LOGGER.setIndent(provenance.length());

    PROFILER_MEAS("nPoints", partitionPoints.size());

    LOG("triangulateDAC called on level " << provenance << " with "
        << partitionPoints.size() << " points" << std::endl);


    if (provenance.length() - 1 < recursionDepth && partitionPoints.size() > BASE_CUTOFF) {
        VTUNE_TASK(TriangulateRecursive);

        LOG("Recursive case" << std::endl);

        // partition input
        LOG("Partioning" << std::endl);
        INDENT

        VTUNE_TASK(Partitioning);
        const auto partioning = splitter->partition(partitionPoints, this->points, provenance);
        VTUNE_END_TASK(Partitioning);

        DEDENT

        std::vector<dSimplices<D, Precision>> partialDTs;
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
                    partialDTs[i] = _triangulate(partioning[i].points, partioning[i].bounds,
                                                 provenance + std::to_string(i));
                    DEDENT
                }
        );
        VTUNE_END_TASK(TriangulatePartitions);

        VTUNE_TASK(ExtractEdge);

        Concurrent_Growing_Point_Ids edgePointIds(partitionPoints.size() / 4);
        Concurrent_Growing_Simplex_Ids edgeSimplexIds(partitionPoints.size());

        tbb::parallel_for(
                std::size_t(0), partioning.size(),
                [&](const uint i) {
                    LOGGER.setIndent(
                            provenance.length()); // new thread, initialize Logger indent
                    INDENT
                    VLOG("Partition " << i << ": extracting edge" << std::endl);
                    //TODO how to handle edge detection
                    getEdge(partialDTs[i], partioning, i, edgePointIds, edgeSimplexIds);
                    DEDENT
                }
        );
        VTUNE_END_TASK(ExtractEdge);

        PROFILER_MEAS("edgePoints", edgePointIds.handle().size());
        PROFILER_MEAS("edgePointsPerc", edgePointIds.handle().size() / partitionPoints.size());


        PROFILER_MEAS("edgeSimplices", edgeSimplexIds.handle().size());
        PROFILER_MEAS("edgeSimplicesPerc", edgeSimplexIds.handle().size() / ([&]() -> std::size_t {
            std::size_t size = 0;
            for (uint i = 0; i < partialDTs.size(); ++i)
                size += partialDTs[i].simplices.size();
            return size;
        })());

        LOG("Edge has " << edgeSimplexIds.handle().size() << " simplices with "
            << edgePointIds.handle().size() << " points" << std::endl);

        INDENT

        VTUNE_TASK(TriangulateEdge);
        //TODO how to identify the edge triangulation

        dSimplices<D, Precision> edgeDT = parallelEdgeTria ?
                                          _triangulate(Point_Ids(std::move(edgePointIds.data())), bounds,
                                                       provenance + "e")
                                                           :
                                          _triangulateBase(Point_Ids(std::move(edgePointIds.data())), bounds,
                                                           provenance + "e");
        VTUNE_END_TASK(TriangulateEdge);

        PROFILER_MEAS("edgeDT", edgeDT.simplices.size());
        PROFILER_MEAS("edgeDTPerc", edgeDT.simplices.size() / ([&]() -> std::size_t {
            std::size_t size = 0;
            for (uint i = 0; i < partialDTs.size(); ++i)
                size += partialDTs[i].simplices.size();
            return size;
        })());

        DEDENT

        return mergeTriangulation(std::move(partialDTs), Simplex_Ids(std::move(edgeSimplexIds.data())), std::move(edgeDT), partioning,
                                  provenance);

    } else { // base case
        return _triangulateBase(partitionPoints, bounds, provenance);
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

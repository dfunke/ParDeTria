#pragma once

#pragma once

#include <cassert>
#include <random>
#include <algorithm>
#include <set>
#include "../Partitioner.h"
#include "steps.h"

namespace LoadBalancing
{
    std::vector<std::tuple<size_t, size_t>> findPartitionCenterEdges(const Graph& graph, const std::vector<int>& partition) {
        std::vector<std::tuple<size_t, size_t>> result;
        for(size_t n = 0; n+1 < graph.nodeRecords.size(); ++n) {
            auto begin = graph.nodeRecords[n];
            auto end = graph.nodeRecords[n+1];
            for(int i = begin; i < end; ++i) {
                auto adj = graph.adjacency[i];
                if(0 == partition[n] && 0 != partition[adj]) {
                    result.push_back(std::make_tuple(n, adj));
                }
            }
        }
        return result;
    }
    
    template <uint D, typename Precision>
    std::vector<dPoint<D, Precision>> makePartitionCenterPoints(const std::vector<std::tuple<size_t, size_t>>& centerEdges,
                                                                const dPoints<D, Precision>& samplePoints) {
        std::vector<dPoint<D, Precision>> result;
        for(const auto& edge : centerEdges) {
            const auto& point1 = samplePoints[std::get<0>(edge)];
            const auto& point2 = samplePoints[std::get<1>(edge)];
            result.push_back((point1.coords + point2.coords) * 0.5);
        }
        return result;
    }
    
    template <uint D, typename Precision>
    std::tuple<dBox<D, Precision>, dBox<D, Precision>> estimateBoundingBoxes(std::vector<dPoint<D, Precision>> borderPoints,
                                                                             const dBox<D, Precision>& boundingBox) {
        std::tuple<dBox<D, Precision>, dBox<D, Precision>> result;
        if(borderPoints.size() > 0) {
            size_t midIndex = borderPoints.size() / 2;
            std::array<Precision, D> median;
            std::array<Precision, D> borderLength;
            borderLength.fill(0);
            for(size_t d = 0; d < D; ++d) {
                std::sort(borderPoints.begin(), borderPoints.end(),
                                [d](const dPoint<D, Precision>& left, const dPoint<D, Precision>& right) -> bool {
                                    return left.coords[d] < right.coords[d];
                                });
                median[d] = borderPoints[midIndex].coords[d];
                for(size_t i = 0; i+1 < borderPoints.size(); ++i) {
                    borderLength[d] = lenSquared(borderPoints[i].coords - borderPoints[i+1].coords);
                }
                // todo: consider: borderLength[d] /= lenSquared(borderPoints.first() - borderPoints.last());
            }
            
            size_t splitDimension = std::distance(borderLength.begin(), std::min_element(borderLength.begin(), borderLength.end()));
            result = splitBox(boundingBox, splitDimension, median[splitDimension]);
        } else {
            result = splitBox(boundingBox, 0, (boundingBox.low[0] + boundingBox.high[0])/2);
        }
        return result;
    }
    
    template <uint D, typename Precision>
    std::tuple<Point_Ids, Point_Ids> seperatePointIds(const dPoints<D, Precision>& points,
                                                                                  const Point_Ids& ids,
                                                                                  const dBox<D, Precision>& boundingBox) {
        std::tuple<Point_Ids, Point_Ids> result;
        for(auto id : ids) {
            if(dPoint<D, Precision>::isFinite(id)) {
                if(boundingBox.contains(points[id].coords)) {
                    std::get<0>(result).insert(id);
                } else {
                    std::get<1>(result).insert(id);
                }
            }
        }
        
        for (tIdType k = dPoint<D, Precision>::cINF; k != 0; ++k) {
            std::get<0>(result).insert(k);
            std::get<1>(result).insert(k);
        }
        
        return result;
    }
    
    template <uint D, typename Precision>
    struct BinaryBoxEstimatingSamplePartitioner : public Partitioner<D, Precision>
    {
        BinaryBoxEstimatingSamplePartitioner(size_t sampleSize, size_t sampleSeed, size_t minNumOfPoints)
            : mSampleSize(sampleSize), mRand(sampleSeed), mPointsCutoff(minNumOfPoints)
        { }
        PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        const dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds) override
        {
            return buildTreeRecursively(bounds, points, pointIds, 10);
        }
        
        std::string info() const override
        {
            return "binary box-estimating sample partitioner";
        }
        
    private:
        size_t mSampleSize;
        std::mt19937 mRand;
        size_t mPointsCutoff;
        
        PartitionTree<D, Precision> buildTreeRecursively(const dBox<D, Precision>& bounds,
                                        const dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds,
                                        size_t remainingRecursions) {
            PartitionTree<D, Precision> tree;
                
            if(remainingRecursions > 0 && pointIds.size() >= std::max(mPointsCutoff, mSampleSize)) {
                auto sample = generateSample<D, Precision>(mSampleSize, pointIds, mRand);
                dPoints<D, Precision> samplePoints;
                for(auto id : sample) {
                    assert((dPoint<D, Precision>::isFinite(id)));
                    samplePoints.emplace_back(points[id]);
                }
                auto simplices = triangulateSample(bounds, samplePoints);
                auto graph = makeGraph(simplices);
                auto partitioning = partitionGraph(graph, 2, mRand);
                auto centerEdges = findPartitionCenterEdges(graph, partitioning);
                auto centerPoints = makePartitionCenterPoints(centerEdges, samplePoints);
                auto boundingBoxes = estimateBoundingBoxes(centerPoints, bounds);
                auto pointIdsPair = seperatePointIds(points, pointIds, std::get<0>(boundingBoxes));
                
                auto leftSubtree = buildTreeRecursively(std::get<0>(boundingBoxes), points, std::get<0>(pointIdsPair), remainingRecursions - 1);
                auto rightSubtree = buildTreeRecursively(std::get<1>(boundingBoxes), points, std::get<1>(pointIdsPair), remainingRecursions - 1);
            
                tree.bounds = bounds;
                typename PartitionTree<D, Precision>::ChildContainer children{std::move(leftSubtree), std::move(rightSubtree)};
                tree.attachment = std::move(children);
            } else {
                tree.bounds = bounds;
                tree.attachment = pointIds;
            }
            
            return tree;                                    
        }
    };
}

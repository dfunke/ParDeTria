#pragma once

#include <cassert>
#include <random>
#include <algorithm>
#include <set>
#include "SamplePartitioner.h"
#include "Sampler.h"

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
                                                                const std::vector<dPoint<D, Precision>>& samplePoints) {
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
    struct BinaryBoxEstimatingSamplePartitioner : public SamplePartitioner<D, Precision>
    {
        BinaryBoxEstimatingSamplePartitioner(size_t minNumOfPoints, Sampler<D, Precision> sampler)
            : mPointsCutoff(minNumOfPoints), mSampler(std::move(sampler))
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
        
        const Sampling<D, Precision>& sampling() const {
            return mSampling;
        }
        
    private:
        Sampling<D, Precision> mSampling;
        size_t mPointsCutoff;
        Sampler<D, Precision> mSampler;
        
        PartitionTree<D, Precision> buildTreeRecursively(const dBox<D, Precision>& bounds,
                                        const dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds,
                                        size_t remainingRecursions) {
            PartitionTree<D, Precision> tree;
                
            if(remainingRecursions > 0 && pointIds.size() >= mPointsCutoff) {
                mSampling = mSampler(bounds, points, pointIds, 2);
                auto centerEdges = findPartitionCenterEdges(mSampling.graph, mSampling.partition);
                auto centerPoints = makePartitionCenterPoints(centerEdges, mSampling.points);
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

#pragma once
#include <cassert>
#include <memory>
#include <set>
#include <limits>
#include <algorithm>
#include <tuple>
#include <functional>
#include <random>
#include <kaHIP_interface.h>
#include "CGALTriangulator.h"
#include "load_balancing/VectorOperations.h"
#include "load_balancing/BoxUtils.h"

namespace LoadBalancing
{
    template <typename Generator_t>
    std::vector<tIdType> generateSample(size_t sampleSize, size_t pointCloudSize, size_t offset, Generator_t& rand) {
        std::uniform_int_distribution<size_t> dist(0, pointCloudSize - 1);
        
        std::set<tIdType> sample;
        while(sample.size() < sampleSize) {
            sample.insert(offset + dist(rand));
        }
        
        std::vector<tIdType> result;
        std::copy(sample.begin(), sample.end(), std::back_inserter(result));
        return result;
    }
    template <uint D, typename Precision, typename Generator_t>
    std::vector<tIdType> generateSample(size_t sampleSize, const Point_Ids& pointIds, Generator_t& rand) {
        
        std::vector<tIdType> result;
        std::sample(pointIds.begin(), pointIds.end(), std::back_inserter(result), sampleSize, rand);
        result.erase(std::remove_if(result.begin(), result.end(), [](tIdType id){
            return !dPoint<D, Precision>::isFinite(id);
            
        }),
        result.end());
        return result;
    }

    template <uint D, typename Precision>
    dSimplices<D, Precision> triangulateSample(const dBox<D, Precision>& bounds,
                                                const dPoints<D, Precision>& samplePoints) {
        
        auto sp = samplePoints;
        CGALTriangulator<D, Precision> triangulator(bounds, sp);
        return triangulator.triangulate();
    }

    struct Graph {
        std::vector<int> nodeRecords;
        std::vector<int> adjacency;
    };

    Graph adjacencyListToGraph(const std::vector<std::vector<size_t>>& adjacencyList) {
        Graph result;
        
        int currentRecordBegin = 0;
        for(const auto& adjacentPoints : adjacencyList) {
                result.nodeRecords.push_back(currentRecordBegin);
                currentRecordBegin += adjacentPoints.size();
                
            for(auto adjacentPoint : adjacentPoints) {
                result.adjacency.push_back(adjacentPoint);
            }
        }
        result.nodeRecords.push_back(result.adjacency.size());
        
        return result;
    }

    template <uint D, typename Precision>
    Graph makeGraph(const dSimplices<D, Precision>& simplices) {
        std::vector<std::vector<size_t>> adjacencyList;
        for(auto simplex : simplices) {
            for(auto pointId : simplex.vertices) {
                if(dPoint<D, Precision>::isFinite(pointId)) {
                    if(pointId >= adjacencyList.size()) {
                        adjacencyList.resize(pointId + 1);
                    }
                    std::copy_if(simplex.vertices.begin(), simplex.vertices.end(),
                                std::back_inserter(adjacencyList[pointId]), 
                                [](tIdType id) -> bool {
                                    return dPoint<D, Precision>::isFinite(id);
                                });
                }
            }
        }
        return adjacencyListToGraph(adjacencyList);
    }

    template <typename Generator_t>
    std::vector<bool> partitionGraph(Graph& graph, Generator_t& rand) {
        assert(graph.nodeRecords.size() - 1 <= std::numeric_limits<int>::max());
        std::uniform_int_distribution<int> dist(std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
        int n = graph.nodeRecords.size() - 1;
        int edgecut;
        std::vector<int> part(n);
        int nparts = 2;
        double imbalance = 0.05;
        kaffpa(&n,
                nullptr,
                graph.nodeRecords.data(),
                nullptr,
                graph.adjacency.data(),
                &nparts,
                &imbalance,
                true,
                dist(rand),
                FAST,
                &edgecut,
                part.data());
        
        std::vector<bool> result(part.size());
        std::transform(part.begin(), part.end(), result.begin(), [](int x) -> bool {
            return x > 0;
        });
        return result;
    }

    std::vector<std::tuple<size_t, size_t>> findPartitionCenterEdges(const Graph& graph, const std::vector<bool>& partition) {
        std::vector<std::tuple<size_t, size_t>> result;
        for(size_t n = 0; n+1 < graph.nodeRecords.size(); ++n) {
            auto begin = graph.nodeRecords[n];
            auto end = graph.nodeRecords[n+1];
            for(int i = begin; i < end; ++i) {
                auto adj = graph.adjacency[i];
                if(!partition[n] && partition[adj]) {
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
        return splitBox(boundingBox, splitDimension, median[splitDimension]);
    }
    
    template <uint D, typename Precision>
    std::tuple<Point_Ids, Point_Ids> seperatePointIds(const dPoints<D, Precision>& points,
                                                                                  const Point_Ids& ids,
                                                                                  const dBox<D, Precision>& boundingBox) {
        std::tuple<Point_Ids, Point_Ids> result;
        for(auto id : ids) {
            if(boundingBox.contains(points[id].coords)) {
                std::get<0>(result).insert(id);
            } else {
                std::get<1>(result).insert(id);
            }
        }
        return result;
    }

    template <uint D, typename Precision>
    struct Hyperplane {
        dVector<D, Precision> normal;
        dVector<D, Precision> position;
    };
}

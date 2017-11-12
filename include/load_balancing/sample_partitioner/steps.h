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
    std::vector<int> partitionGraph(Graph& graph, size_t numPartitions, Generator_t& rand) {
        assert(numPartitions <= std::numeric_limits<int>::max());
        assert(graph.nodeRecords.size() - 1 <= std::numeric_limits<int>::max());
        
        std::uniform_int_distribution<int> dist(std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
        int n = graph.nodeRecords.size() - 1;
        int edgecut;
        std::vector<int> part(n);
        int nparts = numPartitions;
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
        
        return part;
    }

    template <uint D, typename Precision>
    struct Hyperplane {
        dVector<D, Precision> normal;
        dVector<D, Precision> position;
    };
}

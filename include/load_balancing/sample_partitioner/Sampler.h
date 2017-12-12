#pragma once

#include <cassert>
#include <random>
#include <algorithm>
#include <kaHIP_interface.h>
#include "CGALTriangulator.h"
#include "load_balancing/VectorOperations.h"
#include "load_balancing/BoxUtils.h"

namespace LoadBalancing
{
    struct Graph {
        std::vector<int> nodeRecords;
        std::vector<int> adjacency;
    };
    
    template <uint D, typename Precision>
    struct Sampling
    {
        Graph graph;
        std::vector<int> partition;
        dPoints<D, Precision> points;
    };
    
    template <uint D, typename Precision>
    struct Sampler
    {
        Sampler(size_t sampleSeed, std::function<size_t(size_t)> sampleSize)
            : mRand(sampleSeed), mSampleSize(std::move(sampleSize))
            {}
        
        auto operator()(const dBox<D, Precision>& bounds,
                                    const dPoints<D, Precision>& points,
                                    const Point_Ids& pointIds,
                                    size_t partitionSize)
        {
            auto sample = generateSample(mSampleSize(points.finite_size()), pointIds, mRand);
            dPoints<D, Precision> samplePoints;
            for(auto id : sample) {
                assert((dPoint<D, Precision>::isFinite(id)));
                samplePoints.emplace_back(points[id]);
            }
            auto simplices = triangulateSample(bounds, samplePoints);
            auto graph = makeGraph(simplices);
            auto partition = partitionGraph(graph, partitionSize, mRand);
            return Sampling<D, Precision>{std::move(graph), std::move(partition), std::move(samplePoints)};
        }
        
    private:
        std::mt19937 mRand;
        std::function<size_t(size_t)>  mSampleSize;
        
        template <typename Generator_t>
        static std::vector<tIdType> generateSample(size_t sampleSize, size_t pointCloudSize, size_t offset, Generator_t& rand) {
            std::uniform_int_distribution<size_t> dist(0, pointCloudSize - 1);
            
            std::set<tIdType> sample;
            while(sample.size() < sampleSize) {
                sample.insert(offset + dist(rand));
            }
            
            std::vector<tIdType> result;
            std::copy(sample.begin(), sample.end(), std::back_inserter(result));
            return result;
        }
        
        template <typename Generator_t>
        static std::vector<tIdType> generateSample(size_t sampleSize, const Point_Ids& pointIds, Generator_t& rand) {
            
            std::vector<tIdType> result;
            std::sample(pointIds.begin(), pointIds.end(), std::back_inserter(result), sampleSize, rand);
            result.erase(std::remove_if(result.begin(), result.end(), [](tIdType id){
                return !dPoint<D, Precision>::isFinite(id);
                
            }),
            result.end());
            return result;
        }

        static dSimplices<D, Precision> triangulateSample(const dBox<D, Precision>& bounds,
                                                    const dPoints<D, Precision>& samplePoints) {
            
            auto sp = samplePoints;
            CGALTriangulator<D, Precision> triangulator(bounds, sp);
            return triangulator.triangulate();
        }
        
        static void sanitizeAdjacencyList(std::vector<std::vector<size_t>>& adjacencyList) {
            for(size_t i = 0; i < adjacencyList.size(); ++i) {
                auto& adjacentNodes = adjacencyList[i];
                
                std::sort(begin(adjacentNodes), end(adjacentNodes));
                adjacentNodes.erase(std::unique(begin(adjacentNodes), end(adjacentNodes)), end(adjacentNodes));
                
                auto lowerBound = std::lower_bound(begin(adjacentNodes), end(adjacentNodes), i);
                auto upperBound = std::upper_bound(begin(adjacentNodes), end(adjacentNodes), i);
                assert(std::distance(lowerBound, upperBound) <= 1);
                adjacentNodes.erase(lowerBound, upperBound);
                
            }
        }

        static Graph adjacencyListToGraph(const std::vector<std::vector<size_t>>& adjacencyList) {
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

        static Graph makeGraph(const dSimplices<D, Precision>& simplices) {
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
            sanitizeAdjacencyList(adjacencyList);
            return adjacencyListToGraph(adjacencyList);
        }

        template <typename Generator_t>
        static std::vector<int> partitionGraph(Graph& graph, size_t numPartitions, Generator_t& rand) {
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
    };
}

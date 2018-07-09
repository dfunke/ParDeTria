#pragma once

#include <cassert>
#include <random>
#include <algorithm>
#include <limits>
#include <vector>
#include <kaHIP_interface.h>
#include "CGALTriangulator.h"
#include "load_balancing/VectorOperations.h"
#include "load_balancing/BoxUtils.h"
#include "load_balancing/VectorOperations.h"

namespace LoadBalancing
{
	template <typename In, typename Out>
	struct Mapper
	{
		Mapper(In minIn, In maxIn, Out minOut, Out maxOut)
			: minIn(minIn), maxIn(maxIn), minOut(minOut), maxOut(maxOut)
		{
			assert(minIn < maxIn);
		}

		Out operator()(In x) {
			x = std::min(std::max(x, minIn), maxIn);
			return (x - minIn) * (maxOut - minOut)/(maxIn - minIn) + minOut;
		}
		
	private:
		In minIn, maxIn;
		Out minOut, maxOut;
	};
	
    struct Graph {
        std::vector<int> nodeRecords;
        std::vector<int> adjacency;
	std::vector<int> edgeWeights;
    };
    
    template <uint D, typename Precision>
    struct Sampling
    {
        Graph graph;
        std::vector<int> partition;
        std::vector<dVector<D, Precision>> points;
        dSimplices<D, Precision> simplices;

        Sampling clone() const {
            Sampling res;
            res.graph = graph;
            res.partition = partition;
            res. points = points;

            return res;
        }
    };
    
    template <uint D, typename Precision>
    struct Sampler
    {
	    /**
	     * Precondition for edgeWeight:
	     *   - For all (Precision) x, x >= 0: edgeWeight(x) is defined
	     *   - edgeWeight is strictly monotonic
	     */
        Sampler(size_t sampleSeed, std::function<size_t(size_t)> sampleSize, std::function<int(Precision)> edgeWeight)
            : mRand(sampleSeed), mSampleSize(std::move(sampleSize)), mUniformEdges(false), mEdgeWeight(std::move(edgeWeight))
            {}
        
        Sampler(size_t sampleSeed, std::function<size_t(size_t)> sampleSize)
            : mRand(sampleSeed), mSampleSize(std::move(sampleSize)), mUniformEdges(true)
        {}
        
        auto operator()(const dBox<D, Precision>& bounds,
                                    dPoints<D, Precision>& points,
                                    const Point_Ids& pointIds,
                                    size_t partitionSize)
        {
            auto sample = generateSample(mSampleSize(points.finite_size()), pointIds, mRand);
            auto simplices = triangulateSample(bounds, points, sample);
            auto graph = makeGraph(simplices, points, lenSquared(bounds.low - bounds.high));
            auto partition = partitionGraph(graph, partitionSize, mRand);
	        std::vector<dVector<D, Precision>> samplePointVector;
            for(auto id : sample) {
                samplePointVector.push_back(points[id].coords);
            }
            return Sampling<D, Precision>{std::move(graph), std::move(partition), std::move(samplePointVector), std::move(simplices)};
        }
        
    private:
        std::mt19937 mRand;
        std::function<size_t(size_t)>  mSampleSize;
        bool mUniformEdges;
        std::function<Precision(int)> mEdgeWeight;
        
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
        
        template <typename Generator_t>
        std::vector<tIdType> generateSample(size_t sampleSize, const Point_Ids& pointIds, Generator_t& rand) {
            
            std::vector<tIdType> result;
            std::sample(pointIds.begin(), pointIds.end(), std::back_inserter(result), sampleSize, rand);
            result.erase(std::remove_if(result.begin(), result.end(), [](tIdType id){
                return !dPoint<D, Precision>::isFinite(id);
                
            }),
            result.end());
            return result;
        }

        dSimplices<D, Precision> triangulateSample(const dBox<D, Precision>& bounds,
                                                   dPoints<D, Precision>& points,
                                                   const std::vector<tIdType>& samplePoints) {
            
            CGALTriangulator<D, Precision> triangulator(bounds, points);
            return triangulator.triangulate(samplePoints);
        }
        
        void sanitizeAdjacencyList(std::vector<std::vector<std::pair<size_t, int>>>& adjacencyList) {
	    	assert(adjacencyList.size() > 0);
			assert(adjacencyList[0].size() == 0);

            for(size_t i = 0; i < adjacencyList.size(); ++i) {
                auto& adjacentNodes = adjacencyList[i];
                
                std::sort(begin(adjacentNodes), end(adjacentNodes));
                adjacentNodes.erase(std::unique(begin(adjacentNodes), end(adjacentNodes)), end(adjacentNodes));

		auto compare = [](const std::pair<size_t, int>& left, const std::pair<size_t, int>& right) -> bool {
		   return left.first < right.first;
		};
                
                auto lowerBound = std::lower_bound(begin(adjacentNodes), end(adjacentNodes), std::make_pair(i, 0), compare);
                auto upperBound = std::upper_bound(begin(adjacentNodes), end(adjacentNodes), std::make_pair(i, 0), compare);
                assert(std::distance(lowerBound, upperBound) <= 1);
                adjacentNodes.erase(lowerBound, upperBound);
                
		std::transform(begin(adjacentNodes), end(adjacentNodes), begin(adjacentNodes), [](const auto& node) -> auto { return std::make_pair(node.first - 1, node.second); });
	    }
	    adjacencyList.erase(adjacencyList.begin());

    	    assert((std::all_of(begin(adjacencyList), end(adjacencyList), [&adjacencyList](const auto& adjacentNodes) -> bool {
				    return std::all_of(begin(adjacentNodes), end(adjacentNodes), [&adjacencyList](const auto& node) -> bool {
						    return node.first < adjacencyList.size();
						    });
				    })));
            
        }

        Graph adjacencyListToGraph(const std::vector<std::vector<std::pair<size_t, int>>>& adjacencyList) {
            Graph result;
            
            int currentRecordBegin = 0;
            for(const auto& adjacentPoints : adjacencyList) {
                    result.nodeRecords.push_back(currentRecordBegin);
                    currentRecordBegin += adjacentPoints.size();
                    
                for(auto adjacentPoint : adjacentPoints) {
                    result.adjacency.push_back(adjacentPoint.first);
				    result.edgeWeights.push_back(adjacentPoint.second);
                }
            }
            result.nodeRecords.push_back(result.adjacency.size());
            
            return result;
        }

        Graph makeGraph(const dSimplices<D, Precision>& simplices, const dPoints<D, Precision>& samplePoints,
			Precision maxDistSquared) {
	        Mapper<Precision, int> map(0.0, maxDistSquared, 1, std::numeric_limits<int>::max());

            std::vector<std::vector<std::pair<size_t, int>>> adjacencyList;
            for(auto simplex : simplices) {
                for(auto pointId : simplex.vertices) {
                    if(dPoint<D, Precision>::isFinite(pointId)) {
                        if(pointId >= adjacencyList.size()) {
                            adjacencyList.resize(pointId + 1);
                        }
						for(auto id : simplex.vertices){
							if(dPoint<D, Precision>::isFinite(id)) {
								const Precision distSquared = lenSquared(samplePoints[pointId].coords - samplePoints[id].coords);
								std::pair<tIdType, int> edge;
								edge.first = id;
								edge.second = mUniformEdges ? 1 : map(mEdgeWeight(distSquared));
								adjacencyList[pointId].push_back(edge);
							}
						}
                    }
                }
            }
            sanitizeAdjacencyList(adjacencyList);
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
                    mUniformEdges ? nullptr : graph.edgeWeights.data(),
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

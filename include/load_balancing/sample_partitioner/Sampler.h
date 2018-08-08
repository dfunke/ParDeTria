#pragma once

#include <cassert>
#include <random>
#include <algorithm>
#include <limits>
#include <vector>
#include <unordered_map>
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
			In xClamped = std::clamp(x, minIn, maxIn);
			return (xClamped - minIn) * (maxOut - minOut)/(maxIn - minIn) + minOut;
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
        using DistanceToEdgeWeightFunction = std::function<Precision(Precision)>;

	    /**
	     * Precondition for edgeWeight:
	     *   - For all (Precision) x, x >= 0: edgeWeight(x) is defined
	     *   - edgeWeight is strictly monotonic
	     */
        Sampler(size_t sampleSeed, std::function<size_t(size_t)> sampleSize, int kaffpaMode,
                DistanceToEdgeWeightFunction edgeWeight)
            : mRand(sampleSeed), mSampleSize(std::move(sampleSize)), mUniformEdges(false),
			  mKaffpaMode(kaffpaMode), mEdgeWeight(std::move(edgeWeight))
        {
	        assert(mKaffpaMode == FAST || mKaffpaMode == ECO || mKaffpaMode == STRONG);
        }
        
        Sampler(size_t sampleSeed, std::function<size_t(size_t)> sampleSize, int kaffpaMode)
            : mRand(sampleSeed), mSampleSize(std::move(sampleSize)), mUniformEdges(true),
			  mKaffpaMode(kaffpaMode)
        {
	        assert(mKaffpaMode == FAST || mKaffpaMode == ECO || mKaffpaMode == STRONG);
        }
        
        auto operator()(const dBox<D, Precision>& bounds,
                                    dPoints<D, Precision>& points,
                                    const Point_Ids& pointIds,
                                    size_t partitionSize)
        {
            auto sample = generateSample(mSampleSize(points.finite_size()), pointIds, mRand);
            auto simplices = triangulateSample(bounds, points, sample);
            auto idTranslation = generateIdTranslation(sample);
            auto graph = makeGraph(simplices, points, idTranslation,
                                   lenSquared(bounds.low - bounds.high));
            auto partition = partitionGraph(graph, partitionSize, mRand);

            VTUNE_TASK(CollectPartititioning);
            std::vector<dVector<D, Precision>> samplePointVector(idTranslation.size());
            for(auto id : sample) {
                samplePointVector[idTranslation[id]] = points[id].coords;
            }
            return Sampling<D, Precision>{std::move(graph), std::move(partition),
	            std::move(samplePointVector), std::move(simplices)};
        }
        
    private:
        std::mt19937 mRand;
        std::function<size_t(size_t)>  mSampleSize;
        bool mUniformEdges;
        int mKaffpaMode;
        DistanceToEdgeWeightFunction mEdgeWeight;
        
        template <typename Generator_t>
        Point_Ids generateSample(size_t sampleSize, const Point_Ids& pointIds, Generator_t& rand) {
			VTUNE_TASK(GenerateSample);
            std::vector<tIdType> resultIds;
            std::sample(pointIds.begin(), pointIds.end(), std::back_inserter(resultIds),
                        sampleSize, rand);
            resultIds.erase(std::remove_if(resultIds.begin(), resultIds.end(), [](tIdType id){
                return !dPoint<D, Precision>::isFinite(id);
                
            }),
            resultIds.end());
            Point_Ids result;
            for(auto id : resultIds) {
	            result.insert(id);
			}
            return result;
        }

        dSimplices<D, Precision> triangulateSample(const dBox<D, Precision>& bounds,
                                                   dPoints<D, Precision>& points,
                                                   const Point_Ids& samplePoints) {
            
            VTUNE_TASK(TriangulateSample);
            CGALTriangulator<D, Precision> triangulator(bounds, points);
            return triangulator.triangulateSome(samplePoints, bounds);
        }

		std::unordered_map<tIdType, size_t> generateIdTranslation(const Point_Ids& ids) {
			VTUNE_TASK(SampleIDTranslation);
            std::unordered_map<tIdType, size_t> result;
			size_t count = 0;
			for(auto id : ids) {
				result[id] = count++;
			}
			assert(count == result.size());
			return result;
		}
        
        void sanitizeAdjacencyList(std::vector<std::vector<std::pair<size_t, int>>>& adjacencyList) {
	    	assert(adjacencyList.size() > 0);
    	    assert((std::all_of(begin(adjacencyList), end(adjacencyList),
								[&adjacencyList](const auto& adjacentNodes) -> bool {
								return std::all_of(begin(adjacentNodes), end(adjacentNodes),
				                       [&adjacencyList](const auto& node) -> bool {
				                       return node.first < adjacencyList.size();
				                       });
								})));

            for(size_t i = 0; i < adjacencyList.size(); ++i) {
                auto& adjacentNodes = adjacencyList[i];
                
				auto compare = [](const std::pair<size_t, int>& left,
				                  const std::pair<size_t, int>& right) -> bool {
				   return left.first < right.first;
				};
				
                std::sort(begin(adjacentNodes), end(adjacentNodes), compare);
                adjacentNodes.erase(std::unique(begin(adjacentNodes),
                                                end(adjacentNodes)),
                                    end(adjacentNodes));
	    }

    	    assert((std::all_of(begin(adjacencyList), end(adjacencyList),
								[&adjacencyList](const auto& adjacentNodes) -> bool {
								return std::all_of(begin(adjacentNodes), end(adjacentNodes),
				                       [&adjacencyList](const auto& node) -> bool {
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

        Graph makeGraph(const dSimplices<D, Precision>& simplices,
                        const dPoints<D, Precision>& samplePoints,
                        std::unordered_map<tIdType, size_t> idTranslation,
                        Precision maxDistSquared) {

            VTUNE_TASK(SampleMakeGraph);

            auto lowerWeight = mUniformEdges ? 1 : mEdgeWeight(0.0);
	        auto upperWeight = mUniformEdges ? 2 : mEdgeWeight(1.0);
	        auto [minWeight, maxWeight] = std::minmax(lowerWeight, upperWeight);
	        Mapper<Precision, int> map(minWeight, maxWeight, 1, 99);

            std::vector<std::vector<std::pair<size_t, int>>> adjacencyList(idTranslation.size());
            for(auto simplex : simplices) {
                for(auto pointId : simplex.vertices) {
	                auto i = idTranslation[pointId];
                    if(dPoint<D, Precision>::isFinite(pointId)) {
						for(auto id : simplex.vertices){
							auto j = idTranslation[id];
							if(dPoint<D, Precision>::isFinite(id) && i != j) {
								const Precision distSquared =
									lenSquared(samplePoints[pointId].coords - samplePoints[id].coords);
								const Precision normalizedDistSquared = distSquared / maxDistSquared;
								std::pair<tIdType, int> edge;
								edge.first = j;
								edge.second = mUniformEdges ?
									1 : map(mEdgeWeight(normalizedDistSquared));
								adjacencyList[i].push_back(edge);
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

            VTUNE_TASK(PartitionGraph);

            std::uniform_int_distribution<int> dist(std::numeric_limits<int>::min(),
													std::numeric_limits<int>::max());
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
                    mKaffpaMode,
                    &edgecut,
                    part.data());
            
            return part;
        }
    };
}

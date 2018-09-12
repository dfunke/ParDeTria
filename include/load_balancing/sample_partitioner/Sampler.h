#pragma once

#include <cassert>
#include <random>
#include <algorithm>
#include <numeric>
#include <limits>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <kaHIP_interface.h>
#include <tbb/concurrent_vector.h>
#include <tbb/enumerable_thread_specific.h>
#include "CGALTriangulator.h"
#include "load_balancing/VectorOperations.h"
#include "load_balancing/BoxUtils.h"
#include "load_balancing/VectorOperations.h"
#include "datastructures/LP_Map.hxx"

namespace LoadBalancing
{
	template <typename In, typename Out>
	struct Mapper
	{
		Mapper(In minIn, In maxIn, Out minOut, Out maxOut)
			: minIn(minIn), maxIn(maxIn), minOut(minOut), maxOut(maxOut), scale((maxOut - minOut)/(maxIn - minIn))
		{
			assert(minIn < maxIn);
		}

		Out operator()(In x) {
			In xClamped = std::clamp(x, minIn, maxIn);
			return (xClamped - minIn) * scale + minOut;
		}
		
	private:
		const In minIn, maxIn;
		const Out minOut, maxOut;
		const Out scale;
	};
	
    struct Graph {

        Graph(const std::size_t & n = 1, const std::size_t & m = 1){
            nodeRecords.resize(n+1);
            adjacency.resize(m);
            edgeWeights.resize(m);
        }

        void resizeEdges(const std::size_t & m){
            adjacency.resize(m);
            edgeWeights.resize(m);
        }

        std::vector<int> nodeRecords;
        std::vector<int> adjacency;
	    std::vector<int> edgeWeights;

	    void toFile(const std::string & filename){

	        std::ofstream o(filename);
	        const char SEP = ' ';

	        o << nodeRecords.size() - 1 << SEP << adjacency.size() / 2 << SEP << "1" << std::endl;
	        for(int n = 0; n < static_cast<int>(nodeRecords.size()) - 1; ++n){
	            for(int e = nodeRecords[n]; e < nodeRecords[n+1]; ++e){
	                if(e > nodeRecords[n])
	                    o << SEP;
	                o << adjacency[e] + 1 << SEP << edgeWeights[e];
	            }
	            o << std::endl;
	        }
	    }
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
	        assert(mKaffpaMode == FAST || mKaffpaMode == ECO || mKaffpaMode == STRONG || mKaffpaMode == FASTSOCIALMULTITRY_PARALLEL);
        }
        
        Sampler(size_t sampleSeed, std::function<size_t(size_t)> sampleSize, int kaffpaMode)
            : mRand(sampleSeed), mSampleSize(std::move(sampleSize)), mUniformEdges(true),
			  mKaffpaMode(kaffpaMode)
        {
	        assert(mKaffpaMode == FAST || mKaffpaMode == ECO || mKaffpaMode == STRONG || mKaffpaMode == FASTSOCIALMULTITRY_PARALLEL);
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

        Graph adjacencyListToGraph(const std::vector<std::vector<std::pair<size_t, int>>>& adjacencyList, const std::size_t & n, const std::size_t & m) {
            Graph result(n, m);

            result.nodeRecords.clear();
            result.adjacency.clear();
            result.edgeWeights.clear();

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

        Graph makeGraph_sequential(const dSimplices<D, Precision>& simplices,
                        const dPoints<D, Precision>& samplePoints,
                        std::unordered_map<tIdType, size_t> idTranslation,
                        Precision maxDistSquared) {

            VTUNE_TASK(SampleMakeGraph);

            auto lowerWeight = mUniformEdges ? 1 : mEdgeWeight(0.0);
            auto upperWeight = mUniformEdges ? 2 : mEdgeWeight(1.0);
            auto [minWeight, maxWeight] = std::minmax(lowerWeight, upperWeight);
            Mapper<Precision, int> map(minWeight, maxWeight, 1, 99);

            std::size_t n = idTranslation.size();
            std::size_t m = 0;

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
                                ++m;
                            }
                        }
                    }
                }
            }
            sanitizeAdjacencyList(adjacencyList);
            return adjacencyListToGraph(adjacencyList, n, m);
        }

        Graph makeGraph(const dSimplices<D, Precision>& simplices,
                        const dPoints<D, Precision>& samplePoints,
                        std::unordered_map<tIdType, size_t> idTranslation,
                        const Precision maxDistSquared) {

            VTUNE_TASK(SampleMakeGraph);

            const auto lowerWeight = mUniformEdges ? 1 : mEdgeWeight(0.0);
	        const auto upperWeight = mUniformEdges ? 2 : mEdgeWeight(1.0);
	        const auto [minWeight, maxWeight] = std::minmax(lowerWeight, upperWeight);
	        Mapper<Precision, int> map(minWeight, maxWeight, 1, 99);


            VTUNE_TASK(SampleMakeGraph_COLLECT);

//            tbb::concurrent_unordered_map<std::pair<tIdType, tIdType>, int> edgeMap(D*simplices.upperBound() - simplices.lowerBound());

            typedef Concurrent_LP_Map<tIdType, int> tEdgeMap;
            typedef GrowingHashTable<tEdgeMap> tgEdgeMap;
            typedef GrowingHashTableHandle<tEdgeMap> hEdgeMap;

            tgEdgeMap edgeMap(D * simplices.upperBound() - simplices.lowerBound());
            tbb::enumerable_thread_specific<hEdgeMap,
                    tbb::cache_aligned_allocator<hEdgeMap>,
                    tbb::ets_key_usage_type::ets_key_per_instance> tsEdgeMapHandle(std::ref(edgeMap));

            Graph graph(idTranslation.size());

            tbb::parallel_for(simplices.range(), [&](auto &r) {

                auto &edgeMapHandle = tsEdgeMapHandle.local();

                // we need an explicit iterator loop here for the it < r.end() comparision
                // range-based for loop uses it != r.end() which doesn't work

                for (auto it = r.begin(); it < r.end(); ++it) {
                    auto &simplex = *it;
                    for(auto pointId : simplex.vertices) {
                        auto i = idTranslation[pointId];
                        if(dPoint<D, Precision>::isFinite(pointId)) {
                            for(auto id : simplex.vertices){
                                auto j = idTranslation[id];
                                if(dPoint<D, Precision>::isFinite(id) && i != j) {
                                    const Precision distSquared =
                                            lenSquared(samplePoints[pointId].coords - samplePoints[id].coords);
                                    const Precision normalizedDistSquared = distSquared / maxDistSquared;
                                    std::size_t edge = (i << (sizeof(tIdType) * CHAR_BIT / 2)) + j;
//                                    auto edge = std::make_pair(i, j);
                                    int weight = mUniformEdges ?
                                                  1 : map(mEdgeWeight(normalizedDistSquared));
                                    if(edgeMapHandle.insert(edge, weight)) {
//                                    if(edgeMap.insert(std::make_pair(edge, weight)).second) {
                                        __atomic_fetch_add(&graph.nodeRecords[i+1], 1, __ATOMIC_RELAXED);
                                    }

                                }
                            }
                        }
                    }
                }
            });
            VTUNE_END_TASK(SampleMakeGraph_COLLECT);

            VTUNE_TASK(SampleMakeGraph_PREFIX);
            std::partial_sum(graph.nodeRecords.begin(), graph.nodeRecords.end(), graph.nodeRecords.begin());
            graph.resizeEdges(graph.nodeRecords.back());
            VTUNE_END_TASK(SampleMakeGraph_PREFIX);

            VTUNE_TASK(SampleMakeGraph_ASSIGN);
            std::vector<int> index(idTranslation.size());
            tbb::parallel_for(edgeMap.handle().range(), [&] (const auto & r) {
//            tbb::parallel_for(edgeMap.range(), [&] (const auto & r) {

                for(auto v=r.begin(); v !=r.end(); ++v ){
                    const auto & edge = *v;
                    tIdType i = edge.first >> (sizeof(tIdType) * CHAR_BIT / 2);
                    tIdType j = edge.first - (i << (sizeof(tIdType) * CHAR_BIT / 2));
//                    tIdType i = edge.first.first;
//                    tIdType j = edge.first.second;

                    assert(0 <= i && i < idTranslation.size());
                    assert(0 <= j && j < idTranslation.size());
                    assert(i != j);

                    int base = graph.nodeRecords[i];
                    int offset = __atomic_fetch_add(&index[i], 1, __ATOMIC_RELAXED);

                    graph.adjacency[base+offset] = j;
                    graph.edgeWeights[base+offset] = edge.second;
                }
            });

//            tbb::parallel_for(tbb::blocked_range(std::size_t(0), idTranslation.size()), [&] (const auto &r){
//                for (auto n = r.begin(); n < r.end(); ++n) {
//                    for (int i = graph.nodeRecords[n]; i < graph.nodeRecords[n+1]; ++i) {
//                        auto value = graph.adjacency[i];
//                        auto data = graph.edgeWeights[i];
//                        int hole = i;
//
//                        for (; hole > graph.nodeRecords[n] && value < graph.adjacency[hole - 1]; --hole) {
//                            graph.adjacency[hole] = graph.adjacency[hole - 1];
//                            graph.edgeWeights[hole] = graph.edgeWeights[hole - 1];
//                        }
//
//                        graph.adjacency[hole] = value;
//                        graph.edgeWeights[hole] = data;
//                    }
//                }
//            });

            return graph;
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
                    mKaffpaMode, numPartitions,
                    &edgecut,
                    part.data());
            
            return part;
        }
    };
}

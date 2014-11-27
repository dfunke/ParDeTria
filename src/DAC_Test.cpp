
//std library
#include <fstream>
#include <iostream>
#include <vector>
#include <functional>
#include <set>
#include <stdexcept>
#include <queue>

//debug
#ifndef NDEBUG
#include <csignal>
#endif

//own
#include "Timings.h"
#include "Random.h"
#include "Geometry.h"
#include "Painter.h"
#include "Logger.h"
#include "CGAL_Interface.h"

//**************************
dBox bounds;
dSimplices realDT;
//**************************

dPoints genPoints(const uint n, const dBox & bounds, std::function<tCoordinate()> & dice){

	dPoints points;

	for(uint i = 0; i < n; ++i){
		points[i].id = i;
		for(uint d = 0; d < D; ++d){
			points[i].coords[d] = bounds.coords[d] + bounds.dim[d] * dice();
		}
	}

	//TODO checks for colliding points, straight lines, etc.

	return points;

}

Partition partition(dPoints & points){

	// do mid-point based partitioning for now
	auto stats = getPointStats(points);

	LOG << "Midpoint is " <<  stats.mid << std::endl;

	Partition partition(pow(2,D));

	for(auto & p : points){
		uint part = 0;
		for(uint dim = 0; dim < D; ++dim){
			part |= (p.coords[dim] > stats.mid.coords[dim]) << dim;
		}

		//LOG << "Adding " << p << " to " << part << std::endl;
		partition[part].push_back(p.id);
	}

	//add points at the extermities
	for(uint i = 0; i < partition.size(); ++i){
		auto stats = getPointStats(points, &partition[i]);
		VLOG << "Partition " << i << ": " << stats.min << " - " << stats.mid << " - " << stats.max << std::endl;
		for(uint k = 0; k < pow(2,D); ++k){

			//low point
			dPoint p = stats.mid;
			p.id = dPoint::cINF | i << D | k;

			for(uint d = 0; d < D; ++d)
				p.coords[d] += (k & (1 << d) ? 1 : -1) * 2 * (stats.max.coords[d] - stats.min.coords[d]);

			points.push_back(p);
			partition[i].push_back(p.id);
		}
	}

	return partition;

}

dSimplices getEdge(const dSimplices & simplices){

	dSimplices edgeSimplices;
	std::set<uint> edgeIdx;

	for(const auto & s : simplices){
		/*if(!s.isFinite()){
			//we have an infinte vertex, at least one of its neighbors must be finite
			for(uint i = 0; i < D+1; ++i){
				if(dSimplex::isFinite(s.neighbors[i])
						&& simplices[norm(s.neighbors[i])].isFinite()
						&& edgeIdx.insert(s.neighbors[i]).second)
					//the neighbor is finite and has not been added to the edge simplices yet
					edgeSimplices.push_back(simplices[norm(s.neighbors[i])]);
			}
		}*/
		if(!s.isFinite())
			edgeSimplices.push_back(s);
	}

	return edgeSimplices;

}

dPoints extractPoints(const dSimplices & simplices, const dPoints & inPoints){

	dPoints outPoints;
	std::set<uint> idx;

	for(const auto & s : simplices){
		for(uint i = 0; i < D+1; ++i){
			if(dPoint::isFinite(s.vertices[i]) && idx.insert(s.vertices[i]).second)
				outPoints.push_back(inPoints[s.vertices[i]]);
		}
	}

	return outPoints;

}

void updateNeighbors(dSimplices & simplices, dPoints & points){

	INDENT
	for(dSimplex & simplex : simplices){

		PLOG << "Updating neighbors of " << simplex << std::endl;

#ifndef NDEBUG
		dSimplex saveSimplex = simplex;
#endif

		simplex.neighbors.clear();

		INDENT
		for(uint v = 0; v < D+1; ++v){

			//for every point, look where else its used
			//if the other simplex shares a second point -> it is a neighbor

			const dPoint & vertex = points[simplex.vertices[v]];

			if(IS_PROLIX){
				PLOG << "Vertex " << vertex << " used in ";
				LOGGER.logContainer(vertex.simplices, Logger::Verbosity::PROLIX);
			}

			INDENT
			for(const uint u : vertex.simplices){

				if(u != simplex.id && simplices.contains(u) && simplex.isNeighbor(simplices[u])){

					PLOG << "Neighbor with " << simplices[u] << std::endl;

					simplex.neighbors.insert(u);

					LOGGER.logContainer(simplex.neighbors, Logger::Verbosity::PROLIX);

					assert(simplex.neighbors.size() <= D+1);
					//u will be updated in its own round;
				}

			}
			DEDENT

		}
		DEDENT

		if(!(simplex.neighbors.size() > 0 && simplex.neighbors.size() <= D+1)) {
			LOG << "Error: not enough neighbors for simplex " << simplex << std::endl;
			Painter errorPainter(bounds);
			errorPainter.draw(points);
			errorPainter.drawPartition(points);

			errorPainter.setColor(0,0,0,0.4);
			errorPainter.draw(simplices, points, true);

			errorPainter.setColor(1,0,0); //simplex in red
			errorPainter.draw(simplex, points, true);

			errorPainter.setColor(1,1,0,0.4);//neighbors in yellow
			errorPainter.drawNeighbors(simplex, simplices, points, true);
#ifndef NDEBUG
			if(!simplex.equals(saveSimplex)) {
				LOG << "Error: was before " << saveSimplex << std::endl;

				errorPainter.setColor(0,0,1); //old simplex in blue
				errorPainter.draw(saveSimplex, points, true);

				errorPainter.setColor(0,1,1,0.4);//old simplex neighbors in cyan
				errorPainter.drawNeighbors(saveSimplex, simplices, points, true);
			}
#endif
			errorPainter.savePNG("errors/neighbors_" + std::to_string(simplex.id) + ".png");
		}


		//assert(simplex.neighbors.size() > 0 && simplex.neighbors.size() <= D+1);

		//TODO it might be better NOT to save the infinite simplex as neighbor
		if(simplex.neighbors.size() < D+1){ //if it is a simplex at the border, it might have one infinite simplex as neighbor
			simplex.neighbors.insert(uint(dSimplex::cINF));
		}
	}
	DEDENT

}

void eliminateDuplicates(dSimplices & DT, dPoints & points) {

	Ids saveSimplices(DT.begin_keys(), DT.end_keys());

	for(const auto & s : saveSimplices){
		if(!DT.contains(s))
			continue;

		auto duplicates = DT.findSimplices(DT[s].vertices, true);
		if(duplicates.size() > 1){
			//keep the one with the minimum index
			auto minSimplex = *std::min_element(duplicates.begin_keys(), duplicates.end_keys());
			for(const auto & del : duplicates){
				if(del.id == minSimplex)
					continue; //keep the lowest one

				//delete simplex from where-used list
				for(const auto & p : del.vertices)
					points[p].simplices.erase(del.id);
				DT.erase(del.id);
			}
		}
	}

}

void mergeTriangulation(dSimplices & DT, const dSimplices & edgeDT,
		const Partition & partitioning, dPoints & points) {

	auto partitionPoint =
			[&] (const uint & point) -> uint {

				for(uint p = 0; p < partitioning.size(); ++p){
					if(std::find(partitioning[p].begin(), partitioning[p].end(), point) != partitioning[p].end())
						return p;
				}

				throw std::out_of_range(std::to_string(point) + "not found");

			};

	auto partition =
			[&] (const dSimplex & simplex) -> uint {

				for(uint p = 0; p < partitioning.size(); ++p){
					if(std::all_of(simplex.vertices.begin(), simplex.vertices.end(),
							[&] (const uint v){
								return std::find(partitioning[p].begin(), partitioning[p].end(), v) != partitioning[p].end();
							}))
						return p;
				}

				uint p0 = partitionPoint(simplex.vertices[0]);
				uint p1 = partitionPoint(simplex.vertices[1]);
				uint p2 = partitionPoint(simplex.vertices[2]);

				throw std::out_of_range("Partition of simplex " + std::to_string(simplex.id) + " ambiguous: "
						+ std::to_string(p0) + std::to_string(p1) + std::to_string(p2));

			};

	auto partitionContains =
			[&] (const dSimplex & simplex, const uint partition) -> bool {
				for(uint i = 0; i < D+1; ++i) {
					if(std::find(partitioning[partition].begin(),
							     partitioning[partition].end(),
								 simplex.vertices[i]) == partitioning[partition].end())
					return false;
				}

				return true;
			};

	typedef std::pair<tCoordinate, uint> tCandidate;
	auto cmpCandidates = [] (const tCandidate & a, const tCandidate & b) { return a.first > b.first; };
	typedef std::priority_queue<tCandidate, std::vector<tCandidate>, decltype(cmpCandidates)> tCandidateQueue;

	//we use the overflow of the uint to zero to abort the loop
	for (uint infVertex = dPoint::cINF; infVertex != 0; ++infVertex) {
		//loop over all infinite vertices contained in DT
		assert(points.contains(infVertex));

		if(IS_PROLIX){
			PLOG << "Vertex " << points[infVertex] << " used in ";
			LOGGER.logContainer(points[infVertex].simplices, Logger::Verbosity::PROLIX);
		}

		auto saveSimplices = points[infVertex].simplices; //needs to be copied since we modify the original set

		INDENT
		for (auto & i : saveSimplices) {

			if(!DT.contains(i))
				continue; //TODO is this correct?

			dSimplex & s = DT[i];
			VLOG<< "Edge Simplex " << s << std::endl;

			Painter detailPainter(bounds);
			detailPainter.draw(points);
			detailPainter.drawPartition(points);

			detailPainter.setColor(1,1,0,.4); //simplex neighbors in yellow
			detailPainter.drawNeighbors(s, DT, points, true);

			detailPainter.setColor(1,0,0); // to be replaced simplex in red
			detailPainter.draw(s, points, true);

			detailPainter.savePNG("details/merging_" + std::to_string(infVertex) + "_" + std::to_string(s.id) + ".png");

			uint sp = partition(s); //partition number of s

			assert(!s.isFinite());

			//the simplex should have D finite points
			dPoints finitePoints;
			for(uint v = 0; v < D+1; ++v) {
				if(dPoint::isFinite(s.vertices[v]))
				finitePoints.push_back(points[s.vertices[v]]);
			}

			auto realSimplices = realDT.findSimplices(finitePoints);
			detailPainter.setColor(0,0,0,0.4); //real simplices in grey
			detailPainter.draw(realSimplices, points);

			detailPainter.savePNG("details/merging_" + std::to_string(infVertex) + "_" + std::to_string(s.id) + ".png");


			if(finitePoints.size() < D){
				VLOG << "Found only " << finitePoints.size() << " finite points, skipping" << std::endl;
				continue;
				//TODO: handle less than D finite points
			}

			tCandidateQueue candPQ(cmpCandidates);

			VLOG << "Candidate mergers" << std::endl;
			INDENT
			for(const auto & x : edgeDT) {
				if(x.containsAll(finitePoints) && !partitionContains(x, sp)) {
					auto cc = x.circumcircle(points);
					candPQ.emplace(cc.radius, x.id);
					VLOG << x << " - r = " << cc.radius << std::endl;
				}
			}
			DEDENT

			if(candPQ.empty()){
				INDENT
				VLOG << "NO candidates found" << std::endl;
				DEDENT
				continue;
			}

			//take the canddiate with the smallest circumcircle
			const auto & mergeSimplex = edgeDT[candPQ.top().second];


			if(IS_PROLIX){

				dSimplices candidates;
				//gather candidates
				while(!candPQ.empty()){
					candidates.push_back(edgeDT[candPQ.top().second]);
					candPQ.pop();
				}

				detailPainter.setColor(0,1,1,.4); //candidate neighbors in cyan
				detailPainter.drawNeighbors(candidates, edgeDT, points, true);

				detailPainter.setColor(0, 0, 1,0.4); //candidate in blue
				detailPainter.draw(candidates, points, true);

				detailPainter.setColor(0,1,0); //accepted candidate in green
				detailPainter.draw(mergeSimplex, points,true);

				detailPainter.savePNG("details/merging_" + std::to_string(infVertex) + "_" + std::to_string(s.id) + ".png");
			}

			//find infinite point in original simplex
			uint infIdx = 0;
			while(dPoint::isFinite(s.vertices[infIdx])) ++infIdx;
			assert(infIdx < D+1);

			//find differing point in mergeSimplex
			uint diffIdx = 0;
			while(s.contains(mergeSimplex.vertices[diffIdx])) ++diffIdx;
			assert(diffIdx < D+1);

			//do the change
			VLOG << "Changing " << infIdx << ": " << points[s.vertices[infIdx]] << " to "
			<< diffIdx << ": " << points[mergeSimplex.vertices[diffIdx]] << std::endl;;

			points[s.vertices[infIdx]].simplices.erase(s.id);
			points[mergeSimplex.vertices[diffIdx]].simplices.erase(mergeSimplex.id);

			s.vertices[infIdx] = mergeSimplex.vertices[diffIdx];

			points[s.vertices[infIdx]].simplices.insert(s.id);


		}
		DEDENT
	}

	updateNeighbors(DT, points);

}

int main(int argc, char* argv[]) {

	if(argc == 2)
		LOGGER.setLogLevel(static_cast<Logger::Verbosity>(std::stoi(argv[1])));
	else
		LOGGER.setLogLevel(Logger::Verbosity::LIVE);

	uint N = 1e2;

	for(uint i = 0; i < D; ++i){
		bounds.coords[i] = 0;
		bounds.dim[i] = 100;
	}

	Painter basePainter(bounds);

	std::uniform_real_distribution<tCoordinate> distribution(0,1);
	std::function<tCoordinate()> dice = std::bind (distribution, generator);

	auto points = genPoints(N, bounds, dice);

	LOG << "Partioning" << std::endl;
	INDENT
	auto partioning = partition(points);
	DEDENT

	basePainter.draw(points);
	basePainter.drawPartition(points);
	basePainter.savePNG("01_points.png");

	for(uint i = 0; i < partioning.size(); ++i){
		PLOG << "Partition " << i << ": ";
		for(auto & p : partioning[i])
			CONT << p << " ";
		CONT << std::endl;
	}

	LOG << "Real triangulation" << std::endl;
	INDENT
	realDT = delaunayCgal(points,nullptr, true);
	LOG << "Real triangulation contains " << realDT.size() << " tetrahedra" << std::endl << std::endl;
	DEDENT

	Painter paintRealDT = basePainter;

	paintRealDT.setColor(0, 0, 1);
	paintRealDT.draw(realDT, points);
	paintRealDT.setColor(0, 0, 0);

	paintRealDT.savePNG("00_realDT.png");

	std::vector<dSimplices> partialDTs;
	partialDTs.resize(partioning.size());
	Painter paintPartialDTs = basePainter;

	INDENT
	for(uint i = 0; i < partioning.size(); ++i){
		LOG << "Partition " << i << std::endl;
		INDENT
		auto dt = delaunayCgal(points, &partioning[i]);
		LOG << "Triangulation " << i << " contains " << dt.size() << " tetrahedra" << std::endl << std::endl;
		DEDENT

		LOG << "Verifying CGAL triangulation" << std::endl;
		dt.verify(points.project(partioning[i]));

#ifndef NDEBUG
		auto saveDT = dt;
#endif

		LOG << "Updating neighbors" << std::endl;
		updateNeighbors(dt, points);

		LOG << "Verifying updated triangulation" << std::endl;
		dt.verify(points.project(partioning[i]));

		assert(saveDT == dt); //only performed if not NDEBUG

		partialDTs[i] = dt;

		paintPartialDTs.setColor(Painter::tetradicColor(i) ,0.4);
		paintPartialDTs.draw(dt, points);

		Painter paintPartialDT = basePainter;
		paintPartialDT.draw(dt, points, true);
		paintPartialDT.savePNG("02_partialDT_" + std::to_string(i) + ".png");

	}
	DEDENT

	paintPartialDTs.savePNG("02_partialDTs.png");

	LOG << "Extracting edges" << std::endl;
	INDENT

	dPoints edgePoints;
	Painter paintEdges = paintPartialDTs;
	paintEdges.setColor(1, 0, 0);

	for(uint i = 0; i < partioning.size(); ++i){

		auto edge = getEdge(partialDTs[i]);
		auto ep = extractPoints(edge, points);
		//points are in different partitions, there can be no overlap
		edgePoints.insert(ep.begin(), ep.end());

		VLOG << "Edge has " << edge.size() << " simplices with " << ep.size() << " points" << std::endl;

		paintEdges.draw(edge, points, false);
	}
	DEDENT

	paintEdges.draw(edgePoints);
	paintEdges.savePNG("03_edgeMarked.png");

	LOG << "Triangulating edges" << std::endl;
	INDENT
	auto edgeDT = delaunayCgal(edgePoints);

#ifndef NDEBUG
	auto saveEdgeDT = edgeDT;
#endif
	updateNeighbors(edgeDT, edgePoints);
	edgeDT.verify(edgePoints);

	assert(saveEdgeDT == edgeDT); //only performed if not NDEBUG
	LOG << "Edge triangulation contains " << edgeDT.size() << " tetrahedra" << std::endl << std::endl;
	DEDENT

	Painter paintEdgeDT = basePainter;

	paintEdgeDT.setColor(0, 1, 0);
	paintEdgeDT.draw(edgeDT, points);
	paintEdgeDT.setColor(0, 0, 0);

	paintEdgeDT.savePNG("04_edgeDT.png");

	Painter paintMerged = basePainter;
	paintMerged.setColor(0,1,0,0.4);

	LOG << "Merging partial DTs into one triangulation" << std::endl;
	dSimplices mergedDT;
	for(uint i = 0; i < partialDTs.size(); ++i)
		mergedDT.insert(partialDTs[i].begin(), partialDTs[i].end());

	mergeTriangulation(mergedDT, edgeDT, partioning, points);

	mergedDT.verify(points); //this should NOT fail;

	paintMerged.draw(mergedDT, points);

	paintMerged.setColor(1,0,0);
	paintMerged.draw(edgePoints);
	paintMerged.savePNG("05_mergedDT.png");

	paintMerged.setColor(0,0,0, 0.1);
	paintMerged.draw(realDT, points);

	paintMerged.savePNG("05_mergedDT_overlay.png");


	LOG << "Finished" << std::endl;
}

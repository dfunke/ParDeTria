
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

Partitioning partition(dPoints & points){

	// do mid-point based partitioning for now
	auto stats = getPointStats(points);

	LOG << "Midpoint is " <<  stats.mid << std::endl;

	Partitioning partitioning;
	for(uint i = 0; i < pow(2,D); ++i){
		partitioning[i].id = i;

		for(uint d = 0; d < D; ++d){
			partitioning[i].bounds.coords[d] = i & (1 << d) ? stats.mid.coords[d]
														    : stats.min.coords[d];
			partitioning[i].bounds.dim[d]    = i & (1 << d) ? stats.max.coords[d] - stats.mid.coords[d]
														    : stats.mid.coords[d] - stats.min.coords[d];
		}
	}

#ifndef NDEBUG
	auto inPartition = [&] (const dPoint & p, const uint partition) -> bool {
		for(uint i = 0; i < D; ++i)
			if(!(partitioning[partition].bounds.coords[i] <= p.coords[i] && p.coords[i] <= partitioning[partition].bounds.coords[i] + partitioning[partition].bounds.dim[i]))
				return false;
		return true;
	};
#endif

	for(auto & p : points){
		uint part = 0;
		for(uint dim = 0; dim < D; ++dim){
			part |= (p.coords[dim] > stats.mid.coords[dim]) << dim;
		}

		assert(inPartition(p, part));

		//LOG << "Adding " << p << " to " << part << std::endl;
		partitioning[part].points.insert(p.id);
	}

	//add points at the extermities
	for(uint i = 0; i < partitioning.size(); ++i){
		auto stats = getPointStats(points, &partitioning[i].points);
		VLOG << "Partition " << i << ": " << stats.min << " - " << stats.mid << " - " << stats.max << std::endl;
		for(uint k = 0; k < pow(2,D); ++k){

			//low point
			dPoint p = stats.mid;
			p.id = dPoint::cINF | i << D | k;

			for(uint d = 0; d < D; ++d)
				p.coords[d] += (k & (1 << d) ? 1 : -1) * 2 * (stats.max.coords[d] - stats.min.coords[d]);

			points.push_back(p);
			partitioning[i].points.insert(p.id);
		}
	}

	return partitioning;

}

Ids getEdge(const dSimplices & simplices, const Partition & partition, const dPoints & points){

	Ids edgeSimplices;

	//we use the overflow of the uint to zero to abort the loop
	for (uint infVertex = (dPoint::cINF | (partition.id << D));
			  infVertex != (dPoint::cINF | (partition.id << D)) + (1 << D); ++infVertex) {

		assert(std::find(partition.points.begin(), partition.points.end(), infVertex) != partition.points.end()); //the infinite point must be in the partition
		assert(points.contains(infVertex));

		PLOG << "Vertex " << points[infVertex] << std::endl;

		INDENT
		for(const auto & s : points[infVertex].simplices){
			if(!simplices.contains(s))
				continue;

			PLOG << "Adding " << simplices[s] << " to edge" << std::endl;
			edgeSimplices.insert(simplices[s].id);

			/* Walk along the neighbors,
			 * testing for each neighbor whether its circumsphere is within the partition or not
			 */
			std::deque<uint> wq;
			std::set<uint> wqa;

			//work queue of simplices to check for circum circle criterian
			wq.insert(wq.end(), simplices[s].neighbors.begin(), simplices[s].neighbors.end());

			//keep track of already inspected simplices
			wqa.insert(simplices[s].neighbors.begin(), simplices[s].neighbors.end());
			INDENT
			while(!wq.empty()){
				PLOG << "Inspecting " << wq.size() << " simplices" << std::endl;
				uint x = wq.front(); wq.pop_front();
				if(simplices.contains(x) && !partition.bounds.contains(simplices[x].circumcircle(points))){
					PLOG << "Adding " << simplices[x] << " to edge -> circumcircle criterion" << std::endl;
					edgeSimplices.insert(simplices[x].id);

					for(const auto & n : simplices[x].neighbors){
						if(wqa.insert(n).second){
							//n was not yet inspected
							wq.push_back(n);
						}
					}
				}
			}
			DEDENT
		}
		DEDENT
	}

	return edgeSimplices;

}

Ids extractPoints(const Ids & edgeSimplices, const dSimplices & simplices){

	Ids outPoints;
	std::set<uint> idx;

	for(const auto & id : edgeSimplices){

		assert(simplices.contains(id));

		for(uint i = 0; i < D+1; ++i){
			if(dPoint::isFinite(simplices[id].vertices[i]) && idx.insert(simplices[id].vertices[i]).second)
				outPoints.insert(simplices[id].vertices[i]);
		}
	}

	//add the extreme infinite points to the set
	for(uint i = 0; i < std::pow(2,D); ++i)
		outPoints.insert(dPoint::cINF | i << D | i);

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

					//assert(simplex.neighbors.size() <= D+1);
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

	uint inSimplices = DT.size();
	LOG << "Eliminating duplicates: " << inSimplices << " simplices" << std::endl;

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

	LOG << inSimplices - DT.size() << " duplicates found" << std::endl;

}

void mergeTriangulation(dSimplices & DT, const Ids & edgeSimplices, const dSimplices & edgeDT,
		const Partitioning & partitioning, dPoints & points) {

	auto partitionPoint =
			[&] (const uint & point) -> uint {

				for(uint p = 0; p < partitioning.size(); ++p){
					if(std::find(partitioning[p].points.begin(), partitioning[p].points.end(), point) != partitioning[p].points.end())
						return p;
				}

				throw std::out_of_range(std::to_string(point) + "not found");

			};

	auto partitionContains =
			[&] (const dSimplex & simplex) -> bool {
					uint p0 = partitionPoint(simplex.vertices[0]);
					uint p1 = partitionPoint(simplex.vertices[1]);
					uint p2 = partitionPoint(simplex.vertices[2]);

					return p0 == p1 && p1 == p2;
			};

	typedef std::pair<tCoordinate, uint> tCandidate;
	auto cmpCandidates = [] (const tCandidate & a, const tCandidate & b) { return a.first > b.first; };
	typedef std::priority_queue<tCandidate, std::vector<tCandidate>, decltype(cmpCandidates)> tCandidateQueue;

	PainterBulkWriter paintWriter;

	for (auto & i : edgeSimplices) {

		if(!DT.contains(i))
			continue; //TODO is this correct?

		dSimplex & s = DT[i];
		VLOG<< "Edge Simplex " << s << std::endl;

		std::string png = "details/merging_" + std::to_string(s.id) + ".png";

		if(IS_PROLIX){
			paintWriter[png] = Painter(bounds);
			paintWriter[png].draw(points);
			paintWriter[png].drawPartition(points);

			paintWriter[png].setColor(1,1,0,.4); //simplex neighbors in yellow
			paintWriter[png].drawNeighbors(s, DT, points, true);

			paintWriter[png].setColor(1,0,0); // to be replaced simplex in red
			paintWriter[png].draw(s, points, true);
		}

		//assert(!s.isFinite());

		//the simplex should have D finite points
		dPoints finitePoints;
		for(uint v = 0; v < D+1; ++v) {
			if(dPoint::isFinite(s.vertices[v]))
			finitePoints.push_back(points[s.vertices[v]]);
		}

		if(IS_PROLIX){
			auto realSimplices = realDT.findSimplices(finitePoints);
			paintWriter[png].setColor(0,0,0,0.4); //real simplices in grey
			paintWriter[png].draw(realSimplices, points);
		}


		/*if(finitePoints.size() < D){
			VLOG << "Found only " << finitePoints.size() << " finite points, skipping" << std::endl;
			continue;
			//TODO: handle less than D finite points
		}*/

		tCandidateQueue candPQ(cmpCandidates);

		VLOG << "Candidate mergers" << std::endl;
		INDENT
		if(s.isFinite()){
			auto cc = s.circumcircle(points);
			candPQ.emplace(cc.radius, s.id);
			VLOG << s << " - r = " << cc.radius << std::endl;
		}

		for(const auto & x : edgeDT) {
			if(x.containsAll(finitePoints) && !partitionContains(x)) {
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

		auto drawCandidates = [&] () {
			dSimplices candidates;
			//gather candidates
			while(!candPQ.empty()){

				if(edgeDT.contains(candPQ.top().second))
					candidates.push_back(edgeDT[candPQ.top().second]);

				candPQ.pop();
			}

			paintWriter[png].setColor(0,1,1,.4); //candidate neighbors in cyan
			paintWriter[png].drawNeighbors(candidates, edgeDT, points, true);

			paintWriter[png].setColor(0, 0, 1,0.4); //candidate in blue
			paintWriter[png].draw(candidates, points, true);
		};

		//take the canddiate with the smallest circumcircle
		if(!edgeDT.contains(candPQ.top().second)){
			//the topmost simplex must be the one from DT itself -> nothing to be merged
			assert(DT.contains(candPQ.top().second) && s == DT[candPQ.top().second]);

			if(IS_PROLIX)
				drawCandidates();

			continue;
		}

		const dSimplex & mergeSimplex = edgeDT[candPQ.top().second];

		if(IS_PROLIX){
			drawCandidates();
			paintWriter[png].setColor(0,1,0); //accepted candidate in green
			paintWriter[png].draw(mergeSimplex, points,true);
		}

		VLOG << "Merging " << s << " and " << mergeSimplex << std::endl;

		//find the points they are differing in
		INDENT
		for(uint i = 0; i < D+1; ++i){
			if(!mergeSimplex.contains(s.vertices[i])){
				//vertex i of S is NOT in mergeSimplex
				//find partner to exchange it with

				uint j = 0;
				while(s.contains(mergeSimplex.vertices[j])) ++j;
				assert(j < D+1);

				//do the change
				VLOG << "Changing " << i << ": " << points[s.vertices[i]] << " to "
				<< j << ": " << points[mergeSimplex.vertices[j]] << std::endl;;

				//remove the simplices from the where-used list of the points
				points[s.vertices[i]].simplices.erase(s.id);
				points[mergeSimplex.vertices[j]].simplices.erase(mergeSimplex.id);

				//perform the change
				s.vertices[i] = mergeSimplex.vertices[j];

				//update where-used list
				points[s.vertices[i]].simplices.insert(s.id);
			}
		}
		DEDENT
	}
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

	for(auto & p : partioning)
		PLOG << "Partition " << p << std::endl;;

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
		auto dt = delaunayCgal(points, &partioning[i].points);
		LOG << "Triangulation " << i << " contains " << dt.size() << " tetrahedra" << std::endl << std::endl;
		DEDENT

		LOG << "Verifying CGAL triangulation" << std::endl;
		dt.verify(points.project(partioning[i].points));

#ifndef NDEBUG
		auto saveDT = dt;
#endif

		LOG << "Updating neighbors" << std::endl;
		updateNeighbors(dt, points);

		LOG << "Verifying updated triangulation" << std::endl;
		dt.verify(points.project(partioning[i].points));

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

	Ids edgePointIds;
	Ids edgeSimplexIds;
	Painter paintEdges = paintPartialDTs;
	paintEdges.setColor(1, 0, 0);

	for(uint i = 0; i < partioning.size(); ++i){

		//points are in different partitions, there can be no overlap

		auto edge = getEdge(partialDTs[i], partioning[i], points);
		edgeSimplexIds.insert(edge.begin(), edge.end());

		auto ep = extractPoints(edge, partialDTs[i]);
		edgePointIds.insert(ep.begin(), ep.end());

		VLOG << "Edge has " << edge.size() << " simplices with " << ep.size() << " points" << std::endl;

		paintEdges.draw(partialDTs[i].project(edge), points, false);
	}
	DEDENT

	paintEdges.draw(points.project(edgePointIds));
	paintEdges.savePNG("03_edgeMarked.png");

	LOG << "Triangulating edges" << std::endl;
	INDENT
	auto edgeDT = delaunayCgal(points, &edgePointIds);

#ifndef NDEBUG
	auto saveEdgeDT = edgeDT;
#endif
	updateNeighbors(edgeDT, points);
	edgeDT.verify(points.project(edgePointIds));

	assert(saveEdgeDT == edgeDT); //only performed if not NDEBUG
	LOG << "Edge triangulation contains " << edgeDT.size() << " tetrahedra" << std::endl << std::endl;
	DEDENT

	Painter paintEdgeDT = basePainter;

	paintEdgeDT.setColor(0, 1, 0);
	paintEdgeDT.draw(edgeDT, points, true);
	paintEdgeDT.setColor(1, 0, 0);
	paintEdgeDT.draw(points.project(edgePointIds));

	paintEdgeDT.savePNG("04_edgeDT.png");

	auto paintMerging = [&] (const dSimplices & dt, const std::string & name){
		Painter painter = basePainter;

		painter.setColor(0,1,0,0.4);
		painter.draw(dt, points);

		painter.setColor(1,0,0);
		painter.draw(points.project(edgePointIds));
		painter.savePNG(name + ".png");

		painter.setColor(0,0,0, 0.1);
		painter.draw(realDT, points);

		painter.savePNG(name + "_overlay.png");
	};

	LOG << "Merging partial DTs into one triangulation" << std::endl;
	dSimplices mergedDT;
	for(uint i = 0; i < partialDTs.size(); ++i)
		mergedDT.insert(partialDTs[i].begin(), partialDTs[i].end());

	paintMerging(mergedDT, "05a_mergedDT_plain");

	mergeTriangulation(mergedDT, edgeSimplexIds, edgeDT, partioning, points);

	paintMerging(mergedDT, "05b_mergedDT_merged");

	LOG << "Eliminating duplicates" << std::endl;
	eliminateDuplicates(mergedDT,points);

	paintMerging(mergedDT, "05c_mergedDT_duplicates");

	LOG << "Updating neighbors" << std::endl;
	updateNeighbors(mergedDT, points);

	mergedDT.verify(points);

	paintMerging(mergedDT, "05d_mergedDT_finished");

	LOG << "Finished" << std::endl;
}


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

			points.insert(p);
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

	PainterBulkWriter paintWriter;

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
			LOG << "Error: wrong number of neighbors for simplex " << simplex << std::endl;

			if(IS_PROLIX){
				std::string png = "img/neighbors_" + std::to_string(simplex.id) + ".png";

				paintWriter[png] = Painter(bounds);
				paintWriter[png].draw(points);
				paintWriter[png].drawPartition(points);

				paintWriter[png].setColor(0,0,0,0.4);
				paintWriter[png].draw(simplices, points, true);

				paintWriter[png].setColor(1,0,0); //simplex in red
				paintWriter[png].draw(simplex, points, true);

				paintWriter[png].setColor(1,1,0,0.4);//neighbors in yellow
				paintWriter[png].drawNeighbors(simplex, simplices, points, true);
#ifndef NDEBUG
				if(!(simplex.equalVertices(saveSimplex) && simplex.equalNeighbors(saveSimplex))) {
					LOG << "Error: was before " << saveSimplex << std::endl;

					paintWriter[png].setColor(0,0,1); //old simplex in blue
					paintWriter[png].draw(saveSimplex, points, true);

					paintWriter[png].setColor(0,1,1,0.4);//old simplex neighbors in cyan
					paintWriter[png].drawNeighbors(saveSimplex, simplices, points, true);
				}
#endif
			}
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

dSimplices mergeTriangulation(std::vector<dSimplices> & partialDTs, const Ids & edgeSimplices, const dSimplices & edgeDT,
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

	LOG << "Merging partial DTs into one triangulation" << std::endl;
	dSimplices DT;
	for(uint i = 0; i < partialDTs.size(); ++i)
		DT.insert(partialDTs[i].begin(), partialDTs[i].end());

	auto edgePointIds = extractPoints(edgeSimplices, DT);

	auto paintMerging = [&] (const dSimplices & dt, const std::string & name) -> Painter {
		Painter painter(bounds);

		painter.draw(points);
		painter.drawPartition(points);

		painter.setColor(0,1,0,0.4);
		painter.draw(dt, points);

		painter.setColor(1,0,0);
		painter.draw(points.project(edgePointIds));
		painter.savePNG(name + ".png");

		painter.setColor(0,0,0, 0.1);
		painter.draw(realDT, points);

		painter.savePNG(name + "_overlay.png");

		return painter;
	};

	paintMerging(DT, "05a_merging_merged");

	auto deletedSimplices = DT.project(edgeSimplices);

	//delete all simplices belonging to the edge from DT
	LOG << "Striping triangulation from edge" << std::endl;
	for(const uint id : edgeSimplices){
		assert(DT.contains(id));

		//remove simplex from where used list
		for(uint p = 0; p < D+1; ++p){
			points[DT[id].vertices[p]].simplices.erase(id);
		}

		DT.erase(id);
	}

	auto painter = paintMerging(DT, "05b_merging_stripped");

	painter.setColor(0,0,1,0.4);
	painter.setLineDash({2,4});
	painter.draw(edgeDT, points);

	painter.savePNG("05b_merging_stripped+edge_overlay.png");

	painter.setColor(1,0,0,0.4);
	painter.setLineDash({2,4});
	painter.draw(deletedSimplices, points);
	painter.unsetLineDash();

	painter.savePNG("05b_merging_stripped+edge+deleted_overlay.png");

	//merge partial DTs and edge DT
	for(const auto & edgeSimplex : edgeDT){
		if(!partitionContains(edgeSimplex)){
			DT.insert(edgeSimplex);
		} else {
			auto it = std::find_if(deletedSimplices.begin(), deletedSimplices.end(),
					[&] (const dSimplex & s) { return s.equalVertices(edgeSimplex); } );
			if(it != deletedSimplices.end())
				DT.insert(edgeSimplex);
		}
	}

	paintMerging(DT, "05c_merging_edge");

	eliminateDuplicates(DT, points);
	updateNeighbors(DT, points);

	paintMerging(DT, "05d_merging_finished");

	return DT;

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

	auto mergedDT = mergeTriangulation(partialDTs, edgeSimplexIds, edgeDT, partioning, points);

	LOG << "Consistency check of triangulation" << std::endl;
	mergedDT.verify(points);

	LOG << "Cross check with real triangulation" << std::endl;
	auto report = mergedDT.verify(realDT);

	PainterBulkWriter writer;
	for(const auto & missingSimplex : report.missing){
		std::string img = "img/missing_" + std::to_string(missingSimplex.id) + ".png";

		writer[img] = basePainter;

		writer[img].setColor(1,0,0);
		writer[img].draw(missingSimplex, points, true);

		auto mySimplices = mergedDT.findSimplices(missingSimplex.vertices);

		writer[img].setColor(0,1,0);
		writer[img].draw(mySimplices, points, true);
	}

	LOG << "Finished" << std::endl;
}

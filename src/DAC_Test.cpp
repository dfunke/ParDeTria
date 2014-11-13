
//std library
#include <fstream>
#include <iostream>
#include <vector>
#include <functional>
#include <set>

//own
#include "Timings.h"
#include "Random.h"
#include "Geometry.h"
#include "Painter.h"
#include "Logger.h"
#include "CGAL_Interface.h"

dPoints genPoints(const uint n, const dBox & bounds, std::function<tCoordinate()> & dice){

	dPoints points(n);


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
		partition[part].push_back(p);
	}

	//add points at the extermities
	for(uint i = 0; i < partition.size(); ++i){
		auto stats = getPointStats(partition[i]);
		VLOG << "Partition " << i << ": " << stats.min << " - " << stats.mid << " - " << stats.max << std::endl;
		for(uint d = 0; d < D; ++d){

			//low point
			dPoint low = stats.mid;
			low.id = dPoint::cINF | i << D | ((1 << d) & ~(1));
			low.coords[d] = stats.min.coords[d] - (stats.max.coords[d] - stats.min.coords[d]);
			partition[i].push_back(low);
			points.push_back(low);

			//high point
			dPoint high = stats.mid;
			high.id = dPoint::cINF | i << D | ((1 << d) | 1);
			high.coords[d] = stats.max.coords[d] + (stats.max.coords[d] - stats.min.coords[d]);
			partition[i].push_back(high);
			points.push_back(high);

			VLOG << "Partition " << i << " - dimension " << d << ": " << low << " - " << high << std::endl;
		}
	}

	return partition;

}

dSimplices getEdge(const dSimplices & simplices){

	dSimplices edgeSimplices;
	std::set<uint> edgeIdx;

	auto norm = [&] (uint i) {
		return i - simplices[0].id; //is only called if simplices contains at least one element
	};

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

void mergeTriangulation(dSimplices & mergeDT, const dSimplices & otherDT, const dPoints & partition, const dPoints & points){

	auto edgeSimplices = getEdge(mergeDT);
	auto edgePoints = extractPoints(edgeSimplices, points);

	auto partitionContains = [&] (const dSimplex & simplex){
		for(uint i = 0; i < D+1; ++i){
			if(std::find(partition.begin(), partition.end(), simplex.vertices[i]) == partition.end())
				return false;
		}

		return true;
	};

	for(const auto & s : edgeSimplices){

		VLOG << "Edge Simplex " << s << std::endl;

		if(s.isFinite())
			continue; //this simplex has only finite points -> skip it //TODO alter edge detection

		//the simplex should have D finite points
		dPoints finitePoints;
		for(uint v = 0; v < D+1; ++v){
			if(dPoint::isFinite(s.vertices[v]))
				finitePoints.push_back(points[s.vertices[v]]);
		}
		VLOG << "Found " << finitePoints.size() << " finite points" << std::endl;

		VLOG << "Candidate mergers" << std::endl;
		INDENT
		dSimplices simplices;
		for(const auto & x : otherDT){
			if(x.contains(finitePoints) && !partitionContains(x)){
				simplices.push_back(x);
				VLOG << x << std::endl;
			}
		}
		DEDENT
	}

}

int main(int argc, char* argv[]) {

	Logger::getInstance().setLogLevel(Logger::Verbosity::LIVEVERBOSE);

	uint N = 1e2;

	dBox bounds;
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
	auto part = partition(points);
	DEDENT

	basePainter.draw(points);
	basePainter.savePNG("01_points.png");

	for(uint i = 0; i < part.size(); ++i){
		PLOG << "Partition " << i << ": ";
		for(auto & p : part[i])
			CONT << p << " ";
		CONT << std::endl;
	}

	std::vector<dSimplices> partialDT;
	Painter paintPartialDTs = basePainter;

	for(uint i = 0; i < part.size(); ++i){
		LOG << "Partition " << i << std::endl;
		INDENT
		partialDT.push_back(delaunayCgal(part[i]));
		LOG << "Triangulation " << i << " contains " << partialDT[i].size() << " tetrahedra" << std::endl << std::endl;
		DEDENT

		paintPartialDTs.draw(partialDT[i], points);

		Painter paintPartialDT = basePainter;
		paintPartialDT.draw(partialDT[i], points);
		paintPartialDT.savePNG("02_partialDT_" + std::to_string(i) + ".png");
	}
	paintPartialDTs.savePNG("02_partialDTs.png");

	LOG << "Extracting edges" << std::endl;
	INDENT

	dPoints edgePoints;
	Painter paintEdges = paintPartialDTs;

	for(uint i = 0; i < part.size(); ++i){

		auto edge = getEdge(partialDT[i]);
		auto ep = extractPoints(edge, points);
		//points are in different partitions, there can be no overlap
		edgePoints.insert(edgePoints.end(), ep.begin(), ep.end());

		VLOG << "Edge has " << edge.size() << " simplices with " << ep.size() << " points" << std::endl;

		paintEdges.setColor(1, 0, 0);
		paintEdges.draw(edge, points, false);
		paintEdges.setColor(0, 0, 0);

	}
	DEDENT

	paintEdges.setColor(1, 0, 0);
	paintEdges.draw(edgePoints, false);
	paintEdges.setColor(0, 0, 0);

	paintEdges.savePNG("03_edgeMarked.png");

	LOG << "Triangulating edges" << std::endl;
	INDENT
	auto edgeDT = delaunayCgal(edgePoints);
	LOG << "Edge triangulation contains " << edgeDT.size() << " tetrahedra" << std::endl << std::endl;
	DEDENT

	Painter paintEdgeDT = basePainter;

	paintEdgeDT.setColor(0, 1, 0);
	paintEdgeDT.draw(edgeDT, points);
	paintEdgeDT.setColor(0, 0, 0);

	paintEdgeDT.savePNG("04_edgeDT.png");

	for(uint i = 0; i < partialDT.size(); ++i){
		LOG << "Merge partition " << i << std::endl;
		INDENT
		mergeTriangulation(partialDT[i], edgeDT, part[i], points);
		DEDENT
	}

	LOG << "Real triangulation" << std::endl;
	INDENT
	auto realDT = delaunayCgal(points,true);
	LOG << "Real triangulation contains " << realDT.size() << " tetrahedra" << std::endl << std::endl;
	DEDENT

	Painter paintRealDT = basePainter;

	paintRealDT.setColor(0, 0, 1);
	paintRealDT.draw(realDT, points);
	paintRealDT.setColor(0, 0, 0);

	paintRealDT.savePNG("05_realDT.png");

}

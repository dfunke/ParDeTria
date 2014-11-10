
//std library
#include <fstream>
#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <functional>
#include <limits>
#include <sstream>
#include <stdexcept>

//CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<uint, K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> CT; //CGAL triangulation

//own
#include "Timings.h"
#include "Random.h"
#include "Geometry.h"
#include "Painter.h"

//#############################################################################
uint indentLevel = 0;
std::string indent(){
	std::stringstream ss;
	for(uint i = 0; i < indentLevel; ++i)
		ss << "\t";
	return ss.str();
}

#define COUT std::cout << indent()

//#############################################################################

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

Partition partition(const dPoints & points){

	// do mid-point based partitioning for now
	dPoint midpoint;
	for(uint dim = 0; dim < D; ++dim){
		auto minmax = std::minmax_element(points.begin(), points.end(),
				[dim] (const dPoint & a, const dPoint & b) {
					return a.coords[dim] < b.coords[dim];
				});
		midpoint.coords[dim] = ((*minmax.second).coords[dim] - (*minmax.first).coords[dim]) / 2;
	}

	COUT << "Midpoint is " <<  midpoint << std::endl;

	Partition partition(pow(2,D));

	for(auto & p : points){
		uint part = 0;
		for(uint dim = 0; dim < D; ++dim){
			part |= (p.coords[dim] > midpoint.coords[dim]) << dim;
		}

		//COUT << "Adding " << p << " to " << part << std::endl;
		partition[part].push_back(p);
	}

	return partition;

}

uint tetrahedronID = 0;
dSimplices delaunayCgal(const dPoints & points){

	//copy points into CGAL structure
	std::vector<std::pair<CT::Point, uint> > cPoints;
	cPoints.reserve(points.size());
	for(const auto & p : points){
		CT::Point cp(p.coords[0], p.coords[1]);
		cPoints.push_back(std::make_pair(cp, p.id));
	}
	
	CT t;

	//auto start = Clock::now();
	t.insert(cPoints.begin(), cPoints.end());
	//auto end = Clock::now();

	COUT << "Triangulation is " << (t.is_valid() ? "" : "NOT ") << "valid" << std::endl;
	COUT << "finite faces/vertices "
			  << t.number_of_faces() << "/" << t.number_of_vertices()  << std::endl;
	
	COUT << "Collecting simplices" << std::endl;
	++indentLevel;

	dSimplices tria;
	for(auto it = t.all_faces_begin(); it != t.all_faces_end(); ++it){
		dSimplex a;
		a.id = tetrahedronID++;

		for(uint i = 0; i < D+1; ++i){
			if(!t.is_infinite(it->vertex(i))){
				a.vertices[i] = it->vertex(i)->info();
			} else {
				a.vertices[i] = INF;
			}

		}

		COUT << a << std::endl;
		tria.push_back(a);
	}
	--indentLevel;

	auto cmp = [&] (const CT::Face_handle & f) -> dSimplex {
		dSimplex a;
		a.id = INF;
		for(uint i = 0; i < D+1; ++i){
			if(!t.is_infinite(f->vertex(i))){
				a.vertices[i] = f->vertex(i)->info();
			} else {
				a.vertices[i] = INF;
			}
		}

		return a;
	};

	COUT << "Collecting neighbors" << std::endl;

	++indentLevel;
	for(auto it = t.all_faces_begin(); it != t.all_faces_end(); ++it){
		auto tet = std::find(tria.begin(), tria.end(), cmp(it));
		
		for(uint i = 0; i < D+1; ++i){
			auto n = it->neighbor(i);
			auto nn = std::find(tria.begin(), tria.end(), cmp(n));
			tet->neighbors[i] = nn->id;
		}

		COUT << *tet << std::endl;
	}
	--indentLevel;

	return tria;
}

dSimplices getEdge(const dSimplices & simplices){

	dSimplices edgeSimplices;
	std::set<uint> edgeIdx;

	auto norm = [&] (uint i) {
		return i - simplices[0].id; //is only called if simplices contains at least one element
	};

	for(const auto & s : simplices){
		if(!s.isFinite()){
			//we have an infinte vertex, at least one of its neighbors must be finite
			for(uint i = 0; i < D+1; ++i){
				if(simplices[norm(s.neighbors[i])].isFinite() && edgeIdx.insert(s.neighbors[i]).second)
					//the neighbor is finite and has not been added to the edge simplices yet
					edgeSimplices.push_back(simplices[norm(s.neighbors[i])]);
			}
		}
	}

	return edgeSimplices;

}

dPoints extractPoints(const dSimplices & simplices, const dPoints & inPoints){

	dPoints outPoints;
	std::set<uint> idx;

	for(const auto & s : simplices){
		for(uint i = 0; i < D+1; ++i){
			if(s.vertices[i] != INF && idx.insert(s.vertices[i]).second)
				outPoints.push_back(inPoints[s.vertices[i]]);
		}
	}

	return outPoints;

}

int main(int argc, char* argv[]) {

	uint N = 1e2;

	dBox bounds;
	for(uint i = 0; i < D; ++i){
		bounds.coords[i] = 0;
		bounds.dim[i] = 100;
	}

	dBox img;
	for(uint i = 0; i < D; ++i){
		img.coords[i] = 0;
		img.dim[i] = 2000;
	}

	Painter painter(bounds, img);

	std::uniform_real_distribution<tCoordinate> distribution(0,1);
	std::function<tCoordinate()> dice = std::bind (distribution, generator);

	auto points = genPoints(N, bounds, dice);

	COUT << "Partioning" << std::endl;
	++indentLevel;
	auto part = partition(points);
	--indentLevel;

	painter.draw(points);
	painter.savePNG("01_points.png");

	for(uint i = 0; i < part.size(); ++i){
		COUT << "Partition " << i << ": ";
		for(auto & p : part[i])
			std::cout << p << " ";
		std::cout << std::endl;
	}

	std::vector<dSimplices> partialDT;
	for(uint i = 0; i < part.size(); ++i){
		COUT << "Partition " << i << std::endl;
		++indentLevel;
		partialDT.push_back(delaunayCgal(part[i]));
		COUT << "Triangulation " << i << " contains " << partialDT[i].size() << " tetrahedra" << std::endl << std::endl;
		--indentLevel;

		painter.draw(partialDT[i], points);
	}
	painter.savePNG("02_partialDTs.png");

	dPoints edgePoints;
	for(uint i = 0; i < part.size(); ++i){

		auto edge = getEdge(partialDT[i]);
		auto ep = extractPoints(edge, points);
		//points are in different partitions, there can be no overlap
		edgePoints.insert(edgePoints.end(), ep.begin(), ep.end());

		painter.setColor(1, 0, 0);
		painter.draw(edge, points);
		painter.setColor(0, 0, 0);

	}
	painter.setColor(1, 0, 0);
	painter.draw(edgePoints, false);
	painter.setColor(0, 0, 0);
	painter.savePNG("03_edgeMarked.png");

	COUT << "Triangulating edges" << std::endl;
	++indentLevel;
	auto edgeDT = delaunayCgal(edgePoints);
	COUT << "Edge triangulation contains " << edgeDT.size() << " tetrahedra" << std::endl << std::endl;
	--indentLevel;

	painter.setColor(0, 1, 0);
	painter.draw(edgeDT, points);
	painter.setColor(0, 0, 0);

	painter.savePNG("04_edgeDT.png");

	COUT << "Real triangulation" << std::endl;
	++indentLevel;
	auto realDT = delaunayCgal(points);
	COUT << "Real triangulation contains " << realDT.size() << " tetrahedra" << std::endl << std::endl;
	--indentLevel;

	painter.setColor(0, 0, 1);
	painter.draw(realDT, points);
	painter.setColor(0, 0, 0);

	painter.savePNG("05_realDT.png");

}

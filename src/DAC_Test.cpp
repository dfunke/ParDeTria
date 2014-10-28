
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

//#############################################################################

//dimensionality of our problem
const uint D = 2;
const uint INF = std::numeric_limits<uint>::max();

typedef float tCoordinate;

struct Point {
	uint id;
	tCoordinate coords[D];

	bool operator==(const Point & a) const {
		if(!(this->id == INF ^ a.id == INF))
			return this->id == a.id;
		else{
			bool eq = true;
			for(uint i = 0; i < D; ++i){
				if(this->coords[i] != a.coords[i])
					eq = false;
			}
			return eq;
		}
	}
};


struct dSimplex {
	uint id;
	uint vertices[D+1];
	uint neighbors[D+1];

	bool operator==(const dSimplex & a) const {
		if(!(this->id == INF ^ a.id == INF))
			return this->id == a.id;
		else {
			bool eq = true;
			for(uint i = 0; i < D+1; ++i){
				bool found = false;
				for(uint j = 0; j < D+1; ++j){
					if(vertices[i] == vertices[j])
						found = true;
				}
				if(!found)
					eq = false;
			}
			return eq;
		}
	}
};

std::ostream & operator<<(std::ostream & o, const Point & p){
	o << "[" << p.coords[0];
	for(uint i = 1; i < D; ++i)
		o << ", " << p.coords[i];
	o << "]";
	return o;
}

typedef std::vector<Point> Points;
typedef std::vector<Points> Partition;
typedef std::vector<dSimplex> Triangulation;

struct Box {
	tCoordinate coords[D];
	tCoordinate dim[D];
};

//#############################################################################

Points genPoints(const uint n, const Box & bounds, std::function<tCoordinate()> & dice){

	Points points(n);


	for(uint i = 0; i < n; ++i){
		points[i].id = i;
		for(uint d = 0; d < D; ++d){
			points[i].coords[d] = bounds.coords[d] + bounds.dim[d] * dice();
		}
	}

	//TODO checks for colliding points, straight lines, etc.

	return points;

}

Partition partition(const Points & points){

	// do mid-point based partitioning for now
	Point midpoint;
	for(uint dim = 0; dim < D; ++dim){
		auto minmax = std::minmax_element(points.begin(), points.end(),
				[dim] (const Point & a, const Point & b) {
					return a.coords[dim] < b.coords[dim];
				});
		midpoint.coords[dim] = ((*minmax.second).coords[dim] - (*minmax.first).coords[dim]) / 2;
	}

	std::cout << "Midpoint is " <<  midpoint << std::endl;

	Partition partition(8);

	for(auto & p : points){
		uint part = 0;
		for(uint dim = 0; dim < D; ++dim){
			part |= (p.coords[dim] > midpoint.coords[dim]) << dim;
		}

		//std::cout << "Adding " << p << " to " << part << std::endl;
		partition[part].push_back(p);
	}

	return partition;

}

uint tetrahedronID = 0;
Triangulation delaunayCgal(const Points & points){

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

	std::cout << "Triangulation is " << (t.is_valid() ? "" : "NOT ") << "valid" << std::endl;
	std::cout << "finite faces/vertices " 
						<< t.number_of_faces() << "/" << t.number_of_vertices()  << std::endl;
	std::cout << std::endl;
	
	Triangulation tria;
	std::map<uint, std::set<uint> > pointToSimplex;
	for(auto it = t.all_faces_begin(); it != t.all_faces_end(); ++it){
		dSimplex a;
		a.id = tetrahedronID++;

		std::cout << "Simplex " << a.id << " vertices: [";

		for(uint i = 0; i < D+1; ++i){
			if(!t.is_infinite(it->vertex(i))){
				a.vertices[i] = it->vertex(i)->info();
				pointToSimplex[it->vertex(i)->info()].insert(a.id);
			} else {
				a.vertices[i] = INF;
				//pointToSimplex[INF].insert(a.id);
			}

			std::cout << a.vertices[i] << " ";
		}

		std::cout << "]" << std::endl;
		tria.push_back(a);
	}

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

	for(auto it = t.finite_faces_begin(); it != t.finite_faces_end(); ++it){
		auto tet = std::find(tria.begin(), tria.end(), cmp(it));
		std::cout << "Simplex " << tet->id << " neighbors: [";
		
		for(uint i = 0; i < D+1; ++i){
			auto n = it->neighbor(i);
			auto nn = std::find(tria.begin(), tria.end(), cmp(n));
			tet->neighbors[i] = nn->id;
			std::cout << nn->id << " ";
		}
		std::cout << "]" << std::endl;
	}

	return tria;
}

int main(int argc, char* argv[]) {

	uint N = 1e1;

	Box bounds;
	for(uint i = 0; i < D; ++i){
		bounds.coords[i] = 0;
		bounds.dim[i] = 100;
	}

	std::uniform_real_distribution<tCoordinate> distribution(0,1);
	std::function<tCoordinate()> dice = std::bind (distribution, generator);

	auto points = genPoints(N, bounds, dice);
	auto part = partition(points);

	for(uint i = 0; i < part.size(); ++i){
		std::cout << "Partition " << i << ": ";
		for(auto & p : part[i])
			std::cout << p << " ";
		std::cout << std::endl;
		std::cout << std::endl;
	}

	std::vector<Triangulation> partialDT;
	for(uint i = 0; i < part.size(); ++i){
		partialDT.push_back(delaunayCgal(part[i]));
		std::cout << "Triangulation " << i << " contains " << partialDT[i].size() << " tetrahedra" << std::endl;
	}
}

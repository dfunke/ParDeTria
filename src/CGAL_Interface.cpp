#include "CGAL_Interface.h"
#include "Logger.h"

#include <algorithm>

//CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<uint, K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> CT; //CGAL triangulation

uint tetrahedronID = 0;
dSimplices delaunayCgal(dPoints & points, const Ids * ids, bool filterInfinite){

	//copy points into CGAL structure
	std::vector<std::pair<CT::Point, uint> > cPoints;
	cPoints.reserve(points.size());
	for(const auto & p : points){
		if((ids != nullptr && std::find(ids->begin(), ids->end(), p.id) == ids->end())
				|| (filterInfinite && !p.isFinite()))
			continue;

		CT::Point cp(p.coords[0], p.coords[1]);
		cPoints.push_back(std::make_pair(cp, p.id));
	}

	CT t;

	//auto start = Clock::now();
	t.insert(cPoints.begin(), cPoints.end());
	//auto end = Clock::now();

	LOG << "CGAL triangulation is " << (t.is_valid() ? "" : "NOT ") << "valid" << std::endl;
	LOG << "finite faces/vertices "
			  << t.number_of_faces() << "/" << t.number_of_vertices()  << std::endl;

	LOG << "Collecting simplices" << std::endl;
	INDENT

	dSimplices tria;
	for(auto it = t.finite_faces_begin(); it != t.finite_faces_end(); ++it){
		dSimplex a;
		a.id = tetrahedronID++;

		for(uint i = 0; i < D+1; ++i){
			a.vertices[i] = it->vertex(i)->info();
			points[a.vertices[i]].simplices.insert(a.id);
		}

		PLOG << a << std::endl;
		tria.insert(a);
	}
	DEDENT

	auto cmp = [&] (const CT::Face_handle & f) -> dSimplex {
		dSimplex a;
		a.id = dSimplex::cINF;
		for(uint i = 0; i < D+1; ++i){
			//if(!t.is_infinite(f->vertex(i))){
				a.vertices[i] = f->vertex(i)->info();
			//} else {
			//	a.vertices[i] = dPoint::cINF;
			//}
		}

		return a;
	};

	LOG << "Collecting neighbors" << std::endl;

	INDENT
	for(auto it = t.finite_faces_begin(); it != t.finite_faces_end(); ++it){
		auto tet = std::find(tria.begin(), tria.end(), cmp(it));

		for(uint i = 0; i < D+1; ++i){
			auto n = it->neighbor(i);
			auto nn = std::find(tria.begin(), tria.end(), cmp(n));

			tet->neighbors.insert((nn != tria.end()) ? nn->id : dSimplex::cINF);

				/*LOG << "Did not find neighbor " << i << " of " << *tet << " - handle "
				<< n->vertex(0)->info() << " " << t.is_infinite(n->vertex(0)) << "; "
				<< n->vertex(1)->info() << " " << t.is_infinite(n->vertex(1)) << "; "
				<< n->vertex(2)->info() << " " << t.is_infinite(n->vertex(2)) << "; "
				<< std::endl;*/
		}

		PLOG << *tet << std::endl;
	}
	DEDENT

	return tria;
}

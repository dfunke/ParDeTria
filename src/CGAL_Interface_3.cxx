#include "CGAL_Interface.h"
#include "utils/Logger.h"

#include <algorithm>

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

template <>
dSimplices<3> delaunayCgal<3>(dPoints<3> &points, const Ids *ids,
                              bool filterInfinite) {

  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Triangulation_vertex_base_with_info_3<uint, K> Vb;
  typedef CGAL::Triangulation_data_structure_3<Vb> Tds;
  typedef CGAL::Delaunay_triangulation_3<K, Tds> CT;

  const uint D = 3;

  // copy points into CGAL structure
  std::vector<std::pair<CT::Point, uint>> cPoints;
  cPoints.reserve(points.size());
  for (const auto &p : points) {
    if ((ids != nullptr &&
         std::find(ids->begin(), ids->end(), p.id) == ids->end()) ||
        (filterInfinite && !p.isFinite()))
      continue;

    CT::Point cp(p.coords[0], p.coords[1], p.coords[2]);
    cPoints.push_back(std::make_pair(cp, p.id));
  }

  CT t;

  // auto start = Clock::now();
  t.insert(cPoints.begin(), cPoints.end());
  // auto end = Clock::now();

  VLOG << "CGAL triangulation is " << (t.is_valid() ? "" : "NOT ") << "valid"
       << std::endl;
  VLOG << "finite cells/vertices " << t.number_of_finite_cells() << "/"
       << t.number_of_vertices() << std::endl;

  PLOG << "Collecting simplices" << std::endl;
  INDENT

  dSimplices<3> tria;
  for (auto it = t.finite_cells_begin(); it != t.finite_cells_end(); ++it) {
    dSimplex<3> a;
    a.id = tetrahedronID++;

    for (uint i = 0; i < D + 1; ++i) {
      a.vertices[i] = it->vertex(i)->info();
      points[a.vertices[i]].simplices.insert(a.id);
    }

    PLOG << a << std::endl;
    tria.insert(a);
  }
  DEDENT

  auto cmp = [&](const CT::Cell_handle &f) -> dSimplex<3> {
    dSimplex<3> a;
    a.id = dSimplex<3>::cINF;
    for (uint i = 0; i < D + 1; ++i) {
      // if(!t.is_infinite(f->vertex(i))){
      a.vertices[i] = f->vertex(i)->info();
      //} else {
      //	a.vertices[i] = dPoint::cINF;
      //}
    }

    return a;
  };

  PLOG << "Collecting neighbors" << std::endl;

  INDENT
  for (auto it = t.finite_cells_begin(); it != t.finite_cells_end(); ++it) {
    auto tet = std::find(tria.begin(), tria.end(), cmp(it));

    for (uint i = 0; i < D + 1; ++i) {
      auto n = it->neighbor(i);
      auto nn = std::find(tria.begin(), tria.end(), cmp(n));

      tet->neighbors.insert((nn != tria.end()) ? nn->id : dSimplex<3>::cINF);

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

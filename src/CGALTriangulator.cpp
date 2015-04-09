#include "CGALTriangulator.h"

#include <atomic>
#include <type_traits>

#include "utils/ASSERT.h"

// define a static counter for the tetrahedronID
std::atomic<uint> tetrahedronID(0);

// CGAL
#define CGAL_LINKED_WITH_TBB
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <CGAL/Unique_hash_map.h>

template <uint D, typename Precision, class Tria, bool Parallel = true> class CGALHelper;

template <typename Precision, class Tria> class CGALHelper<2, Precision, Tria> {

public:
  CGALHelper(__attribute__((unused)) const dBox<2, Precision> &bounds) {}

  typename Tria::Finite_faces_iterator begin(Tria &t) {
    return t.finite_faces_begin();
  }

  typename Tria::Finite_faces_iterator end(Tria &t) {
    return t.finite_faces_end();
  }

  uint size(Tria &t) { return t.number_of_faces(); }

  typename Tria::Point make_point(const dPoint<2, Precision> &p) {
    return typename Tria::Point(p.coords[0], p.coords[1]);
  }

  Tria make_tria() {
    Tria t;
    return t;
  }

  using Handle = typename Tria::Face_handle;
};

template <typename Precision, class Tria> class CGALHelper<3, Precision, Tria, false> {

public:
  CGALHelper(__attribute__((unused)) const dBox<3, Precision> &bounds) { }

  typename Tria::Finite_cells_iterator begin(Tria &t) {
    return t.finite_cells_begin();
  }

  typename Tria::Finite_cells_iterator end(Tria &t) {
    return t.finite_cells_end();
  }

  uint size(Tria &t) { return t.number_of_finite_cells(); }

  typename Tria::Point make_point(const dPoint<3, Precision> &p) {
    return typename Tria::Point(p.coords[0], p.coords[1], p.coords[2]);
  }

  Tria make_tria() {
    Tria t;
    return t;
  }

  using Handle = typename Tria::Cell_handle;
};

template <typename Precision, class Tria> class CGALHelper<3, Precision, Tria, true> {

public:
  CGALHelper(const dBox<3, Precision> &bounds)
      : lockingDS(CGAL::Bbox_3(bounds.low[0], bounds.low[1], bounds.low[2],
                               bounds.high[0], bounds.high[1], bounds.high[2]),
                  50) {}

  typename Tria::Finite_cells_iterator begin(Tria &t) {
    return t.finite_cells_begin();
  }

  typename Tria::Finite_cells_iterator end(Tria &t) {
    return t.finite_cells_end();
  }

  uint size(Tria &t) { return t.number_of_finite_cells(); }

  typename Tria::Point make_point(const dPoint<3, Precision> &p) {
    return typename Tria::Point(p.coords[0], p.coords[1], p.coords[2]);
  }

  Tria make_tria() {
    Tria t(&lockingDS);
    return t;
  }

  using Handle = typename Tria::Cell_handle;

private:
  typename Tria::Lock_data_structure lockingDS;
};

template <uint D, typename Precision, class Tria, bool Parallel>
dSimplices<D, Precision>
_delaunayCgal(const Ids &ids, dPoints<D, Precision> &points,
              const dBox<D, Precision> &bounds
              /*, bool filterInfinite */) {

  CGALHelper<D, Precision, Tria, Parallel> helper(bounds);

  // copy points into CGAL structure
  std::vector<std::pair<typename Tria::Point, uint>> cPoints;
  cPoints.reserve(ids.size());
  for (const auto &id : ids) {

    const auto &p = points[id];

    /* if (filterInfinite && !p.isFinite())
      continue; */

    auto cp = helper.make_point(p);
    cPoints.push_back(std::make_pair(cp, p.id));
  }

  Tria t = helper.make_tria();
  t.insert(cPoints.begin(), cPoints.end());

  ASSERT(t.is_valid());

  // VLOG("CGAL triangulation is " << (t.is_valid() ? "" : "NOT ") << "valid"
  //     << std::endl);
  // VLOG("finite cells/vertices " << t.number_of_finite_cells() << "/"
  //     << t.number_of_vertices() << std::endl);

  PLOG("Collecting simplices" << std::endl);
  INDENT

  dSimplices<D, Precision> tria;
  tria.reserve(helper.size(t));

  CGAL::Unique_hash_map<typename CGALHelper<D, Precision, Tria, Parallel>::Handle, uint>
      simplexLookup(0, helper.size(t));

  for (auto it = helper.begin(t); it != helper.end(t); ++it) {
    dSimplex<D, Precision> a;
    a.id = tetrahedronID++;

    for (uint i = 0; i < D + 1; ++i) {
      dPoint<D, Precision> &point = points[it->vertex(i)->info()];
      ASSERT(point.id == it->vertex(i)->info());

      a.vertices[i] = point.id;

      // update where-used data structure
      tria.whereUsed[point.id].emplace_back(a.id);
    }
    // sort vertices by ascending point id
    std::sort(a.vertices.begin(), a.vertices.end());
    a.fingerprint();

    PLOG(a << std::endl);
    tria.insert(a);
    simplexLookup[it] = a.id;
  }
  DEDENT

  PLOG("Collecting neighbors" << std::endl);

  INDENT
  for (auto it = helper.begin(t); it != helper.end(t); ++it) {
    auto &tet = tria[simplexLookup[it]];
    tet.neighbors.reserve(D + 1);

    for (uint i = 0; i < D + 1; ++i) {
      if (simplexLookup.is_defined(it->neighbor(i))) {
        const auto &n = tria[simplexLookup[it->neighbor(i)]];
        if (n.id != tet.id) {
          tet.neighbors.insert(n.id);
        }
      } else {
        tet.neighbors.insert(dSimplex<D, Precision>::cINF);
      }
    }

    PLOG(tet << std::endl);
  }
  DEDENT

  return tria;
}

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

// partial specializations

template <typename Precision, bool Parallel> class CGALTriangulator<2, Precision, Parallel>
        : public Triangulator<2, Precision> {

public:
  CGALTriangulator(const dBox<2, Precision> &_bounds, dPoints<2, Precision> &_points)
  : Triangulator<2, Precision>(_bounds, _points) {};

protected:
  dSimplices<2, Precision> _triangulate(const Ids &ids,
                                              const dBox<2, Precision> &bounds,
                                              __attribute__((unused)) const std::string provenance
                                              /*, bool filterInfinite */) {

    typedef CGAL::Triangulation_vertex_base_with_info_2<uint, K> Vb;
    typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
    typedef CGAL::Delaunay_triangulation_2<K, Tds> CT;

    return _delaunayCgal<2, Precision, CT, Parallel>(ids, this->points, bounds /*, filterInfinite */);
  }
};

template <typename Precision> class CGALTriangulator<3, Precision, true>
        : public Triangulator<3, Precision> {

public:
  CGALTriangulator(const dBox<3, Precision> &_bounds, dPoints<3, Precision> &_points)
          : Triangulator<3, Precision>(_bounds, _points) {};

protected:

  dSimplices<3, Precision> _triangulate(const Ids &ids,
                                              const dBox<3, Precision> &bounds,
                                              __attribute__((unused)) const std::string provenance
                                              /*, bool filterInfinite */) {

    typedef CGAL::Triangulation_vertex_base_with_info_3<uint, K> Vb;
    typedef CGAL::Triangulation_data_structure_3<
        Vb, CGAL::Triangulation_cell_base_3<K>, CGAL::Parallel_tag> Tds;
    typedef CGAL::Delaunay_triangulation_3<K, Tds> CT;

    return _delaunayCgal<3, Precision, CT, true>(ids, this->points, bounds /*, filterInfinite */);
  }
};

template <typename Precision> class CGALTriangulator<3, Precision, false>
        : public Triangulator<3, Precision> {

public:
  CGALTriangulator(const dBox<3, Precision> &_bounds, dPoints<3, Precision> &_points)
          : Triangulator<3, Precision>(_bounds, _points) {};

protected:

  dSimplices<3, Precision> _triangulate(const Ids &ids,
                                        const dBox<3, Precision> &bounds,
                                        __attribute__((unused)) const std::string provenance
          /*, bool filterInfinite */) {

    typedef CGAL::Triangulation_vertex_base_with_info_3<uint, K> Vb;
    typedef CGAL::Triangulation_data_structure_3<
            Vb, CGAL::Triangulation_cell_base_3<K>, CGAL::Sequential_tag> Tds;
    typedef CGAL::Delaunay_triangulation_3<K, Tds> CT;

    return _delaunayCgal<3, Precision, CT, false>(ids, this->points, bounds /*, filterInfinite */);
  }
};

// specializations
template class CGALTriangulator<2, float>;

template class CGALTriangulator<3, float, true>;
template class CGALTriangulator<3, float, false>;

template class CGALTriangulator<2, double>;

template class CGALTriangulator<3, double, true>;
template class CGALTriangulator<3, double, false>;

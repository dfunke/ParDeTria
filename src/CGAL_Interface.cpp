#include "CGAL_Interface.h"

#include <atomic>

#include "utils/ASSERT.h"

// define a static counter for the tetrahedronID
std::atomic<uint> tetrahedronID(0);

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <CGAL/Unique_hash_map.h>

template <uint D, typename Precision> class CGALHelper {

public:
  template <class IT, class Tria> static IT begin(Tria &t);

  template <class IT, class Tria> static IT end(Tria &t);

  template <class Tria> static uint size(Tria &t);

  template <class Tria>
  static typename Tria::Point make_point(const dPoint<D, Precision> &p);
};

template <typename Precision> class CGALHelper<2, Precision> {

public:
  template <class Tria>
  static typename Tria::Finite_faces_iterator begin(Tria &t) {
    return t.finite_faces_begin();
  }

  template <class Tria>
  static typename Tria::Finite_faces_iterator end(Tria &t) {
    return t.finite_faces_end();
  }

  template <class Tria> static uint size(Tria &t) {
    return t.number_of_faces();
  }

  template <class Tria>
  static typename Tria::Point make_point(const dPoint<2, Precision> &p) {
    return typename Tria::Point(p.coords[0], p.coords[1]);
  }

  template <class Tria> using Handle = typename Tria::Face_handle;
};

template <typename Precision> class CGALHelper<3, Precision> {

public:
  template <class Tria>
  static typename Tria::Finite_cells_iterator begin(Tria &t) {
    return t.finite_cells_begin();
  }

  template <class Tria>
  static typename Tria::Finite_cells_iterator end(Tria &t) {
    return t.finite_cells_end();
  }

  template <class Tria> static uint size(Tria &t) {
    return t.number_of_finite_cells();
  }

  template <class Tria>
  static typename Tria::Point make_point(const dPoint<3, Precision> &p) {
    return typename Tria::Point(p.coords[0], p.coords[1], p.coords[2]);
  }

  template <class Tria> using Handle = typename Tria::Cell_handle;
};

template <uint D, typename Precision, class Tria>
dSimplices<D, Precision> _delaunayCgal(dPoints<D, Precision> &points,
                                       const Ids *ids, bool filterInfinite) {

  // copy points into CGAL structure
  std::vector<std::pair<typename Tria::Point, uint>> cPoints;
  cPoints.reserve(points.size());
  for (const auto &p : points) {
    if ((ids != nullptr &&
         std::find(ids->begin(), ids->end(), p.id) == ids->end()) ||
        (filterInfinite && !p.isFinite()))
      continue;

    auto cp = CGALHelper<D, Precision>::template make_point<Tria>(p);
    cPoints.push_back(std::make_pair(cp, p.id));
  }

  Tria t;

  // auto start = Clock::now();
  t.insert(cPoints.begin(), cPoints.end());
  // auto end = Clock::now();

  ASSERT(t.is_valid());

  // VLOG("CGAL triangulation is " << (t.is_valid() ? "" : "NOT ") << "valid"
  //     << std::endl);
  // VLOG("finite cells/vertices " << t.number_of_finite_cells() << "/"
  //     << t.number_of_vertices() << std::endl);

  PLOG("Collecting simplices" << std::endl);
  INDENT

  dSimplices<D, Precision> tria;
  tria.reserve(CGALHelper<D, Precision>::size(t));

  CGAL::Unique_hash_map<
      typename CGALHelper<D, Precision>::template Handle<Tria>, uint>
      simplexLookup(0, CGALHelper<D, Precision>::size(t));

  for (auto it = CGALHelper<D, Precision>::begin(t);
       it != CGALHelper<D, Precision>::end(t); ++it) {
    dSimplex<D, Precision> a;
    a.id = tetrahedronID++;

    for (uint i = 0; i < D + 1; ++i) {
      dPoint<D, Precision> &point = points[it->vertex(i)->info()];

      a.vertices[i] = point.id;

      std::lock_guard<SpinMutex> lock(point.mtx);
      point.simplices.insert(a.id);
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
  for (auto it = CGALHelper<D, Precision>::begin(t);
       it != CGALHelper<D, Precision>::end(t); ++it) {
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

    PLOG(*tet << std::endl);
  }
  DEDENT

  return tria;
}

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

// partial specializations

template <typename Precision> class CGALInterface<2, Precision> {
public:
  static dSimplices<2, Precision> triangulate(dPoints<2, Precision> &points,
                                              const Ids *ids,
                                              bool filterInfinite) {

    typedef CGAL::Triangulation_vertex_base_with_info_2<uint, K> Vb;
    typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
    typedef CGAL::Delaunay_triangulation_2<K, Tds> CT;

    return _delaunayCgal<2, Precision, CT>(points, ids, filterInfinite);
  }
};

template <typename Precision> class CGALInterface<3, Precision> {
public:
  static dSimplices<3, Precision> triangulate(dPoints<3, Precision> &points,
                                              const Ids *ids,
                                              bool filterInfinite) {

    typedef CGAL::Triangulation_vertex_base_with_info_3<uint, K> Vb;
    typedef CGAL::Triangulation_data_structure_3<Vb> Tds;
    typedef CGAL::Delaunay_triangulation_3<K, Tds> CT;

    return _delaunayCgal<3, Precision, CT>(points, ids, filterInfinite);
  }
};

// specializations
template class CGALInterface<2, float>;
template class CGALInterface<3, float>;

template class CGALInterface<2, double>;
template class CGALInterface<3, double>;

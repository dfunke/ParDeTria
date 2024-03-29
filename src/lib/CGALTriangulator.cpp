#include "CGALTriangulator.h"

#ifndef NDEBUG

#include <csignal>

#endif

#include <atomic>
#include <type_traits>

#include "utils/ASSERT.h"
#include "utils/VTuneAdapter.h"

// CGAL
#define CGAL_LINKED_WITH_TBB
#define CGAL_CONCURRENT_TRIANGULATION_3_ADD_TEMPORARY_POINTS_ON_FAR_SPHERE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "mods/Indexed_Triangulation_data_structure_2.h"
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include "mods/Indexed_Triangulation_data_structure_3.h"
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

// boost
#include <boost/iterator/transform_iterator.hpp>

#include <CGAL/Unique_hash_map.h>

template<typename GT,
        class Ih = _detail::IndexHandler,
        typename Fb_ = CGAL::Triangulation_face_base_2<GT>>
class Triangulation_dSimplexAdapter_2
        : public Fb_ {
public:
    typedef typename Fb_::Vertex_handle Vertex_handle;
    typedef typename Fb_::Face_handle Face_handle;

    template<typename TDS2>
    struct Rebind_TDS {
        typedef typename Fb_::template Rebind_TDS<TDS2>::Other Fb2;
        typedef Triangulation_dSimplexAdapter_2<GT, Ih, Fb2> Other;
    };

    Triangulation_dSimplexAdapter_2(const Ih &idHandler)
            : Fb_(), m_idHandler(idHandler) {

        m_id = m_idHandler.getId();
    }

    Triangulation_dSimplexAdapter_2(const Triangulation_dSimplexAdapter_2 &o) // do not copy the ID, just the faces and neighbors
            : Fb_(o), m_idHandler(o.m_idHandler) {

        m_id = m_idHandler.getId();
    }

    Triangulation_dSimplexAdapter_2(Vertex_handle v0,
                                    Vertex_handle v1,
                                    Vertex_handle v2,
                                    const Ih &idHandler)
            : Fb_(v0, v1, v2), m_idHandler(idHandler) {

        m_id = m_idHandler.getId();
    }

    Triangulation_dSimplexAdapter_2(Vertex_handle v0,
                                    Vertex_handle v1,
                                    Vertex_handle v2,
                                    Face_handle n0,
                                    Face_handle n1,
                                    Face_handle n2,
                                    const Ih &idHandler)
            : Fb_(v0, v1, v2, n0, n1, n2), m_idHandler(idHandler) {

        m_id = m_idHandler.getId();
    }

    ~Triangulation_dSimplexAdapter_2() {
        m_idHandler.releaseId(m_id);
    }

public:
    tIdType m_id;
    const Ih &m_idHandler;
};

template<typename GT,
        class Ih = _detail::IndexHandler,
        typename Cb = CGAL::Triangulation_cell_base_3<GT>>
class Triangulation_dSimplexAdapter_3
        : public Cb {
public:
    typedef typename Cb::Vertex_handle Vertex_handle;
    typedef typename Cb::Cell_handle Cell_handle;

    template<typename TDS2>
    struct Rebind_TDS {
        typedef typename Cb::template Rebind_TDS<TDS2>::Other Cb2;
        typedef Triangulation_dSimplexAdapter_3<GT, Ih, Cb2> Other;
    };

    Triangulation_dSimplexAdapter_3(const Ih &idHandler)
            : Cb(), m_idHandler(idHandler) {

        m_id = m_idHandler.getId();
    }

    Triangulation_dSimplexAdapter_3(const Triangulation_dSimplexAdapter_3 &o) // do not copy the ID, just the faces and neighbors
            : Cb(o), m_idHandler(o.m_idHandler) {

        m_id = m_idHandler.getId();
    }

    Triangulation_dSimplexAdapter_3(Vertex_handle v0, Vertex_handle v1,
                                    Vertex_handle v2, Vertex_handle v3,
                                    const Ih &idHandler)
            : Cb(v0, v1, v2, v3), m_idHandler(idHandler) {

        m_id = m_idHandler.getId();
    }

    Triangulation_dSimplexAdapter_3(Vertex_handle v0, Vertex_handle v1,
                                    Vertex_handle v2, Vertex_handle v3,
                                    Cell_handle n0, Cell_handle n1,
                                    Cell_handle n2, Cell_handle n3,
                                    const Ih &idHandler)
            : Cb(v0, v1, v2, v3, n0, n1, n2, n3),
              m_idHandler(idHandler) {

        m_id = m_idHandler.getId();
    }

    ~Triangulation_dSimplexAdapter_3() {
        m_idHandler.releaseId(m_id);
    }

public:
    tIdType m_id;
    const Ih &m_idHandler;
};

template<uint D, typename Precision, class Tria, bool Parallel = false>
class CGALHelper;

template<typename Precision, class Tria, bool Parallel>
class CGALHelper<2, Precision, Tria, Parallel> {

public:
    using Handle = typename Tria::Face_handle;
    using Iterator = typename Tria::Finite_faces_iterator;

public:
    CGALHelper(__attribute__((unused)) const dBox<2, Precision> &bounds,
               __attribute__((unused)) const uint N) { }

    Iterator begin(const Tria &t) {
        return t.finite_faces_begin();
    }

    Iterator end(const Tria &t) {
        return t.finite_faces_end();
    }

    std::size_t finite_size(Tria &t) { return t.number_of_faces(); }
    std::size_t infinite_size(Tria &t) {

        auto size = t.number_of_faces();
        auto fc = t.incident_faces(t.infinite_vertex()), done(fc);
        if (fc != 0) {
            do { ++size;
            }while(++fc != done);
        }

        return size;
    }

    typename Tria::Point make_point(const dPoint<2, Precision> &p) {
        return typename Tria::Point(p.coords[0], p.coords[1]);
    }

    Tria make_tria() {
        Tria t;
        return t;
    }

    void insert(Tria &tria, const Point_Ids &ids, dPoints<2, Precision> &points) {
        // transform points into CGAL points with info
        auto transform = [&points, this](const tIdType i) -> std::pair<typename Tria::Point, tIdType> {
            const auto &p = points[i];
            return std::make_pair(make_point(p), i);
        };

        tria.insert(boost::make_transform_iterator(ids.begin(), transform),
                    boost::make_transform_iterator(ids.end(), transform));
    }
};

template<typename Precision, class Tria>
class CGALHelper<3, Precision, Tria, false> {

public:
    using Handle = typename Tria::Cell_handle;
    using Iterator = typename Tria::Finite_cells_iterator;

public:
    CGALHelper(__attribute__((unused)) const dBox<3, Precision> &bounds,
               __attribute__((unused)) const uint N) { }

    Iterator begin(const Tria &t) {
        return t.finite_cells_begin();
    }

    Iterator end(const Tria &t) {
        return t.finite_cells_end();
    }

    std::size_t finite_size(Tria &t) { return t.number_of_finite_cells(); }
    std::size_t infinite_size(Tria &t) { return t.number_of_cells(); }

    typename Tria::Point make_point(const dPoint<3, Precision> &p) {
        return typename Tria::Point(p.coords[0], p.coords[1], p.coords[2]);
    }

    Tria make_tria() {
        Tria t;
        return t;
    }

    /*void insert(Tria & tria, const Point_Ids &ids, dPoints<3, Precision> &points){

        tria.custom_insert_with_info(ids, points);
    }*/

    void insert(Tria &tria, const Point_Ids &ids, dPoints<3, Precision> &points) {
        // transform points into CGAL points with info
        auto transform = [&points, this](const tIdType i) -> std::pair<typename Tria::Point, tIdType> {
            const auto &p = points[i];
            return std::make_pair(make_point(p), i);
        };

        tria.insert(boost::make_transform_iterator(ids.begin(), transform),
                    boost::make_transform_iterator(ids.end(), transform));
    }
};

template<typename Precision, class Tria>
class CGALHelper<3, Precision, Tria, true> {

public:
    using Handle = typename Tria::Cell_handle;
    using Iterator = typename Tria::Finite_cells_iterator;

public:
    CGALHelper(const dBox<3, Precision> &bounds, const uint N)
            : lockingDS(CGAL::Bbox_3(bounds.low[0], bounds.low[1], bounds.low[2],
                                     bounds.high[0], bounds.high[1], bounds.high[2]),
                        std::min(
                                std::max(1u, N),
                                (uint) std::floor(std::cbrt(std::numeric_limits<int>::max()))
                        )) { }

    Iterator begin(const Tria &t) {
        return t.finite_cells_begin();
    }

    Iterator end(const Tria &t) {
        return t.finite_cells_end();
    }

    std::size_t finite_size(Tria &t) { return t.number_of_finite_cells(); }
    std::size_t infinite_size(Tria &t) { return t.number_of_cells(); }

    typename Tria::Point make_point(const dPoint<3, Precision> &p) {
        return typename Tria::Point(p.coords[0], p.coords[1], p.coords[2]);
    }

    Tria make_tria() {
        Tria t(&lockingDS);
        return t;
    }

    void insert(Tria &tria, const Point_Ids &ids, dPoints<3, Precision> &points) {
        // transform points into CGAL points with info
        auto transform = [&points, this](const tIdType i) -> std::pair<typename Tria::Point, tIdType> {
            const auto &p = points[i];
            return std::make_pair(make_point(p), i);
        };

        tria.insert(boost::make_transform_iterator(ids.begin(), transform),
                    boost::make_transform_iterator(ids.end(), transform));
    }

private:
    typename Tria::Lock_data_structure lockingDS;
};

template<uint D, typename Precision, class Tria, bool Parallel>
dSimplices<D, Precision> _delaunayCgal(const Point_Ids &ids, dPoints<D, Precision> &points,
                                       const dBox<D, Precision> &bounds,
                                       const uint gridOccupancy
        /*, bool filterInfinite */) {

    PLOG("delaunayCGAL called with " << ids.size() << " points and bounds " << bounds << std::endl);

    CGALHelper<D, Precision, Tria, Parallel> helper(bounds, gridOccupancy);

    VTUNE_TASK(CgalTriangulation);
    Tria t = helper.make_tria();
    helper.insert(t, ids, points);
    VTUNE_END_TASK(CgalTriangulation);

    ASSERT(t.is_valid());

    PLOG("CGAL triangulation is " << (t.is_valid() ? "" : "NOT ") << "valid" << std::endl);

    PLOG("Collecting simplices" << std::endl);
    INDENT

    VTUNE_TASK(CollectCgal);
    //auto triaSize = helper.size(t);
    auto lastId = t.tds().maxId();
    t.tds().disableId();

    tIdType startId = dSimplices<D, Precision>::simplexID.fetch_add(lastId);

    dSimplices<D, Precision> DT(startId, startId + lastId, 2*(helper.infinite_size(t) - helper.finite_size(t)));

#ifndef NDEBUG
    std::set<tIdType> idCheck;
#endif

    for (auto it = helper.begin(t); it != helper.end(t); ++it) {
        dSimplex<D, Precision> a;
        a.id = startId + it->m_id;

        ASSERT(idCheck.insert(a.id).second);

        for (uint d = 0; d < D + 1; ++d) {
            a.vertices[d] = it->vertex(d)->info();
            a.neighbors[d] = t.is_infinite(it->neighbor(d)) ? dSimplex<D, Precision>::cINF
                                                            : startId + it->neighbor(d)->m_id;
        }

        // sort vertices by ascending point id
        a.sortVertices();
        a.genFingerprint();

        //check whether vertex belongs to the convex hull
        if (!a.isFinite())
            DT.convexHull.push_back(a.id);

        ASSERT((a.id != dSimplex<D, Precision>::cINF));

        PLOG(a << std::endl);

        DT.unsafe_at(a.id) = std::move(a);
    }
    DEDENT

    return DT;
}

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

// partial specializations

template<typename Precision, bool Parallel>
class CGALTriangulator<2, Precision, Parallel>
        : public Triangulator<2, Precision> {

public:
    CGALTriangulator(const dBox<2, Precision> &_bounds, dPoints<2, Precision> &_points, const uint _gridOccupancy)
            : Triangulator<2, Precision>(_bounds, _points), gridOccupancy(_gridOccupancy) { };

protected:
    dSimplices<2, Precision> _triangulate(const Point_Ids &ids,
                                          const dBox<2, Precision> &bounds,
                                          __attribute__((unused)) const std::string provenance
            /*, bool filterInfinite */) {

        typedef CGAL::Triangulation_vertex_base_with_info_2<tIdType, K> Vb;
        typedef Triangulation_dSimplexAdapter_2<K> Cb;
        typedef CGAL::Indexed_Triangulation_data_structure_2<Vb, Cb> Tds;
        typedef CGAL::Delaunay_triangulation_2<K, Tds> CT;

        return _delaunayCgal<2, Precision, CT, Parallel>(ids, this->points, bounds,
                                                         this->gridOccupancy /*, filterInfinite */);
    }

protected:
    const uint gridOccupancy;
};

template<typename Precision>
class CGALTriangulator<3, Precision, true>
        : public Triangulator<3, Precision> {

public:
    CGALTriangulator(const dBox<3, Precision> &_bounds, dPoints<3, Precision> &_points, const uint _gridOccupancy)
            : Triangulator<3, Precision>(_bounds, _points), gridOccupancy(_gridOccupancy) { };

protected:

    dSimplices<3, Precision> _triangulate(const Point_Ids &ids,
                                          const dBox<3, Precision> &bounds,
                                          __attribute__((unused)) const std::string provenance
            /*, bool filterInfinite */) {

        typedef CGAL::Triangulation_vertex_base_with_info_3<tIdType, K> Vb;
        typedef Triangulation_dSimplexAdapter_3<K, _detail::Concurrent_IndexHandler> Cb;
        typedef CGAL::Indexed_Triangulation_data_structure_3<Vb, Cb, CGAL::Parallel_tag, _detail::Concurrent_IndexHandler> Tds;
        typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Default, CGAL::Spatial_lock_grid_3<CGAL::Tag_non_blocking>> CT;

        return _delaunayCgal<3, Precision, CT, true>(ids, this->points, bounds,
                                                     this->gridOccupancy /*, filterInfinite */);
    }

protected:
    const uint gridOccupancy;
};

template<typename Precision>
class CGALTriangulator<3, Precision, false>
        : public Triangulator<3, Precision> {

public:
    CGALTriangulator(const dBox<3, Precision> &_bounds, dPoints<3, Precision> &_points, const uint _gridOccupancy)
            : Triangulator<3, Precision>(_bounds, _points), gridOccupancy(_gridOccupancy) { };

protected:

    dSimplices<3, Precision> _triangulate(const Point_Ids &ids,
                                          const dBox<3, Precision> &bounds,
                                          __attribute__((unused)) const std::string provenance
            /*, bool filterInfinite */) {

        typedef CGAL::Triangulation_vertex_base_with_info_3<tIdType, K> Vb;
        typedef Triangulation_dSimplexAdapter_3<K> Cb;
        typedef CGAL::Indexed_Triangulation_data_structure_3<Vb, Cb, CGAL::Sequential_tag> Tds;
        typedef CGAL::Delaunay_triangulation_3<K, Tds> CT;

        return _delaunayCgal<3, Precision, CT, false>(ids, this->points, bounds,
                                                      this->gridOccupancy /*, filterInfinite */);
    }

protected:
    const uint gridOccupancy;
};

// specializations
template
class CGALTriangulator<2, float, true>;

template
class CGALTriangulator<2, float, false>;

template
class CGALTriangulator<3, float, true>;

template
class CGALTriangulator<3, float, false>;

template
class CGALTriangulator<2, double, true>;

template
class CGALTriangulator<2, double, false>;

template
class CGALTriangulator<3, double, true>;

template
class CGALTriangulator<3, double, false>;

// pure CGAL triangulators for comparision study

template<uint D, typename Precision, class Tria, bool Parallel>
dSimplices<D, Precision> _pureCgal(const Point_Ids &ids, dPoints<D, Precision> &points,
                                   const dBox<D, Precision> &bounds,
                                   const uint gridOccupancy
        /*, bool filterInfinite */) {

    CGALHelper<D, Precision, Tria, Parallel> helper(bounds, gridOccupancy);

    // transform points into CGAL points with info
    auto transform = [&points, &helper](const tIdType i) -> std::pair<typename Tria::Point, tIdType> {
        const auto &p = points[i];
        return std::make_pair(helper.make_point(p), i);
    };

    Tria t = helper.make_tria();
    t.insert(boost::make_transform_iterator(ids.begin(), transform),
             boost::make_transform_iterator(ids.end(), transform));

    dSimplices<D, Precision> dummy(0, 1, 1);
    return dummy;
}

// partial specializations

template<typename Precision, bool Parallel>
class PureCGALTriangulator<2, Precision, Parallel>
        : public Triangulator<2, Precision> {

public:
    PureCGALTriangulator(const dBox<2, Precision> &_bounds, dPoints<2, Precision> &_points, const uint _gridOccupancy)
            : Triangulator<2, Precision>(_bounds, _points), gridOccupancy(_gridOccupancy) { };

protected:
    dSimplices<2, Precision> _triangulate(const Point_Ids &ids,
                                          const dBox<2, Precision> &bounds,
                                          __attribute__((unused)) const std::string provenance
            /*, bool filterInfinite */) {

        typedef CGAL::Triangulation_vertex_base_with_info_2<tIdType, K> Vb;
        typedef CGAL::Triangulation_face_base_2<K> Cb;
        typedef CGAL::Triangulation_data_structure_2<Vb, Cb> Tds;
        typedef CGAL::Delaunay_triangulation_2<K, Tds> CT;

        return _pureCgal<2, Precision, CT, Parallel>(ids, this->points, bounds,
                                                     this->gridOccupancy /*, filterInfinite */);
    }

protected:
    const uint gridOccupancy;
};

template<typename Precision>
class PureCGALTriangulator<3, Precision, true>
        : public Triangulator<3, Precision> {

public:
    PureCGALTriangulator(const dBox<3, Precision> &_bounds, dPoints<3, Precision> &_points, const uint _gridOccupancy)
            : Triangulator<3, Precision>(_bounds, _points), gridOccupancy(_gridOccupancy) { };

protected:

    dSimplices<3, Precision> _triangulate(const Point_Ids &ids,
                                          const dBox<3, Precision> &bounds,
                                          __attribute__((unused)) const std::string provenance
            /*, bool filterInfinite */) {

        typedef CGAL::Triangulation_vertex_base_with_info_3<tIdType, K> Vb;
        typedef CGAL::Triangulation_cell_base_3<K> Cb;
        typedef CGAL::Triangulation_data_structure_3<Vb, Cb, CGAL::Parallel_tag> Tds;
        typedef CGAL::Delaunay_triangulation_3<K, Tds> CT;

        return _pureCgal<3, Precision, CT, true>(ids, this->points, bounds,
                                                 this->gridOccupancy /*, filterInfinite */);
    }

protected:
    const uint gridOccupancy;
};

template<typename Precision>
class PureCGALTriangulator<3, Precision, false>
        : public Triangulator<3, Precision> {

public:
    PureCGALTriangulator(const dBox<3, Precision> &_bounds, dPoints<3, Precision> &_points, const uint _gridOccupancy)
            : Triangulator<3, Precision>(_bounds, _points), gridOccupancy(_gridOccupancy) { };

protected:

    dSimplices<3, Precision> _triangulate(const Point_Ids &ids,
                                          const dBox<3, Precision> &bounds,
                                          __attribute__((unused)) const std::string provenance
            /*, bool filterInfinite */) {

        typedef CGAL::Triangulation_vertex_base_with_info_3<tIdType, K> Vb;
        typedef CGAL::Triangulation_cell_base_3<K> Cb;
        typedef CGAL::Triangulation_data_structure_3<Vb, Cb, CGAL::Sequential_tag> Tds;
        typedef CGAL::Delaunay_triangulation_3<K, Tds> CT;

        return _pureCgal<3, Precision, CT, false>(ids, this->points, bounds,
                                                  this->gridOccupancy /*, filterInfinite */);
    }

protected:
    const uint gridOccupancy;
};

// specializations
template
class PureCGALTriangulator<2, float, true>;

template
class PureCGALTriangulator<2, float, false>;

template
class PureCGALTriangulator<3, float, true>;

template
class PureCGALTriangulator<3, float, false>;

template
class PureCGALTriangulator<2, double, true>;

template
class PureCGALTriangulator<2, double, false>;

template
class PureCGALTriangulator<3, double, true>;

template
class PureCGALTriangulator<3, double, false>;

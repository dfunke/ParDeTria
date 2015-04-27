#include "CGALTriangulator.h"

#include <atomic>
#include <type_traits>

#include "utils/ASSERT.h"

// define a static counter for the tetrahedronID
std::atomic<uint> gAtomicTetrahedronID(0);
std::atomic<uint> gAtomicCgalID(0);

// CGAL
#define CGAL_LINKED_WITH_TBB
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Delaunay_triangulation_3.h> use modified version
#include "mods/Delaunay_triangulation_3.h"
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

// boost
#include <boost/iterator/transform_iterator.hpp>

#include <CGAL/Unique_hash_map.h>



template <typename Precision, typename GT,
        typename Cb = CGAL::Triangulation_cell_base_3<GT> >
class Triangulation_dSimplexAdapter_3
        : public Cb
{
public:
  typedef typename Cb::Vertex_handle                   Vertex_handle;
  typedef typename Cb::Cell_handle                     Cell_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other       Cb2;
    typedef Triangulation_dSimplexAdapter_3<Precision, GT, Cb2>  Other;
  };

  Triangulation_dSimplexAdapter_3()
          : Cb() {

    m_simplex.id = gAtomicCgalID++;
  }

  Triangulation_dSimplexAdapter_3(Vertex_handle v0, Vertex_handle v1,
                                      Vertex_handle v2, Vertex_handle v3)
          : Cb(v0, v1, v2, v3) {

    m_simplex.id = gAtomicCgalID++;
  }

  Triangulation_dSimplexAdapter_3(Vertex_handle v0, Vertex_handle v1,
                                      Vertex_handle v2, Vertex_handle v3,
                                      Cell_handle   n0, Cell_handle   n1,
                                      Cell_handle   n2, Cell_handle   n3)
          : Cb(v0, v1, v2, v3, n0, n1, n2, n3) {

    m_simplex.id = gAtomicCgalID++;
  }

  // SETTING

  void set_vertex(int i, Vertex_handle v)
  {
    Cb::set_vertex(i, v);
    m_simplex.vertices[i] = v->info();
  }

  void set_neighbor(int i, Cell_handle n)
  {
    Cb::set_neighbor(i, n);
    if(n != Cell_handle())
      m_simplex.neighbors[i] = n->m_simplex.id;
  }

  void set_vertices()
  {
    Cb::set_vertices();
    m_simplex.vertices[0] = dPoint<3, Precision>::cINF;
    m_simplex.vertices[1] = dPoint<3, Precision>::cINF;
    m_simplex.vertices[2] = dPoint<3, Precision>::cINF;
    m_simplex.vertices[3] = dPoint<3, Precision>::cINF;
  }

  void set_vertices(Vertex_handle v0, Vertex_handle v1,
                    Vertex_handle v2, Vertex_handle v3)
  {
    Cb::set_vertices(v0, v1, v2, v3);
    m_simplex.vertices[0] = v0->info();
    m_simplex.vertices[1] = v1->info();
    m_simplex.vertices[2] = v2->info();
    m_simplex.vertices[3] = v3->info();
  }

  void set_neighbors()
  {
    Cb::set_neighbors();
    m_simplex.neighbors[0] = dSimplex<3, Precision>::cINF;
    m_simplex.neighbors[1] = dSimplex<3, Precision>::cINF;
    m_simplex.neighbors[2] = dSimplex<3, Precision>::cINF;
    m_simplex.neighbors[3] = dSimplex<3, Precision>::cINF;
  }

  void set_neighbors(Cell_handle n0, Cell_handle n1,
                     Cell_handle n2, Cell_handle n3)
  {
    Cb::set_neighbors(n0, n1,n2,n3);
    if(n0 != Cell_handle())
      m_simplex.neighbors[0] = n0->m_simplex.id;
    if(n1 != Cell_handle())
      m_simplex.neighbors[1] = n1->m_simplex.id;
    if(n2 != Cell_handle())
      m_simplex.neighbors[2] = n2->m_simplex.id;
    if(n3 != Cell_handle())
      m_simplex.neighbors[3] = n3->m_simplex.id;
  }

public:
  dSimplex<3, Precision> m_simplex;
};

template <uint D, typename Precision, class Tria, bool Parallel = false> class CGALHelper;

template <typename Precision, class Tria, bool Parallel> class CGALHelper<2, Precision, Tria, Parallel> {

public:
  CGALHelper(__attribute__((unused)) const dBox<2, Precision> &bounds,
             __attribute__((unused)) const uint N) {}

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
  CGALHelper(__attribute__((unused)) const dBox<3, Precision> &bounds,
             __attribute__((unused)) const uint N) { }

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
  CGALHelper(const dBox<3, Precision> &bounds, const uint N)
      : lockingDS(CGAL::Bbox_3(bounds.low[0], bounds.low[1], bounds.low[2],
                               bounds.high[0], bounds.high[1], bounds.high[2]),
                               std::min(
                                       std::max(1u, N),
                                       (uint) std::floor(std::cbrt(std::numeric_limits<int>::max()))
                               )) {}

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
              const dBox<D, Precision> &bounds,
              const uint gridOccupancy
              /*, bool filterInfinite */) {

  CGALHelper<D, Precision, Tria, Parallel> helper(bounds, std::cbrt(ids.size() / gridOccupancy));

  // transform points into CGAL points with info
  auto transform = [&points, &helper](const uint i) -> std::pair<typename Tria::Point, uint> {
      const auto &p = points[i];
      return std::make_pair(helper.make_point(p), p.id);
  };

  Tria t = helper.make_tria();
  t.insert(boost::make_transform_iterator(ids.begin(), transform),
           boost::make_transform_iterator(ids.end(), transform));

  ASSERT(t.is_valid());

  // VLOG("CGAL triangulation is " << (t.is_valid() ? "" : "NOT ") << "valid"
  //     << std::endl);
  // VLOG("finite cells/vertices " << t.number_of_finite_cells() << "/"
  //     << t.number_of_vertices() << std::endl);

  PLOG("Collecting simplices" << std::endl);
  INDENT

  dSimplices<D, Precision> tria;
  tria.reserve(helper.size(t));
  tria.wuFaces.reserve(3*helper.size(t));

  CGAL::Unique_hash_map<typename CGALHelper<D, Precision, Tria, Parallel>::Handle, uint>
      simplexLookup(0, helper.size(t));

  uint tetrahedronID = gAtomicTetrahedronID.fetch_add(helper.size(t), std::memory_order::memory_order_relaxed);
#ifndef NDEBUG
  uint saveTetrahedronID = tetrahedronID;
#endif

  dSimplex<D, Precision> a;
  for (auto it = helper.begin(t); it != helper.end(t); ++it) {
    a.id = tetrahedronID++;

    for (uint i = 0; i < D + 1; ++i) {
      dPoint<D, Precision> &point = points[it->vertex(i)->info()];
      ASSERT(point.id == it->vertex(i)->info());

      a.vertices[i] = point.id;
    }

    // sort vertices by ascending point id
    std::sort(a.vertices.begin(), a.vertices.end());
    a.fingerprint();

    //check whether vertex belongs to the convex hull
    if(!a.isFinite())
      tria.convexHull.insert(a.id);

    //update where-used data structure for faces
    for (uint i = 0; i < D + 1; ++i) {
      uint facetteHash = a.vertexFingerprint ^ a.vertices[i];

      tria.wuFaces.emplace(facetteHash, a.id);
    }

    PLOG(a << std::endl);
    tria.insert(a);
    simplexLookup[it] = a.id;
  }
  DEDENT

  ASSERT(tetrahedronID == saveTetrahedronID + helper.size(t));

  PLOG("Collecting neighbors" << std::endl);

  INDENT
  for (auto it = helper.begin(t); it != helper.end(t); ++it) {
    auto &tet = tria[simplexLookup[it]];
    uint neighborIdx = 0;
    for (uint i = 0; i < D + 1; ++i) {
      if (simplexLookup.is_defined(it->neighbor(i))) {
        const auto &n = tria[simplexLookup[it->neighbor(i)]];
        if (n.id != tet.id) {
          tet.neighbors[neighborIdx++] = n.id;
        }
      }
    }

    while(neighborIdx < D + 1)
      tet.neighbors[neighborIdx++] = dSimplex<D, Precision>::cINF;

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
  CGALTriangulator(const dBox<2, Precision> &_bounds, dPoints<2, Precision> &_points, const uint _gridOccupancy)
  : Triangulator<2, Precision>(_bounds, _points), gridOccupancy(_gridOccupancy) {};

protected:
  dSimplices<2, Precision> _triangulate(const Ids &ids,
                                              const dBox<2, Precision> &bounds,
                                              __attribute__((unused)) const std::string provenance
                                              /*, bool filterInfinite */) {

    typedef CGAL::Triangulation_vertex_base_with_info_2<uint, K> Vb;
    typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
    typedef CGAL::Delaunay_triangulation_2<K, Tds> CT;

    return _delaunayCgal<2, Precision, CT, Parallel>(ids, this->points, bounds, this->gridOccupancy /*, filterInfinite */);
  }

protected:
  const uint gridOccupancy;
};

template <typename Precision> class CGALTriangulator<3, Precision, true>
        : public Triangulator<3, Precision> {

public:
  CGALTriangulator(const dBox<3, Precision> &_bounds, dPoints<3, Precision> &_points, const uint _gridOccupancy)
          : Triangulator<3, Precision>(_bounds, _points), gridOccupancy(_gridOccupancy) {};

protected:

  dSimplices<3, Precision> _triangulate(const Ids &ids,
                                              const dBox<3, Precision> &bounds,
                                              __attribute__((unused)) const std::string provenance
                                              /*, bool filterInfinite */) {

    typedef CGAL::Triangulation_vertex_base_with_info_3<uint, K> Vb;
    typedef Triangulation_dSimplexAdapter_3<Precision, K> Cb;
    typedef CGAL::Triangulation_data_structure_3<Vb, Cb, CGAL::Parallel_tag> Tds;
    typedef CGAL::Delaunay_triangulation_3<K, Tds> CT;

    return _delaunayCgal<3, Precision, CT, true>(ids, this->points, bounds, this->gridOccupancy /*, filterInfinite */);
  }

protected:
  const uint gridOccupancy;
};

template <typename Precision> class CGALTriangulator<3, Precision, false>
        : public Triangulator<3, Precision> {

public:
  CGALTriangulator(const dBox<3, Precision> &_bounds, dPoints<3, Precision> &_points, const uint _gridOccupancy)
          : Triangulator<3, Precision>(_bounds, _points), gridOccupancy(_gridOccupancy) {};

protected:

  dSimplices<3, Precision> _triangulate(const Ids &ids,
                                        const dBox<3, Precision> &bounds,
                                        __attribute__((unused)) const std::string provenance
          /*, bool filterInfinite */) {

    typedef CGAL::Triangulation_vertex_base_with_info_3<uint, K> Vb;
    typedef Triangulation_dSimplexAdapter_3<Precision, K> Cb;
    typedef CGAL::Triangulation_data_structure_3<Vb, Cb, CGAL::Sequential_tag> Tds;
    typedef CGAL::Delaunay_triangulation_3<K, Tds> CT;

    return _delaunayCgal<3, Precision, CT, false>(ids, this->points, bounds, this->gridOccupancy /*, filterInfinite */);
  }

protected:
  const uint gridOccupancy;
};

// specializations
template class CGALTriangulator<2, float, true>;
template class CGALTriangulator<2, float, false>;

template class CGALTriangulator<3, float, true>;
template class CGALTriangulator<3, float, false>;

template class CGALTriangulator<2, double, true>;
template class CGALTriangulator<2, double, false>;

template class CGALTriangulator<3, double, true>;
template class CGALTriangulator<3, double, false>;

//#define CGAL_BENCHMARK
#ifdef CGAL_BENCHMARK

#include <boost/program_options.hpp>

#include "utils/Random.h"
#include "utils/Timings.h"
#include "utils/System.h"
#include "utils/DBConnection.h"

#define D 3
#define Precision double
#define REPS 10

int main(int argc, char *argv[]) {

  // define program options

  namespace po = boost::program_options;

  uint maxN, minN = 10;

  uint maxThreads = tbb::task_scheduler_init::default_num_threads();
  uint minThreads = 1;
  uint occ = 0;
  std::string runFile;
  std::string run;

  po::options_description cCommandLine("Command Line Options");
  cCommandLine.add_options()("n", po::value<uint>(&maxN),
                             "maximum number of points");
  cCommandLine.add_options()("minN", po::value<uint>(&minN),
                             "minimum number of points");
  cCommandLine.add_options()("threads", po::value<uint>(&maxThreads),
                             "maximum number of threads");
  cCommandLine.add_options()("minThreads", po::value<uint>(&minThreads),
                             "minimum number of threads");
  cCommandLine.add_options()("occ", po::value<uint>(&occ),
                             "occupancy of grid lock data structure");
  cCommandLine.add_options()("help", "produce help message");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cCommandLine), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << cCommandLine << std::endl;
    return EXIT_SUCCESS;
  }

  //***************************************************************************

  typedef CGAL::Triangulation_vertex_base_with_info_3<uint, K> Vb;
  typedef CGAL::Triangulation_data_structure_3<
          Vb, CGAL::Triangulation_cell_base_3<K>, CGAL::Parallel_tag> Tds;
  typedef CGAL::Delaunay_triangulation_3<K, Tds> CT;

  //setup
  std::vector<uint> occupancies;
  if(vm.count("occ"))
    occupancies = { occ };
  else
    occupancies = { 10, 50, 100, 1000 };

  dBox<D, Precision> bounds(dVector<D, Precision>({{0,0,0}}), dVector<D, Precision>({{100,100,100}}));
  LOGGER.setLogLevel(Logger::Verbosity::SILENT);

  DBConnection db("db_" + getHostname() + ".dat", "pureCGAL");

  //determine the latest run number
  uint runNumber = db.getMaximum<uint>("run-number") + 1;

  //loop over number of points
  for (uint nPoints = minN; nPoints <= maxN; nPoints += pow(10, floor(log10(nPoints)))) {

    tGenerator gen(START_SEED);
    auto dice = RandomFactory<Precision>::make('u', gen);
    auto points = genPoints(nPoints, bounds, dice);

    // copy points into CGAL structure
    std::vector<std::pair<typename CT::Point, uint>> cPoints;
    cPoints.reserve(points.size());
    for (const auto &p : points) {
      auto cp = CT::Point(p.coords[0], p.coords[1], p.coords[2]);
      cPoints.push_back(std::make_pair(cp, p.id));
    }

    //loop over number of threads, if algo is multi-threaded
    for (uint threads = minThreads; threads <= maxThreads; threads <<= 1) {

      // load scheduler with specified number of threads
      tbb::task_scheduler_init init(threads);

      //loop over occupancies
      for(uint occ : occupancies) {

        ExperimentRun run;
        run.addTrait("run-number", runNumber);
        run.addTrait("dist", "u");
        run.addTrait("nP", nPoints);
        run.addTrait("alg", "m");
        run.addTrait("threads", threads);
        run.addTrait("occupancy", occ);

        run.addTrait("start-time", getDatetime());

        std::cout << run.str() << std::endl;

        try {

          for(uint i = 0; i < REPS; ++i) {

            int numGridCells = std::min(
                    std::max(1, (int) std::floor(std::cbrt(points.size() / occ))),
                    (int) std::floor(std::cbrt(std::numeric_limits<int>::max())));

            CT::Lock_data_structure lds(
                    CGAL::Bbox_3(bounds.low[0], bounds.low[1], bounds.low[2],
                         bounds.high[0], bounds.high[1], bounds.high[2]),
                    numGridCells);

            CT tria(&lds);

            auto t1 = Clock::now();
            tria.insert(cPoints.begin(), cPoints.end());
            auto t2 = Clock::now();

            run.addMemory(getCurrentRSS());
            run.addTime(std::chrono::duration_cast<tDuration>(t2 - t1));
          }

        } catch (std::exception &e) {
          std::cerr << "Experiment: " << run.str() << std::endl;
          std::cerr << "\tException raised: " << e.what() << std::endl;
        }

        run.addTrait("end-time", getDatetime());
        db.save(run);

        std::cout << "\tAverage time: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(run.avgTime()).count()
        << " ms\tAverage mem: "
        << run.avgMem() / 1e6 << " MB" << std::endl;

        }

      }
    }

  return EXIT_SUCCESS;
}

#endif
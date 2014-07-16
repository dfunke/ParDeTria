cmake_minimum_required (VERSION 2.6)

project (Delaunay)

find_package(CGAL)
find_package(root)

include_directories(lib/delaunay lib/fortune lib/triangle)

add_executable (Delaunay Delaunay.cpp)
target_link_libraries(Delaunay delaunay voronoi triangle Hist CGAL)

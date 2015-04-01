# CMake script for src directory

# build utils
add_subdirectory(utils)

add_library (Geometry
  Geometry.cpp
)

target_link_libraries(Geometry Logging)

add_library(InterfaceCGAL
  CGAL_Interface.cpp
)

target_link_libraries(InterfaceCGAL
  Logging
  Geometry
  ${CGAL_LIBRARY}
  ${GMP_LIBRARIES}
  ${TBB_MALLOC_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
)

add_library(Painting2D Painter_2.cxx)
target_link_libraries(Painting2D ${Cairomm_LIBRARIES})

add_library(Painting3D Painter_3.cxx)
target_link_libraries(Painting3D ${VTK_LIBRARIES})

add_library(Painting
  Painter.cpp
)

target_link_libraries(Painting
  Logging
  Geometry
  ${Boost_LIBRARIES}
  debug Painting2D
  debug Painting3D
)

add_library(Triangulate
  Triangulator.cpp
  Partitioner.cpp
)

target_link_libraries(Triangulate Painting Logging Geometry)

add_executable (DACTest
  DAC_Test.cpp
)

target_link_libraries(DACTest
  Logging
  Measurements
  Geometry
  Triangulate
  Painting
  InterfaceCGAL
  ${Boost_LIBRARIES}
  ${TBB_MALLOC_PROXY_LIBRARIES}
)

add_executable (Benchmark
  Benchmark.cpp
)

target_link_libraries(Benchmark
  Logging
  Measurements
  Geometry
  Triangulate
  Painting
  InterfaceCGAL
  ${Boost_LIBRARIES}
  ${TBB_MALLOC_PROXY_LIBRARIES}
)

add_executable(EdgeStudy
  EdgeStudy.cpp
)

target_link_libraries(EdgeStudy
  Logging
  Geometry
  Triangulate
  InterfaceCGAL
)

add_executable(WhereUsedStudy
  WhereUsedStudy.cpp
)

target_link_libraries(WhereUsedStudy
  Logging
  Geometry
  Triangulate
  InterfaceCGAL
)


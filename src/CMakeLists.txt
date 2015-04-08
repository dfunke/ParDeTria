# CMake script for src directory

# build utils
add_subdirectory(utils)

add_library (Geometry
  Geometry.cpp
)

target_link_libraries(Geometry Logging)

add_library(CGALTriangulator
  CGALTriangulator.cpp
)

add_library(Triangulator
  Triangulator.cpp
)

target_link_libraries(CGALTriangulator
  Triangulator
  Logging
  Geometry
  ${CGAL_LIBRARY}
  ${GMP_LIBRARIES}
  ${TBB_MALLOC_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
)

add_library(DCTriangulator
  DCTriangulator.cpp
  Partitioner.cpp
)

target_link_libraries(DCTriangulator
  Triangulator
  CGALTriangulator
  Painting
  Logging
  Geometry
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

add_executable (DACTest
  DAC_Test.cpp
)

target_link_libraries(DACTest
  Logging
  Measurements
  Geometry
  DCTriangulator
  Painting
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
  DCTriangulator
  Painting
  ${Boost_LIBRARIES}
  ${TBB_MALLOC_PROXY_LIBRARIES}
)

add_executable(EdgeStudy
  EdgeStudy.cpp
)

target_link_libraries(EdgeStudy
  Logging
  Geometry
  DCTriangulator
)

add_executable(WhereUsedStudy
  WhereUsedStudy.cpp
)

target_link_libraries(WhereUsedStudy
  Logging
  Geometry
  DCTriangulator
)


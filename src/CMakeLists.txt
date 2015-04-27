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
  Logging
  Geometry
)

add_executable (DACTest
  DAC_Test.cpp
)

target_link_libraries(DACTest
  Logging
  Measurements
  Geometry
  DCTriangulator
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
  DBConn
  ${Boost_LIBRARIES}
  ${TBB_MALLOC_PROXY_LIBRARIES}
  ${EJDB_LIBRARIES}
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

add_executable(CGALBenchmark
  CGALTriangulator.cpp
)

target_link_libraries(CGALBenchmark
  Triangulator
  Logging
  Geometry
  Measurements
  DBConn
  ${CGAL_LIBRARY}
  ${GMP_LIBRARIES}
  ${TBB_MALLOC_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  ${Boost_LIBRARIES}
  ${EJDB_LIBRARIES}
)

set_target_properties(CGALBenchmark PROPERTIES COMPILE_DEFINITIONS "CGAL_BENCHMARK")


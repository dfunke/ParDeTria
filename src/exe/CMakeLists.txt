# CMake script for src directory

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


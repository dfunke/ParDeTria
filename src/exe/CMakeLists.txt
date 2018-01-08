# CMake script for src directory
add_executable (LBExperiments
  LB_Experiments.cpp
)

target_link_libraries(LBExperiments
  Logging
  Measurements
  Geometry
  Triangulator
  CGALTriangulator
  SpatialSort
  Partitioner
  kahip
  ${Boost_LIBRARIES}
  ${TBB_MALLOC_PROXY_LIBRARIES}
)

if(PAINT)
add_executable(LBPainting
  LB_Painting.cpp
)

target_link_libraries(LBPainting
  Painter
  Logging
  Measurements
  Geometry
  Triangulator
  CGALTriangulator
  SpatialSort
  Partitioner
  kahip
  ${Boost_LIBRARIES}
  ${TBB_MALLOC_PROXY_LIBRARIES}
)
endif()  

add_executable (DACTest
  DAC_Test.cpp
)

target_link_libraries(DACTest
  Logging
  Measurements
  Geometry
  Triangulator
  DCTriangulator
  kahip
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

add_executable (Studies
  Studies.cpp
)

target_link_libraries(Studies
  Logging
  Measurements
  Geometry
  DCTriangulator
  DBConn
  ${Boost_LIBRARIES}
  ${TBB_MALLOC_PROXY_LIBRARIES}
  ${EJDB_LIBRARIES}
)

add_executable (BenchmarkDS
  BenchmarkDS.cpp
)

target_link_libraries(BenchmarkDS
  Logging
  Measurements
  ${TBB_MALLOC_PROXY_LIBRARIES}
)

add_executable (DS_Study
        DS_Studies.cpp
)

target_link_libraries(DS_Study
        Logging
        Measurements
        Geometry
        DCTriangulator
        ${Boost_LIBRARIES}
        ${TBB_MALLOC_PROXY_LIBRARIES}
)

add_executable (getRun
        GetRun.cpp
        )

target_link_libraries(getRun
        DBConn
        Measurements
        ${EJDB_LIBRARIES}
        )

file(INSTALL runExperiments.sh DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} USE_SOURCE_PERMISSIONS)

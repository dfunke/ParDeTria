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
  ${KaHIP_LIBRARIES}
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${TBB_LIBRARIES}
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
  ${KaHIP_LIBRARIES}
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${TBB_LIBRARIES}
)
endif()  

add_executable (DACTest EXCLUDE_FROM_ALL
  DAC_Test.cpp
)

target_link_libraries(DACTest
  Logging
  Measurements
  Geometry
  Triangulator
  DCTriangulator
  ${KaHIP_LIBRARIES}
  ${Boost_LIBRARIES}
  ${TBB_LIBRARIES}
)

add_executable (Benchmark EXCLUDE_FROM_ALL
  Benchmark.cpp
)

target_link_libraries(Benchmark
  Logging
  Measurements
  Geometry
  DCTriangulator
  DBConn
  ${Boost_LIBRARIES}
  ${TBB_LIBRARIES}
  ${EJDB_LIBRARES}
)

add_executable (Studies EXCLUDE_FROM_ALL
  Studies.cpp
)

target_link_libraries(Studies
  Logging
  Measurements
  Geometry
  DCTriangulator
  DBConn
  ${Boost_LIBRARIES}
  ${TBB_LIBRARIES}
  ${EJDB_LIBRARIES}
)

add_executable (BenchmarkDS EXCLUDE_FROM_ALL
  BenchmarkDS.cpp
)

target_link_libraries(BenchmarkDS
  Logging
  Measurements
  ${TBB_LIBRARIES}
)

add_executable (DS_Study EXCLUDE_FROM_ALL
        DS_Studies.cpp
)

target_link_libraries(DS_Study
        Logging
        Measurements
        Geometry
        DCTriangulator
        ${Boost_LIBRARIES}
        ${TBB_LIBRARIES}
)

add_executable (getRun EXCLUDE_FROM_ALL
        GetRun.cpp
        )

target_link_libraries(getRun
        DBConn
        Measurements
        ${EJDB_LIBRARIES}
        )

file(INSTALL runExperiments.sh DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} USE_SOURCE_PERMISSIONS)

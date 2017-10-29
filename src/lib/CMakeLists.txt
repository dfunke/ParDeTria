# CMake script for src directory

add_library (Geometry
  Geometry.cpp
)

target_link_libraries(Geometry Predicates Logging Measurements)

add_library(Triangulator
  Triangulator.cpp
)

add_library(LoadBalancedDCTriangulator
  LoadBalancedDCTriangulator.cpp
  Partitioner.cpp
)

target_link_libraries(LoadBalancedDCTriangulator
  Triangulator
  CGALTriangulator
  SpatialSort
  Logging
  Geometry
)

add_library(CGALTriangulator
  CGALTriangulator.cpp
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

add_library(SpatialSort
  SpatialSort.cpp
)

target_link_libraries(SpatialSort
        ${CGAL_LIBRARY}
)

target_link_libraries(CGALTriangulator
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
  SpatialSort
  Logging
  Geometry
)

if(VTUNE)
  target_link_libraries(DCTriangulator VTune)
  target_link_libraries(CGALTriangulator VTune)
  target_link_libraries(SpatialSort VTune)
endif()

# CMake script for src directory

if(PAINT)
add_library (Painter
  Painter.cpp
  Painter_2.cxx
)

target_link_libraries(Painter Geometry Logging ${Cairomm_LIBRARIES})
endif()

add_library (Geometry EXCLUDE_FROM_ALL
  Geometry.cpp
)

target_link_libraries(Geometry Predicates Logging Measurements)

add_library(Triangulator EXCLUDE_FROM_ALL
  Triangulator.cpp
)

add_library(CGALTriangulator EXCLUDE_FROM_ALL
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

add_library(SpatialSort EXCLUDE_FROM_ALL
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

add_library(DCTriangulator EXCLUDE_FROM_ALL
  DCTriangulator.cpp
)

target_link_libraries(DCTriangulator
  Triangulator
  CGALTriangulator
  SpatialSort
  Logging
  Geometry
  Partitioner
)

add_library(Partitioner EXCLUDE_FROM_ALL
  Partitioner.cpp
)

if(VTUNE)
  target_link_libraries(DCTriangulator VTune)
  target_link_libraries(CGALTriangulator VTune)
  target_link_libraries(SpatialSort VTune)
endif()

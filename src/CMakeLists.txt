cmake_minimum_required (VERSION 3.0)
cmake_policy(VERSION 3.0)

project (Delaunay)

# adds the project-specific cmake module directory cmake/Modules to the cmake
# search path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(CGAL REQUIRED)
find_package(ROOT REQUIRED COMPONENTS Hist)
find_package(Cairomm REQUIRED)
find_package(TBB REQUIRED)

add_subdirectory(lib)

set(BUILTIN_LIBS lib/delaunay lib/fortune lib/triangle)

include_directories(${BUILTIN_LIBS} ${ROOT_INCLUDE_DIRS} ${CGAL_INCLUDE_DIRS} ${Cairomm_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -std=c++1y -Wall")

add_executable (Delaunay Delaunay.cpp)
target_link_libraries(Delaunay delaunay fortune triangle pthread ${ROOT_LIBRARIES} ${CGAL_LIBRARY} ${GMP_LIBRARIES})

add_executable (DelaunayRoot DelaunayRoot.cpp)
target_link_libraries(DelaunayRoot ${ROOT_LIBRARIES} ${CGAL_LIBRARY})

add_executable (DACTest DAC_Test.cpp Geometry.cpp Logger.cpp Painter.cpp CGAL_Interface.cpp)
target_link_libraries(DACTest ${CGAL_LIBRARY} ${GMP_LIBRARIES} ${Cairomm_LIBRARIES} ${TBB_LIBRARIES})

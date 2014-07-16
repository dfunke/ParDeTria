cmake_minimum_required (VERSION 3.0)
cmake_policy(VERSION 3.0)

project (Delaunay)

find_package(CGAL REQUIRED)
find_package(ROOT REQUIRED COMPONENTS Hist)

if(ROOT_FOUND)
	message(STATUS "Found ROOT")
endif(ROOT_FOUND)

set(BUILTIN_LIBS lib/delaunay lib/fortune lib/triangle)

include_directories(${BUILTIN_LIBS} ${ROOT_INCLUDE_DIRS} ${CGAL_INCLUDE_DIRS})
link_directories(${BUILTIN_LIBS})

set(CMAKE_CXX_FLAGS "-std=c++1y -DLINUX -DNO_TIMER -DANSI_DECLARATORS -O3 -g -Wall -fmessage-length=0 -lm")
set(CMAKE_EXE_LINKER_FLAGS "-g -DJUHU")

add_custom_target(
	libDelaunay
	COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/lib/delaunay
)

add_custom_target(
	libTriangle
	COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/lib/triangle trilibrary
)

add_custom_target(
	libFortune
	COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/lib/fortune libvoronoi.a
)
add_executable (Delaunay Delaunay.cpp)

target_link_libraries(Delaunay delaunay voronoi triangle pthread ${ROOT_LIBRARIES} ${CGAL_LIBRARY} ${GMP_LIBRARIES})
add_dependencies(Delaunay libDelaunay libTriangle libFortune)

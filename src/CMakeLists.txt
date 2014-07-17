cmake_minimum_required (VERSION 3.0)
cmake_policy(VERSION 3.0)

project (Delaunay)

find_package(CGAL REQUIRED)
find_package(ROOT REQUIRED COMPONENTS Hist)

add_subdirectory(lib)

set(BUILTIN_LIBS lib/delaunay lib/fortune lib/triangle)

include_directories(${BUILTIN_LIBS} ${ROOT_INCLUDE_DIRS} ${CGAL_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "-std=c++1y -DLINUX -DNO_TIMER -DANSI_DECLARATORS -O3 -g -Wall -fmessage-length=0 -lm")

add_executable (Delaunay Delaunay.cpp)

target_link_libraries(Delaunay delaunay fortune triangle pthread ${ROOT_LIBRARIES} ${CGAL_LIBRARY} ${GMP_LIBRARIES})

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

include_directories(${ROOT_INCLUDE_DIRS} ${CGAL_INCLUDE_DIRS} ${Cairomm_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -std=c++1y -pedantic -Wall -Wextra -Werror")

if (CMAKE_CXX_COMPILER MATCHES ".*clang")
  add_definitions(-D_LIBCPP_VERSION)
endif ()

if (VALGRIND)
  message(STATUS "Valgrind mode enabled")
  add_definitions(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)
endif ()

if (STUDY)
  message(STATUS "Study mode enabled")
  add_definitions(-DSTUDY)

  if  (NOT FORCE_OUTPUT)
    add_definitions(-DNO_OUTPUT)
  else ()
    message(STATUS "Forcing PNG output")
  endif ()
endif ()

add_executable (DACTest DAC_Test.cpp Geometry.cpp Logger.cpp Painter.cpp CGAL_Interface.cpp Triangulator.cpp Partitioner.cpp Timings.cpp RSS.c)
target_link_libraries(DACTest ${CGAL_LIBRARY} ${GMP_LIBRARIES} ${Cairomm_LIBRARIES} ${TBB_LIBRARIES})

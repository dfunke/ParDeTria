cmake_minimum_required (VERSION 3.0)
cmake_policy(VERSION 3.0)

project (Delaunay)

enable_testing()
find_package(GTest)

# adds the project-specific cmake module directory cmake/Modules to the cmake
# search path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(CGAL REQUIRED)
find_package(Cairomm REQUIRED)
find_package(TBB REQUIRED)
find_package(VTK REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem serialization program_options)

include_directories(${CMAKE_SOURCE_DIR}
  ${ROOT_INCLUDE_DIRS}
  ${CGAL_INCLUDE_DIRS}
  ${Cairomm_INCLUDE_DIRS}
)

include(${VTK_USE_FILE})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -std=c++1y -pedantic -Wall -Wextra -Werror -Wno-extra-semi")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -std=c11 -pedantic -Wall -Wextra -Werror -Wno-extra-semi")

if (VALGRIND)
  message(STATUS "Valgrind mode enabled")
  add_definitions(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)
endif ()

if (SANITIZE)
  message(WARNING "ThreadSanitizer is enabled")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=thread")
endif()

# add git-cmake-format
add_subdirectory(${CMAKE_SOURCE_DIR}/cmake/clang-format)

add_library (Logging
  utils/Logger.cpp
  utils/ccolor.cpp
  utils/StringUtils.cxx
)

target_link_libraries(Logging ${TBB_LIBRARIES})

add_library (Geometry
  Geometry.cpp
)

target_link_libraries(Geometry Logging)

add_library(InterfaceCGAL
  CGAL_Interface.cpp
)

target_link_libraries(InterfaceCGAL Logging Geometry ${CGAL_LIBRARY} ${GMP_LIBRARIES})

add_library(Painting
  Painter.cpp
)

target_link_libraries(Painting Logging Geometry ${Cairomm_LIBRARIES} ${VTK_LIBRARIES} ${Boost_LIBRARIES})

add_library(Triangulate
  Triangulator.cpp
  Partitioner.cpp
)

target_link_libraries(Triangulate Painting Logging Geometry)

add_executable (DACTest
  DAC_Test.cpp
  utils/Timings.cpp
  utils/RSS.c
)

target_link_libraries(DACTest
  Logging
  Geometry
  Triangulate
  Painting
  InterfaceCGAL
  ${Boost_LIBRARIES}
  ${TBB_MALLOC_PROXY_LIBRARIES}
)

add_executable (Benchmark
  Benchmark.cpp
  utils/Timings.cpp
  utils/RSS.c
)

target_link_libraries(Benchmark
  Logging
  Geometry
  Triangulate
  Painting
  InterfaceCGAL
  ${Boost_LIBRARIES}
)

add_executable (VTUShow
  utils/VTU_Show.cpp
)

target_link_libraries(VTUShow ${VTK_LIBRARIES} ${Boost_LIBRARIES})

add_executable (archiveConverter
  utils/ArchiveConverter.cpp
)

target_link_libraries(archiveConverter ${Boost_LIBRARIES})

add_executable(EdgeStudy
  EdgeStudy.cpp
)

target_link_libraries(EdgeStudy
  Logging
  Geometry
  Triangulate
  InterfaceCGAL
)

add_executable(WhereUsedStudy
  WhereUsedStudy.cpp
)

target_link_libraries(WhereUsedStudy
  Logging
  Geometry
  Triangulate
  InterfaceCGAL
)

# add testing subdirectory
add_subdirectory(${CMAKE_SOURCE_DIR}/test)

cmake_minimum_required (VERSION 3.0)
cmake_policy(VERSION 3.0)

project (Delaunay C CXX)

# prohibit in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif()


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

find_program(CCACHE_PROG ccache)

if (NOT CCACHE_PROG MATCHES NOTFOUND)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif ()

enable_testing()
find_package(GTest)

# adds the project-specific cmake module directory cmake/Modules to the cmake
# search path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")

find_package(CGAL REQUIRED)
find_package(Cairomm REQUIRED)
find_package(TBB REQUIRED)
find_package(VTK REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem serialization program_options)
find_package(Threads REQUIRED)

include_directories("${PROJECT_SOURCE_DIR}/include")

include_directories(
  SYSTEM
  ${CGAL_INCLUDE_DIRS}
  ${Cairomm_INCLUDE_DIRS}
)

include(${VTK_USE_FILE})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -std=c++1y -pedantic -Wall -Wextra -Werror")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -std=c11 -pedantic -Wall -Wextra -Werror")

option(VALGRIND "Enable Valgrind compatibility options" OFF)
if (VALGRIND)
  message(STATUS "Valgrind mode enabled")
  add_definitions(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)
endif ()

option(SANITIZE "Enable thread sanitizer" OFF)
if (SANITIZE)
  message(WARNING "ThreadSanitizer is enabled")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=thread")
endif()

# add src subdirectory
add_subdirectory(src)

# add testing subdirectory
add_subdirectory(tests)

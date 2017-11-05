cmake_minimum_required (VERSION 3.0)
cmake_policy(VERSION 3.0)

project (Delaunay C CXX)

# prohibit in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif()


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")

find_program(CCACHE_PROG ccache)

if (NOT CCACHE_PROG MATCHES NOTFOUND)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif ()

find_package(GTest)
if(GTest_FOUND)
	enable_testing()
endif()

# adds the project-specific cmake module directory cmake/Modules to the cmake
# search path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")

set(Boost_NO_BOOST_CMAKE ON)
set(Boost_REALPATH ON)
find_package(Boost 1.58 REQUIRED COMPONENTS system filesystem serialization program_options)

find_package(CGAL REQUIRED)
find_package(TBB REQUIRED)

find_package(Threads REQUIRED)
find_package(EJDB REQUIRED)
find_package(VTune)
find_package(KaHIP REQUIRED)

include_directories("${PROJECT_SOURCE_DIR}/include")

include_directories(
  SYSTEM
  ${CGAL_INCLUDE_DIRS}
)

include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -std=c++17 -pedantic -Wall -Wextra -Werror")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -std=c11 -pedantic -Wall -Wextra -Werror")

option(VALGRIND "Enable Valgrind compatibility options" OFF)
if (VALGRIND)
  message(STATUS "Valgrind mode enabled")
  add_definitions(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)
endif ()

option(TIMER "Enable fine grained TaskTimer" OFF)
if (TIMER)
    message(STATUS "Timing enabled")
    add_definitions(-DENABLE_TASK_TIMER)
endif ()

option(SANITIZE_THREAD "Enable thread sanitizer" OFF)
if (SANITIZE_THREAD)
  message(WARNING "ThreadSanitizer is enabled")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=thread")
endif()

option(SANITIZE_ADDRESS "Enable address sanitizer" OFF)
if (SANITIZE_ADDRESS)
  message(WARNING "AdressSanitizer is enabled")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address")
endif()

option(PROFILING "Enable fine-grained profiling" OFF)
if(PROFILING)
  message(STATUS "Profiling enabled")
  add_definitions(-DENABLE_PROFILING)
endif()

option(VTUNE "Enable instrumentation for VTUNE" OFF)
if(VTUNE)
    if(VTUNE_FOUND)
        message(STATUS "VTUNE instrumentation enabled")
        include_directories(
          SYSTEM
          ${VTUNE_INCLUDE_DIRS}
        )
        add_definitions(-DENABLE_VTUNE)
    else()
        message(ERROR "VTUNE requested but not found")
    endif()
endif()

# add src subdirectory
add_subdirectory(src)

if(GTest_FOUND)
	# add testing subdirectory
	add_subdirectory(tests)
endif()

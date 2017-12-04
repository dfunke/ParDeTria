# - Try to find libkdtree++
# Once done, this will define
#
#  LibKDTree_FOUND - system has KaHIP
#  --KaHIP_INCLUDE_DIRS - the KaHIP include directories
#  EJDB_LIBRARIES - link these to use KaHIP

include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(LibKDTree_PKGCONF libkahip)

# Include dir
find_path(LibKDTree_INCLUDE_DIR
  NAMES kdtree++/kdtree.hpp
  PATHS ${LibKDTree_PKGCONF_INCLUDE_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(KaHIP_PROCESS_INCLUDES LibKDTree_INCLUDE_DIR)
libfind_process(LibKDTree)


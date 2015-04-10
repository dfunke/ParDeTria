# CMake script for src directory
# also deals with utils subdir

include_directories("${PROJECT_SOURCE_DIR}/include/utils")
configure_file("version.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/version.cpp" @ONLY)

add_library (Logging
  Logger.cpp
  ccolor.cpp
  StringUtils.cxx
)

target_link_libraries(Logging ${TBB_LIBRARIES})

add_library (Measurements
  Timings.cpp
  System.c
  System.cpp
  ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
)

add_executable (VTUShow
  VTU_Show.cpp
)

target_link_libraries(VTUShow ${VTK_LIBRARIES} ${Boost_LIBRARIES})

add_executable (archiveConverter
  ArchiveConverter.cpp
)

target_link_libraries(archiveConverter ${Boost_LIBRARIES} ${TBB_LIBRARIES})

add_library (DBConn
  DBConnection.cpp
  ${EJDB_LIBRARIES}
)


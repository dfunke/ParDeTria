# CMake script for src directory
# also deals with utils subdir

include_directories("${PROJECT_SOURCE_DIR}/include/utils")
configure_file("version.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/version.cpp" @ONLY)

if (VTUNE)
    add_library(VTune
            VTuneAdapter.cpp
            )
    target_link_libraries(VTune
            ${VTUNE_LIBRARIES}
            ${CMAKE_THREAD_LIBS_INIT}
            )
endif ()

add_library(Logging EXCLUDE_FROM_ALL
        Logger.cpp
        ccolor.cpp
        StringUtils.cxx
        )

target_link_libraries(Logging ${TBB_LIBRARIES})

add_library(Measurements EXCLUDE_FROM_ALL
        Timings.cpp
        System.c
        System.cpp
        ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
        )

if (VTUNE)
    target_link_libraries(Measurements
            VTune
            )
endif ()

target_link_libraries(Measurements ${TBB_LIBRARIES})

add_library(DBConn EXCLUDE_FROM_ALL
        DBConnection.cpp
        ${EJDB_LIBRARIES}
        )

add_library(Predicates EXCLUDE_FROM_ALL
        Predicates.cpp
        )


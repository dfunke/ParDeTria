file(GLOB tests_SRC
  "Test*.cpp"
)

add_executable(runTests
  ${tests_SRC}
)

target_link_libraries(runTests
  Geometry Logging Measurements SpatialSort
  ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

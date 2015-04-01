add_executable(runTests
  TestGeometry3D.cpp
)

target_link_libraries(runTests
  Geometry Logging
  ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES})

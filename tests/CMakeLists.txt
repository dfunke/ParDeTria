add_executable(runTests
  TestGeometry3D.cpp
  TestExperimentRun.cpp
)

target_link_libraries(runTests
  Geometry Logging Measurements
  ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES})

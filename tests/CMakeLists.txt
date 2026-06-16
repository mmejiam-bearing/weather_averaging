add_executable(weather_unit_tests
  test_averaging.cpp
  test_io.cpp
  ../src/averaging.cpp
  ../src/io.cpp
)

target_include_directories(weather_unit_tests PRIVATE
  ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(weather_unit_tests PRIVATE
  weather_common
  GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(weather_unit_tests)

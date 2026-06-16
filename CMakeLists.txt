cmake_minimum_required(VERSION 3.20)

project(
  weather_daily_avg
  VERSION 1.0.0
  DESCRIPTION "Energy-based daily averaging of GFS-Wave fields for maritime routing"
  LANGUAGES CXX
)

# ─────────────────────────────────────────────
# C++ Standard
# ─────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # No GNU extensions; strict ISO C++23 only

# ─────────────────────────────────────────────
# Build type defaults
# ─────────────────────────────────────────────
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type specified — defaulting to Release")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
endif()

# ─────────────────────────────────────────────
# Compiler flags
# ─────────────────────────────────────────────
# Optimization and vectorization
set(RELEASE_FLAGS "-O3 -funroll-loops")
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # Enable native ISA tuning only if the user hasn't already set it
  if(NOT CMAKE_CXX_FLAGS MATCHES "march")
    string(APPEND RELEASE_FLAGS " -march=native")
  endif()
endif()

# Strict warnings — all new code must compile cleanly under these
set(WARNING_FLAGS "-Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wnon-virtual-dtor")

# Clang-specific: disable -Wambiguous-reversed-operator warning from Eigen internals
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  string(APPEND WARNING_FLAGS " -Wno-ambiguous-reversed-operator")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${RELEASE_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   -g -fsanitize=address,undefined")
add_compile_options("$<$<CONFIG:Debug>:-DDEBUG>")

# ─────────────────────────────────────────────
# Dependencies
# ─────────────────────────────────────────────

# --- OpenMP ---
find_package(OpenMP REQUIRED)
if(NOT OpenMP_CXX_FOUND)
  message(FATAL_ERROR "OpenMP not found. Install libomp-dev (apt) or llvm-openmp (brew).")
endif()

# --- Eigen3 ---
# Allow user to override via -DEIGEN3_INCLUDE_DIR=/path/to/eigen
if(DEFINED EIGEN3_INCLUDE_DIR)
  message(STATUS "Using user-provided Eigen3 include dir: ${EIGEN3_INCLUDE_DIR}")
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}"
  )
else()
  find_package(Eigen3 3.4 REQUIRED NO_MODULE)
  message(STATUS "Found Eigen3: ${EIGEN3_VERSION} at ${EIGEN3_INCLUDE_DIRS}")
endif()

# --- zlib (required by cnpy) ---
find_package(ZLIB REQUIRED)

# --- cnpy (git submodule at ./cnpy) ---
if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cnpy/cnpy.h")
  message(FATAL_ERROR
    "cnpy submodule not found at ${CMAKE_CURRENT_SOURCE_DIR}/cnpy/cnpy.h\n"
    "Run: git submodule update --init --recursive"
  )
endif()

# --- GoogleTest (fetched once, used by tests/ only) ---
option(BUILD_TESTS "Build the test suite (requires network access to fetch GoogleTest)" ON)

if(BUILD_TESTS)
  include(FetchContent)
  set(GTEST_COMMIT "v1.14.0")
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        ${GTEST_COMMIT}
    GIT_SHALLOW    TRUE
  )
  # Prevent GTest from installing alongside the project
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
endif()

add_library(cnpy STATIC cnpy/cnpy.cpp)
target_include_directories(cnpy PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/cnpy")
target_link_libraries(cnpy PUBLIC ZLIB::ZLIB)
set_target_properties(cnpy PROPERTIES
  CXX_STANDARD 23
  CXX_STANDARD_REQUIRED ON
)

# ─────────────────────────────────────────────
# Common interface library (shared flags + deps)
# ─────────────────────────────────────────────
add_library(weather_common INTERFACE)
target_include_directories(weather_common INTERFACE
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
)
target_compile_options(weather_common INTERFACE
  "$<$<CXX_COMPILER_ID:GNU,Clang>:${WARNING_FLAGS}>"
)
target_link_libraries(weather_common INTERFACE
  Eigen3::Eigen
  OpenMP::OpenMP_CXX
  cnpy
  ZLIB::ZLIB
)
# Suppress Eigen internal deprecation warnings in release builds
target_compile_definitions(weather_common INTERFACE
  EIGEN_NO_DEBUG              # disables Eigen assertions in Release
  $<$<CONFIG:Release>:NDEBUG>
)

# ─────────────────────────────────────────────
# Main executable
# ─────────────────────────────────────────────
add_executable(weather_daily_avg
  src/main.cpp
  src/averaging.cpp
  src/io.cpp
)
target_link_libraries(weather_daily_avg PRIVATE weather_common)

# Export compile commands for clangd / IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ─────────────────────────────────────────────
# Tests
# ─────────────────────────────────────────────
if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
  message(STATUS "Tests enabled (GoogleTest ${GTEST_COMMIT})")
endif()

# ─────────────────────────────────────────────
# Install rules
# ─────────────────────────────────────────────
include(GNUInstallDirs)
install(TARGETS weather_daily_avg
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

# ─────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────
message(STATUS "--------------------------------------------")
message(STATUS "  Project      : ${PROJECT_NAME} v${PROJECT_VERSION}")
message(STATUS "  Build type   : ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ standard : C++${CMAKE_CXX_STANDARD}")
message(STATUS "  Compiler     : ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  OpenMP       : ${OpenMP_CXX_VERSION}")
message(STATUS "  Build tests  : ${BUILD_TESTS}")
message(STATUS "--------------------------------------------")

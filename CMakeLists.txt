cmake_minimum_required(VERSION 3.20)
project(chebyshev_solver LANGUAGES CXX VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(Sanitizers)

# ---- Núcleo header-only ----
add_library(numa INTERFACE)
target_include_directories(numa INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_compile_features(numa INTERFACE cxx_std_20)

find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
  target_link_libraries(numa INTERFACE OpenMP::OpenMP_CXX)
  target_compile_definitions(numa INTERFACE NUMA_HAS_OPENMP=1)
  message(STATUS "OpenMP habilitado")
endif()

option(NUMA_BUILD_TESTS "Build tests" ON)
option(NUMA_BUILD_BENCH "Build benchmarks" ON)
option(NUMA_BUILD_APPS  "Build demos" ON)
option(NUMA_BUILD_CUDA  "Build CUDA demos" OFF)
option(NUMA_ENABLE_SANITIZERS "ASan/UBSan" OFF)

if(NUMA_ENABLE_SANITIZERS)
  enable_sanitizers(numa)
endif()

# ---- Tests ----
if(NUMA_BUILD_TESTS)
  enable_testing()
  include(FetchContent)
  FetchContent_Declare(catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.2
    GIT_SHALLOW TRUE)
  FetchContent_MakeAvailable(catch2)

  add_executable(numa_tests
    tests/test_core.cpp
    tests/test_linalg.cpp
    tests/test_recurrence.cpp
    tests/test_methods.cpp
    tests/test_golden.cpp)
  target_link_libraries(numa_tests PRIVATE numa Catch2::Catch2WithMain)
  list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
  include(Catch)
  catch_discover_tests(numa_tests)
endif()

# ---- Benchmark ----
if(NUMA_BUILD_BENCH)
  include(FetchContent)
  FetchContent_Declare(benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.8.3
    GIT_SHALLOW TRUE)
  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(benchmark)

  add_executable(numa_bench apps/bench_all.cpp)
  target_link_libraries(numa_bench PRIVATE numa benchmark::benchmark)
endif()

# ---- Demos ----
if(NUMA_BUILD_APPS)
  foreach(app demo_newton demo_broyden demo_kantorovich
              demo_chebyshev demo_chebyshev_omp)
    add_executable(${app} apps/${app}.cpp)
    target_link_libraries(${app} PRIVATE numa)
  endforeach()
endif()

# ---- CUDA ----
if(NUMA_BUILD_CUDA)
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  if(NOT CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 75 80 86 89)
  endif()

  add_executable(demo_chebyshev_gpu apps/demo_chebyshev_gpu.cu)
  target_include_directories(demo_chebyshev_gpu PRIVATE include)
  set_target_properties(demo_chebyshev_gpu PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()

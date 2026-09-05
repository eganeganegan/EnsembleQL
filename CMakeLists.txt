cmake_minimum_required(VERSION 3.20)
project(EnsembleQL VERSION 0.1.0 LANGUAGES CXX)

option(ENSEMBLEQL_BUILD_PYTHON "Build the pybind11 extension" ON)
option(ENSEMBLEQL_BUILD_TESTS "Build C++ tests" ON)
option(ENSEMBLEQL_BUILD_BENCHMARKS "Build microbenchmarks" OFF)
option(ENSEMBLEQL_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

add_library(ensembleql_core
    cpp/src/units.cpp
    cpp/src/topology.cpp
    cpp/src/trajectory.cpp
    cpp/src/selection.cpp
    cpp/src/geometry.cpp
    cpp/src/observable.cpp
    cpp/src/event.cpp
    cpp/src/temporal.cpp
    cpp/src/parser.cpp
    cpp/src/planner.cpp
    cpp/src/engine.cpp)
add_library(EnsembleQL::core ALIAS ensembleql_core)
set_target_properties(ensembleql_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(ensembleql_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cpp/include>
    $<INSTALL_INTERFACE:include>)
target_compile_features(ensembleql_core PUBLIC cxx_std_20)
if(MSVC)
  target_compile_options(ensembleql_core PRIVATE /W4)
  if(ENSEMBLEQL_WARNINGS_AS_ERRORS)
    target_compile_options(ensembleql_core PRIVATE /WX)
  endif()
else()
  target_compile_options(ensembleql_core PRIVATE -Wall -Wextra -Wpedantic)
  if(ENSEMBLEQL_WARNINGS_AS_ERRORS)
    target_compile_options(ensembleql_core PRIVATE -Werror)
  endif()
endif()

if(ENSEMBLEQL_BUILD_PYTHON)
  find_package(pybind11 CONFIG QUIET)
  if(pybind11_FOUND)
    pybind11_add_module(_core cpp/src/bindings.cpp)
    target_link_libraries(_core PRIVATE EnsembleQL::core)
    install(TARGETS _core DESTINATION ensembleql)
  else()
    message(STATUS "pybind11 not found; skipping Python extension")
  endif()
endif()

if(ENSEMBLEQL_BUILD_TESTS)
  enable_testing()
  add_executable(ensembleql_cpp_tests tests/cpp/test_core.cpp)
  target_link_libraries(ensembleql_cpp_tests PRIVATE EnsembleQL::core)
  target_compile_definitions(ensembleql_cpp_tests PRIVATE ENSEMBLEQL_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
  target_compile_features(ensembleql_cpp_tests PRIVATE cxx_std_20)
  add_test(NAME ensembleql_cpp_tests COMMAND ensembleql_cpp_tests)
  find_package(Catch2 3 QUIET)
  if(Catch2_FOUND)
    add_executable(ensembleql_catch2_tests tests/cpp/test_catch2.cpp)
    target_link_libraries(ensembleql_catch2_tests PRIVATE EnsembleQL::core Catch2::Catch2WithMain)
    target_compile_features(ensembleql_catch2_tests PRIVATE cxx_std_20)
    include(Catch)
    catch_discover_tests(ensembleql_catch2_tests)
  else()
    message(STATUS "Catch2 3 not found; comprehensive dependency-free C++ tests remain enabled")
  endif()
endif()

if(ENSEMBLEQL_BUILD_BENCHMARKS)
  add_executable(ensembleql_benchmarks benchmarks/benchmarks.cpp)
  target_link_libraries(ensembleql_benchmarks PRIVATE EnsembleQL::core)
endif()

install(TARGETS ensembleql_core ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(DIRECTORY cpp/include/ DESTINATION include)

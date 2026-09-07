cmake_minimum_required(VERSION 3.20)
project(EnsembleQL VERSION 0.1.0 LANGUAGES CXX)

option(ENSEMBLEQL_BUILD_PYTHON "Build the pybind11 extension" ON)
option(ENSEMBLEQL_BUILD_TESTS "Build C++ tests" ON)
option(ENSEMBLEQL_BUILD_BENCHMARKS "Build microbenchmarks" OFF)
option(ENSEMBLEQL_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(ENSEMBLEQL_ENABLE_OPENMP "Enable OpenMP acceleration when available" ON)
option(ENSEMBLEQL_ENABLE_CHEMFILES "Enable chemfiles trajectory backends when available" ON)
option(ENSEMBLEQL_FETCH_CHEMFILES "Fetch chemfiles 0.10.4 when it is not installed" OFF)
option(ENSEMBLEQL_REQUIRE_CHEMFILES "Fail configuration when chemfiles is unavailable" OFF)

# The native core and any fetched static dependencies can be linked into the
# Python extension. Set this before creating targets (including FetchContent
# targets) so Linux builds compile every contributing object as PIC.
if(ENSEMBLEQL_BUILD_PYTHON)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

add_library(ensembleql_core
    cpp/src/units.cpp
    cpp/src/topology.cpp
    cpp/src/trajectory.cpp
    cpp/src/pdb_reader.cpp
    cpp/src/selection.cpp
    cpp/src/geometry.cpp
    cpp/src/molecular.cpp
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
if(ENSEMBLEQL_ENABLE_OPENMP)
  find_package(OpenMP QUIET COMPONENTS CXX)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(ensembleql_core PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(ensembleql_core PUBLIC ENSEMBLEQL_HAS_OPENMP=1)
    message(STATUS "OpenMP acceleration enabled")
  else()
    message(STATUS "OpenMP not found; using serial kernels")
  endif()
endif()
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

if(ENSEMBLEQL_ENABLE_CHEMFILES)
  find_package(chemfiles 0.10 CONFIG QUIET)
  if(NOT chemfiles_FOUND AND ENSEMBLEQL_FETCH_CHEMFILES)
    include(FetchContent)
    # chemfiles 0.10.4 predates CMake 4's removal of compatibility below 3.5.
    set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
    set(CHFL_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(CHFL_BUILD_DOCUMENTATION OFF CACHE BOOL "" FORCE)
    if(APPLE)
      # chemfiles 0.10.4's bundled zlib assumes large-file APIs that recent
      # Apple SDKs no longer expose. Linux keeps the upstream bundled default.
      set(CHFL_SYSTEM_ZLIB ON CACHE BOOL "" FORCE)
    else()
      set(CHFL_SYSTEM_ZLIB OFF CACHE BOOL "" FORCE)
    endif()
    FetchContent_Declare(chemfiles_dependency
      GIT_REPOSITORY https://github.com/chemfiles/chemfiles.git
      GIT_TAG 0.10.4
      GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(chemfiles_dependency)
    if(APPLE AND TARGET chemfiles_objects)
      # macOS uses 64-bit off_t without Linux's off64_t compatibility type.
      get_target_property(_chemfiles_definitions chemfiles_objects COMPILE_DEFINITIONS)
      list(FILTER _chemfiles_definitions EXCLUDE REGEX "^_LARGEFILE64_SOURCE")
      list(FILTER _chemfiles_definitions EXCLUDE REGEX "^_LFS64_LARGEFILE")
      set_target_properties(chemfiles_objects PROPERTIES COMPILE_DEFINITIONS "${_chemfiles_definitions}")
      target_compile_definitions(chemfiles_objects PRIVATE gzopen64=gzopen gzseek64=gzseek)
    endif()
  endif()

  if(chemfiles_FOUND OR TARGET chemfiles OR TARGET chemfiles::chemfiles)
    set(ENSEMBLEQL_HAS_CHEMFILES ON)
    target_sources(ensembleql_core PRIVATE cpp/src/chemfiles_reader.cpp)
    target_compile_definitions(ensembleql_core PUBLIC ENSEMBLEQL_HAS_CHEMFILES=1)
    if(TARGET chemfiles::chemfiles)
      target_link_libraries(ensembleql_core PRIVATE chemfiles::chemfiles)
    else()
      target_link_libraries(ensembleql_core PRIVATE chemfiles)
    endif()
    message(STATUS "chemfiles backend enabled (extended trajectory formats)")
  elseif(ENSEMBLEQL_REQUIRE_CHEMFILES)
    message(FATAL_ERROR "chemfiles >= 0.10 is required but was not found")
  else()
    message(STATUS "chemfiles not found; extended trajectory formats disabled")
  endif()
elseif(ENSEMBLEQL_REQUIRE_CHEMFILES)
  message(FATAL_ERROR "ENSEMBLEQL_REQUIRE_CHEMFILES requires ENSEMBLEQL_ENABLE_CHEMFILES=ON")
endif()

if(ENSEMBLEQL_BUILD_PYTHON)
  set(PYBIND11_FINDPYTHON ON)
  find_package(pybind11 CONFIG QUIET)
  if(pybind11_FOUND)
    pybind11_add_module(_core cpp/src/bindings.cpp)
    target_link_libraries(_core PRIVATE EnsembleQL::core)
    install(TARGETS _core DESTINATION ensembleql COMPONENT python)
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
  if(ENSEMBLEQL_HAS_CHEMFILES)
    add_executable(ensembleql_chemfiles_tests tests/cpp/test_chemfiles.cpp)
    target_link_libraries(ensembleql_chemfiles_tests PRIVATE EnsembleQL::core)
    if(TARGET chemfiles::chemfiles)
      target_link_libraries(ensembleql_chemfiles_tests PRIVATE chemfiles::chemfiles)
    else()
      target_link_libraries(ensembleql_chemfiles_tests PRIVATE chemfiles)
    endif()
    target_compile_features(ensembleql_chemfiles_tests PRIVATE cxx_std_20)
    target_compile_definitions(ensembleql_chemfiles_tests PRIVATE
      ENSEMBLEQL_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
      ENSEMBLEQL_BINARY_DIR="${CMAKE_CURRENT_BINARY_DIR}")
    add_test(NAME ensembleql_chemfiles_tests COMMAND ensembleql_chemfiles_tests)
  endif()
endif()

if(ENSEMBLEQL_BUILD_BENCHMARKS)
  add_executable(ensembleql_benchmarks benchmarks/benchmarks.cpp)
  target_link_libraries(ensembleql_benchmarks PRIVATE EnsembleQL::core)
endif()

install(TARGETS ensembleql_core ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(DIRECTORY cpp/include/ DESTINATION include)

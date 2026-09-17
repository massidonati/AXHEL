# Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.13)
project(AXHEL 
        VERSION 1.0.0 
        DESCRIPTION "ARM aXeleration for Homomorphic Encryption Library"
        LANGUAGES C CXX)

#------------------------------------
# CMake modules
#------------------------------------

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)
include(CMakePackageConfigHelpers)

set(AXHEL_CMAKE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/axhel")
list(APPEND CMAKE_MODULE_PATH "${AXHEL_CMAKE_PATH}")
include(axhel-util)

#------------------------------------
# Build types
#------------------------------------

if(CMAKE_BUILD_TYPE)
  set(AXHEL_SUPPORTED_BUILD_TYPES
      Debug
      Release)
  list(FIND AXHEL_SUPPORTED_BUILD_TYPES ${CMAKE_BUILD_TYPE} INDEX_FOUND)
  if(${INDEX_FOUND} EQUAL -1)
    message(
      FATAL_ERROR
        "CMAKE_BUILD_TYPE must be one of Debug, Release"
      )
  endif()
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(AXHEL_DEBUG ON)
else()
  set(AXHEL_DEBUG OFF)
endif()

#------------------------------------
# C++ configuration
#------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_INSTALL_MESSAGE LAZY)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Create compilation database compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Allow installed shared libraries to find dependencies in the same directory
set(CMAKE_INSTALL_RPATH "\$ORIGIN")

#------------------------------------
# AXHEL options
#------------------------------------

option(AXHEL_BENCHMARK "Enable benchmarking" OFF)
option(AXHEL_BENCHMARK_FORCE_NATIVE "Force native backend only in AXHEL benchmark builds" OFF)
option(AXHEL_BENCHMARK_FORCE_SCALAR "Force strict scalar backend (no compiler auto-vectorization) in AXHEL benchmark builds" OFF)
#option(AXHEL_DOCS "Enable documentation building" OFF)
option(AXHEL_SHARED_LIB "Generate a shared library" OFF)
option(AXHEL_TESTING "Enable AXHEL correctness and regression tests" OFF)
option(AXHEL_TEST_FETCH_GTEST "Fetch GoogleTest if not installed" OFF)
option(AXHEL_TEST_FORCE_NATIVE "Force native backend only in AXHEL test builds" OFF)
option(AXHEL_TREAT_WARNING_AS_ERROR "Treat all compile-time warnings as errors" OFF)
option(AXHEL_OPT_REPORT "Enable GCC/Clang optimization reports" OFF)
set(AXHEL_CPU "native" CACHE STRING "ARM CPU target passed to the compiler through -mcpu")

check_cxx_compiler_flag("-mcpu=${AXHEL_CPU}" AXHEL_CPU_FLAG_SUPPORTED)

if(NOT AXHEL_CPU_FLAG_SUPPORTED)
    message(
        FATAL_ERROR
        "The compiler does not support -mcpu=${AXHEL_CPU}. "
        "Use a compiler-supported ARM CPU name or AXHEL_CPU=native."
    )
endif()

#------------------------------------
# Configuration summary
#------------------------------------

message(STATUS "CMAKE_BUILD_TYPE:              ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_C_COMPILER:              ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER:            ${CMAKE_CXX_COMPILER}")
message(STATUS "AXHEL_CPU:                     ${AXHEL_CPU}")
message(STATUS "AXHEL_BENCHMARK:               ${AXHEL_BENCHMARK}")
message(STATUS "AXHEL_BENCHMARK_FORCE_NATIVE:  ${AXHEL_BENCHMARK_FORCE_NATIVE}")
message(STATUS "AXHEL_BENCHMARK_FORCE_SCALAR:  ${AXHEL_BENCHMARK_FORCE_SCALAR}")
message(STATUS "AXHEL_DEBUG:                   ${AXHEL_DEBUG}")
#message(STATUS "AXHEL_DOCS:                    ${AXHEL_DOCS}")
message(STATUS "AXHEL_SHARED_LIB:              ${AXHEL_SHARED_LIB}")
message(STATUS "AXHEL_TESTING:                 ${AXHEL_TESTING}")
message(STATUS "AXHEL_TEST_FORCE_NATIVE:       ${AXHEL_TEST_FORCE_NATIVE}")
message(STATUS "AXHEL_TREAT_WARNING_AS_ERROR:  ${AXHEL_TREAT_WARNING_AS_ERROR}")
message(STATUS "AXHEL_OPT_REPORT:              ${AXHEL_OPT_REPORT}")

#------------------------------------
# Compiler detection and configuration
#------------------------------------

axhel_check_compiler_version()
axhel_add_compiler_definition()

if (AXHEL_TREAT_WARNING_AS_ERROR)
  add_compile_options(-Werror)
endif()

#------------------------------------
# Project directories
#------------------------------------

set(AXHEL_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})
set(AXHEL_SRC_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/axhel)
set(AXHEL_INC_ROOT_DIR ${AXHEL_SRC_ROOT_DIR}/include) # Public headers

message(STATUS "CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}")

#------------------------------------
# Detect SVE support
#------------------------------------

# Always perform the normal SVE capability detection. The test/benchmark
# force-Native/strict-scalar switches only change the effective backend in
# their dedicated builds; neither bypasses hardware/compiler detection.
axhel_check_compile_flag("${AXHEL_CMAKE_PATH}/test-sve.cpp" AXHEL_HAS_SVE)
set(AXHEL_SVE_DETECTED ${AXHEL_HAS_SVE})

if(AXHEL_TEST_FORCE_NATIVE)
    if(NOT AXHEL_TESTING)
        message(FATAL_ERROR "AXHEL_TEST_FORCE_NATIVE is valid only with AXHEL_TESTING=ON")
    endif()
    set(AXHEL_HAS_SVE 0)
endif()

if(AXHEL_BENCHMARK_FORCE_NATIVE AND AXHEL_BENCHMARK_FORCE_SCALAR)
    message(FATAL_ERROR "AXHEL_BENCHMARK_FORCE_NATIVE and AXHEL_BENCHMARK_FORCE_SCALAR are mutually exclusive")
endif()

if(AXHEL_BENCHMARK_FORCE_NATIVE)
    if(NOT AXHEL_BENCHMARK)
        message(FATAL_ERROR "AXHEL_BENCHMARK_FORCE_NATIVE is valid only with AXHEL_BENCHMARK=ON")
    endif()
    set(AXHEL_HAS_SVE 0)
endif()

if(AXHEL_BENCHMARK_FORCE_SCALAR)
    if(NOT AXHEL_BENCHMARK)
        message(FATAL_ERROR "AXHEL_BENCHMARK_FORCE_SCALAR is valid only with AXHEL_BENCHMARK=ON")
    endif()
    # Strict scalar is a distinct benchmark configuration: exclude SVE sources
    # and dispatch through the public Native path. Auto-vectorization itself is
    # disabled on the AXHEL target in axhel/CMakeLists.txt.
    set(AXHEL_HAS_SVE 0)
endif()

message(STATUS "AXHEL_SVE_DETECTED: ${AXHEL_SVE_DETECTED}")
message(STATUS "AXHEL_HAS_SVE (effective backend): ${AXHEL_HAS_SVE}")

#------------------------------------
# Installation configuration
#------------------------------------

if(NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR})
endif()

# Defines CMAKE_INSTALL_LIBDIR, CMAKE_INSTALL_INCLUDEDIR, CMAKE_INSTALL_BINDIR, etc.
include(GNUInstallDirs)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_INSTALL_LIBDIR})

message(STATUS "CMAKE_INSTALL_LIBDIR:     ${CMAKE_INSTALL_LIBDIR}")
message(STATUS "CMAKE_INSTALL_INCLUDEDIR: ${CMAKE_INSTALL_INCLUDEDIR}")
message(STATUS "CMAKE_INSTALL_PREFIX:     ${CMAKE_INSTALL_PREFIX}")

#------------------------------------
# AXHEL library
#------------------------------------

add_subdirectory(axhel)

if(AXHEL_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

if(AXHEL_BENCHMARK)
    add_subdirectory(benchmark)
endif()

#------------------------------------
# Optional components
#------------------------------------

#if (AXHEL_DOCS)
#  add_subdirectory(docs)
#endif()

#------------------------------------
# Final check
#------------------------------------

add_custom_target(check COMMAND pre-commit install && pre-commit run --all-files)


# Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.13)
project(AXHEL 
        VERSION 0.1 
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
#option(AXHEL_DOCS "Enable documentation building" OFF)
option(AXHEL_SHARED_LIB "Generate a shared library" OFF)
option(AXHEL_TESTING "Enables unit-tests" OFF)
option(AXHEL_TREAT_WARNING_AS_ERROR "Treat all compile-time warnings as errors" OFF)

#------------------------------------
# Configuration summary
#------------------------------------

message(STATUS "CMAKE_BUILD_TYPE:              ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_C_COMPILER:              ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER:            ${CMAKE_CXX_COMPILER}")
message(STATUS "AXHEL_BENCHMARK:               ${AXHEL_BENCHMARK}")
message(STATUS "AXHEL_DEBUG:                   ${AXHEL_DEBUG}")
#message(STATUS "AXHEL_DOCS:                    ${AXHEL_DOCS}")
message(STATUS "AXHEL_SHARED_LIB:              ${AXHEL_SHARED_LIB}")
message(STATUS "AXHEL_TESTING:                 ${AXHEL_TESTING}")
message(STATUS "AXHEL_TREAT_WARNING_AS_ERROR:  ${AXHEL_TREAT_WARNING_AS_ERROR}")

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

axhel_check_compile_flag("${AXHEL_CMAKE_PATH}/test-sve.cpp" AXHEL_HAS_SVE)

message(STATUS "AXHEL_HAS_SVE: ${AXHEL_HAS_SVE}")

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


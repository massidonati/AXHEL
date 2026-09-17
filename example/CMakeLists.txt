# Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
# SPDX-License-Identifier: Apache-2.0


cmake_minimum_required(VERSION 3.13)
project(axhel_example LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build AXHEL directly from the repository source tree.
add_subdirectory(
    "${CMAKE_CURRENT_LIST_DIR}/.."
    "${CMAKE_BINARY_DIR}/axhel"
)

if (NOT TARGET AXHEL::axhel)
    message(FATAL_ERROR "TARGET AXHEL::axhel not found")
endif()

add_executable(example example.cpp)
target_link_libraries(example PRIVATE AXHEL::axhel)
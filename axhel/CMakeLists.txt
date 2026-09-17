# Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
# SPDX-License-Identifier: Apache-2.0

#---------------------------------
# Source files
#---------------------------------

set(NATIVE_SRC
    eltwise/eltwise-add-mod.cpp
    eltwise/eltwise-sub-mod.cpp
    eltwise/eltwise-mul-mod.cpp
    eltwise/eltwise-fma-mod.cpp
    eltwise/eltwise-reduce-mod.cpp
    number-theory/multiply-factor.cpp
    ntt/ntt.cpp
    ntt/ntt-native.cpp
    ntt/ntt-tables.cpp
)

set(SVE_SRC)
if (AXHEL_HAS_SVE)
    set(SVE_SRC
        eltwise/eltwise-add-mod-sve.cpp
        eltwise/eltwise-sub-mod-sve.cpp
        eltwise/eltwise-mul-mod-sve.cpp
        eltwise/eltwise-fma-mod-sve.cpp
        eltwise/eltwise-reduce-mod-sve.cpp
        ntt/ntt-sve.cpp
    )
endif()

set(AXHEL_SRC 
    ${NATIVE_SRC}
    ${SVE_SRC}
)

#---------------------------------
# Library target
#---------------------------------

if (AXHEL_SHARED_LIB)
    add_library(axhel SHARED ${AXHEL_SRC})
else()
     add_library(axhel STATIC ${AXHEL_SRC})
endif()

add_library(AXHEL::axhel ALIAS axhel)

# Test-only hooks are compiled into the library only for AXHEL test builds.
# PRIVATE keeps AXHEL_TESTING out of the installed/public target interface.
if(AXHEL_TESTING)
    target_compile_definitions(axhel PRIVATE AXHEL_TESTING=1)
endif()

axhel_add_asan_flag(axhel)

# GCC 11 may trigger an internal compiler error when AddressSanitizer
# instruments SVE intrinsic-heavy template kernels.
#
# Keep ASan enabled for the AXHEL target in Debug builds, but disable
# instrumentation only for the SVE translation units affected by the
# compiler bug.
if(
    AXHEL_DEBUG
    AND AXHEL_HAS_SVE
    AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
    AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12.0
)
    message(
        STATUS
        "Disabling AddressSanitizer instrumentation for SVE sources "
        "due to a GCC ${CMAKE_CXX_COMPILER_VERSION} compiler limitation"
    )

    set_source_files_properties(
        ${SVE_SRC}
        PROPERTIES
            COMPILE_OPTIONS "-fno-sanitize=address"
    )
endif()

#---------------------------------
# Generated configuration header
#---------------------------------

set(AXHEL_GENERATED_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
set(AXHEL_DEFINES_IN_FILENAME "${CMAKE_CURRENT_SOURCE_DIR}/include/axhel/util/defines.hpp.in")
set(AXHEL_DEFINES_FILENAME "${AXHEL_GENERATED_INCLUDE_DIR}/axhel/util/defines.hpp")
configure_file("${AXHEL_DEFINES_IN_FILENAME}" "${AXHEL_DEFINES_FILENAME}")

#---------------------------------
# Target properties
#---------------------------------

set_target_properties(axhel PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(axhel PROPERTIES VERSION ${AXHEL_VERSION})

if (AXHEL_DEBUG)
    set_target_properties(axhel PROPERTIES OUTPUT_NAME "axhel_debug")
else()
    set_target_properties(axhel PROPERTIES OUTPUT_NAME "axhel")
endif()

#---------------------------------
# Include directories
#---------------------------------

target_include_directories(axhel
    PRIVATE ${AXHEL_SRC_ROOT_DIR}                              # Private headers
    PUBLIC  $<BUILD_INTERFACE:${AXHEL_GENERATED_INCLUDE_DIR}>
    PUBLIC  $<BUILD_INTERFACE:${AXHEL_INC_ROOT_DIR}>           # Public headers
    PUBLIC  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>   # Public headers
)

#---------------------------------
# Compiler options
#---------------------------------

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(
        axhel 
        PRIVATE 
            -Wall 
            -Wextra
            -Wconversion 
            -Wshadow 
            -Wno-sign-conversion
            -Wno-unknown-pragmas 
            "-mcpu=${AXHEL_CPU}"
            $<$<CONFIG:Debug>:-fno-omit-frame-pointer>
            $<$<CONFIG:Release>:-fomit-frame-pointer>
            -Wno-implicit-int-conversion
    )
    # Avoid 3rd-party dependency warnings when including AXHEL as a dependency
    target_compile_options(
        axhel 
        PUBLIC
            -Wno-unknown-warning
            -Wno-unknown-warning-option
    )
    # Third benchmark baseline: preserve the normal Release optimization level
    # and CPU tuning, but explicitly disable both forms of compiler
    # auto-vectorization on the AXHEL library itself. This is intentionally
    # separate from AXHEL_BENCHMARK_FORCE_NATIVE.
    if (AXHEL_BENCHMARK_FORCE_SCALAR)
        if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            target_compile_options(
                axhel
                PRIVATE
                    -fno-tree-vectorize
                    -fno-tree-loop-vectorize
                    -fno-tree-slp-vectorize
            )
            message(STATUS "AXHEL strict-scalar benchmark: GCC loop/SLP vectorization disabled")
        elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            target_compile_options(
                axhel
                PRIVATE
                    -fno-vectorize
                    -fno-slp-vectorize
            )
            message(STATUS "AXHEL strict-scalar benchmark: Clang loop/SLP vectorization disabled")
        endif()
    endif()

    if (AXHEL_OPT_REPORT)
         target_compile_options(
            axhel
            PRIVATE
                -fopt-info-loop-optimized
                -fopt-info-loop-missed
                -fopt-info-vec-optimized
                -fopt-info-vec-missed
        )       
    endif()
else()
    message(FATAL_ERROR "Unsupported compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()

#---------------------------------
# Install public headers
#---------------------------------

install(
    DIRECTORY 
        "${AXHEL_INC_ROOT_DIR}/"
    DESTINATION 
        "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING
        PATTERN "*.hpp"
        PATTERN "*.h")

# defines.hpp is generated in the build tree and must be installed separately.
install(
    FILES
        "${AXHEL_DEFINES_FILENAME}"

    DESTINATION
        "${CMAKE_INSTALL_INCLUDEDIR}/axhel/util"
)

#---------------------------------
# Install library and export target
#---------------------------------

set(
    AXHEL_CONFIG_INSTALL_DIR
    "${CMAKE_INSTALL_LIBDIR}/cmake/axhel-${AXHEL_VERSION}"
)

install(
    TARGETS
        axhel
    EXPORT
        AXHELTargets
    ARCHIVE DESTINATION
        "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION
        "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION
        "${CMAKE_INSTALL_BINDIR}"
    INCLUDES DESTINATION
        "${CMAKE_INSTALL_INCLUDEDIR}"
)

install(
    EXPORT
        AXHELTargets
    FILE
        AXHELTargets.cmake
    NAMESPACE
        AXHEL::
    DESTINATION
        "${AXHEL_CONFIG_INSTALL_DIR}"
)

#---------------------------------
# Package configuration files
#---------------------------------

set(
    AXHEL_TARGET_FILENAME
    "${CMAKE_CURRENT_BINARY_DIR}/AXHELTargets.cmake"
)

set(
    AXHEL_CONFIG_IN_FILENAME
    "${AXHEL_CMAKE_PATH}/AXHELConfig.cmake.in"
)

set(
    AXHEL_CONFIG_FILENAME
    "${CMAKE_CURRENT_BINARY_DIR}/AXHELConfig.cmake"
)

set(
    AXHEL_CONFIG_VERSION_FILENAME
    "${CMAKE_CURRENT_BINARY_DIR}/AXHELConfigVersion.cmake"
)

configure_package_config_file(
    "${AXHEL_CONFIG_IN_FILENAME}"
    "${AXHEL_CONFIG_FILENAME}"

    INSTALL_DESTINATION
        "${AXHEL_CONFIG_INSTALL_DIR}"
)

write_basic_package_version_file(
    "${AXHEL_CONFIG_VERSION_FILENAME}"
    VERSION "${AXHEL_VERSION}"
    COMPATIBILITY ExactVersion
)

install(
    FILES
        "${AXHEL_CONFIG_FILENAME}"
        "${AXHEL_CONFIG_VERSION_FILENAME}"
    DESTINATION
        "${AXHEL_CONFIG_INSTALL_DIR}"
)

#---------------------------------
# Export from build tree
#---------------------------------

export(
    EXPORT
        AXHELTargets
    FILE
        "${AXHEL_TARGET_FILENAME}"
    NAMESPACE
        AXHEL::
)

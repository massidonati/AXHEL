# Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
# SPDX-License-Identifier: Apache-2.0

set(NATIVE_SRC
    eltwise/eltwise-add-mod.cpp
    eltwise/eltwise-sub-mod.cpp
    eltwise/eltwise-mul-mod.cpp
    eltwise/eltwise-fma-mod.cpp
    eltwise/eltwise-reduce-mod.cpp
    number-theory/multiply-factor.cpp
)

if (AXHEL_HAS_SVE)
    set(SVE_SRC
        eltwise/eltwise-add-mod-sve.cpp
        eltwise/eltwise-sub-mod-sve.cpp
        eltwise/eltwise-mul-mod-sve.cpp
        eltwise/eltwise-fma-mod-sve.cpp
        eltwise/eltwise-reduce-mod-sve.cpp
    )
endif()

if (AXHEL_HAS_SVE2)
    set(SVE2_SRC
        eltwise/eltwise-mul-mod-sve2.cpp
        eltwise/eltwise-fma-mod-sve2.cpp
    )
endif()

set(AXHEL_SRC "${NATIVE_SRC};${SVE_SRC};${SVE2_SRC}")

if (AXHEL_SHARED_LIB)
    add_library(axhel SHARED ${AXHEL_SRC})
else()
     add_library(axhel STATIC ${AXHEL_SRC})
endif()
add_library(AXHEL::axhel ALIAS axhel)

axhel_add_asan_flag(axhel)

set(AXHEL_DEFINES_IN_FILENAME ${CMAKE_CURRENT_SOURCE_DIR}/include/axhel/util/defines.hpp.in)
set(AXHEL_DEFINES_FILENAME ${CMAKE_CURRENT_SOURCE_DIR}/include/axhel/util/defines.hpp)
configure_file(${AXHEL_DEFINES_IN_FILENAME} ${AXHEL_DEFINES_FILENAME})

set_target_properties(axhel PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(axhel PROPERTIES VERSION ${AXHEL_VERSION})
if (AXHEL_DEBUG)
    set_target_properties(axhel PROPERTIES OUTPUT_NAME "axhel_debug")
else()
    set_target_properties(axhel PROPERTIES OUTPUT_NAME "axhel")
endif()

target_include_directories(axhel
    PRIVATE ${AXHEL_SRC_ROOT_DIR}                              # Private headers
    PUBLIC  $<BUILD_INTERFACE:${AXHEL_INC_ROOT_DIR}>           # Public headers
    PUBLIC  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>   # Public headers
)

#if(CpuFeatures_FOUND)
#    target_include_directories(axhel PUBLIC ${CpuFeatures_INCLUDE_DIR}) # Public headers
#endif()


#TODO: specializzare i flag di compilazione
#-march=armv8-a+sve
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(axhel PRIVATE -Wall -Wconversion -Wshadow -Wextra
        -Wno-unknown-pragmas -mcpu=native -O3 -fomit-frame-pointer
        -Wno-sign-conversion
        -Wno-implicit-int-conversion
    )
    # Avoid 3rd-party dependency warnings when including AXHEL as a dependency
    target_compile_options(axhel PUBLIC
        -Wno-unknown-warning
        -Wno-unknown-warning-option
    )

else()
    message(FATAL_ERROR "Unsupported compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()

install(DIRECTORY ${AXHEL_INC_ROOT_DIR}/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/
        FILES_MATCHING
        PATTERN "*.hpp"
        PATTERN "*.h")

install(TARGETS axhel DESTINATION ${CMAKE_INSTALL_LIBDIR})

set(AXHEL_TARGET_FILENAME ${CMAKE_CURRENT_BINARY_DIR}/cmake/axhel-${AXHEL_VERSION}/AXHELTargets.cmake)
set(AXHEL_CONFIG_IN_FILENAME ${AXHEL_CMAKE_PATH}/AXHELConfig.cmake.in)
set(AXHEL_CONFIG_FILENAME ${AXHEL_ROOT_DIR}/cmake/axhel-${AXHEL_VERSION}/AXHELConfig.cmake)
set(AXHEL_CONFIG_VERSION_FILENAME ${CMAKE_CURRENT_BINARY_DIR}/cmake/axhel-${AXHEL_VERSION}/AXHELConfigVersion.cmake)
set(AXHEL_CONFIG_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/axhel-${AXHEL_VERSION}/)

install(
    EXPORT AXHELTargets
    NAMESPACE AXHEL::
    DESTINATION ${AXHEL_CONFIG_INSTALL_DIR}
)

write_basic_package_version_file(
  ${AXHEL_CONFIG_VERSION_FILENAME}
  VERSION ${AXHEL_VERSION}
  COMPATIBILITY ExactVersion)

include(CMakePackageConfigHelpers)
    configure_package_config_file(
        ${AXHEL_CONFIG_IN_FILENAME} ${AXHEL_CONFIG_FILENAME}
        INSTALL_DESTINATION ${AXHEL_CONFIG_INSTALL_DIR}
    )

install(
    TARGETS axhel
    EXPORT AXHELTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

install(FILES ${AXHEL_CONFIG_FILENAME}
              ${AXHEL_CONFIG_VERSION_FILENAME}
        DESTINATION ${AXHEL_CONFIG_INSTALL_DIR})

export(EXPORT AXHELTargets
        FILE ${AXHEL_TARGET_FILENAME})

# Pkgconfig
#get_target_property(AXHEL_TARGET_NAME axhel OUTPUT_NAME)

#configure_file(${AXHEL_ROOT_DIR}/pkgconfig/axhel.pc.in
#    ${AXHEL_ROOT_DIR}/pkgconfig/axhel.pc @ONLY)

#if(EXISTS ${AXHEL_ROOT_DIR}/pkgconfig/axhel.pc)
#  install(
#      FILES ${AXHEL_ROOT_DIR}/pkgconfig/axhel.pc
#      DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
#endif()
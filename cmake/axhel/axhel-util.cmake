# Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
# SPDX-License-Identifier: Apache-2.0

# Check if source can be compiled and run 
function(axhel_check_compile_flag SOURCE_FILE OUTPUT_FLAG)
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        set(NATIVE_COMPILE_OPTIONS "-mcpu=native")
    else()
        message(FATAL_ERROR "Compiler not supported: ${CMAKE_CXX_COMPILER_ID}")
    endif()

    try_run(CAN_RUN CAN_COMPILE ${CMAKE_BINARY_DIR}
        "${SOURCE_FILE}"
        COMPILE_DEFINITIONS ${NATIVE_COMPILE_OPTIONS}
        OUTPUT_VARIABLE TRY_COMPILE_OUTPUT
    )
    if (CAN_COMPILE AND CAN_RUN EQUAL 0)
        message(STATUS "Setting ${OUTPUT_FLAG}")
        add_compile_definitions(${OUTPUT_FLAG})
        set(${OUTPUT_FLAG} 1 PARENT_SCOPE)
    else()
        message(STATUS "Compile flag not found: ${OUTPUT_FLAG}")
        set(${OUTPUT_FLAG} 0 PARENT_SCOPE)
    endif()    
endfunction()

# Check compiler versions
function(axhel_check_compiler_version)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7.0)
            message(FATAL_ERROR "Requires GCC >= 7.0")
        endif()
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
            message(FATAL_ERROR "Requires Clang >= 5.0")
        endif()
    endif()
endfunction()

# If the input variable is set, stores its value in a _CACHE variable
function(axhel_cache_variable variable)
  if (DEFINED ${variable})
    set(${variable}_CACHE ${${variable}} PARENT_SCOPE)
  endif()
endfunction()

# If the input variable is cached, restores its value from the cache
function(axhel_uncache_variable variable)
  if (DEFINED ${variable}_CACHE)
    set(${variable} ${${variable}_CACHE} CACHE BOOL "" FORCE )
  endif()
endfunction()

# Defines compiler-specific CMake variables
function(axhel_add_compiler_definition)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(AXHEL_USE_GNU ON PARENT_SCOPE)
  elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(AXHEL_USE_CLANG ON PARENT_SCOPE)
  else()
      message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}")
  endif()
endfunction()

# Link with AddressSanitizer in Debug mode on Mac/Linux
function(axhel_add_asan_flag target)
  if(AXHEL_DEBUG AND UNIX)
    target_compile_options(${target} PUBLIC -fsanitize=address)
    target_link_options(${target} PUBLIC -fsanitize=address)
    set(AXHEL_ASAN_LINK "-fsanitize=address" PARENT_SCOPE)
  else()
    set(AXHEL_ASAN_LINK "" PARENT_SCOPE)
  endif()
endfunction()

# Add dependency to the target archive
function(axhel_create_archive target dependency)
  # For proper export of AXHELConfig.cmake / AXHELTargets.cmake,
  # we avoid explicitly linking dependencies via target_link_libraries, since
  # this would add dependencies to the exported axhel target.
  add_dependencies(${target} ${dependency})

  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_custom_command(TARGET ${target} POST_BUILD
                      COMMAND ar -x $<TARGET_FILE:${target}>
                      COMMAND ar -x $<TARGET_FILE:${dependency}>
                      COMMAND ar -qcs $<TARGET_FILE:${target}> *.o
                      COMMAND rm -f *.o
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      DEPENDS ${target} ${dependency})
  else()
    message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}")
  endif()
endfunction()
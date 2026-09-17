# Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
# SPDX-License-Identifier: Apache-2.0

add_executable(axhel-benchmark
    src/axhel_benchmark.cpp
)

target_link_libraries(axhel-benchmark PRIVATE AXHEL::axhel)
target_compile_features(axhel-benchmark PRIVATE cxx_std_17)

target_compile_definitions(axhel-benchmark PRIVATE
    AXHEL_BENCHMARK_CPU_TARGET="${AXHEL_CPU}"
    AXHEL_BENCHMARK_SVE_DETECTED=${AXHEL_SVE_DETECTED}
)

if(AXHEL_BENCHMARK_FORCE_SCALAR)
    target_compile_definitions(axhel-benchmark PRIVATE AXHEL_BENCHMARK_STRICT_SCALAR=1)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(axhel-benchmark PRIVATE
        -Wall
        -Wextra
        -Wconversion
        -Wshadow
        -Wno-sign-conversion
        "-mcpu=${AXHEL_CPU}"
    )

    # Keep the entire strict-scalar benchmark binary free from compiler
    # auto-vectorization as an additional guard. The critical flags are also
    # applied independently to the AXHEL library target itself.
    if(AXHEL_BENCHMARK_FORCE_SCALAR)
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            target_compile_options(axhel-benchmark PRIVATE
                -fno-tree-vectorize
                -fno-tree-loop-vectorize
                -fno-tree-slp-vectorize
            )
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            target_compile_options(axhel-benchmark PRIVATE
                -fno-vectorize
                -fno-slp-vectorize
            )
        endif()
    endif()
endif()

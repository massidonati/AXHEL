# AXHEL deterministic extended regression suite (V5)
include(CTest)

find_package(GTest QUIET)
if(NOT GTest_FOUND)
    if(AXHEL_TEST_FETCH_GTEST)
        include(FetchContent)
        FetchContent_Declare(
            googletest
            URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.zip
        )
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(googletest)
    else()
        message(FATAL_ERROR
            "GoogleTest not found. Install libgtest-dev or configure "
            "-DAXHEL_TEST_FETCH_GTEST=ON")
    endif()
endif()

include(GoogleTest)

set(AXHEL_PUBLIC_TEST_SOURCES
    public/test_eltwise.cpp
    public/test_reduce.cpp
    public/test_ntt.cpp
    public/test_polymul.cpp
)

set(AXHEL_INTERNAL_TEST_SOURCES
    internal/test_internal_native.cpp
)

# SVE-only white-box tests are compiled only when the effective backend is SVE.
if(AXHEL_HAS_SVE)
    list(APPEND AXHEL_INTERNAL_TEST_SOURCES
        internal/test_internal_sve.cpp
    )
endif()

add_executable(axhel-tests-regression
    ${AXHEL_PUBLIC_TEST_SOURCES}
    ${AXHEL_INTERNAL_TEST_SOURCES}
)

target_link_libraries(axhel-tests-regression PRIVATE AXHEL::axhel GTest::gtest_main)
target_compile_features(axhel-tests-regression PRIVATE cxx_std_17)

# White-box tests intentionally include AXHEL private source-tree headers.
target_include_directories(axhel-tests-regression PRIVATE
    ${PROJECT_SOURCE_DIR}/axhel
    ${PROJECT_SOURCE_DIR}/axhel/include
    ${CMAKE_CURRENT_SOURCE_DIR}/common
    ${CMAKE_CURRENT_SOURCE_DIR}/internal
)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Compile tests for exactly the same CPU target as the AXHEL library.
    target_compile_options(axhel-tests-regression PRIVATE
        -Wall
        -Wextra
        "-mcpu=${AXHEL_CPU}"
    )
endif()

gtest_discover_tests(axhel-tests-regression
    TEST_PREFIX "regression::"
    PROPERTIES LABELS "regression"
    DISCOVERY_TIMEOUT 60)

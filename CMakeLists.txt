cmake_minimum_required(VERSION 3.16)
project(NBodvGravitationalIntegrator CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler flags: /O2 /W4 for MSVC; -O2 -Wall -Wextra for GCC/Clang
add_compile_options(
    $<$<CXX_COMPILER_ID:MSVC>:/O2>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-O2>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wextra>
)

# ── Main executable ──────────────────────────────────────────────────────────
add_executable(nbody src/main.cpp)
target_include_directories(nbody PRIVATE src)

# ── Doctest (single-header, fetched at configure time) ───────────────────────
include(FetchContent)
FetchContent_Declare(
    doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG        v2.4.11
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(doctest)

# ── Test executable ──────────────────────────────────────────────────────────
add_executable(nbody_tests tests/test_main.cpp)
target_include_directories(nbody_tests PRIVATE src)
target_link_libraries(nbody_tests PRIVATE doctest::doctest)
# Suppress third-party warnings in the test target under MSVC
target_compile_options(nbody_tests PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/W3>)

enable_testing()
add_test(NAME nbody_unit_tests COMMAND nbody_tests)

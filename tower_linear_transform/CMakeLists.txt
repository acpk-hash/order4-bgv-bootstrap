cmake_minimum_required(VERSION 3.16)
project(tower_bgv LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -march=native")

find_library(NTL_LIB ntl PATHS /usr/local/lib)
find_library(GMP_LIB gmp)
find_library(PTHREAD_LIB pthread)
include_directories(/usr/local/include)
include_directories(${CMAKE_SOURCE_DIR}/src)

# Core library
add_library(tower_bgv_lib STATIC
    src/poly_ring.cpp
    src/ciphertext.cpp
    src/keys.cpp
    src/automorph.cpp
    src/tower_transform.cpp
    src/poly_eval.cpp
    src/bootstrap.cpp
)
target_link_libraries(tower_bgv_lib ${NTL_LIB} ${GMP_LIB} ${PTHREAD_LIB})

# Tests
add_executable(test_basic test/test_basic.cpp)
target_link_libraries(test_basic tower_bgv_lib)

add_executable(test_tower test/test_tower.cpp)
target_link_libraries(test_tower tower_bgv_lib)

# Benchmark
add_executable(bench_tower bench/bench_tower_vs_bsgs.cpp)
target_link_libraries(bench_tower tower_bgv_lib)

add_executable(bench_fair bench/bench_fair_comparison.cpp)
target_link_libraries(bench_fair tower_bgv_lib)

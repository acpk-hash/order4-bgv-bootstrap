cmake_minimum_required(VERSION 3.16)
project(order4_bench CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Point this at the HElib install that carries the order-four evaluator. The
# benchmark is self-contained apart from this: the parameter table, the
# evaluator selection and the driver all live in this directory. Bootstrapping
# itself cannot be vendored, so HElib remains an external dependency.
#
#   cmake -S . -B build -Dhelib_DIR=/path/to/helib/share/cmake/helib
if(NOT DEFINED helib_DIR)
  message(WARNING
    "helib_DIR is not set. Configure with\n"
    "  cmake -S . -B build -Dhelib_DIR=<prefix>/share/cmake/helib\n"
    "where <prefix> is the HElib install carrying the order-four evaluator.")
endif()

find_package(helib REQUIRED)

add_executable(bench_boot bench_boot.cpp)
target_include_directories(bench_boot PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(bench_boot helib ntl gmp)

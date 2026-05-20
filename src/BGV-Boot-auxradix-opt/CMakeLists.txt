cmake_minimum_required(VERSION 3.16)
project(bgv_bts)

SET(CMAKE_CXX_STANDARD 17)

# Isolated artifact build. Do not link against a baseline install tree.
# The top-level build script passes helib_DIR explicitly.  The fallback below
# keeps direct CMake invocation portable inside this repository.
if(NOT helib_DIR)
  get_filename_component(ARTIFACT_ROOT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
  set(helib_DIR "${ARTIFACT_ROOT}/local/helib_auxradix_opt/share/cmake/helib")
endif()
find_package(helib REQUIRED)

add_executable(fatboot fatboot.cpp)
target_link_libraries(fatboot helib ntl gmp)
add_executable(test_batch test_batch_amortization.cpp)
target_link_libraries(test_batch helib ntl gmp)
add_executable(pcmm_bgv pcmm_bgv.cpp)
target_link_libraries(pcmm_bgv helib ntl gmp)
add_executable(batchBGVBoot batchBGVBoot.cpp)
target_link_libraries(batchBGVBoot helib ntl gmp)
add_executable(batchAutomorph batchAutomorph.cpp)
target_link_libraries(batchAutomorph helib ntl gmp)
add_executable(batchFatboot batchFatboot.cpp)
target_link_libraries(batchFatboot helib ntl gmp)
add_executable(batchE2E batchE2E.cpp)
target_link_libraries(batchE2E helib ntl gmp)

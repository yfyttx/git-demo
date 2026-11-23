cmake_minimum_required(VERSION 3.10)
project(QuantumLDPC LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)

set(QLDPC_COMMON_SOURCES
    src/gf2p/gf2p.cpp
    src/ldpc_codes/binary_codes.cpp
    src/ldpc_codes/nonbinary_codes.cpp
    src/utils/matrix_utils.cpp
    src/decoder/NB_BP_Decoder.cpp
    src/decoder/NB_BP_Factorized.cpp
    peg-unige/BigGirth.C
    peg-unige/Random.C
    peg-unige/CyclesOfGraph.C
)

add_library(qldpc_common STATIC ${QLDPC_COMMON_SOURCES})
target_include_directories(qldpc_common PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/peg-unige
)
target_link_libraries(qldpc_common PUBLIC Eigen3::Eigen)

add_executable(construct src/main.cpp)
target_link_libraries(construct PRIVATE qldpc_common)

add_executable(ldpc_constructor src/main.cpp)
target_link_libraries(ldpc_constructor PRIVATE qldpc_common)

add_executable(simulate_qec src/simulate_qec.cpp)
target_link_libraries(simulate_qec PRIVATE qldpc_common)

add_executable(kasai_paper src/kasai_paper.cpp)
target_link_libraries(kasai_paper PRIVATE qldpc_common)

add_executable(gf256_cascade_construct experiments/gf256_cascade/construct.cpp)
target_link_libraries(gf256_cascade_construct PRIVATE qldpc_common)

add_executable(gf256_cascade_simulate experiments/gf256_cascade/simulate.cpp)
target_link_libraries(gf256_cascade_simulate PRIVATE qldpc_common)
add_executable(audit_small_logical_map experiments/gf256_cascade/audit_small_logical_map.cpp)
target_link_libraries(audit_small_logical_map PRIVATE qldpc_common)

add_executable(test_gf256_soft_x_prior
    experiments/gf256_cascade/test_soft_x_prior.cpp
)
target_link_libraries(test_gf256_soft_x_prior PRIVATE qldpc_common)

add_executable(test_gf256_cascade_coset_invariants
    experiments/gf256_cascade/test_cascade_coset_invariants.cpp
)
target_link_libraries(test_gf256_cascade_coset_invariants PRIVATE qldpc_common)

add_executable(test_generic_g2_turbo_exact
    experiments/gf256_cascade/test_generic_g2_turbo_exact.cpp
)
target_link_libraries(test_generic_g2_turbo_exact PRIVATE qldpc_common)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(qldpc_common PRIVATE -O3 -Wall)
    target_compile_options(construct PRIVATE -O3 -Wall)
    target_compile_options(ldpc_constructor PRIVATE -O3 -Wall)
    target_compile_options(simulate_qec PRIVATE -O3 -Wall)
    target_compile_options(kasai_paper PRIVATE -O3 -Wall)
    target_compile_options(gf256_cascade_construct PRIVATE -O3 -Wall)
    target_compile_options(gf256_cascade_simulate PRIVATE -O3 -Wall)
    target_compile_options(test_gf256_soft_x_prior PRIVATE -O3 -Wall)
    target_compile_options(test_gf256_cascade_coset_invariants PRIVATE -O3 -Wall)
    target_compile_options(test_generic_g2_turbo_exact PRIVATE -O3 -Wall)
endif()

add_executable(audit_bp_bridge experiments/gf256_cascade/audit_bp_bridge.cpp)
target_link_libraries(audit_bp_bridge PRIVATE qldpc_common)

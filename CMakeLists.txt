cmake_minimum_required(VERSION 3.10)
project(QuantumLDPC)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)

include_directories(
    include
    peg-unige
    ${EIGEN3_INCLUDE_DIR}
)

# 公共逻辑源文件（不含 main 函数）
set(COMMON_SOURCES
    src/gf2p/gf2p.cpp
    src/ldpc_codes/binary_codes.cpp
    src/ldpc_codes/nonbinary_codes.cpp
    src/utils/matrix_utils.cpp
    src/decoder/NB_BP_Decoder.cpp   # 新增译码器实现
    peg-unige/BigGirth.C
    peg-unige/Random.C
    peg-unige/CyclesOfGraph.C
)

# 目标 1: 矩阵构造器
add_executable(ldpc_constructor src/main.cpp ${COMMON_SOURCES})
target_link_libraries(ldpc_constructor Eigen3::Eigen)

# 目标 2: 量子纠错仿真器 (新)
add_executable(simulate_qec src/simulate_qec.cpp ${COMMON_SOURCES})
target_link_libraries(simulate_qec Eigen3::Eigen)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(ldpc_constructor PRIVATE -O3 -Wall)
    target_compile_options(simulate_qec PRIVATE -O3 -Wall)
endif()
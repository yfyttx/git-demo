cmake_minimum_required(VERSION 3.10)

# 项目名称
project(QuantumLDPC)

# 设置 C++ 标准 (推荐 C++17 或更高)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---------------------------------------------------------
# 1. 查找依赖库
# ---------------------------------------------------------
# 查找 Eigen3 (矩阵计算库)
find_package(Eigen3 REQUIRED)

# ---------------------------------------------------------
# 2. 包含头文件目录
# ---------------------------------------------------------
# include_directories 会告诉编译器去哪里找 .h 文件
include_directories(
    include       # 你的头文件 (如 gf2p/gf2p.h, ldpc_codes/binary_codes.h)
    peg-unige     # PEG 库的头文件 (如 BigGirth.h)
    ${EIGEN3_INCLUDE_DIR}
)

# ---------------------------------------------------------
# 3. 定义源文件
# ---------------------------------------------------------
set(SOURCES
    # --- 你的核心代码 ---
    src/main.cpp
    src/gf2p/gf2p.cpp
    src/ldpc_codes/binary_codes.cpp
    src/ldpc_codes/nonbinary_codes.cpp
    # 如果 matrix_utils.cpp 包含实现，也需要加上
    src/utils/matrix_utils.cpp 

    # --- PEG 官方库代码 ---
    # 注意：不要包含 MainPEG.C，因为它里面有 main 函数，会冲突
    peg-unige/BigGirth.C
    peg-unige/Random.C
    peg-unige/CyclesOfGraph.C
)

# ---------------------------------------------------------
# 4. 生成可执行文件
# ---------------------------------------------------------
add_executable(ldpc_constructor ${SOURCES})

# ---------------------------------------------------------
# 5. 链接库
# ---------------------------------------------------------
target_link_libraries(ldpc_constructor Eigen3::Eigen)

# ---------------------------------------------------------
# (可选) 设置编译选项，例如开启优化
# ---------------------------------------------------------
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(ldpc_constructor PRIVATE -O3 -Wall)
endif()

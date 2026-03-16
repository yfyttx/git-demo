cmake_minimum_required(VERSION 3.10)

# 项目名称
project(QuantumLDPC)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---------------------------------------------------------
# 1. 查找依赖库
# ---------------------------------------------------------
find_package(Eigen3 REQUIRED)

# ---------------------------------------------------------
# 2. 定义源文件 (建议使用相对路径确保稳健性)
# ---------------------------------------------------------
set(SOURCES
    src/main.cpp
    # 关键点：添加下面这个文件，因为 constructBinaryHMatrices 在这里定义
    src/ldpc_codes/binary_codes.cpp  
    
    # 如果你还有其他的工具类，也一并加上
    # src/gf2p/gf2p.cpp 
    # src/utils/matrix_utils.cpp 

    # PEG 库文件
    peg-unige/BigGirth.C
    peg-unige/Random.C
    peg-unige/CyclesOfGraph.C
)

# ---------------------------------------------------------
# 3. 生成可执行文件
# ---------------------------------------------------------
# 这里把生成的文件名改为 construct，方便刚才的 Python 脚本调用
add_executable(construct ${SOURCES})

# ---------------------------------------------------------
# 4. 包含目录与链接库
# ---------------------------------------------------------
target_include_directories(construct PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/peg-unige
    ${EIGEN3_INCLUDE_DIRS}
)

target_link_libraries(construct PRIVATE Eigen3::Eigen)

# ---------------------------------------------------------
# 5. 编译优化 (这对大规模矩阵生成至关重要)
# ---------------------------------------------------------
if(MSVC)
    target_compile_options(construct PRIVATE /O2 /W3)
else()
    # 开启 -O3 优化和高级向量化，加快矩阵乘法和 PEG 运行速度
    target_compile_options(construct PRIVATE -O3 -Wall -march=native)
endif()

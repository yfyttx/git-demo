cmake_minimum_required(VERSION 3.10)
project(QuantumLDPC)

set(CMAKE_CXX_STANDARD 17)

# 查找 Eigen3 (通常需要安装: sudo apt install libeigen3-dev 或 brew install eigen)
find_package(Eigen3 REQUIRED)

# 包含头文件目录
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(include)

# 添加源文件
file(GLOB_RECURSE SOURCES "src/*.cpp")

# 创建可执行文件
add_executable(ldpc_constructor ${SOURCES})

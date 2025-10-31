cmake_minimum_required(VERSION 3.10)
project(QuantumLDPCConstructor CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include_directories(
    "C:/msys64/mingw64/include/eigen3"
    ${PROJECT_SOURCE_DIR}/include
)

set(SOURCES
    src/main.cpp
    src/gf2p/gf2p.cpp
    src/ldpc_codes/binary_codes.cpp
    src/ldpc_codes/nonbinary_codes.cpp
    src/utils/matrix_utils.cpp
)

add_executable(ldpc_constructor ${SOURCES})



message(STATUS "Using hardcoded Eigen path: C:/msys64/mingw64/include/eigen3")
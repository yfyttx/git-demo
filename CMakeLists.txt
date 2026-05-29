cmake_minimum_required(VERSION 3.10)
project(QuantumLDPC LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED)

set(COMMON_SOURCES
    src/gf2p/gf2p.cpp
    src/ldpc_codes/binary_codes.cpp
    src/ldpc_codes/nonbinary_codes.cpp
    src/utils/matrix_utils.cpp
    src/decoder/NB_BP_Decoder.cpp
    peg-unige/BigGirth.C
    peg-unige/Random.C
    peg-unige/CyclesOfGraph.C
)

add_executable(ldpc_constructor
    src/main.cpp
    ${COMMON_SOURCES}
)

add_executable(simulate_qec
    src/simulate_qec.cpp
    ${COMMON_SOURCES}
)

foreach(target_name IN ITEMS ldpc_constructor simulate_qec)
    target_include_directories(${target_name} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/peg-unige
    )
    target_link_libraries(${target_name} PRIVATE Eigen3::Eigen)
endforeach()

add_custom_target(construct DEPENDS ldpc_constructor)

if(MSVC)
    target_compile_options(ldpc_constructor PRIVATE /O2 /W3)
    target_compile_options(simulate_qec PRIVATE /O2 /W3)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(ldpc_constructor PRIVATE -O3 -Wall)
    target_compile_options(simulate_qec PRIVATE -O3 -Wall)
endif()

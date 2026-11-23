#ifndef NONBINARY_CODES_H
#define NONBINARY_CODES_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <utility>
#include <random>

// 前向声明，避免循环引用
class GF2p; 

// --- 新增的函数声明 ---

// 1. 根据二进制骨架，构建在 Z_(2^p-1) 上的线性方程组
Eigen::MatrixXi buildLinearSystem(
    const std::vector<std::pair<int, int>>& Hc_nonzero,
    const Eigen::MatrixXi& Hc,
    const Eigen::MatrixXi& Hd,
    int mod);

// 2. 在 Z_mod 环上求解 Ax=0；随机自由变量由调用方提供的 RNG 决定
std::vector<int> gaussElimination(
    const Eigen::MatrixXi& system, 
    int mod,
    std::mt19937& gen);

// 3. 为 H_delta 的每一行，求解一个在 GF(2^p) 上的局部方程组
std::vector<int> solveDeltaRow(
    int m_prime,
    const Eigen::MatrixXi& Hc,
    const Eigen::MatrixXi& Hd, 
    const Eigen::MatrixXi& H_gamma,
    const GF2p& gf);

// Kasai et al., Sec. II-B: choose Gamma labels from the modular cycle
// constraints, then solve each Delta row in the nullspace of Gamma.
std::pair<Eigen::SparseMatrix<int>, Eigen::SparseMatrix<int>>
constructKasaiNonbinaryPair(const Eigen::MatrixXi& hc,
                           const Eigen::MatrixXi& hd,
                           GF2p& gf, std::mt19937& rng);

#endif // NONBINARY_CODES_H

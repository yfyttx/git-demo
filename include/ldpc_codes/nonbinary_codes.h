#ifndef NONBINARY_CODES_H
#define NONBINARY_CODES_H

#include <Eigen/Dense>
#include <vector>
#include <utility>

// 前向声明，避免循环引用
class GF2p; 

// --- 新增的函数声明 ---

// 1. 根据二进制骨架，构建在 Z_(2^p-1) 上的线性方程组
Eigen::MatrixXi buildLinearSystem(
    const std::vector<std::pair<int, int>>& Hc_nonzero,
    const Eigen::MatrixXi& Hc,
    const Eigen::MatrixXi& Hd,
    int mod);

// 2. 在 Z_mod 域上用高斯消元法求解 Ax=0
std::vector<int> gaussElimination(
    const Eigen::MatrixXi& system, 
    int mod);

// 3. 为 H_delta 的每一行，求解一个在 GF(2^p) 上的局部方程组
std::vector<int> solveDeltaRow(
    int m_prime,
    const Eigen::MatrixXi& Hc,
    const Eigen::MatrixXi& Hd, 
    const Eigen::MatrixXi& H_gamma,
    const GF2p& gf);

#endif // NONBINARY_CODES_H

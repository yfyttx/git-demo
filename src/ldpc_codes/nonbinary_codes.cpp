#include "ldpc_codes/nonbinary_codes.h"
#include "gf2p/gf2p.h"
#include <algorithm>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <vector>

using namespace Eigen;

// ==========================================
// 辅助函数 (内部使用)
// ==========================================

// 获取与 H_delta 第 m' 行相关的 H_gamma 中的非零位置 E(m')
std::pair<std::vector<int>, std::vector<std::pair<int, int>>>
getNandE(int m_prime, const Eigen::MatrixXi &Hc, const Eigen::MatrixXi &Hd) {
    std::vector<int> N;
    for (int n = 0; n < Hd.cols(); ++n) {
        if (Hd(m_prime, n) == 1)
            N.push_back(n);
    }
    std::vector<std::pair<int, int>> E;
    for (int n : N) {
        for (int m = 0; m < Hc.rows(); ++m) {
            if (Hc(m, n) == 1)
                E.emplace_back(m, n);
        }
    }
    return {N, E};
}

// 追踪环路
std::vector<std::pair<int, int>> trace_cycle(const std::vector<std::pair<int, int>> &E) {
    if (E.size() < 2)
        return {};

    std::map<int, std::vector<int>> row_to_cols;
    std::map<int, std::vector<int>> col_to_rows;

    for (const auto &edge : E) {
        row_to_cols[edge.first].push_back(edge.second);
        col_to_rows[edge.second].push_back(edge.first);
    }

    std::vector<std::pair<int, int>> path;
    auto start_edge = E[0];
    path.push_back(start_edge);

    int current_row = start_edge.first;
    int current_col = start_edge.second;

    for (size_t i = 0; i < E.size() - 1; ++i) {
        int next_row = -1;
        for (int r : col_to_rows[current_col]) {
            if (r != current_row) {
                next_row = r;
                break;
            }
        }
        if (next_row == -1)
            return {};
        path.emplace_back(next_row, current_col);
        current_row = next_row;
        if (path.size() >= E.size())
            break;

        int next_col = -1;
        for (int c : row_to_cols[current_row]) {
            if (c != current_col) {
                next_col = c;
                break;
            }
        }
        if (next_col == -1)
            return {};
        path.emplace_back(current_row, next_col);
        current_col = next_col;
    }
    return path;
}

// ==========================================
// 核心函数实现
// ==========================================

// [实现 1] 建立线性方程组
MatrixXi buildLinearSystem(const std::vector<std::pair<int, int>> &Hc_nonzero,
                           const Eigen::MatrixXi &Hc, const Eigen::MatrixXi &Hd, int mod) {

    int var_count = Hc_nonzero.size();
    if (var_count == 0)
        return MatrixXi::Zero(0, 1);

    std::map<std::pair<int, int>, int> pos_to_idx;
    for (int i = 0; i < var_count; ++i)
        pos_to_idx[Hc_nonzero[i]] = i;

    std::vector<std::vector<int>> equations;

    for (int m_prime = 0; m_prime < Hd.rows(); ++m_prime) {
        auto [N, E] = getNandE(m_prime, Hc, Hd);
        if (E.size() < 2)
            continue;

        std::vector<std::pair<int, int>> cycle_path = trace_cycle(E);
        if (cycle_path.size() != E.size())
            continue;

        std::vector<int> eq(var_count, 0);
        for (size_t i = 0; i < cycle_path.size(); ++i) {
            int var_idx = pos_to_idx.at(cycle_path[i]);
            // 论文 Eq(14): 沿环路交替求和
            eq[var_idx] = (i % 2 == 0) ? 1 : -1;
        }
        equations.push_back(eq);
    }

    if (equations.empty())
        return MatrixXi::Zero(0, var_count + 1);

    MatrixXi system(equations.size(), var_count + 1);
    for (size_t i = 0; i < equations.size(); ++i) {
        for (int j = 0; j < var_count; ++j) {
            system(i, j) = (equations[i][j] % mod + mod) % mod;
        }
        system(i, var_count) = 0;
    }
    return system;
}

// [实现 2] 高斯消元 (解决 Linker Error)
std::vector<int> gaussElimination(const Eigen::MatrixXi &system, int mod) {
    int var_count = system.cols() - 1;
    if (var_count <= 0)
        return {};

    Eigen::MatrixXi A = system;
    int rows = A.rows();

    if (rows == 0) {
        std::vector<int> x(var_count);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(1, mod - 1);
        for (int i = 0; i < var_count; ++i)
            x[i] = distrib(gen);
        return x;
    }

    int rank = 0;
    std::vector<int> pivot_col(rows, -1);
    for (int j = 0; j < var_count && rank < rows; ++j) {
        int pivot_row_idx = rank;
        while (pivot_row_idx < rows && A(pivot_row_idx, j) == 0)
            pivot_row_idx++;
        if (pivot_row_idx < rows) {
            A.row(rank).swap(A.row(pivot_row_idx));
            pivot_col[rank] = j;
            int inv = -1;
            for (int i = 1; i < mod; ++i)
                if ((A(rank, j) * i) % mod == 1)
                    inv = i;
            if (inv == -1)
                throw std::runtime_error("Modular inverse does not exist");
            for (int k = j; k < A.cols(); ++k)
                A(rank, k) = (A(rank, k) * inv) % mod;
            for (int i = 0; i < rows; ++i) {
                if (i != rank) {
                    int factor = A(i, j);
                    for (int k = j; k < A.cols(); ++k) {
                        long long temp = A(i, k);
                        temp -= (long long)factor * A(rank, k);
                        A(i, k) = (int)((temp % mod + mod) % mod);
                    }
                }
            }
            rank++;
        }
    }

    std::vector<int> x(var_count);
    std::vector<bool> is_pivot_var(var_count, false);
    for (int i = 0; i < rank; ++i)
        is_pivot_var[pivot_col[i]] = true;

    // 随机化自由变量
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(1, mod - 1);
    for (int j = 0; j < var_count; ++j)
        if (!is_pivot_var[j])
            x[j] = distrib(gen);

    // 回代
    for (int i = rank - 1; i >= 0; --i) {
        int p_col = pivot_col[i];
        long long sum = 0;
        for (int j = p_col + 1; j < var_count; ++j)
            sum += (long long)A(i, j) * x[j];
        x[p_col] = (int)(((long long)A(i, var_count) - (sum % mod) + mod) % mod);
    }
    return x;
}

// [实现 3] 求解 H_delta 的行 (解决 Linker Error)
std::vector<int> solveDeltaRow(int m_prime, const Eigen::MatrixXi &Hc, const Eigen::MatrixXi &Hd,
                               const Eigen::MatrixXi &H_gamma, const GF2p &gf) {

    auto [N, E] = getNandE(m_prime, Hc, Hd);
    int L_size = N.size();
    if (L_size == 0)
        return std::vector<int>(Hd.cols(), 0);

    // 映射到局部小矩阵
    std::map<int, int> m_to_local_row;
    std::map<int, int> n_to_local_col;
    {
        std::set<int> unique_rows;
        for (const auto &p : E)
            unique_rows.insert(p.first);
        int r_idx = 0;
        for (int row : unique_rows)
            m_to_local_row[row] = r_idx++;
        int c_idx = 0;
        for (int col : N)
            n_to_local_col[col] = c_idx++;
    }

    int num_rows = m_to_local_row.size();
    int num_cols = n_to_local_col.size();

    Eigen::MatrixXi A = Eigen::MatrixXi::Zero(num_rows, num_cols);
    for (const auto &p : E) {
        if (m_to_local_row.count(p.first) && n_to_local_col.count(p.second)) {
            A(m_to_local_row.at(p.first), n_to_local_col.at(p.second)) = H_gamma(p.first, p.second);
        }
    }

    // 在 GF(2^p) 域上高斯消元
    int rank = 0;
    std::vector<int> pivot_of_row(num_rows, -1);
    for (int j = 0; j < num_cols && rank < num_rows; ++j) {
        int pivot_row = rank;
        while (pivot_row < num_rows && A(pivot_row, j) == 0)
            pivot_row++;
        if (pivot_row < num_rows) {
            A.row(rank).swap(A.row(pivot_row));
            pivot_of_row[rank] = j;
            int inv = gf.inv(A(rank, j));
            for (int k = j; k < num_cols; ++k)
                A(rank, k) = gf.mul(A(rank, k), inv);
            for (int i = 0; i < num_rows; ++i) {
                if (i != rank && A(i, j) != 0) {
                    int factor = A(i, j);
                    for (int k = j; k < num_cols; ++k) {
                        A(i, k) = gf.add(A(i, k), gf.mul(factor, A(rank, k)));
                    }
                }
            }
            rank++;
        }
    }

    // 提取非零解
    std::vector<int> delta(num_cols, 0);
    std::vector<bool> is_pivot_var(num_cols, false);
    for (int i = 0; i < rank; ++i)
        if (pivot_of_row[i] != -1)
            is_pivot_var[pivot_of_row[i]] = true;

    std::vector<int> free_vars;
    for (int j = 0; j < num_cols; ++j)
        if (!is_pivot_var[j])
            free_vars.push_back(j);

    if (!free_vars.empty()) {
        delta[free_vars[0]] = 1; // 设第一个自由变量为1
        for (int i = rank - 1; i >= 0; --i) {
            int p_col = pivot_of_row[i];
            int sum = 0;
            for (int f_col : free_vars) {
                sum = gf.add(sum, gf.mul(A(i, f_col), delta[f_col]));
            }
            delta[p_col] = sum;
        }
    }

    // 映射回全局向量
    std::vector<int> delta_row(Hd.cols(), 0);
    for (int i = 0; i < num_cols; ++i) {
        delta_row[N[i]] = delta[i];
    }
    return delta_row;
}

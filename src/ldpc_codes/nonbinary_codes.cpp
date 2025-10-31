#include "ldpc_codes/nonbinary_codes.h"
#include "gf2p/gf2p.h"
#include <map>
#include <vector>
#include <random>
#include <algorithm>
#include <set>
#include <iostream> // 用于调试输出

using namespace Eigen;

// 辅助函数: 根据 H_delta 的某一行 m'，找到 Hc 中相关的非零位置 E(m')
std::pair<std::vector<int>, std::vector<std::pair<int, int>>> getNandE(
    int m_prime, const Eigen::MatrixXi& Hc, const Eigen::MatrixXi& Hd) {
    std::vector<int> N;
    for (int n = 0; n < Hd.cols(); ++n) {
        if (Hd(m_prime, n) == 1) {
            N.push_back(n);
        }
    }
    std::vector<std::pair<int, int>> E;
    for (int n : N) {
        for (int m = 0; m < Hc.rows(); ++m) {
            if (Hc(m, n) == 1) {
                E.emplace_back(m, n);
            }
        }
    }
    return {N, E};
}

// [新增] 辅助函数: 追踪由 E(m') 形成的 2L 环路
// Tanner图中的环由行节点(check node)和列节点(variable node)交替连接构成
std::vector<std::pair<int, int>> trace_cycle(const std::vector<std::pair<int, int>>& E) {
    if (E.size() < 2) return {};

    std::map<int, std::vector<int>> row_to_cols; // 检查节点的邻接表
    std::map<int, std::vector<int>> col_to_rows; // 变量节点的邻接表

    for (const auto& edge : E) {
        row_to_cols[edge.first].push_back(edge.second);
        col_to_rows[edge.second].push_back(edge.first);
    }

    std::vector<std::pair<int, int>> path;
    std::set<std::pair<int, int>> visited_edges;

    // 从任意一条边开始遍历
    auto start_edge = E[0];
    path.push_back(start_edge);
    visited_edges.insert(start_edge);

    int current_row = start_edge.first;
    int current_col = start_edge.second;

    // 循环 E.size() - 1 次来找到剩下的所有边
    for (size_t i = 0; i < E.size() - 1; ++i) {
        // 当前在(current_row, current_col)。先沿着列走，找到下一个行。
        int next_row = -1;
        // 一个列节点应该连接两个行节点
        for (int r : col_to_rows.at(current_col)) {
            if (r != current_row) {
                next_row = r;
                break;
            }
        }
        
        // 如果找不到下一个节点，说明图结构有问题，不是一个简单的环
        if (next_row == -1) return {}; 
        
        path.emplace_back(next_row, current_col);
        current_row = next_row;

        if (path.size() >= E.size()) break;

        // 现在在(current_row, current_col)。再沿着行走，找到下一个列。
        int next_col = -1;
        // 一个行节点应该连接两个列节点
        for (int c : row_to_cols.at(current_row)) {
            if (c != current_col) {
                next_col = c;
                break;
            }
        }
        
        if (next_col == -1) return {};

        path.emplace_back(current_row, next_col);
        current_col = next_col;
    }
    return path;
}


// [重写] 核心函数: 根据论文公式(14)建立线性方程组
// 每个方程对应 H_delta 的一行，确保其相关的行列式为零
MatrixXi buildLinearSystem(
    const std::vector<std::pair<int, int>>& Hc_nonzero,
    const Eigen::MatrixXi& Hc,
    const Eigen::MatrixXi& Hd,
    int mod) {

    int var_count = Hc_nonzero.size();
    if (var_count == 0) return MatrixXi::Zero(0, 1);
    
    // 建立从矩阵位置 (row, col) 到变量索引的映射
    std::map<std::pair<int, int>, int> pos_to_idx;
    for (int i = 0; i < var_count; ++i) {
        pos_to_idx[Hc_nonzero[i]] = i;
    }

    std::vector<std::vector<int>> equations;

    // 遍历 H_delta 的每一行，为每一行建立一个约束方程
    for (int m_prime = 0; m_prime < Hd.rows(); ++m_prime) {
        auto [N, E] = getNandE(m_prime, Hc, Hd);

        if (E.size() < 2) continue; // 不是环，无法形成约束

        // 追踪环路，得到有序的边列表
        std::vector<std::pair<int, int>> cycle_path = trace_cycle(E);
        
        // 健壮性检查：如果追踪到的路径长度和边的数量不符，说明不是一个简单的环，跳过
        if (cycle_path.size() != E.size()) {
             std::cerr << "Warning: Failed to trace a simple cycle for m_prime = " << m_prime << ". Skipping constraint." << std::endl;
             continue;
        }

        std::vector<int> eq(var_count, 0);
        for (size_t i = 0; i < cycle_path.size(); ++i) {
            auto pos = cycle_path[i];
            int var_idx = pos_to_idx.at(pos);

            // 根据论文公式(14)，构造交替和
            // Σ log(γ_E1) - Σ log(γ_E2) = 0
            if (i % 2 == 0) {
                eq[var_idx] = 1;
            } else {
                eq[var_idx] = -1; 
            }
        }
        equations.push_back(eq);
    }

    if (equations.empty()) {
        std::cerr << "Warning: No linear system constraints were generated." << std::endl;
        return MatrixXi::Zero(0, var_count + 1);
    }

    // 将方程列表转换为 Eigen 矩阵
    MatrixXi system(equations.size(), var_count + 1);
    for (size_t i = 0; i < equations.size(); ++i) {
        for (int j = 0; j < var_count; ++j) {
            // 处理模运算中的负数
            system(i, j) = (equations[i][j] % mod + mod) % mod;
        }
        system(i, var_count) = 0; // 方程右侧恒为0
    }

    return system;
}

// [保留] 高斯消元求解函数 (无需修改)
std::vector<int> gaussElimination(const Eigen::MatrixXi& system, int mod) {
    int var_count = system.cols() - 1;
    if (var_count <= 0) return {};

    Eigen::MatrixXi A = system;
    int rows = A.rows();

    // 如果没有方程（例如，没有找到任何环），返回一个随机解
    if (rows == 0) {
        std::vector<int> x(var_count);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(1, mod - 1); // 避免0
        for (int i = 0; i < var_count; ++i) x[i] = distrib(gen);
        return x;
    }

    int rank = 0;
    std::vector<int> pivot_col(rows, -1);
    for (int j = 0; j < var_count && rank < rows; ++j) {
        int pivot_row_idx = rank;
        while (pivot_row_idx < rows && A(pivot_row_idx, j) == 0) pivot_row_idx++;
        if (pivot_row_idx < rows) {
            A.row(rank).swap(A.row(pivot_row_idx));
            pivot_col[rank] = j;
            int inv = -1;
            for (int i = 1; i < mod; ++i) if ((A(rank, j) * i) % mod == 1) inv = i;
            if (inv == -1) throw std::runtime_error("Modular inverse does not exist");
            for (int k = j; k < A.cols(); ++k) A(rank, k) = (A(rank, k) * inv) % mod;
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
    for(int i = 0; i < rank; ++i) is_pivot_var[pivot_col[i]] = true;

    // 为自由变量赋随机值
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(1, mod - 1); // 避免0
    for (int j = 0; j < var_count; ++j) if (!is_pivot_var[j]) x[j] = distrib(gen);

    // 回代求解主元变量
    for (int i = rank - 1; i >= 0; --i) {
        int p_col = pivot_col[i];
        long long sum = 0;
        for (int j = p_col + 1; j < var_count; ++j) sum += (long long)A(i, j) * x[j];
        x[p_col] = (int)(((long long)A(i, var_count) - (sum % mod) + mod) % mod);
    }
    return x;
}

// [保留] 为 H_delta 的每一行求解非零解 (无需修改)
std::vector<int> solveDeltaRow(
    int m_prime,
    const Eigen::MatrixXi& Hc,
    const Eigen::MatrixXi& Hd,
    const Eigen::MatrixXi& H_gamma,
    const GF2p& gf) {
    auto [N, E] = getNandE(m_prime, Hc, Hd);
    int L_size = N.size();
    if (L_size == 0) return std::vector<int>(Hd.cols(), 0);

    // 建立局部坐标映射
    std::map<int, int> m_to_local_row;
    std::map<int, int> n_to_local_col;
    {
        std::set<int> unique_rows;
        for(const auto& p : E) unique_rows.insert(p.first);
        int r_idx = 0;
        for(int row : unique_rows) m_to_local_row[row] = r_idx++;
        int c_idx = 0;
        for(int col : N) n_to_local_col[col] = c_idx++;
    }
    
    int num_rows = m_to_local_row.size();
    int num_cols = n_to_local_col.size();
    
    Eigen::MatrixXi A = Eigen::MatrixXi::Zero(num_rows, num_cols);
    for (const auto& p : E) {
        if (m_to_local_row.count(p.first) && n_to_local_col.count(p.second)) {
            A(m_to_local_row.at(p.first), n_to_local_col.at(p.second)) = H_gamma(p.first, p.second);
        }
    }

    // 在 GF(2^p) 上进行高斯消元
    int rank = 0;
    std::vector<int> pivot_of_row(num_rows, -1);
    for (int j = 0; j < num_cols && rank < num_rows; ++j) {
        int pivot_row = rank;
        while (pivot_row < num_rows && A(pivot_row, j) == 0) pivot_row++;
        if (pivot_row < num_rows) {
            A.row(rank).swap(A.row(pivot_row));
            pivot_of_row[rank] = j;
            int inv = gf.inv(A(rank, j));
            for (int k = j; k < num_cols; ++k) A(rank, k) = gf.mul(A(rank, k), inv);
            for (int i = 0; i < num_rows; ++i) {
                if (i != rank) {
                    int factor = A(i, j);
                    if (factor != 0) {
                        for (int k = j; k < num_cols; ++k) {
                            int term = gf.mul(factor, A(rank, k));
                            A(i, k) = gf.add(A(i, k), term); // 在 GF(2^p) 中，加法就是异或
                        }
                    }
                }
            }
            rank++;
        }
    }

    // 寻找一个非零特解
    std::vector<int> delta(num_cols, 0);
    std::vector<bool> is_pivot_var(num_cols, false);
    for(int i=0; i < rank; ++i) if(pivot_of_row[i] != -1) is_pivot_var[pivot_of_row[i]] = true;
    
    std::vector<int> free_vars_indices;
    for(int j = 0; j < num_cols; ++j) if(!is_pivot_var[j]) free_vars_indices.push_back(j);

    // 因为 H_gamma 的构造保证了系统是奇异的，所以一定存在自由变量
    if (!free_vars_indices.empty()) {
        // 将第一个自由变量设为1，以获得一个非零解
        delta[free_vars_indices[0]] = 1; 

        // 回代求解主元变量
        for (int i = rank - 1; i >= 0; --i) {
            int p_col = pivot_of_row[i];
            int sum = 0;
            // 只对自由变量求和
            for (int f_col : free_vars_indices) {
                sum = gf.add(sum, gf.mul(A(i, f_col), delta[f_col]));
            }
            delta[p_col] = sum; // 在 GF(2^p) 中, x = -sum 就是 x = sum
        }
    }

    // 将局部解映射回完整的行向量
    std::vector<int> delta_row(Hd.cols(), 0);
    for (int i = 0; i < num_cols; ++i) {
        delta_row[N[i]] = delta[i];
    }
    return delta_row;
}
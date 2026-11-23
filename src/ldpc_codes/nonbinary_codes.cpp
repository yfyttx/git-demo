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
std::vector<int> gaussElimination(const Eigen::MatrixXi &system, int mod,
                                  std::mt19937 &gen) {
    int var_count = system.cols() - 1;
    if (var_count <= 0)
        return {};

    // GF(1024) labels are exponents modulo 1023 = 3 * 11 * 31. Reduce
    // independently over the three prime fields: Gaussian elimination over
    // Z/1023Z cannot divide by a non-unit pivot.
    if (mod == 1023) {
        constexpr int primes[] = {3, 11, 31};
        std::vector<std::vector<int>> components;
        for (int prime : primes) {
            Eigen::MatrixXi A = system;
            for (int row = 0; row < A.rows(); ++row)
                for (int col = 0; col < A.cols(); ++col)
                    A(row, col) = (A(row, col) % prime + prime) % prime;
            std::vector<int> pivots;
            int rank = 0;
            for (int col = 0; col < var_count && rank < A.rows(); ++col) {
                int pivot = rank;
                while (pivot < A.rows() && A(pivot, col) == 0) ++pivot;
                if (pivot == A.rows()) continue;
                A.row(rank).swap(A.row(pivot));
                int inverse = 1;
                while ((A(rank, col) * inverse) % prime != 1) ++inverse;
                for (int j = col; j <= var_count; ++j)
                    A(rank, j) = (A(rank, j) * inverse) % prime;
                for (int row = 0; row < A.rows(); ++row) {
                    if (row == rank || A(row, col) == 0) continue;
                    const int factor = A(row, col);
                    for (int j = col; j <= var_count; ++j)
                        A(row, j) = (A(row, j) - factor * A(rank, j) % prime
                                     + prime) % prime;
                }
                pivots.push_back(col);
                ++rank;
            }
            for (int row = rank; row < A.rows(); ++row)
                if (A(row, var_count) != 0)
                    throw std::runtime_error("QC label system is inconsistent modulo a CRT prime");
            std::vector<int> solution(var_count, 0);
            std::vector<char> is_pivot(var_count, 0);
            for (int col : pivots) is_pivot[col] = 1;
            std::uniform_int_distribution<> distribution(0, prime - 1);
            for (int col = 0; col < var_count; ++col)
                if (!is_pivot[col]) solution[col] = distribution(gen);
            for (int row = rank - 1; row >= 0; --row) {
                const int pivot_col = pivots[row];
                int value = A(row, var_count);
                for (int col = pivot_col + 1; col < var_count; ++col)
                    value = (value - A(row, col) * solution[col] % prime
                             + prime) % prime;
                solution[pivot_col] = value;
            }
            components.push_back(std::move(solution));
        }
        std::vector<int> combined(var_count);
        for (int col = 0; col < var_count; ++col)
            combined[col] = (682 * components[0][col]
                           + 837 * components[1][col]
                           + 528 * components[2][col]) % mod;
        for (int row = 0; row < system.rows(); ++row) {
            long long residual = -system(row, var_count);
            for (int col = 0; col < var_count; ++col)
                residual += static_cast<long long>(system(row, col)) * combined[col];
            if ((residual % mod + mod) % mod != 0)
                throw std::runtime_error("CRT solution failed the modulo-1023 constraints");
        }
        return combined;
    }

    // The GF(256) label constraints are integer exponent congruences modulo
    // 255, which is composite. Ordinary Gaussian elimination modulo 255 is
    // invalid whenever a pivot is not a unit. Preserve the historical result
    // whenever unit-pivot elimination succeeds; use CRT only for systems that
    // the legacy routine cannot reduce.
    if (mod == 255) {
        Eigen::MatrixXi legacy = system;
        const int rows = legacy.rows();
        int rank = 0;
        std::vector<int> pivots(rows, -1);
        bool nonunit_pivot = false;
        for (int col = 0; col < var_count && rank < rows; ++col) {
            int pivot = rank;
            while (pivot < rows && legacy(pivot, col) == 0) ++pivot;
            if (pivot == rows) continue;
            legacy.row(rank).swap(legacy.row(pivot));
            int inverse = -1;
            for (int value = 1; value < mod; ++value)
                if ((legacy(rank, col) * value) % mod == 1) inverse = value;
            if (inverse == -1) {
                nonunit_pivot = true;
                break;
            }
            pivots[rank] = col;
            for (int j = col; j < legacy.cols(); ++j)
                legacy(rank, j) = (legacy(rank, j) * inverse) % mod;
            for (int row = 0; row < rows; ++row) {
                if (row == rank) continue;
                const int factor = legacy(row, col);
                for (int j = col; j < legacy.cols(); ++j) {
                    const long long value = legacy(row, j) -
                        static_cast<long long>(factor) * legacy(rank, j);
                    legacy(row, j) = static_cast<int>((value % mod + mod) % mod);
                }
            }
            ++rank;
        }
        if (!nonunit_pivot) {
            std::vector<int> solution(var_count, 0);
            std::vector<char> is_pivot(var_count, 0);
            for (int row = 0; row < rank; ++row) is_pivot[pivots[row]] = 1;
            std::uniform_int_distribution<> distrib(0, mod - 1);
            for (int col = 0; col < var_count; ++col)
                if (!is_pivot[col]) solution[col] = distrib(gen);
            for (int row = rank - 1; row >= 0; --row) {
                const int pivot_col = pivots[row];
                long long sum = 0;
                for (int col = pivot_col + 1; col < var_count; ++col)
                    sum += static_cast<long long>(legacy(row, col)) * solution[col];
                solution[pivot_col] = static_cast<int>(
                    (legacy(row, var_count) - sum % mod + mod) % mod);
            }
            bool valid = true;
            for (int row = 0; row < system.rows(); ++row) {
                long long residual = -system(row, var_count);
                for (int col = 0; col < var_count; ++col)
                    residual += static_cast<long long>(system(row, col)) * solution[col];
                valid &= (residual % mod + mod) % mod == 0;
            }
            if (valid) return solution;
            throw std::runtime_error("Legacy modulo-255 solution failed verification");
        }

        // Solve the square-free CRT components separately, then combine each
        // variable's residues when legacy unit-pivot elimination gets stuck.
        constexpr int primes[] = {3, 5, 17};
        std::vector<std::vector<int>> component_solutions;
        component_solutions.reserve(3);
        for (int prime : primes) {
            Eigen::MatrixXi A = system;
            for (int row = 0; row < A.rows(); ++row)
                for (int col = 0; col < A.cols(); ++col)
                    A(row, col) = (A(row, col) % prime + prime) % prime;
            std::vector<int> pivots;
            int rank = 0;
            for (int col = 0; col < var_count && rank < A.rows(); ++col) {
                int pivot = rank;
                while (pivot < A.rows() && A(pivot, col) == 0) ++pivot;
                if (pivot == A.rows()) continue;
                A.row(rank).swap(A.row(pivot));
                int inverse = 1;
                while ((A(rank, col) * inverse) % prime != 1) ++inverse;
                for (int j = col; j <= var_count; ++j)
                    A(rank, j) = (A(rank, j) * inverse) % prime;
                for (int row = 0; row < A.rows(); ++row) {
                    if (row == rank || A(row, col) == 0) continue;
                    const int factor = A(row, col);
                    for (int j = col; j <= var_count; ++j)
                        A(row, j) = (A(row, j) - factor * A(rank, j) % prime
                                      + prime) % prime;
                }
                pivots.push_back(col);
                ++rank;
            }
            for (int row = rank; row < A.rows(); ++row)
                if (A(row, var_count) != 0)
                    throw std::runtime_error("QC label system is inconsistent modulo a CRT prime");
            std::vector<int> solution(var_count, 0);
            std::vector<char> is_pivot(var_count, 0);
            for (int col : pivots) is_pivot[col] = 1;
            std::uniform_int_distribution<> distrib(0, prime - 1);
            for (int col = 0; col < var_count; ++col)
                if (!is_pivot[col]) solution[col] = distrib(gen);
            for (int row = rank - 1; row >= 0; --row) {
                const int pivot_col = pivots[row];
                int value = A(row, var_count);
                for (int col = pivot_col + 1; col < var_count; ++col)
                    value = (value - A(row, col) * solution[col] % prime + prime) % prime;
                solution[pivot_col] = value;
            }
            component_solutions.push_back(std::move(solution));
        }
        std::vector<int> combined(var_count, 0);
        for (int col = 0; col < var_count; ++col) {
            // CRT idempotents for moduli 3, 5, and 17 are 85, 51, and 120.
            combined[col] = (85 * component_solutions[0][col]
                           + 51 * component_solutions[1][col]
                           + 120 * component_solutions[2][col]) % 255;
        }
        for (int row = 0; row < system.rows(); ++row) {
            long long residual = -system(row, var_count);
            for (int col = 0; col < var_count; ++col)
                residual += static_cast<long long>(system(row, col)) * combined[col];
            if ((residual % mod + mod) % mod != 0)
                throw std::runtime_error("CRT solution failed the modulo-255 constraints");
        }
        return combined;
    }

    Eigen::MatrixXi A = system;
    int rows = A.rows();

    if (rows == 0) {
        std::vector<int> x(var_count);
        std::uniform_int_distribution<> distrib(0, mod - 1);
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
    std::uniform_int_distribution<> distrib(0, mod - 1);
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
    for (int i = 0; i < system.rows(); ++i) {
        long long residual = -system(i, var_count);
        for (int j = 0; j < var_count; ++j) {
            residual += static_cast<long long>(system(i, j)) * x[j];
        }
        if ((residual % mod + mod) % mod != 0) {
            throw std::runtime_error("QC label system has no solution for the sampled variables");
        }
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

std::pair<Eigen::SparseMatrix<int>, Eigen::SparseMatrix<int>>
constructKasaiNonbinaryPair(const Eigen::MatrixXi& hc, const Eigen::MatrixXi& hd,
                           GF2p& gf, std::mt19937& rng) {
    if (hc.cols() != hd.cols()) throw std::invalid_argument("QC support widths differ");
    std::vector<std::pair<int, int>> positions;
    for (int row = 0; row < hc.rows(); ++row)
        for (int col = 0; col < hc.cols(); ++col)
            if (hc(row, col)) positions.emplace_back(row, col);
    const auto system = buildLinearSystem(positions, hc, hd, gf.getMod());
    if (system.rows() != hd.rows())
        throw std::runtime_error("Kasai cycle constraints were not all constructed");
    constexpr int max_label_attempts = 256;
    for (int attempt = 0; attempt < max_label_attempts; ++attempt) {
        const auto exponents = gaussElimination(system, gf.getMod(), rng);
        Eigen::MatrixXi gamma_dense = Eigen::MatrixXi::Zero(hc.rows(), hc.cols());
        std::vector<Eigen::Triplet<int>> gamma_entries, delta_entries;
        for (size_t i = 0; i < positions.size(); ++i) {
            const auto [row, col] = positions[i];
            const int value = gf.exp(exponents[i]);
            gamma_dense(row, col) = value;
            gamma_entries.emplace_back(row, col, value);
        }
        bool support_matches = true;
        for (int row = 0; row < hd.rows() && support_matches; ++row) {
            const auto delta = solveDeltaRow(row, hc, hd, gamma_dense, gf);
            for (int col = 0; col < hd.cols(); ++col) {
                if ((delta[col] != 0) != (hd(row, col) != 0)) {
                    support_matches = false;
                    break;
                }
                if (delta[col]) delta_entries.emplace_back(row, col, delta[col]);
            }
        }
        if (!support_matches) continue;

        Eigen::SparseMatrix<int> gamma(hc.rows(), hc.cols());
        Eigen::SparseMatrix<int> delta(hd.rows(), hd.cols());
        gamma.setFromTriplets(gamma_entries.begin(), gamma_entries.end());
        delta.setFromTriplets(delta_entries.begin(), delta_entries.end());
        return {gamma, delta};
    }
    throw std::runtime_error("Could not find GF labels preserving the requested QC support after "
                             + std::to_string(max_label_attempts) + " attempts");
}

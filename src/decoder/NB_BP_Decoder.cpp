#include "NB_BP_Decoder.h"
#include <iostream>
#include <numeric>
#include <algorithm>
#include <cmath>

NB_BP_Decoder::NB_BP_Decoder(GF2p& _gf, int _max_iters)
    : gf(_gf), q(_gf.getSize()), max_iters(_max_iters) {
    // 修复：通过 getSize 计算扩展度 e (例如 64 对应 6)
    e = static_cast<int>(std::round(std::log2(q)));
    initialize_x_basis_maps();
}

int NB_BP_Decoder::apply_binary_matrix(const std::vector<std::vector<int>>& matrix, int value) const {
    int result = 0;
    for (int i = 0; i < e; ++i) {
        int bit = 0;
        for (int j = 0; j < e; ++j) {
            if (matrix[i][j] && ((value >> j) & 1))
                bit ^= 1;
        }
        if (bit)
            result |= (1 << i);
    }
    return result;
}

void NB_BP_Decoder::initialize_x_basis_maps() {
    std::vector<std::vector<int>> gram(e, std::vector<int>(e, 0));

    auto trace = [&](int value) {
        int acc = 0;
        int cur = value;
        for (int i = 0; i < e; ++i) {
            acc = gf.add(acc, cur);
            cur = gf.mul(cur, cur);
        }
        return acc & 1;
    };

    for (int i = 0; i < e; ++i) {
        for (int j = 0; j < e; ++j)
            gram[i][j] = trace(gf.mul(1 << i, 1 << j));
    }

    std::vector<std::vector<int>> augmented(e, std::vector<int>(2 * e, 0));
    for (int i = 0; i < e; ++i) {
        for (int j = 0; j < e; ++j)
            augmented[i][j] = gram[i][j];
        augmented[i][e + i] = 1;
    }

    int rank = 0;
    for (int col = 0; col < e && rank < e; ++col) {
        int pivot = -1;
        for (int row = rank; row < e; ++row) {
            if (augmented[row][col]) {
                pivot = row;
                break;
            }
        }
        if (pivot == -1)
            continue;

        if (pivot != rank)
            std::swap(augmented[pivot], augmented[rank]);

        for (int row = 0; row < e; ++row) {
            if (row != rank && augmented[row][col]) {
                for (int k = 0; k < 2 * e; ++k)
                    augmented[row][k] ^= augmented[rank][k];
            }
        }
        rank++;
    }

    std::vector<std::vector<int>> dual_to_primal(e, std::vector<int>(e, 0));
    for (int i = 0; i < e; ++i) {
        for (int j = 0; j < e; ++j)
            dual_to_primal[i][j] = augmented[i][e + j];
    }

    x_to_internal_basis.assign(q, 0);
    x_from_internal_basis.assign(q, 0);
    for (int value = 0; value < q; ++value) {
        x_to_internal_basis[value] = apply_binary_matrix(dual_to_primal, value);
        x_from_internal_basis[value] = apply_binary_matrix(gram, value);
    }
}

// ========================
// 物理-逻辑映射辅助函数 (对应论文 Sec 5.4)
// ========================

int NB_BP_Decoder::binToGf(const std::vector<int>& bits, int start_idx) {
    int val = 0;
    for (int i = 0; i < e; ++i) {
        if (bits[start_idx + i]) val |= (1 << i);
    }
    return val;
}

std::vector<int> NB_BP_Decoder::gfToBin(int val) {
    std::vector<int> bits(e);
    for (int i = 0; i < e; ++i) {
        bits[i] = (val >> i) & 1;
    }
    return bits;
}

// ========================
// FWHT 加速算法 (线性复杂度核心)
// ========================

void NB_BP_Decoder::fwht(std::vector<double>& a, bool inverse) {
    int n = a.size();
    for (int len = 1; len < n; len <<= 1) {
        for (int i = 0; i < n; i += 2 * len) {
            for (int j = 0; j < len; j++) {
                double u = a[i + j];
                double v = a[i + len + j];
                a[i + j] = u + v;
                a[i + len + j] = u - v;
            }
        }
    }
    if (inverse) {
        for (double& x : a) x /= n;
    }
}

void NB_BP_Decoder::fwht_2d(JointProb& a, bool inverse) {
    int n = a.size();
    for (int i = 0; i < n; i++) fwht(a[i], inverse);
    for (int j = 0; j < n; j++) {
        std::vector<double> col(n);
        for (int i = 0; i < n; i++) col[i] = a[i][j];
        fwht(col, inverse);
        for (int i = 0; i < n; i++) a[i][j] = col[i];
    }
}

// ========================
// 联合先验概率 (去极化信道采样, 对应论文 Eq 35)
// ========================

JointProb NB_BP_Decoder::create_joint_prior(double p_depolar) {
    JointProb prior(q, std::vector<double>(q, 0.0));
    double total = 0.0;

    for (int x_internal = 0; x_internal < q; ++x_internal) {
        int x_physical = x_from_internal_basis[x_internal];
        for (int z = 0; z < q; ++z) {
            double prob = 1.0;
            for (int bit = 0; bit < e; ++bit) {
                bool has_x = ((x_physical >> bit) & 1) != 0;
                bool has_z = ((z >> bit) & 1) != 0;
                prob *= (!has_x && !has_z) ? (1.0 - p_depolar) : (p_depolar / 3.0);
            }
            prior[x_internal][z] = prob;
            total += prob;
        }
    }

    if (total > 0.0) {
        for (auto& row : prior) {
            for (double& value : row)
                value /= total;
        }
    }
    return prior;
}

// ========================
// 水平步 (校验节点更新, 对应论文 Eq 39-40)
// ========================

void NB_BP_Decoder::horizontal_step_joint(
    const Eigen::SparseMatrix<int, Eigen::RowMajor>& H,
    const std::vector<int>& logical_syn,
    const std::vector<std::map<int, JointProb>>& v2c,
    std::vector<std::map<int, JointProb>>& c2v,
    bool is_Z_channel
) {
    for (int i = 0; i < H.rows(); i++) {
        std::vector<int> row_cols;
        for (Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(H, i); it; ++it)
            row_cols.push_back(it.col());

        std::vector<std::vector<double>> freq_msgs;
        for (int j : row_cols) {
            int hij = H.coeff(i, j);
            std::vector<double> marginal(q, 0.0);
            for (int x = 0; x < q; x++) {
                for (int z = 0; z < q; z++) {
                    int constrained_value = is_Z_channel ? z : x;
                    marginal[constrained_value] += v2c[j].at(i)[x][z];
                }
            }

            std::vector<double> permuted(q, 0.0);
            for (int value = 0; value < q; value++)
                permuted[gf.mul(value, hij)] = marginal[value];

            fwht(permuted, false);
            freq_msgs.push_back(permuted);
        }

        int s_logical = logical_syn[i];
        for (size_t t = 0; t < row_cols.size(); t++) {
            std::vector<double> res_freq(q, 1.0);
            for (size_t s = 0; s < row_cols.size(); s++) {
                if (s == t)
                    continue;
                for (int value = 0; value < q; value++)
                    res_freq[value] *= freq_msgs[s][value];
            }

            fwht(res_freq, true);

            JointProb res(q, std::vector<double>(q, 0.0));
            int ht = H.coeff(i, row_cols[t]);
            std::vector<double> check_msg(q, 0.0);
            for (int value = 0; value < q; value++) {
                int needed_sum = gf.sub(s_logical, gf.mul(ht, value));
                check_msg[value] = std::max(res_freq[needed_sum], 0.0);
            }

            double msg_sum = std::accumulate(check_msg.begin(), check_msg.end(), 0.0);
            if (msg_sum > 0.0) {
                for (double& value : check_msg)
                    value /= msg_sum;
            } else {
                for (double& value : check_msg)
                    value = 1.0 / q;
            }

            for (int x = 0; x < q; x++) {
                for (int z = 0; z < q; z++) {
                    int constrained_value = is_Z_channel ? z : x;
                    res[x][z] = check_msg[constrained_value];
                }
            }

            double sum = 0;
            for(auto& r : res) sum += std::accumulate(r.begin(), r.end(), 0.0);
            for(auto& r : res) for(auto& v : r) v /= (sum + 1e-30);
            c2v[i][row_cols[t]] = res;
        }
    }
}

// ========================
// 垂直步 (变量节点更新, 对应论文 Eq 41-46)
// ========================

void NB_BP_Decoder::vertical_step_joint(
    const Eigen::SparseMatrix<int, Eigen::RowMajor>& HX,
    const Eigen::SparseMatrix<int, Eigen::RowMajor>& HZ,
    const JointProb& joint_prior,
    const std::vector<std::map<int, JointProb>>& c2v_x,
    const std::vector<std::map<int, JointProb>>& c2v_z,
    std::vector<std::map<int, JointProb>>& v2c_x,
    std::vector<std::map<int, JointProb>>& v2c_z,
    std::vector<int>& est_X,
    std::vector<int>& est_Z
) {
    int N = HX.cols();
    for (int j = 0; j < N; j++) {
        JointProb p = joint_prior;
        for (int i = 0; i < HZ.rows(); i++) {
            auto it = c2v_x[i].find(j);
            if (it != c2v_x[i].end()) {
                for (int x = 0; x < q; x++)
                    for (int z = 0; z < q; z++) p[x][z] *= it->second[x][z];
            }
        }
        for (int i = 0; i < HX.rows(); i++) {
            auto it = c2v_z[i].find(j);
            if (it != c2v_z[i].end()) {
                for (int x = 0; x < q; x++)
                    for (int z = 0; z < q; z++) p[x][z] *= it->second[x][z];
            }
        }
        double sum = 0;
        for (auto& row : p) sum += std::accumulate(row.begin(), row.end(), 0.0);
        for (auto& row : p) for (auto& val : row) val /= (sum + 1e-30);

        double max_p = -1;
        for (int x = 0; x < q; x++) {
            for (int z = 0; z < q; z++) {
                if (p[x][z] > max_p) { max_p = p[x][z]; est_X[j] = x; est_Z[j] = z; }
            }
        }
        for (auto& [i, msg] : v2c_x[j]) {
            for (int x = 0; x < q; x++)
                for (int z = 0; z < q; z++) msg[x][z] = p[x][z] / (c2v_x[i].at(j)[x][z] + 1e-30);
        }
        for (auto& [i, msg] : v2c_z[j]) {
            for (int x = 0; x < q; x++)
                for (int z = 0; z < q; z++) msg[x][z] = p[x][z] / (c2v_z[i].at(j)[x][z] + 1e-30);
        }
    }
}

// ========================
// 联合译码主入口
// ========================

std::pair<std::vector<int>, std::vector<int>>
NB_BP_Decoder::joint_decode(
    const Eigen::SparseMatrix<int, Eigen::RowMajor>& HX_nb,
    const Eigen::SparseMatrix<int, Eigen::RowMajor>& HZ_nb,
    const std::vector<int>& syn_x_bin,
    const std::vector<int>& syn_z_bin,
    double p_depolar
) {
    int N = HX_nb.cols();
    auto prior = create_joint_prior(p_depolar);

    // 1. 将二进制综合征映射为逻辑元素
    std::vector<int> logical_syn_x(HZ_nb.rows());
    std::vector<int> logical_syn_z(HX_nb.rows());
    for(int i=0; i<HZ_nb.rows(); i++) logical_syn_x[i] = x_to_internal_basis[binToGf(syn_x_bin, i * e)];
    for(int i=0; i<HX_nb.rows(); i++) logical_syn_z[i] = binToGf(syn_z_bin, i * e);

    // 2. 初始化消息容器 (修复未声明报错)
    std::vector<std::map<int, JointProb>> v2c_x(N), c2v_x(HZ_nb.rows());
    std::vector<std::map<int, JointProb>> v2c_z(N), c2v_z(HX_nb.rows());

    for (int i = 0; i < HZ_nb.rows(); i++) {
        for (Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(HZ_nb, i); it; ++it) {
            v2c_x[it.col()][i] = prior;
            c2v_x[i][it.col()] = JointProb(q, std::vector<double>(q, 1.0/(q*q)));
        }
    }
    for (int i = 0; i < HX_nb.rows(); i++) {
        for (Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(HX_nb, i); it; ++it) {
            v2c_z[it.col()][i] = prior;
            c2v_z[i][it.col()] = JointProb(q, std::vector<double>(q, 1.0/(q*q)));
        }
    }

    std::vector<int> est_X(N, 0), est_Z(N, 0);
    for (int iter = 0; iter < max_iters; iter++) {
        horizontal_step_joint(HZ_nb, logical_syn_x, v2c_x, c2v_x, false);
        horizontal_step_joint(HX_nb, logical_syn_z, v2c_z, c2v_z, true);
        vertical_step_joint(HX_nb, HZ_nb, prior, c2v_x, c2v_z, v2c_x, v2c_z, est_X, est_Z);

        if (check_joint_syndrome(HX_nb, HZ_nb, est_X, est_Z, logical_syn_x, logical_syn_z)) break;
    }
    std::vector<int> est_X_physical(N, 0);
    for (int j = 0; j < N; ++j)
        est_X_physical[j] = x_from_internal_basis[est_X[j]];

    return {est_X_physical, est_Z};
}

bool NB_BP_Decoder::check_joint_syndrome(
    const Eigen::SparseMatrix<int, Eigen::RowMajor>& HX_nb,
    const Eigen::SparseMatrix<int, Eigen::RowMajor>& HZ_nb,
    const std::vector<int>& eX_gf,
    const std::vector<int>& eZ_gf,
    const std::vector<int>& l_syn_x,
    const std::vector<int>& l_syn_z
) {
    for (int i = 0; i < HZ_nb.rows(); i++) {
        int sum = 0;
        for (Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(HZ_nb, i); it; ++it)
            sum = gf.add(sum, gf.mul(it.value(), eX_gf[it.col()]));
        if (sum != l_syn_x[i]) return false;
    }
    for (int i = 0; i < HX_nb.rows(); i++) {
        int sum = 0;
        for (Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(HX_nb, i); it; ++it)
            sum = gf.add(sum, gf.mul(it.value(), eZ_gf[it.col()]));
        if (sum != l_syn_z[i]) return false;
    }
    return true;
}

std::vector<int> NB_BP_Decoder::syndrome_decode(
    const Eigen::SparseMatrix<int, Eigen::RowMajor>& H_nb,
    const std::vector<int>& syn_bin,
    double p_depolar,
    bool is_x_error
) {
    int N = H_nb.cols();
    int M = H_nb.rows();

    std::vector<int> logical_syn(M, 0);
    for (int i = 0; i < M; ++i) {
        int raw = binToGf(syn_bin, i * e);
        logical_syn[i] = is_x_error ? x_to_internal_basis[raw] : raw;
    }

    auto bit_weight = [](int value) {
        int weight = 0;
        while (value) {
            weight += value & 1;
            value >>= 1;
        }
        return weight;
    };

    double bit_error = std::min(std::max(2.0 * p_depolar / 3.0, 1e-12), 1.0 - 1e-12);
    std::vector<double> prior(q, 0.0);
    double prior_sum = 0.0;
    for (int value = 0; value < q; ++value) {
        int physical_value = is_x_error ? x_from_internal_basis[value] : value;
        int weight = bit_weight(physical_value);
        double prob = std::pow(bit_error, weight) * std::pow(1.0 - bit_error, e - weight);
        prior[value] = prob;
        prior_sum += prob;
    }
    for (double& value : prior)
        value /= (prior_sum + 1e-30);

    std::vector<std::map<int, std::vector<double>>> v2c(N), c2v(M);
    for (int i = 0; i < M; ++i) {
        for (Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(H_nb, i); it; ++it) {
            v2c[it.col()][i] = prior;
            c2v[i][it.col()] = std::vector<double>(q, 1.0 / q);
        }
    }

    std::vector<int> est(N, 0);
    auto syndrome_matches = [&]() {
        for (int i = 0; i < M; ++i) {
            int sum = 0;
            for (Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(H_nb, i); it; ++it)
                sum = gf.add(sum, gf.mul(it.value(), est[it.col()]));
            if (sum != logical_syn[i])
                return false;
        }
        return true;
    };

    for (int iter = 0; iter < max_iters; ++iter) {
        for (int i = 0; i < M; ++i) {
            std::vector<int> row_cols;
            for (Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(H_nb, i); it; ++it)
                row_cols.push_back(it.col());

            std::vector<std::vector<double>> freq_msgs;
            freq_msgs.reserve(row_cols.size());
            for (int j : row_cols) {
                int hij = H_nb.coeff(i, j);
                std::vector<double> permuted(q, 0.0);
                for (int value = 0; value < q; ++value)
                    permuted[gf.mul(value, hij)] = v2c[j].at(i)[value];
                fwht(permuted, false);
                freq_msgs.push_back(permuted);
            }

            for (size_t t = 0; t < row_cols.size(); ++t) {
                std::vector<double> res_freq(q, 1.0);
                for (size_t s = 0; s < row_cols.size(); ++s) {
                    if (s == t)
                        continue;
                    for (int value = 0; value < q; ++value)
                        res_freq[value] *= freq_msgs[s][value];
                }
                fwht(res_freq, true);

                int ht = H_nb.coeff(i, row_cols[t]);
                std::vector<double> check_msg(q, 0.0);
                for (int value = 0; value < q; ++value) {
                    int needed_sum = gf.sub(logical_syn[i], gf.mul(ht, value));
                    check_msg[value] = std::max(res_freq[needed_sum], 0.0);
                }

                double msg_sum = std::accumulate(check_msg.begin(), check_msg.end(), 0.0);
                if (msg_sum > 0.0) {
                    for (double& value : check_msg)
                        value /= msg_sum;
                } else {
                    for (double& value : check_msg)
                        value = 1.0 / q;
                }
                c2v[i][row_cols[t]] = check_msg;
            }
        }

        for (int j = 0; j < N; ++j) {
            std::vector<double> posterior = prior;
            for (const auto& edge : v2c[j]) {
                int i = edge.first;
                auto msg_it = c2v[i].find(j);
                if (msg_it == c2v[i].end())
                    continue;
                for (int value = 0; value < q; ++value)
                    posterior[value] *= msg_it->second[value];
            }

            double post_sum = std::accumulate(posterior.begin(), posterior.end(), 0.0);
            for (double& value : posterior)
                value /= (post_sum + 1e-30);

            est[j] = static_cast<int>(
                std::max_element(posterior.begin(), posterior.end()) - posterior.begin());

            for (auto& edge : v2c[j]) {
                int i = edge.first;
                std::vector<double> extrinsic(q, 0.0);
                for (int value = 0; value < q; ++value)
                    extrinsic[value] = posterior[value] / (c2v[i].at(j)[value] + 1e-30);

                double ext_sum = std::accumulate(extrinsic.begin(), extrinsic.end(), 0.0);
                if (ext_sum > 0.0) {
                    for (double& value : extrinsic)
                        value /= ext_sum;
                } else {
                    for (double& value : extrinsic)
                        value = 1.0 / q;
                }
                edge.second = extrinsic;
            }
        }

        if (syndrome_matches())
            break;
    }

    if (is_x_error) {
        for (int& value : est)
            value = x_from_internal_basis[value];
    }
    return est;
}

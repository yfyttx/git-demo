#include "NB_BP_Decoder.h"
#include <iostream>
#include <numeric>
#include <algorithm>
#include <cmath>

NB_BP_Decoder::NB_BP_Decoder(GF2p& _gf, int _max_iters)
    : gf(_gf), q(_gf.getSize()), max_iters(_max_iters) {
    // 修复：通过 getSize 计算扩展度 e (例如 64 对应 6)
    e = static_cast<int>(std::round(std::log2(q)));
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
    double p = p_depolar / 3.0;
    prior[0][0] = 1.0 - p_depolar; 
    for (int i = 1; i < q; i++) {
        prior[i][0] = p; // 纯 X 错误
        prior[0][i] = p; // 纯 Z 错误
    }
    for (int i = 1; i < q; i++)
        for (int j = 1; j < q; j++) prior[i][j] = p; // Y 错误
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

        std::vector<JointProb> freq_msgs;
        for (int j : row_cols) {
            int hij = H.coeff(i, j);
            JointProb permuted(q, std::vector<double>(q, 0.0));
            for(int x=0; x<q; x++) {
                for(int z=0; z<q; z++) {
                    permuted[gf.mul(x, hij)][gf.mul(z, hij)] = v2c[j].at(i)[x][z];
                }
            }
            fwht_2d(permuted, false);
            freq_msgs.push_back(permuted);
        }

        int s_logical = logical_syn[i];
        for (size_t t = 0; t < row_cols.size(); t++) {
            JointProb res_freq(q, std::vector<double>(q, 1.0));
            for (size_t s = 0; s < row_cols.size(); s++) {
                if (s == t) continue;
                for(int x=0; x<q; x++)
                    for(int z=0; z<q; z++) res_freq[x][z] *= freq_msgs[s][x][z];
            }
            fwht_2d(res_freq, true);

            JointProb res(q, std::vector<double>(q, 0.0));
            int ht = H.coeff(i, row_cols[t]);
            for(int x=0; x<q; x++) {
                for(int z=0; z<q; z++) {
                    int val = is_Z_channel ? z : x;
                    int shifted_val = gf.sub(val, s_logical);
                    int ox = is_Z_channel ? x : shifted_val;
                    int oz = is_Z_channel ? shifted_val : z;
                    res[x][z] = res_freq[gf.mul(ox, ht)][gf.mul(oz, ht)];
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
    for(int i=0; i<HZ_nb.rows(); i++) logical_syn_x[i] = binToGf(syn_x_bin, i * e);
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
    return {est_X, est_Z};
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
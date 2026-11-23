#include "NB_BP_Decoder.h"
#include <iostream>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <stdexcept>

NB_BP_Decoder::NB_BP_Decoder(GF2p& _gf, int _max_iters)
    : gf(_gf), q(_gf.getSize()), max_iters(_max_iters), basis(_gf) {
    e = static_cast<int>(std::round(std::log2(q)));
    if ((1 << e) != q) {
        throw std::invalid_argument("The decoder requires a power-of-two field size");
    }
}

// ========================
// 物理-逻辑映射辅助函数 (对应论文 Sec 5.4)
// ========================

int NB_BP_Decoder::binToGf(const std::vector<int>& bits, int start_idx) {
    if (start_idx < 0 || start_idx + e > static_cast<int>(bits.size())) {
        throw std::invalid_argument("Binary syndrome does not contain a complete GF symbol");
    }
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

// ========================
// 联合先验概率：每个物理位独立采样 I/X/Z/Y
// ========================

JointProb NB_BP_Decoder::create_joint_prior(double p_depolar) {
    if (p_depolar < 0.0 || p_depolar > 1.0) {
        throw std::invalid_argument("Depolarizing probability must be between 0 and 1");
    }
    JointProb prior(q, std::vector<double>(q, 0.0));
    for (int x = 0; x < q; ++x) {
        const int physical_x = basis.toDual(x);
        for (int z = 0; z < q; ++z) {
            double probability = 1.0;
            for (int bit = 0; bit < e; ++bit) {
                const bool has_x = (physical_x >> bit) & 1;
                const bool has_z = (z >> bit) & 1;
                probability *= (has_x || has_z) ? p_depolar / 3.0 : 1.0 - p_depolar;
            }
            prior[x][z] = probability;
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
    if (H.rows() != static_cast<int>(logical_syn.size())) {
        throw std::invalid_argument("Syndrome length does not match parity-check rows");
    }
    for (int i = 0; i < H.rows(); i++) {
        std::vector<int> row_cols;
        for (Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(H, i); it; ++it)
            row_cols.push_back(it.col());

        // A check constrains only X or only Z. Marginalize the other component
        // before convolution, then use a one-dimensional Walsh transform.
        std::vector<std::vector<double>> freq_msgs;
        for (int j : row_cols) {
            int hij = H.coeff(i, j);
            std::vector<double> permuted(q, 0.0);
            const JointProb& incoming = v2c[j].at(i);
            for (int x = 0; x < q; ++x) {
                for (int z = 0; z < q; ++z) {
                    const int constrained_value = is_Z_channel ? z : x;
                    permuted[gf.mul(constrained_value, hij)] += incoming[x][z];
                }
            }
            fwht(permuted, false);
            freq_msgs.push_back(permuted);
        }

        const std::vector<double> identity(q, 1.0);
        std::vector<std::vector<double>> prefix(row_cols.size() + 1, identity);
        std::vector<std::vector<double>> suffix(row_cols.size() + 1, identity);
        for (size_t t = 0; t < row_cols.size(); ++t) {
            for (int value = 0; value < q; ++value) {
                prefix[t + 1][value] = prefix[t][value] * freq_msgs[t][value];
            }
        }
        for (size_t t = row_cols.size(); t > 0; --t) {
            for (int value = 0; value < q; ++value) {
                suffix[t - 1][value] = suffix[t][value] * freq_msgs[t - 1][value];
            }
        }

        int s_logical = logical_syn[i];
        for (size_t t = 0; t < row_cols.size(); t++) {
            std::vector<double> res_freq(q);
            for (int value = 0; value < q; ++value) {
                res_freq[value] = prefix[t][value] * suffix[t + 1][value];
            }
            fwht(res_freq, true);

            int ht = H.coeff(i, row_cols[t]);
            std::vector<double> constrained(q);
            double sum = 0.0;
            for (int value = 0; value < q; ++value) {
                const int target = gf.sub(s_logical, gf.mul(value, ht));
                constrained[value] = std::max(0.0, res_freq[target]);
                sum += constrained[value];
            }
            if (sum <= 0.0) {
                throw std::runtime_error("Check-node message has zero probability mass");
            }
            JointProb res(q, std::vector<double>(q));
            for (int x = 0; x < q; ++x) {
                for (int z = 0; z < q; ++z) {
                    res[x][z] = constrained[is_Z_channel ? z : x] / (sum * q);
                }
            }
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
    const auto multiply = [this](JointProb& result, const JointProb& message) {
        // Renormalize after each factor to prevent long products from
        // underflowing to an all-zero message. A tiny floor also preserves
        // support when two approximate BP messages become numerically disjoint.
        constexpr double probability_floor = 1e-300;
        double sum = 0.0;
        for (int x = 0; x < q; ++x) {
            for (int z = 0; z < q; ++z) {
                result[x][z] = std::max(probability_floor, result[x][z] * message[x][z]);
                sum += result[x][z];
            }
        }
        for (int x = 0; x < q; ++x)
            for (int z = 0; z < q; ++z) result[x][z] /= sum;
    };
    const auto normalize = [](JointProb& probabilities) {
        double sum = 0.0;
        for (const auto& row : probabilities)
            sum += std::accumulate(row.begin(), row.end(), 0.0);
        if (sum <= 0.0 || !std::isfinite(sum)) {
            throw std::runtime_error("Variable-node message has zero or invalid probability mass");
        }
        for (auto& row : probabilities)
            for (double& value : row) value /= sum;
    };
    for (int j = 0; j < N; j++) {
        std::vector<std::pair<int, const JointProb*>> x_messages, z_messages;
        for (const auto& [i, outgoing] : v2c_x[j]) {
            (void)outgoing;
            x_messages.emplace_back(i, &c2v_x[i].at(j));
        }
        for (const auto& [i, outgoing] : v2c_z[j]) {
            (void)outgoing;
            z_messages.emplace_back(i, &c2v_z[i].at(j));
        }

        JointProb p = joint_prior;
        for (const auto& [i, message] : x_messages) multiply(p, *message);
        for (const auto& [i, message] : z_messages) multiply(p, *message);
        normalize(p);

        double max_p = -1;
        for (int x = 0; x < q; x++) {
            for (int z = 0; z < q; z++) {
                if (p[x][z] > max_p) { max_p = p[x][z]; est_X[j] = x; est_Z[j] = z; }
            }
        }
        for (auto& [target, outgoing] : v2c_x[j]) {
            outgoing = joint_prior;
            for (const auto& [i, message] : x_messages)
                if (i != target) multiply(outgoing, *message);
            for (const auto& [i, message] : z_messages) multiply(outgoing, *message);
            normalize(outgoing);
        }
        for (auto& [target, outgoing] : v2c_z[j]) {
            outgoing = joint_prior;
            for (const auto& [i, message] : x_messages) multiply(outgoing, *message);
            for (const auto& [i, message] : z_messages)
                if (i != target) multiply(outgoing, *message);
            normalize(outgoing);
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
    if (HX_nb.cols() != HZ_nb.cols()) {
        throw std::invalid_argument("HX and HZ must have the same number of columns");
    }
    if (syn_x_bin.size() != static_cast<size_t>(HZ_nb.rows() * e) ||
        syn_z_bin.size() != static_cast<size_t>(HX_nb.rows() * e)) {
        throw std::invalid_argument("Binary syndrome lengths must equal check rows times field degree");
    }
    int N = HX_nb.cols();
    auto prior = create_joint_prior(p_depolar);

    // 1. 将二进制综合征映射为逻辑元素
    std::vector<int> logical_syn_x(HZ_nb.rows());
    std::vector<int> logical_syn_z(HX_nb.rows());
    for(int i=0; i<HZ_nb.rows(); i++) logical_syn_x[i] = basis.fromDual(binToGf(syn_x_bin, i * e));
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

bool NB_BP_Decoder::recover_cascade_trailing_x(
    const Eigen::SparseMatrix<int, Eigen::RowMajor>& outer_hz_check,
    int inner_length,
    const std::vector<int>& logical_syn_x,
    std::vector<int>& estimate_x
) {
    if (inner_length < 2 || estimate_x.size() !=
            static_cast<size_t>(outer_hz_check.cols() * inner_length)) {
        throw std::invalid_argument("Invalid dimensions for cascade trailing-symbol recovery");
    }
    const int k = inner_length - 1;
    const int upper_rows = outer_hz_check.rows() * k;
    if (logical_syn_x.size() != static_cast<size_t>(upper_rows + outer_hz_check.cols())) {
        throw std::invalid_argument("Cascade syndrome must contain upper and lower blocks");
    }

    // The selector [I_k | 0] makes the upper checks independent across
    // inner coordinates. Confirm those decoded coordinates before solving
    // the final symbol from each all-ones inner check.
    for (int check = 0; check < outer_hz_check.rows(); ++check) {
        for (int coordinate = 0; coordinate < k; ++coordinate) {
            int sum = 0;
            for (Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(outer_hz_check, check);
                 it; ++it) {
                const int variable = it.col() * inner_length + coordinate;
                sum = gf.add(sum, gf.mul(it.value(), estimate_x[variable]));
            }
            if (sum != logical_syn_x[check * k + coordinate]) return false;
        }
    }

    for (int outer_variable = 0; outer_variable < outer_hz_check.cols(); ++outer_variable) {
        int sum = 0;
        for (int coordinate = 0; coordinate < k; ++coordinate) {
            sum = gf.add(sum, estimate_x[outer_variable * inner_length + coordinate]);
        }
        const int lower_syndrome = logical_syn_x[upper_rows + outer_variable];
        estimate_x[outer_variable * inner_length + k] = gf.sub(lower_syndrome, sum);
    }
    return true;
}

std::vector<int> NB_BP_Decoder::decode_nonbinary_symbols(
    const Eigen::SparseMatrix<int, Eigen::RowMajor>& H,
    const std::vector<int>& syndrome,
    double p_bit,
    bool trace_dual_coordinates
) {
    if (H.rows() != static_cast<int>(syndrome.size()) || p_bit < 0.0 || p_bit > 1.0) {
        throw std::invalid_argument("Invalid dimensions or channel probability for NB-BP");
    }
    std::vector<double> channel_prior(q, 0.0);
    for (int symbol = 0; symbol < q; ++symbol) {
        const int bit_pattern = trace_dual_coordinates ? basis.toDual(symbol) : symbol;
        double probability = 1.0;
        for (int bit = 0; bit < e; ++bit) {
            probability *= ((bit_pattern >> bit) & 1) ? p_bit : (1.0 - p_bit);
        }
        channel_prior[symbol] = probability;
    }

    return decode_component_prior(H, syndrome, channel_prior);
}

std::vector<int> NB_BP_Decoder::decode_component_prior(
    const Eigen::SparseMatrix<int, Eigen::RowMajor>& H,
    const std::vector<int>& syndrome,
    const std::vector<double>& channel_prior,
    int* iterations_used
) {
    if (H.rows() != static_cast<int>(syndrome.size()) ||
        channel_prior.size() != static_cast<size_t>(q))
        throw std::invalid_argument("Invalid explicit-prior BP dimensions");
    double total = 0;
    for (double probability : channel_prior) {
        if (!std::isfinite(probability) || probability < 0)
            throw std::invalid_argument("Invalid symbol probability");
        total += probability;
    }
    if (!(total > 0)) throw std::invalid_argument("Empty symbol prior");
    const int variables = H.cols();
    if (iterations_used) *iterations_used = 0;
    std::vector<std::vector<int>> row_cols(H.rows());
    std::vector<std::vector<std::pair<int, int>>> variable_edges(variables);
    std::vector<std::vector<std::vector<double>>> v2c(H.rows()), c2v(H.rows());
    for (int row = 0; row < H.rows(); ++row) {
        for (Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(H, row); it; ++it) {
            const int edge = static_cast<int>(row_cols[row].size());
            row_cols[row].push_back(it.col());
            variable_edges[it.col()].emplace_back(row, edge);
            v2c[row].emplace_back(channel_prior);
            c2v[row].emplace_back(q, 1.0 / q);
        }
    }

    const auto normalize = [](std::vector<double>& probabilities) {
        constexpr double floor = 1e-300;
        double sum = 0.0;
        for (double& value : probabilities) {
            value = std::max(floor, value);
            sum += value;
        }
        for (double& value : probabilities) value /= sum;
    };
    std::vector<int> estimate(variables, 0);
    for (int iteration = 0; iteration < max_iters; ++iteration) {
        for (int row = 0; row < H.rows(); ++row) {
            const size_t degree = row_cols[row].size();
            if (degree == 0) continue;
            std::vector<std::vector<double>> transformed(degree, std::vector<double>(q));
            for (size_t edge = 0; edge < degree; ++edge) {
                const int coefficient = H.coeff(row, row_cols[row][edge]);
                std::vector<double> scaled(q, 0.0);
                for (int symbol = 0; symbol < q; ++symbol) {
                    scaled[gf.mul(coefficient, symbol)] += v2c[row][edge][symbol];
                }
                fwht(scaled, false);
                transformed[edge] = std::move(scaled);
            }

            const std::vector<double> identity(q, 1.0);
            std::vector<std::vector<double>> prefix(degree + 1, identity), suffix(degree + 1, identity);
            for (size_t edge = 0; edge < degree; ++edge)
                for (int value = 0; value < q; ++value)
                    prefix[edge + 1][value] = prefix[edge][value] * transformed[edge][value];
            for (size_t edge = degree; edge > 0; --edge)
                for (int value = 0; value < q; ++value)
                    suffix[edge - 1][value] = suffix[edge][value] * transformed[edge - 1][value];

            for (size_t edge = 0; edge < degree; ++edge) {
                std::vector<double> other_sum(q);
                for (int value = 0; value < q; ++value)
                    other_sum[value] = prefix[edge][value] * suffix[edge + 1][value];
                fwht(other_sum, true);
                const int variable = row_cols[row][edge];
                const int coefficient = H.coeff(row, variable);
                std::vector<double> outgoing(q, 0.0);
                for (int symbol = 0; symbol < q; ++symbol) {
                    const int target = gf.sub(syndrome[row], gf.mul(coefficient, symbol));
                    outgoing[symbol] = std::max(0.0, other_sum[target]);
                }
                normalize(outgoing);
                c2v[row][edge] = std::move(outgoing);
            }
        }

        for (int variable = 0; variable < variables; ++variable) {
            std::vector<double> posterior = channel_prior;
            for (const auto& [row, edge] : variable_edges[variable]) {
                for (int symbol = 0; symbol < q; ++symbol)
                    posterior[symbol] *= c2v[row][edge][symbol];
                normalize(posterior);
            }
            estimate[variable] = static_cast<int>(
                std::distance(posterior.begin(), std::max_element(posterior.begin(), posterior.end())));

            for (const auto& [target_row, target_edge] : variable_edges[variable]) {
                std::vector<double> outgoing = channel_prior;
                for (const auto& [row, edge] : variable_edges[variable]) {
                    if (row == target_row && edge == target_edge) continue;
                    for (int symbol = 0; symbol < q; ++symbol)
                        outgoing[symbol] *= c2v[row][edge][symbol];
                    normalize(outgoing);
                }
                v2c[target_row][target_edge] = std::move(outgoing);
            }
        }

        bool syndrome_ok = true;
        for (int row = 0; row < H.rows() && syndrome_ok; ++row) {
            int sum = 0;
            for (Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(H, row); it; ++it)
                sum = gf.add(sum, gf.mul(it.value(), estimate[it.col()]));
            syndrome_ok = sum == syndrome[row];
        }
        if (iterations_used) *iterations_used = iteration + 1;
        if (syndrome_ok) break;
    }
    return estimate;
}

std::vector<int> NB_BP_Decoder::decode_component(
    const Eigen::SparseMatrix<int, Eigen::RowMajor>& H,
    const std::vector<int>& syndrome,
    double p_bit,
    bool trace_dual_coordinates
) {
    return decode_nonbinary_symbols(H, syndrome, p_bit, trace_dual_coordinates);
}

bool NB_BP_Decoder::decode_cascade_x(
    const Eigen::SparseMatrix<int, Eigen::RowMajor>& outer_hz_check,
    const Eigen::SparseMatrix<int, Eigen::RowMajor>& hz_cascade,
    int inner_length,
    const std::vector<int>& logical_syn_x,
    double p_x_bit,
    std::vector<int>& estimate_x
) {
    if (inner_length < 2 || p_x_bit < 0.0 || p_x_bit > 1.0) {
        throw std::invalid_argument("Invalid inner length or channel probability for cascade decoder");
    }
    const int n1 = outer_hz_check.cols();
    if (n1 <= 0 || hz_cascade.cols() != n1 * inner_length ||
        hz_cascade.rows() < n1) {
        throw std::invalid_argument("Invalid cascade matrix dimensions");
    }
    // HZ has m1*k2 upper rows and n1*(n2-k2) lower rows. Infer the
    // number of inner checks from those dimensions instead of assuming one.
    const int m1 = outer_hz_check.rows();
    const int denominator = n1 - m1;
    const int numerator = hz_cascade.rows() - m1 * inner_length;
    if (denominator <= 0 || numerator <= 0 || numerator % denominator != 0)
        throw std::invalid_argument("Cannot infer inner parity-check count");
    const int inner_checks = numerator / denominator;
    if (inner_checks < 1 || inner_checks >= inner_length)
        throw std::invalid_argument("Invalid inner parity-check count");
    const int k = inner_length - inner_checks;
    const int upper_rows = outer_hz_check.rows() * k;
    if (logical_syn_x.size() != static_cast<size_t>(upper_rows + n1 * inner_checks) ||
        hz_cascade.rows() != upper_rows + n1 * inner_checks) {
        throw std::invalid_argument("Cascade syndrome must contain upper and lower blocks");
    }

    // Decode each nonzero upper-selector coordinate using only the outer
    // checks. The trailing coordinate is absent from [I_(s-1) | 0], so solve
    // it from the corresponding lower all-ones inner-check equation.
    const double x_bit_probability = p_x_bit;
    std::vector<int> staged_estimate(
        static_cast<size_t>(outer_hz_check.cols() * inner_length), 0);
    for (int coordinate = 0; coordinate < k; ++coordinate) {
        std::vector<int> coordinate_syndrome(outer_hz_check.rows());
        for (int check = 0; check < outer_hz_check.rows(); ++check)
            coordinate_syndrome[check] = logical_syn_x[check * k + coordinate];
        const std::vector<int> decoded_coordinate = decode_nonbinary_symbols(
            outer_hz_check, coordinate_syndrome, x_bit_probability, true);
        for (int variable = 0; variable < outer_hz_check.cols(); ++variable)
            staged_estimate[variable * inner_length + coordinate] =
                decoded_coordinate[variable];
    }

    for (int variable = 0; variable < n1; ++variable) {
        for (int check = 0; check < inner_checks; ++check) {
            const int row = upper_rows + variable * inner_checks + check;
            const int parity_col = variable * inner_length + k + check;
            const int parity_coefficient = hz_cascade.coeff(row, parity_col);
            if (parity_coefficient == 0)
                throw std::invalid_argument("H2 parity block is not systematic");
            int known_sum = 0;
            for (Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(
                     hz_cascade, row); it; ++it) {
                if (it.col() == parity_col) continue;
                if (it.col() >= variable * inner_length + k)
                    throw std::invalid_argument("Coupled parity block is unsupported");
                known_sum = gf.add(known_sum,
                    gf.mul(it.value(), staged_estimate[it.col()]));
            }
            staged_estimate[parity_col] = gf.mul(
                gf.sub(logical_syn_x[row], known_sum), gf.inv(parity_coefficient));
        }
    }
    if (!check_x_syndrome(hz_cascade, staged_estimate, logical_syn_x)) return false;
    estimate_x = std::move(staged_estimate);
    return true;
}

bool NB_BP_Decoder::check_joint_syndrome(
    const Eigen::SparseMatrix<int, Eigen::RowMajor>& HX_nb,
    const Eigen::SparseMatrix<int, Eigen::RowMajor>& HZ_nb,
    const std::vector<int>& eX_gf,
    const std::vector<int>& eZ_gf,
    const std::vector<int>& l_syn_x,
    const std::vector<int>& l_syn_z
) {
    return check_x_syndrome(HZ_nb, eX_gf, l_syn_x) &&
           check_z_syndrome(HX_nb, eZ_gf, l_syn_z);
}

bool NB_BP_Decoder::check_x_syndrome(
    const Eigen::SparseMatrix<int, Eigen::RowMajor>& HZ_nb,
    const std::vector<int>& eX_gf,
    const std::vector<int>& l_syn_x
) {
    if (HZ_nb.cols() != static_cast<int>(eX_gf.size()) ||
        HZ_nb.rows() != static_cast<int>(l_syn_x.size())) return false;
    for (int i = 0; i < HZ_nb.rows(); i++) {
        int sum = 0;
        for (Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(HZ_nb, i); it; ++it)
            sum = gf.add(sum, gf.mul(it.value(), eX_gf[it.col()]));
        if (sum != l_syn_x[i]) return false;
    }
    return true;
}

bool NB_BP_Decoder::check_z_syndrome(
    const Eigen::SparseMatrix<int, Eigen::RowMajor>& HX_nb,
    const std::vector<int>& eZ_gf,
    const std::vector<int>& l_syn_z
) {
    if (HX_nb.cols() != static_cast<int>(eZ_gf.size()) ||
        HX_nb.rows() != static_cast<int>(l_syn_z.size())) return false;
    for (int i = 0; i < HX_nb.rows(); i++) {
        int sum = 0;
        for (Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(HX_nb, i); it; ++it)
            sum = gf.add(sum, gf.mul(it.value(), eZ_gf[it.col()]));
        if (sum != l_syn_z[i]) return false;
    }
    return true;
}

#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "NB_BP_Decoder.h"
#include "gf2p/gf2p.h"

using namespace std;
using namespace Eigen;

const string MATRIX_DIR = "matrices";

string matrixPath(const string& filename) {
    return MATRIX_DIR + "/" + filename;
}

bool fileExists(const string& filename) {
    ifstream in(filename);
    return in.good();
}

string legacySuffix(const string& suffix) {
    size_t m_pos = suffix.find("_m");
    size_t n_pos = suffix.find("_n", m_pos == string::npos ? 0 : m_pos);
    if (m_pos == string::npos || n_pos == string::npos)
        return suffix;
    return suffix.substr(0, m_pos) + suffix.substr(n_pos);
}

string matrixPathFor(const string& base, const string& suffix) {
    string current = matrixPath(base + suffix + ".txt");
    if (fileExists(current))
        return current;

    string legacy = legacySuffix(suffix);
    if (legacy != suffix) {
        string old_path = matrixPath(base + legacy + ".txt");
        if (fileExists(old_path)) {
            cerr << "Warning: using legacy matrix filename " << old_path
                 << " for suffix " << suffix << endl;
            return old_path;
        }
    }

    return current;
}

SparseMatrix<int, RowMajor> loadBinaryMatrix(const string& filename) {
    ifstream in(filename);
    if (!in.is_open()) {
        cerr << "Error: Could not open " << filename << endl;
        exit(1);
    }

    int N, M, Q;
    in >> N >> M >> Q;
    string line;
    getline(in, line);
    getline(in, line);
    getline(in, line);

    vector<Triplet<int>> triplets;
    for (int i = 0; i < M; ++i) {
        getline(in, line);
        stringstream ss(line);
        int col_idx;
        int val;
        while (ss >> col_idx >> val) {
            if (val != 0)
                triplets.emplace_back(i, col_idx - 1, 1);
        }
    }

    SparseMatrix<int, RowMajor> H(M, N);
    H.setFromTriplets(triplets.begin(), triplets.end());
    return H;
}

SparseMatrix<int, RowMajor> loadMatrix(const string& filename, GF2p& gf) {
    ifstream in(filename);
    if (!in.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        exit(1);
    }

    int N, M, Q;
    in >> N >> M >> Q;

    string line;
    getline(in, line);
    getline(in, line);
    getline(in, line);

    vector<Triplet<int>> triplets;
    for (int i = 0; i < M; ++i) {
        getline(in, line);
        stringstream ss(line);
        int col_idx;
        int log_val;
        while (ss >> col_idx >> log_val)
            triplets.emplace_back(i, col_idx - 1, gf.exp(log_val));
    }

    SparseMatrix<int, RowMajor> H(M, N);
    H.setFromTriplets(triplets.begin(), triplets.end());
    return H;
}

vector<int> get_binary_syndrome(const SparseMatrix<int, RowMajor>& H,
                                const vector<int>& error_bin) {
    vector<int> syn(H.rows(), 0);
    for (int i = 0; i < H.rows(); ++i) {
        int sum_v = 0;
        for (SparseMatrix<int, RowMajor>::InnerIterator it(H, i); it; ++it)
            sum_v ^= error_bin[it.col()];
        syn[i] = sum_v;
    }
    return syn;
}

using BitRow = vector<uint64_t>;

bool bit_is_set(const BitRow& row, int idx) {
    return ((row[idx >> 6] >> (idx & 63)) & 1ULL) != 0;
}

void bit_flip(BitRow& row, int idx) {
    row[idx >> 6] ^= (1ULL << (idx & 63));
}

bool bit_is_zero(const BitRow& row) {
    for (uint64_t word : row) {
        if (word != 0)
            return false;
    }
    return true;
}

void bit_xor_inplace(BitRow& dst, const BitRow& src) {
    for (size_t i = 0; i < dst.size(); ++i)
        dst[i] ^= src[i];
}

int first_set_bit(const BitRow& row, int n_bits) {
    for (int i = 0; i < n_bits; ++i) {
        if (bit_is_set(row, i))
            return i;
    }
    return -1;
}

BitRow vector_to_bits(const vector<int>& vec) {
    BitRow bits((vec.size() + 63) / 64, 0);
    for (int i = 0; i < static_cast<int>(vec.size()); ++i) {
        if (vec[i] & 1)
            bit_flip(bits, i);
    }
    return bits;
}

vector<BitRow> build_row_basis_gf2(const SparseMatrix<int, RowMajor>& H) {
    int n = H.cols();
    int words = (n + 63) / 64;
    vector<BitRow> basis(n, BitRow(words, 0));

    for (int r = 0; r < H.rows(); ++r) {
        BitRow row(words, 0);
        for (SparseMatrix<int, RowMajor>::InnerIterator it(H, r); it; ++it) {
            if (it.value() & 1)
                bit_flip(row, it.col());
        }

        while (!bit_is_zero(row)) {
            int pivot = first_set_bit(row, n);
            if (bit_is_zero(basis[pivot])) {
                basis[pivot] = row;
                break;
            }
            bit_xor_inplace(row, basis[pivot]);
        }
    }
    return basis;
}

bool in_row_space_gf2(const vector<BitRow>& basis, const vector<int>& vec) {
    BitRow row = vector_to_bits(vec);
    int n = vec.size();

    while (!bit_is_zero(row)) {
        int pivot = first_set_bit(row, n);
        if (pivot < 0 || bit_is_zero(basis[pivot]))
            return false;
        bit_xor_inplace(row, basis[pivot]);
    }
    return true;
}

vector<double> parse_fm_list(int argc, char* argv[], int start_idx) {
    if (argc <= start_idx)
        return {0.001, 0.005, 0.01, 0.02, 0.04};

    vector<double> values;
    for (int i = start_idx; i < argc; ++i) {
        string token = argv[i];
        stringstream ss(token);
        string item;
        while (getline(ss, item, ',')) {
            if (!item.empty())
                values.push_back(atof(item.c_str()));
        }
    }

    if (values.empty())
        values = {0.001, 0.005, 0.01, 0.02, 0.04};
    return values;
}

void sampleDepolarizing(vector<int>& ex_bin, vector<int>& ez_bin, double p_depolar,
                        int N_logical, int e) {
    for (int j = 0; j < N_logical; ++j) {
        for (int b = 0; b < e; ++b) {
            double r = static_cast<double>(rand()) / RAND_MAX;
            int idx = j * e + b;
            if (r < p_depolar / 3.0) {
                ex_bin[idx] = 1;
            } else if (r < 2.0 * p_depolar / 3.0) {
                ez_bin[idx] = 1;
            } else if (r < p_depolar) {
                ex_bin[idx] = 1;
                ez_bin[idx] = 1;
            }
        }
    }
}

vector<int> solve_representative_from_projection(const SparseMatrix<int, RowMajor>& G,
                                                 const vector<int>& target,
                                                 GF2p& gf) {
    int rows = G.rows();
    int cols = G.cols();
    vector<vector<int>> A(rows, vector<int>(cols + 1, 0));

    for (int i = 0; i < rows; ++i) {
        for (SparseMatrix<int, RowMajor>::InnerIterator it(G, i); it; ++it)
            A[i][it.col()] = it.value();
        A[i][cols] = target[i];
    }

    int rank = 0;
    vector<int> pivot_cols;
    for (int col = 0; col < cols && rank < rows; ++col) {
        int pivot = -1;
        for (int row = rank; row < rows; ++row) {
            if (A[row][col] != 0) {
                pivot = row;
                break;
            }
        }
        if (pivot == -1)
            continue;

        if (pivot != rank)
            swap(A[pivot], A[rank]);

        int inv = gf.inv(A[rank][col]);
        for (int c = col; c <= cols; ++c)
            A[rank][c] = gf.mul(A[rank][c], inv);

        for (int row = 0; row < rows; ++row) {
            if (row == rank || A[row][col] == 0)
                continue;
            int factor = A[row][col];
            for (int c = col; c <= cols; ++c)
                A[row][c] = gf.sub(A[row][c], gf.mul(factor, A[rank][c]));
        }

        pivot_cols.push_back(col);
        rank++;
    }

    vector<int> solution(cols, 0);
    for (int r = 0; r < rank; ++r)
        solution[pivot_cols[r]] = A[r][cols];
    return solution;
}

vector<int> decode_hx_split(NB_BP_Decoder& decoder,
                            const SparseMatrix<int, RowMajor>& H3,
                            const SparseMatrix<int, RowMajor>& G2,
                            const vector<int>& syn_z_bin,
                            double p_depolar,
                            GF2p& gf,
                            int e) {
    int k2 = G2.rows();
    int n2 = G2.cols();
    int n1 = H3.cols();

    vector<vector<int>> projected(k2, vector<int>(n1, 0));
    for (int g_row = 0; g_row < k2; ++g_row) {
        vector<int> syn_slice(H3.rows() * e, 0);
        for (int h_row = 0; h_row < H3.rows(); ++h_row) {
            int src_row = h_row * k2 + g_row;
            for (int bit = 0; bit < e; ++bit)
                syn_slice[h_row * e + bit] = syn_z_bin[src_row * e + bit];
        }

        projected[g_row] = decoder.syndrome_decode(H3, syn_slice, p_depolar, false);
    }

    vector<int> est_z(n1 * n2, 0);
    for (int block = 0; block < n1; ++block) {
        vector<int> target(k2, 0);
        for (int g_row = 0; g_row < k2; ++g_row)
            target[g_row] = projected[g_row][block];

        vector<int> representative = solve_representative_from_projection(G2, target, gf);
        for (int col = 0; col < n2; ++col)
            est_z[block * n2 + col] = representative[col];
    }

    return est_z;
}

int bit_weight_int(int value) {
    int weight = 0;
    while (value) {
        weight += value & 1;
        value >>= 1;
    }
    return weight;
}

vector<double> make_symbol_prior(int q, int e, double bit_error) {
    vector<double> prior(q, 0.0);
    double sum = 0.0;
    for (int value = 0; value < q; ++value) {
        int weight = bit_weight_int(value);
        prior[value] = pow(bit_error, weight) * pow(1.0 - bit_error, e - weight);
        sum += prior[value];
    }
    for (double& value : prior)
        value /= (sum + 1e-30);
    return prior;
}

vector<vector<double>> channel_pass_from_extrinsic(
    const vector<vector<double>>& other_extrinsic,
    double p_depolar,
    int e
) {
    int N = static_cast<int>(other_extrinsic.size());
    int q = static_cast<int>(other_extrinsic[0].size());
    vector<vector<double>> priors(N, vector<double>(q, 0.0));

    for (int n = 0; n < N; ++n) {
        for (int current = 0; current < q; ++current) {
            double value = 0.0;
            for (int other = 0; other < q; ++other) {
                double joint = 1.0;
                for (int bit = 0; bit < e; ++bit) {
                    bool has_current = ((current >> bit) & 1) != 0;
                    bool has_other = ((other >> bit) & 1) != 0;
                    joint *= (!has_current && !has_other) ? (1.0 - p_depolar)
                                                          : (p_depolar / 3.0);
                }
                value += joint * other_extrinsic[n][other];
            }
            priors[n][current] = value;
        }
        double sum = accumulate(priors[n].begin(), priors[n].end(), 0.0);
        for (double& value : priors[n])
            value /= (sum + 1e-30);
    }

    return priors;
}

vector<int> hard_decision_from_marginals(const vector<vector<double>>& marginals) {
    vector<int> hard(marginals.size(), 0);
    for (int n = 0; n < static_cast<int>(marginals.size()); ++n) {
        hard[n] = static_cast<int>(
            max_element(marginals[n].begin(), marginals[n].end()) - marginals[n].begin());
    }
    return hard;
}

vector<vector<double>> lift_g2_single_parity_marginals(
    const SparseMatrix<int, RowMajor>& G2,
    const vector<vector<double>>& target_marginals,
    const vector<vector<double>>& z_priors_block,
    GF2p& gf
) {
    int k2 = G2.rows();
    int n2 = G2.cols();
    int q = gf.getSize();
    vector<vector<double>> result(n2, vector<double>(q, 0.0));

    bool simple_form = (n2 == k2 + 1);
    for (int row = 0; row < k2 && simple_form; ++row) {
        vector<pair<int, int>> entries;
        for (SparseMatrix<int, RowMajor>::InnerIterator it(G2, row); it; ++it)
            entries.push_back({it.col(), it.value()});
        simple_form = entries.size() == 2 &&
                      entries[0].first == 0 &&
                      entries[1].first == row + 1 &&
                      entries[0].second == 1 &&
                      entries[1].second == 1;
    }

    if (!simple_form) {
        vector<int> hard_target(k2, 0);
        for (int row = 0; row < k2; ++row) {
            hard_target[row] = static_cast<int>(
                max_element(target_marginals[row].begin(), target_marginals[row].end()) -
                target_marginals[row].begin());
        }
        vector<int> representative = solve_representative_from_projection(G2, hard_target, gf);
        for (int col = 0; col < n2; ++col)
            result[col][representative[col]] = 1.0;
        return result;
    }

    vector<vector<double>> s(k2, vector<double>(q, 0.0));
    for (int row = 0; row < k2; ++row) {
        for (int z0 = 0; z0 < q; ++z0) {
            double sum = 0.0;
            for (int zi = 0; zi < q; ++zi) {
                int target = gf.add(z0, zi);
                sum += z_priors_block[row + 1][zi] * target_marginals[row][target];
            }
            s[row][z0] = sum;
        }
    }

    vector<double> z0_weight(q, 0.0);
    for (int z0 = 0; z0 < q; ++z0) {
        double prod = z_priors_block[0][z0];
        for (int row = 0; row < k2; ++row)
            prod *= s[row][z0];
        z0_weight[z0] = prod;
        result[0][z0] = prod;
    }

    double z0_sum = accumulate(result[0].begin(), result[0].end(), 0.0);
    for (double& value : result[0])
        value /= (z0_sum + 1e-30);

    for (int col = 1; col < n2; ++col) {
        int row_for_col = col - 1;
        for (int zi = 0; zi < q; ++zi) {
            double marginal = 0.0;
            for (int z0 = 0; z0 < q; ++z0) {
                double other_rows = z_priors_block[0][z0];
                for (int row = 0; row < k2; ++row) {
                    if (row != row_for_col)
                        other_rows *= s[row][z0];
                }
                int target = gf.add(z0, zi);
                marginal += other_rows * target_marginals[row_for_col][target];
            }
            result[col][zi] = z_priors_block[col][zi] * marginal;
        }
        double col_sum = accumulate(result[col].begin(), result[col].end(), 0.0);
        for (double& value : result[col])
            value /= (col_sum + 1e-30);
    }

    return result;
}

vector<vector<double>> decode_hx_split_soft(
    NB_BP_Decoder& decoder,
    const SparseMatrix<int, RowMajor>& H3,
    const SparseMatrix<int, RowMajor>& G2,
    const vector<int>& syn_z_bin,
    double p_depolar,
    const vector<vector<double>>& z_priors,
    GF2p& gf,
    int e
) {
    int k2 = G2.rows();
    int n2 = G2.cols();
    int n1 = H3.cols();
    int q = gf.getSize();
    vector<vector<vector<double>>> projected(k2);
    for (int g_row = 0; g_row < k2; ++g_row) {
        vector<vector<double>> target_priors(n1, vector<double>(q, 0.0));
        for (int block = 0; block < n1; ++block) {
            int base = block * n2;
            for (int z0 = 0; z0 < q; ++z0) {
                for (int zi = 0; zi < q; ++zi) {
                    int target = gf.add(z0, zi);
                    target_priors[block][target] +=
                        z_priors[base][z0] * z_priors[base + g_row + 1][zi];
                }
            }
            double target_sum = accumulate(target_priors[block].begin(),
                                           target_priors[block].end(), 0.0);
            for (double& value : target_priors[block])
                value /= (target_sum + 1e-30);
        }

        vector<int> syn_slice(H3.rows() * e, 0);
        for (int h_row = 0; h_row < H3.rows(); ++h_row) {
            int src_row = h_row * k2 + g_row;
            for (int bit = 0; bit < e; ++bit)
                syn_slice[h_row * e + bit] = syn_z_bin[src_row * e + bit];
        }
        projected[g_row] =
            decoder.syndrome_decode_soft_with_prior(H3, syn_slice, target_priors, false);
    }

    vector<vector<double>> z_marginals(n1 * n2, vector<double>(q, 0.0));
    for (int block = 0; block < n1; ++block) {
        vector<vector<double>> target_marginals(k2, vector<double>(q, 0.0));
        for (int g_row = 0; g_row < k2; ++g_row)
            target_marginals[g_row] = projected[g_row][block];

        vector<vector<double>> z_priors_block(n2, vector<double>(q, 0.0));
        for (int col = 0; col < n2; ++col)
            z_priors_block[col] = z_priors[block * n2 + col];

        vector<vector<double>> lifted =
            lift_g2_single_parity_marginals(G2, target_marginals, z_priors_block, gf);
        for (int col = 0; col < n2; ++col)
            z_marginals[block * n2 + col] = lifted[col];
    }

    return z_marginals;
}

pair<vector<vector<double>>, vector<vector<double>>> decode_hx_split_soft_with_extrinsic(
    NB_BP_Decoder& decoder,
    const SparseMatrix<int, RowMajor>& H3,
    const SparseMatrix<int, RowMajor>& G2,
    const vector<int>& syn_z_bin,
    double p_depolar,
    const vector<vector<double>>& z_priors,
    GF2p& gf,
    int e
) {
    vector<vector<double>> marginals =
        decode_hx_split_soft(decoder, H3, G2, syn_z_bin, p_depolar, z_priors, gf, e);
    vector<vector<double>> extrinsics = marginals;
    for (int n = 0; n < static_cast<int>(extrinsics.size()); ++n) {
        for (int value = 0; value < static_cast<int>(extrinsics[n].size()); ++value)
            extrinsics[n][value] /= (z_priors[n][value] + 1e-30);
        double sum = accumulate(extrinsics[n].begin(), extrinsics[n].end(), 0.0);
        for (double& value : extrinsics[n])
            value /= (sum + 1e-30);
    }
    return {marginals, extrinsics};
}

pair<vector<int>, vector<int>> combine_joint_marginals(
    const vector<vector<double>>& x_marginals,
    const vector<vector<double>>& z_marginals,
    const vector<vector<double>>& x_priors,
    const vector<vector<double>>& z_priors,
    double p_depolar,
    int e,
    vector<vector<double>>* next_x_priors = nullptr,
    vector<vector<double>>* next_z_priors = nullptr,
    double damping = 0.5
) {
    int N = static_cast<int>(x_marginals.size());
    int q = static_cast<int>(x_marginals[0].size());
    vector<int> est_x(N, 0), est_z(N, 0);

    if (next_x_priors)
        next_x_priors->assign(N, vector<double>(q, 0.0));
    if (next_z_priors)
        next_z_priors->assign(N, vector<double>(q, 0.0));

    for (int j = 0; j < N; ++j) {
        double best = -1.0;
        double joint_sum = 0.0;
        for (int x = 0; x < q; ++x) {
            double lx = x_marginals[j][x] / (x_priors[j][x] + 1e-30);
            for (int z = 0; z < q; ++z) {
                double joint_prior = 1.0;
                for (int bit = 0; bit < e; ++bit) {
                    bool has_x = ((x >> bit) & 1) != 0;
                    bool has_z = ((z >> bit) & 1) != 0;
                    joint_prior *= (!has_x && !has_z) ? (1.0 - p_depolar)
                                                       : (p_depolar / 3.0);
                }
                double lz = z_marginals[j][z] / (z_priors[j][z] + 1e-30);
                double score = joint_prior * lx * lz;
                joint_sum += score;
                if (next_x_priors)
                    (*next_x_priors)[j][x] += score;
                if (next_z_priors)
                    (*next_z_priors)[j][z] += score;
                if (score > best) {
                    best = score;
                    est_x[j] = x;
                    est_z[j] = z;
                }
            }
        }

        if (next_x_priors) {
            for (int x = 0; x < q; ++x) {
                double posterior = (*next_x_priors)[j][x] / (joint_sum + 1e-30);
                (*next_x_priors)[j][x] =
                    damping * posterior + (1.0 - damping) * x_priors[j][x];
            }
        }
        if (next_z_priors) {
            for (int z = 0; z < q; ++z) {
                double posterior = (*next_z_priors)[j][z] / (joint_sum + 1e-30);
                (*next_z_priors)[j][z] =
                    damping * posterior + (1.0 - damping) * z_priors[j][z];
            }
        }
    }

    return {est_x, est_z};
}

int main(int argc, char* argv[]) {
    srand(static_cast<unsigned int>(time(nullptr)));

    int trials = (argc > 1) ? atoi(argv[1]) : 1000;
    int max_iters = (argc > 2) ? atoi(argv[2]) : 50;
    int target_failures = (argc > 3) ? atoi(argv[3]) : 0;

    string suffix = "";
    if (argc > 4 && string(argv[4]).find(".") == string::npos) {
        suffix = argv[4];
        cout << "Using matrix suffix: " << suffix << endl;
    }

    int fm_start_idx = suffix.empty() ? 4 : 5;
    if (argc > fm_start_idx && string(argv[fm_start_idx]) == "--hybrid-joint") {
        fm_start_idx++;
        cout << "Warning: --hybrid-joint is disabled; using split decoder." << endl;
    } else if (argc > fm_start_idx && string(argv[fm_start_idx]) == "--gd-css-joint") {
        fm_start_idx++;
        cout << "Warning: --gd-css-joint is disabled; using split decoder." << endl;
    }

    vector<double> fm_list;
    if (argc > fm_start_idx) {
        fm_list = parse_fm_list(argc, argv, fm_start_idx);
    } else {
        fm_list = {0.02, 0.015, 0.01, 0.005, 0.001};
    }

    GF2p gf(8, {1, 0, 0, 0, 1, 1, 1, 0, 1});
    int e = gf.getPower();
    NB_BP_Decoder decoder(gf, max_iters);

    SparseMatrix<int, RowMajor> HX_logical = loadMatrix(matrixPathFor("HX", suffix), gf);
    SparseMatrix<int, RowMajor> HZ_logical = loadMatrix(matrixPathFor("HZ", suffix), gf);
    SparseMatrix<int, RowMajor> H3_logical = loadMatrix(matrixPathFor("H3", suffix), gf);
    SparseMatrix<int, RowMajor> G2_logical = loadMatrix(matrixPathFor("G2", suffix), gf);

    SparseMatrix<int, RowMajor> HX_phys = loadBinaryMatrix(matrixPathFor("HX_Phys", suffix));
    SparseMatrix<int, RowMajor> HZ_phys = loadBinaryMatrix(matrixPathFor("HZ_Phys", suffix));
    auto x_stabilizer_basis = build_row_basis_gf2(HX_phys);
    auto z_stabilizer_basis = build_row_basis_gf2(HZ_phys);

    int N_logical = HX_logical.cols();
    int n_phys = N_logical * e;

    cout << "f_m\tFER\tRawFER\tAvgTrueWeight\tTotal\tFailures\tRawFailures\tUpper95IfZero" << endl;

    for (double fm : fm_list) {
        double p_depolar = 1.5 * fm;
        int failures = 0;
        int raw_failures = 0;
        int completed_trials = 0;
        long long total_true_weight = 0;

        for (int t = 0; t < trials; ++t) {
            vector<int> err_x_bin(n_phys, 0);
            vector<int> err_z_bin(n_phys, 0);
            sampleDepolarizing(err_x_bin, err_z_bin, p_depolar, N_logical, e);
            completed_trials++;

            int true_weight = 0;
            for (int i = 0; i < n_phys; ++i) {
                if (err_x_bin[i] != 0 || err_z_bin[i] != 0)
                    true_weight++;
            }
            total_true_weight += true_weight;

            bool raw_error =
                any_of(err_x_bin.begin(), err_x_bin.end(), [](int bit) { return bit != 0; }) ||
                any_of(err_z_bin.begin(), err_z_bin.end(), [](int bit) { return bit != 0; });
            if (raw_error)
                raw_failures++;

            vector<int> sx_bin = get_binary_syndrome(HZ_phys, err_x_bin);
            vector<int> sz_bin = get_binary_syndrome(HX_phys, err_z_bin);

            vector<int> est_x_gf =
                decoder.syndrome_decode(HZ_logical, sx_bin, p_depolar, true);
            vector<int> est_z_gf =
                decode_hx_split(decoder, H3_logical, G2_logical, sz_bin, p_depolar, gf, e);

            vector<int> dec_x_bin(n_phys, 0);
            vector<int> dec_z_bin(n_phys, 0);
            for (int j = 0; j < N_logical; ++j) {
                for (int b = 0; b < e; ++b) {
                    dec_x_bin[j * e + b] = (est_x_gf[j] >> b) & 1;
                    dec_z_bin[j * e + b] = (est_z_gf[j] >> b) & 1;
                }
            }

            vector<int> residual_x(n_phys, 0);
            vector<int> residual_z(n_phys, 0);
            for (int i = 0; i < n_phys; ++i) {
                residual_x[i] = err_x_bin[i] ^ dec_x_bin[i];
                residual_z[i] = err_z_bin[i] ^ dec_z_bin[i];
            }

            bool success =
                in_row_space_gf2(x_stabilizer_basis, residual_x) &&
                in_row_space_gf2(z_stabilizer_basis, residual_z);

            if (!success)
                failures++;

            if (target_failures > 0 && failures >= target_failures)
                break;
        }

        double denom = static_cast<double>(max(completed_trials, 1));
        double upper95_if_zero = (failures == 0 && completed_trials > 0)
                                     ? (3.0 / completed_trials)
                                     : 0.0;
        cout << fm << "\t" << static_cast<double>(failures) / denom << "\t"
             << static_cast<double>(raw_failures) / denom << "\t"
             << static_cast<double>(total_true_weight) / denom << "\t"
             << completed_trials << "\t" << failures << "\t" << raw_failures << "\t"
             << upper95_if_zero << endl;
    }

    return 0;
}

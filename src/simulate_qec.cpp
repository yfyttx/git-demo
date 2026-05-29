#include <Eigen/Sparse>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "NB_BP_Decoder.h"
#include "gf2p/gf2p.h"

using namespace std;
using namespace Eigen;

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

int main(int argc, char* argv[]) {
    srand(static_cast<unsigned int>(time(nullptr)));

    int trials = (argc > 1) ? atoi(argv[1]) : 100;
    int max_iters = (argc > 2) ? atoi(argv[2]) : 50;
    int target_failures = (argc > 3) ? atoi(argv[3]) : 0;
    vector<double> fm_list = parse_fm_list(argc, argv, 4);

    GF2p gf(8, {1, 0, 0, 0, 1, 1, 1, 0, 1});
    int e = gf.getPower();
    NB_BP_Decoder decoder(gf, max_iters);

    SparseMatrix<int, RowMajor> HX_logical = loadMatrix("HX_Safe_GF256.txt", gf);
    SparseMatrix<int, RowMajor> HZ_logical = loadMatrix("HZ_Safe_GF256.txt", gf);
    SparseMatrix<int, RowMajor> H3_logical = loadMatrix("H3_sparse.txt", gf);
    SparseMatrix<int, RowMajor> G2_logical = loadMatrix("G2_sparse.txt", gf);

    SparseMatrix<int, RowMajor> HX_phys = loadBinaryMatrix("HX_Physical_Binary.txt");
    SparseMatrix<int, RowMajor> HZ_phys = loadBinaryMatrix("HZ_Physical_Binary.txt");
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

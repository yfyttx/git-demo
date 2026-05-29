/**
 * main.cpp - Quantum LDPC matrix construction tool.
 *
 * This file resolves the previous merge conflict by keeping the newer GF2p
 * based path and removing the older local SimpleGF implementation.
 */

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../peg-unige/BigGirth.h"
#include "gf2p/gf2p.h"
#include "ldpc_codes/binary_codes.h"
#include "ldpc_codes/nonbinary_codes.h"

using namespace Eigen;
using namespace std;

int calculate_girth(const SparseMatrix<int> &H) {
    int rows = H.rows();
    int cols = H.cols();
    if (H.nonZeros() == 0)
        return 0;

    int min_girth = 999;
    int total_nodes = rows + cols;
    vector<vector<int>> adj(total_nodes);

    for (int k = 0; k < H.outerSize(); ++k) {
        for (SparseMatrix<int>::InnerIterator it(H, k); it; ++it) {
            int r = it.row();
            int c = it.col();
            adj[c].push_back(r + cols);
            adj[r + cols].push_back(c);
        }
    }

    for (int start_node = 0; start_node < total_nodes; ++start_node) {
        if (adj[start_node].empty())
            continue;

        vector<int> dist(total_nodes, -1);
        vector<int> parent(total_nodes, -1);
        queue<int> q;

        dist[start_node] = 0;
        q.push(start_node);

        while (!q.empty()) {
            int u = q.front();
            q.pop();

            for (int v : adj[u]) {
                if (dist[v] == -1) {
                    dist[v] = dist[u] + 1;
                    parent[v] = u;
                    q.push(v);
                } else if (v != parent[u]) {
                    int cycle_len = dist[u] + dist[v] + 1;
                    if (cycle_len < min_girth)
                        min_girth = cycle_len;
                    if (min_girth == 4)
                        return 4;
                }
            }
        }
    }

    return (min_girth == 999) ? 0 : min_girth;
}

pair<SparseMatrix<int>, SparseMatrix<int>>
generateOrthogonalNB_TIT(const MatrixXi &Hc, const MatrixXi &Hd, GF2p &gf) {
    int m2 = Hc.rows();
    int n2 = Hc.cols();
    int m3 = Hd.rows();
    int n3 = Hd.cols();

    vector<pair<int, int>> Hc_nonzero_pos;
    for (int i = 0; i < m2; ++i) {
        for (int j = 0; j < n2; ++j) {
            if (Hc(i, j) == 1)
                Hc_nonzero_pos.emplace_back(i, j);
        }
    }

    MatrixXi system = buildLinearSystem(Hc_nonzero_pos, Hc, Hd, gf.getSize() - 1);
    vector<int> log_gamma = gaussElimination(system, gf.getSize() - 1);

    SparseMatrix<int> H2(m2, n2);
    MatrixXi H2_dense = MatrixXi::Zero(m2, n2);
    vector<Triplet<int>> t2;
    for (size_t i = 0; i < Hc_nonzero_pos.size(); ++i) {
        int val = gf.exp(log_gamma[i]);
        int r = Hc_nonzero_pos[i].first;
        int c = Hc_nonzero_pos[i].second;
        t2.emplace_back(r, c, val);
        H2_dense(r, c) = val;
    }
    H2.setFromTriplets(t2.begin(), t2.end());

    SparseMatrix<int> H3(m3, n3);
    vector<Triplet<int>> t3;
    for (int i = 0; i < m3; ++i) {
        vector<int> delta_row = solveDeltaRow(i, Hc, Hd, H2_dense, gf);
        for (int j = 0; j < n3; ++j) {
            if (delta_row[j] != 0)
                t3.emplace_back(i, j, delta_row[j]);
        }
    }
    H3.setFromTriplets(t3.begin(), t3.end());

    if (t3.empty()) {
        cerr << "CRITICAL ERROR: H3 has no non-zero elements. Check solveDeltaRow." << endl;
    }

    return {H2, H3};
}

void saveToEssaiFormat(const SparseMatrix<int> &H, const string &filename, GF2p &gf) {
    ofstream out(filename);
    if (!out)
        return;

    int M = H.rows();
    int N = H.cols();
    out << N << " " << M << " " << gf.getSize() << "\n";

    vector<int> dv(N, 0);
    vector<int> dc(M, 0);
    for (int k = 0; k < H.outerSize(); ++k) {
        for (SparseMatrix<int>::InnerIterator it(H, k); it; ++it) {
            dv[it.col()]++;
            dc[it.row()]++;
        }
    }

    for (int i = 0; i < N; ++i)
        out << dv[i] << " ";
    out << "\n";

    for (int i = 0; i < M; ++i)
        out << dc[i] << " ";
    out << "\n";

    vector<vector<pair<int, int>>> rows(M);
    for (int k = 0; k < H.outerSize(); ++k) {
        for (SparseMatrix<int>::InnerIterator it(H, k); it; ++it) {
            rows[it.row()].push_back({it.col() + 1, gf.log(it.value())});
        }
    }

    for (int i = 0; i < M; ++i) {
        for (auto &p : rows[i])
            out << p.first << " " << p.second << "    ";
        out << "\n";
    }

    cout << "  -> Saved " << filename << " (" << N << "x" << M << ")\n";
}

SparseMatrix<int> buildH_PEG(int M, int N, int dv, GF2p &gf, const string &name) {
    cout << "  [PEG] Generating " << name << " (M=" << M << ", N=" << N
         << ", dv=" << dv << ")..." << endl;
    vector<int> degSeq(N, dv);
    BigGirth bg(M, N, degSeq.data(), "dummy.dat", 1, 10000, false);
    bg.loadH();

    SparseMatrix<int> H(M, N);
    vector<Triplet<int>> t;
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            if (bg.H[i][j])
                t.emplace_back(i, j, gf.exp(rand() % (gf.getSize() - 1)));
        }
    }
    H.setFromTriplets(t.begin(), t.end());
    return H;
}

SparseMatrix<int> computeGeneratorMatrix(const SparseMatrix<int> &H_sparse, GF2p &gf) {
    int m = H_sparse.rows();
    int n = H_sparse.cols();
    MatrixXi H = MatrixXi::Zero(m, n);

    for (int k = 0; k < H_sparse.outerSize(); ++k) {
        for (SparseMatrix<int>::InnerIterator it(H_sparse, k); it; ++it) {
            H(it.row(), it.col()) = it.value();
        }
    }

    int pivot_row = 0;
    vector<int> pivot_cols;
    vector<bool> is_pivot(n, false);
    vector<bool> col_used(n, false);

    while (pivot_row < m) {
        int best_col = -1;
        int min_degree = n + 1;

        for (int j = 0; j < n; ++j) {
            if (col_used[j])
                continue;

            int current_degree = 0;
            for (int r = pivot_row; r < m; ++r) {
                if (H(r, j) != 0)
                    current_degree++;
            }

            if (current_degree > 0 && current_degree < min_degree) {
                min_degree = current_degree;
                best_col = j;
            }
        }

        if (best_col == -1)
            break;

        int sel_row = -1;
        for (int r = pivot_row; r < m; ++r) {
            if (H(r, best_col) != 0) {
                sel_row = r;
                break;
            }
        }

        if (sel_row != pivot_row)
            H.row(pivot_row).swap(H.row(sel_row));

        int iv = gf.inv(H(pivot_row, best_col));
        for (int k = 0; k < n; ++k) {
            if (H(pivot_row, k) != 0)
                H(pivot_row, k) = gf.mul(H(pivot_row, k), iv);
        }

        for (int r = 0; r < m; ++r) {
            if (r != pivot_row && H(r, best_col) != 0) {
                int factor = H(r, best_col);
                for (int k = 0; k < n; ++k) {
                    if (H(pivot_row, k) != 0 || H(r, k) != 0)
                        H(r, k) = gf.sub(H(r, k), gf.mul(factor, H(pivot_row, k)));
                }
            }
        }

        pivot_cols.push_back(best_col);
        is_pivot[best_col] = true;
        col_used[best_col] = true;
        pivot_row++;
    }

    int rank = static_cast<int>(pivot_cols.size());
    SparseMatrix<int> G(n - rank, n);
    vector<Triplet<int>> tG;
    int g_row_idx = 0;

    for (int j = 0; j < n; ++j) {
        if (is_pivot[j])
            continue;

        tG.emplace_back(g_row_idx, j, 1);
        for (int i = 0; i < rank; ++i) {
            int val = H(i, j);
            if (val != 0)
                tG.emplace_back(g_row_idx, pivot_cols[i], val);
        }
        g_row_idx++;
    }

    G.setFromTriplets(tG.begin(), tG.end());
    return G;
}

SparseMatrix<int> kroneckerNB(const SparseMatrix<int> &A, const SparseMatrix<int> &B, GF2p &gf) {
    SparseMatrix<int> C(A.rows() * B.rows(), A.cols() * B.cols());
    vector<Triplet<int>> t;
    t.reserve(static_cast<size_t>(A.nonZeros()) * static_cast<size_t>(B.nonZeros()));

    for (int k = 0; k < A.outerSize(); ++k) {
        for (SparseMatrix<int>::InnerIterator itA(A, k); itA; ++itA) {
            for (int l = 0; l < B.outerSize(); ++l) {
                for (SparseMatrix<int>::InnerIterator itB(B, l); itB; ++itB) {
                    int val = gf.mul(itA.value(), itB.value());
                    if (val) {
                        t.emplace_back(itA.row() * B.rows() + itB.row(),
                                       itA.col() * B.cols() + itB.col(), val);
                    }
                }
            }
        }
    }

    C.setFromTriplets(t.begin(), t.end());
    return C;
}

SparseMatrix<int> makeIdentityNB(int size) {
    SparseMatrix<int> I(size, size);
    vector<Triplet<int>> t;
    t.reserve(size);

    for (int i = 0; i < size; ++i)
        t.emplace_back(i, i, 1);

    I.setFromTriplets(t.begin(), t.end());
    return I;
}

SparseMatrix<int> makeAllOnesRowNB(int cols) {
    SparseMatrix<int> H(1, cols);
    vector<Triplet<int>> t;
    t.reserve(cols);

    for (int j = 0; j < cols; ++j)
        t.emplace_back(0, j, 1);

    H.setFromTriplets(t.begin(), t.end());
    return H;
}

SparseMatrix<int> verticalStackNB(const SparseMatrix<int> &top, const SparseMatrix<int> &bottom) {
    if (top.cols() != bottom.cols())
        throw runtime_error("verticalStackNB: matrices must have the same number of columns");

    SparseMatrix<int> stacked(top.rows() + bottom.rows(), top.cols());
    vector<Triplet<int>> t;
    t.reserve(static_cast<size_t>(top.nonZeros()) + static_cast<size_t>(bottom.nonZeros()));

    for (int k = 0; k < top.outerSize(); ++k) {
        for (SparseMatrix<int>::InnerIterator it(top, k); it; ++it)
            t.emplace_back(it.row(), it.col(), it.value());
    }

    for (int k = 0; k < bottom.outerSize(); ++k) {
        for (SparseMatrix<int>::InnerIterator it(bottom, k); it; ++it)
            t.emplace_back(top.rows() + it.row(), it.col(), it.value());
    }

    stacked.setFromTriplets(t.begin(), t.end());
    return stacked;
}

int verifyOrthogonality(const SparseMatrix<int> &HX, const SparseMatrix<int> &HZ, const GF2p &gf) {
    if (HX.cols() != HZ.cols())
        return -1;

    SparseMatrix<int, RowMajor> HX_r = HX;
    SparseMatrix<int, RowMajor> HZ_r = HZ;
    int error_count = 0;

    for (int i = 0; i < HX.rows(); ++i) {
        for (int j = 0; j < HZ.rows(); ++j) {
            int dot_product = 0;
            SparseMatrix<int, RowMajor>::InnerIterator itX(HX_r, i);
            SparseMatrix<int, RowMajor>::InnerIterator itZ(HZ_r, j);

            while (itX && itZ) {
                if (itX.col() == itZ.col()) {
                    int term = gf.mul(itX.value(), itZ.value());
                    dot_product = gf.add(dot_product, term);
                    ++itX;
                    ++itZ;
                } else if (itX.col() < itZ.col()) {
                    ++itX;
                } else {
                    ++itZ;
                }
            }

            if (dot_product != 0)
                error_count++;
        }
    }

    return error_count;
}

MatrixXi getCompanionMatrix(int val, GF2p &gf, const vector<int> &irreducible_poly_vec) {
    int e = gf.getPower();
    MatrixXi A = MatrixXi::Zero(e, e);

    if (val == 0)
        return A;

    MatrixXi A_alpha = MatrixXi::Zero(e, e);
    for (int j = 0; j < e - 1; ++j)
        A_alpha(j + 1, j) = 1;

    for (int i = 0; i < e; ++i) {
        if (irreducible_poly_vec[e - i] != 0)
            A_alpha(i, e - 1) = 1;
    }

    int l = gf.log(val);
    A = MatrixXi::Identity(e, e);
    for (int i = 0; i < l; ++i)
        A = (A * A_alpha).unaryExpr([](int x) { return x % 2; });

    return A;
}

SparseMatrix<int> convertToBinary(const SparseMatrix<int> &H_nb, GF2p &gf,
                                  const vector<int> &irreducible_poly_vec,
                                  bool is_transposed) {
    int e = gf.getPower();
    int M_nb = H_nb.rows();
    int N_nb = H_nb.cols();

    SparseMatrix<int> H_bin(M_nb * e, N_nb * e);
    vector<Triplet<int>> triplets;

    for (int k = 0; k < H_nb.outerSize(); ++k) {
        for (SparseMatrix<int>::InnerIterator it(H_nb, k); it; ++it) {
            MatrixXi comp = getCompanionMatrix(it.value(), gf, irreducible_poly_vec);
            if (is_transposed)
                comp = comp.transpose().eval();

            for (int i = 0; i < e; ++i) {
                for (int j = 0; j < e; ++j) {
                    if (comp(i, j) == 1)
                        triplets.emplace_back(it.row() * e + i, it.col() * e + j, 1);
                }
            }
        }
    }

    H_bin.setFromTriplets(triplets.begin(), triplets.end());
    return H_bin;
}

void saveToBinaryFormat(const SparseMatrix<int> &H, const string &filename) {
    ofstream out(filename);
    if (!out)
        return;

    int M = H.rows();
    int N = H.cols();
    out << N << " " << M << " 2\n";
    out << "\n\n";

    SparseMatrix<int, RowMajor> Hr = H;
    for (int i = 0; i < M; ++i) {
        for (SparseMatrix<int, RowMajor>::InnerIterator it(Hr, i); it; ++it)
            out << it.col() + 1 << " 1    ";
        out << "\n";
    }

    cout << "  -> Saved Binary Physical Matrix: " << filename << endl;
}

int main() {
    srand(static_cast<unsigned int>(time(nullptr)));
    vector<int> field_poly = {1, 0, 0, 0, 1, 1, 1, 0, 1};
    GF2p gf(8, field_poly);

    cout << "=== QLDPC: Hybrid (PEG + TIT) Construction ===\n";

    cout << "  [TIT] Generating H1 and H3 with linear constraint solving..." << endl;
    const int J = 2;
    int L = 6;
    int P = 31;
    int sigma = 5;
    int tau = 2;

    pair<MatrixXi, MatrixXi> tit_bin = constructBinaryHMatrices(J, L, P, sigma, tau);
    auto [H1, H3] = generateOrthogonalNB_TIT(tit_bin.first, tit_bin.second, gf);

    int h2_cols = 5;
    h2_cols = max(h2_cols, 2);
    SparseMatrix<int> H2 = makeAllOnesRowNB(h2_cols);
    SparseMatrix<int> G2 = computeGeneratorMatrix(H2, gf);

    cout << "  [Check] H1 non-zeros: " << H1.nonZeros() << endl;
    cout << "  [Check] H2 all-ones row length: " << H2.cols() << endl;
    cout << "  [Check] H3 non-zeros: " << H3.nonZeros() << endl;
    cout << "  [Check] G2 rows: " << G2.rows() << ", cols: " << G2.cols() << endl;
    cout << "  [Check] NB-TIT Orthogonality (H1*H3^T) error count: "
         << verifyOrthogonality(H1, H3, gf) << endl;

    cout << "  [Build] HZ = [H1 x I; I x H2], HX = H3 x G2..." << endl;
    SparseMatrix<int> I1 = makeIdentityNB(H1.cols());
    SparseMatrix<int> I2 = makeIdentityNB(H2.cols());
    SparseMatrix<int> HZ_top = kroneckerNB(H1, I2, gf);
    SparseMatrix<int> HZ_bottom = kroneckerNB(I1, H2, gf);
    SparseMatrix<int> HZ = verticalStackNB(HZ_top, HZ_bottom);
    SparseMatrix<int> HX = kroneckerNB(H3, G2, gf);

    cout << "  [Result] HX: " << HX.rows() << "x" << HX.cols() << endl;
    cout << "  [Result] HZ: " << HZ.rows() << "x" << HZ.cols() << endl;

    cout << "  [H1 Girth]: " << calculate_girth(H1) << endl;
    cout << "  [H2 Girth]: " << calculate_girth(H2) << endl;
    cout << "  [H3 Girth]: " << calculate_girth(H3) << endl;
    cout << "  [Hx Girth]: " << calculate_girth(HX) << endl;
    cout << "  [Hz Girth]: " << calculate_girth(HZ) << endl;

    cout << "\n=== Verifying Quantum Orthogonality ===" << endl;
    cout << "  Checking if HX * HZ^T = 0 in GF(" << gf.getSize() << ")..." << endl;

    int non_orthogonal_elements = verifyOrthogonality(HX, HZ, gf);
    if (non_orthogonal_elements == 0) {
        cout << "  [SUCCESS] -> HX and HZ are perfectly orthogonal in GF(" << gf.getSize()
             << ")!" << endl;
    } else {
        cout << "  [FAILED]  -> Found " << non_orthogonal_elements
             << " non-zero elements in HX * HZ^T." << endl;
    }
    cout << "=======================================\n" << endl;

    cout << "  [Convert] Expanding to binary physical matrices..." << endl;
    SparseMatrix<int> HX_physical = convertToBinary(HX, gf, field_poly, false);
    SparseMatrix<int> HZ_physical = convertToBinary(HZ, gf, field_poly, true);

    cout << "  [Result] Physical HX size: " << HX_physical.rows() << "x"
         << HX_physical.cols() << endl;
    cout << "  [Result] Physical HZ size: " << HZ_physical.rows() << "x"
         << HZ_physical.cols() << endl;

    saveToBinaryFormat(HX_physical, "HX_Physical_Binary.txt");
    saveToBinaryFormat(HZ_physical, "HZ_Physical_Binary.txt");
    saveToEssaiFormat(H1, "H1_sparse.txt", gf);
    saveToEssaiFormat(H2, "H2_sparse.txt", gf);
    saveToEssaiFormat(H3, "H3_sparse.txt", gf);
    saveToEssaiFormat(G2, "G2_sparse.txt", gf);
    saveToEssaiFormat(HX, "HX_Safe_GF256.txt", gf);
    saveToEssaiFormat(HZ, "HZ_Safe_GF256.txt", gf);

    cout << "CHECK: VERSION 3.0 - PRODUCT HZ/HX CONSTRUCTION" << endl;
    return 0;
}

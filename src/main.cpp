/**
 * main.cpp - 量子 LDPC 矩阵构造与导出工具 (H2 增强版)
 * 目标：通过增大 H2 规模确保无 4 环，同时保持 H1/H3 的 TIT 构造逻辑
 */

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>

#include "../peg-unige/BigGirth.h"
#include "../peg-unige/Random.h"
#include "ldpc_codes/binary_codes.h"

using namespace Eigen;
using namespace std;

// ==========================================
// 1. GF(64) 类 - 高元域运算
// ==========================================
class SimpleGF {
    int m, size;
    std::vector<int> exp_table, log_table;

  public:
    SimpleGF(int m_val) : m(m_val) {
        size = 1 << m;
        exp_table.resize(size * 2);
        log_table.resize(size);
        int prim_poly = (m == 6) ? 0x43 : 0x13;
        int x = 1;
        for (int i = 0; i < size - 1; i++) {
            exp_table[i] = x;
            log_table[x] = i;
            x <<= 1;
            if (x & size)
                x ^= prim_poly;
        }
        for (int i = 0; i < size - 1; i++)
            exp_table[i + (size - 1)] = exp_table[i];
        log_table[0] = -1;
    }
    int getSize() const {
        return size;
    }
    int mul(int a, int b) const {
        return (a == 0 || b == 0) ? 0 : exp_table[log_table[a] + log_table[b]];
    }
    int inv(int a) const {
        if (!a)
            return 0;
        return exp_table[(size - 1) - log_table[a]];
    }
    int sub(int a, int b) const {
        return a ^ b;
    }
    int exp(int k) const {
        return exp_table[k % (size - 1)];
    }
    int log(int val) const {
        return (val <= 0 || val >= size) ? -1 : log_table[val];
    }
};

// ==========================================
// 2. 矩阵工具函数
// ==========================================

void saveToEssaiFormat(const SparseMatrix<int> &H, const std::string &filename, SimpleGF &gf) {
    ofstream out(filename);
    if (!out)
        return;
    int M = H.rows(), N = H.cols();
    out << N << " " << M << " " << gf.getSize() << "\n";
    vector<int> dv(N, 0), dc(M, 0);
    for (int k = 0; k < H.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(H, k); it; ++it) {
            dv[it.col()]++;
            dc[it.row()]++;
        }
    for (int i = 0; i < N; ++i)
        out << dv[i] << " ";
    out << "\n";
    for (int i = 0; i < M; ++i)
        out << dc[i] << " ";
    out << "\n";
    vector<vector<pair<int, int>>> rows(M);
    for (int k = 0; k < H.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(H, k); it; ++it)
            rows[it.row()].push_back({it.col() + 1, gf.log(it.value())});
    for (int i = 0; i < M; ++i) {
        for (auto &p : rows[i])
            out << p.first << " " << p.second << "    ";
        out << "\n";
    }
    out.close();
    cout << "  -> Saved " << filename << " (" << N << "x" << M << ")\n";
}

SparseMatrix<int> denseToSparse(const MatrixXi &dense, SimpleGF &gf) {
    SparseMatrix<int> sparse(dense.rows(), dense.cols());
    vector<Triplet<int>> t;
    for (int i = 0; i < dense.rows(); ++i)
        for (int j = 0; j < dense.cols(); ++j)
            if (dense(i, j))
                t.emplace_back(i, j, gf.exp(rand() % (gf.getSize() - 1)));
    sparse.setFromTriplets(t.begin(), t.end());
    return sparse;
}

SparseMatrix<int> buildH_PEG_Irregular(int M, int N, SimpleGF &gf) {
    vector<int> degSeq(N);
    for (int i = 0; i < N; ++i) {
        double pivot = (double)i / N;
        if (pivot < 0.4)
            degSeq[i] = 2;
        else if (pivot < 0.9)
            degSeq[i] = 3;
        else
            degSeq[i] = 4;
    }
    // 目标围长设为 6，尝试 10000 次
    BigGirth bg(M, N, degSeq.data(), "dummy.dat", 6, 10000, false);
    bg.loadH();
    SparseMatrix<int> H(M, N);
    vector<Triplet<int>> t;
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            if (bg.H[i][j])
                t.emplace_back(i, j, gf.exp(rand() % (gf.getSize() - 1)));
    H.setFromTriplets(t.begin(), t.end());
    return H;
}

SparseMatrix<int> computeGeneratorMatrix(const SparseMatrix<int> &H_sparse, SimpleGF &gf) {
    int m = H_sparse.rows(), n = H_sparse.cols();
    MatrixXi H = MatrixXi::Zero(m, n);
    for (int k = 0; k < H_sparse.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(H_sparse, k); it; ++it)
            H(it.row(), it.col()) = it.value();
    int pivot_row = 0;
    vector<int> pivot_cols;
    vector<bool> is_pivot(n, false);
    for (int col = 0; col < n && pivot_row < m; ++col) {
        int sel = -1;
        for (int r = pivot_row; r < m; ++r)
            if (H(r, col)) {
                sel = r;
                break;
            }
        if (sel == -1)
            continue;
        if (sel != pivot_row)
            for (int k = 0; k < n; ++k)
                std::swap(H(pivot_row, k), H(sel, k));
        int iv = gf.inv(H(pivot_row, col));
        for (int k = col; k < n; ++k)
            H(pivot_row, k) = gf.mul(H(pivot_row, k), iv);
        for (int r = 0; r < m; ++r)
            if (r != pivot_row && H(r, col)) {
                int f = H(r, col);
                for (int k = col; k < n; ++k)
                    H(r, k) = gf.sub(H(r, k), gf.mul(f, H(pivot_row, k)));
            }
        pivot_cols.push_back(col);
        is_pivot[col] = true;
        pivot_row++;
    }
    SparseMatrix<int> G(n - (int)pivot_cols.size(), n);
    vector<Triplet<int>> tG;
    int f_idx = 0;
    for (int j = 0; j < n; ++j) {
        if (is_pivot[j])
            continue;
        tG.emplace_back(f_idx, j, 1);
        for (int r = 0; r < (int)pivot_cols.size(); ++r) {
            int val = H(r, j);
            if (val)
                tG.emplace_back(f_idx, pivot_cols[r], val);
        }
        f_idx++;
    }
    G.setFromTriplets(tG.begin(), tG.end());
    return G;
}

SparseMatrix<int> kroneckerNB(const SparseMatrix<int> &A, const SparseMatrix<int> &B,
                              SimpleGF &gf) {
    SparseMatrix<int> C(A.rows() * B.rows(), A.cols() * B.cols());
    vector<Triplet<int>> t;
    for (int k = 0; k < A.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator itA(A, k); itA; ++itA)
            for (int l = 0; l < B.outerSize(); ++l)
                for (SparseMatrix<int>::InnerIterator itB(B, l); itB; ++itB) {
                    int val = gf.mul(itA.value(), itB.value());
                    if (val)
                        t.emplace_back(itA.row() * B.rows() + itB.row(),
                                       itA.col() * B.cols() + itB.col(), val);
                }
    C.setFromTriplets(t.begin(), t.end());
    return C;
}

// ==========================================
// 3. Main 仅生成与保存
// ==========================================
int main() {
    srand(time(0));
    SimpleGF gf(6);
    cout << "=== Quantum LDPC Matrix Generator (H2 30x90 Version) ===\n";

    // H1/H3 使用 TIT 算法 (J=3, L=7, P=31)
    int J = 3, L = 7, P = 31;

    // H2 使用 PEG 算法 (增大到 30x90 以规避 4 环)
    int m2 = 30, n2 = 90;

    cout << "  [1/3] Constructing H1, H3 (TIT) and H2, G2 (PEG)..." << endl;
    pair<MatrixXi, MatrixXi> tit = constructBinaryHMatrices(J, L, P, 2, 3);

    SparseMatrix<int> H1 = denseToSparse(tit.first, gf);
    SparseMatrix<int> H3 = denseToSparse(tit.second, gf);
    SparseMatrix<int> H2 = buildH_PEG_Irregular(m2, n2, gf);
    SparseMatrix<int> G2 = computeGeneratorMatrix(H2, gf);

    // 2. 构造量子校验矩阵 HX, HZ
    cout << "  [2/3] Constructing HX and HZ (Hypergraph Product)..." << endl;
    int n1 = H1.cols();
    SparseMatrix<int> In1(n1, n1);
    In1.setIdentity();
    SparseMatrix<int> In2(n2, n2);
    In2.setIdentity();

    // HX = [ H1 (x) In2 ; In1 (x) H2 ]
    SparseMatrix<int> HX_top = kroneckerNB(H1, In2, gf);
    SparseMatrix<int> HX_bot = kroneckerNB(In1, H2, gf);
    SparseMatrix<int> HX(HX_top.rows() + HX_bot.rows(), HX_top.cols());
    vector<Triplet<int>> tX;
    for (int k = 0; k < HX_top.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(HX_top, k); it; ++it)
            tX.emplace_back(it.row(), it.col(), it.value());
    int r_off = HX_top.rows();
    for (int k = 0; k < HX_bot.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(HX_bot, k); it; ++it)
            tX.emplace_back(it.row() + r_off, it.col(), it.value());
    HX.setFromTriplets(tX.begin(), tX.end());

    // HZ = H3 (x) G2
    SparseMatrix<int> HZ = kroneckerNB(H3, G2, gf);

    // 3. 导出文件
    cout << "  [3/3] Exporting matrices..." << endl;
    saveToEssaiFormat(H1, "H1_sparse.txt", gf);
    saveToEssaiFormat(H3, "H3_sparse.txt", gf);
    saveToEssaiFormat(HX, "HX_Safe_GF64.txt", gf);
    saveToEssaiFormat(HZ, "HZ_Safe_GF64.txt", gf);

    cout << "\n=== Done ===\n";
    cout << "HX Dimension: " << HX.cols() << "x" << HX.rows() << endl;
    cout << "HZ Dimension: " << HZ.cols() << "x" << HZ.rows() << endl;
    return 0;
}

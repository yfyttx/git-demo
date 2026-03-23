/**
 * main.cpp - TIT 架构矩阵构造工具
 * 功能：生成 H1(TIT) 和 H2(PEG)，并计算导出 G2 (H2的生成矩阵)
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
        int prim_poly = (m == 6) ? 0x43 : 0x13; // GF(64) 多项式
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
        if (val <= 0 || val >= size)
            return -1;
        return log_table[val];
    }
};

// ==========================================
// 2. 导出与构造工具
// ==========================================

// 导出为 Essai 格式 (译码器读取格式)
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

// 将 Eigen 稠密矩阵转为非二进制稀疏矩阵
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

// 使用 PEG 算法构造 H2 矩阵
SparseMatrix<int> buildH2_PEG(int M, int N, SimpleGF &gf) {
    vector<int> degSeq(N, 3); // 固定列重为 3
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

// 高斯消元法计算 H 的生成矩阵 G (零空间)
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
            H.row(pivot_row).swap(H.row(sel));

        int iv = gf.inv(H(pivot_row, col));
        for (int k = col; k < n; ++k)
            H(pivot_row, k) = gf.mul(H(pivot_row, k), iv);

        for (int r = 0; r < m; ++r) {
            if (r != pivot_row && H(r, col)) {
                int f = H(r, col);
                for (int k = col; k < n; ++k)
                    H(r, k) = gf.sub(H(r, k), gf.mul(f, H(pivot_row, k)));
            }
        }
        pivot_cols.push_back(col);
        is_pivot[col] = true;
        pivot_row++;
    }

    // 构造 G 矩阵，使得 H * G^T = 0
    int dim_k = n - (int)pivot_cols.size();
    SparseMatrix<int> G(dim_k, n);
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

// ==========================================
// 3. Main 函数
// ==========================================
int main() {
    srand(time(0));
    SimpleGF gf(6); // GF(64)
    cout << "=== TIT Matrix Generator (Syndrome-Splitting Version) ===\n";

    // H1 参数 (TIT 构造: J=3, L=7, P=31)
    int J = 3, L = 7, P = 31;
    // H2 参数 (PEG 构造: 20x60)
    int m2 = 20, n2 = 60;

    cout << "  [1/2] Constructing H1 (TIT) and H2 (PEG)..." << endl;

    // 生成 H1, H3
    pair<MatrixXi, MatrixXi> tit = constructBinaryHMatrices(J, L, P, 2, 3);
    SparseMatrix<int> H1 = denseToSparse(tit.first, gf);

    // 生成 H2 并计算其生成矩阵 G2
    SparseMatrix<int> H2 = buildH2_PEG(m2, n2, gf);
    cout << "  [2/2] Computing G2 (Null space of H2)..." << endl;
    SparseMatrix<int> G2 = computeGeneratorMatrix(H2, gf);

    // 导出文件
    saveToEssaiFormat(H1, "H1_sparse.txt", gf);
    saveToEssaiFormat(H2, "H2_sparse.txt", gf);
    saveToEssaiFormat(G2, "G2_sparse.txt", gf);

    cout << "\n=== Done ===\n";
    cout << "H1 size: " << H1.cols() << "x" << H1.rows() << endl;
    cout << "H2 size: " << H2.cols() << "x" << H2.rows() << endl;
    cout << "G2 (Subblocks) size: " << G2.cols() << "x" << G2.rows() << endl;
    cout << "TIT Split Decoding will run " << G2.rows() << " parallel H1 decoders." << endl;

    return 0;
}

/**
 * main.cpp
 * * QLDPC Constructor: Hybrid (PEG + TIT QC-LDPC)
 * *
 * * Part 1: H1 (PEG) -> G1 (Kernel)
 * * Part 2: H2, H3 (TIT QC-LDPC Pair, where H2*H3^T = 0)
 * *
 * * HX = [ H1 (x) In2 ]
 * * [ In1 (x) H2 ]
 * *
 * * HZ = [ G1 (x) H3 ]
 */

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

// 引入 PEG 和 TIT 构造头文件
#include "../peg-unige/BigGirth.h"
#include "../peg-unige/Random.h"
#include "ldpc_codes/binary_codes.h" // TIT 构造函数在这里

using namespace Eigen;
using namespace std;

// ==========================================
// 1. 通用 GF(64) 类
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
    int div(int a, int b) const {
        if (!b)
            exit(1);
        if (!a)
            return 0;
        int d = log_table[a] - log_table[b];
        while (d < 0)
            d += size - 1;
        return exp_table[d];
    }
    int inv(int a) const {
        return div(1, a);
    }
    int exp(int k) const {
        while (k < 0)
            k += (size - 1);
        return exp_table[k % (size - 1)];
    }
    int log(int val) const {
        return (val <= 0 || val >= size) ? -1 : log_table[val];
    }
    int sub(int a, int b) const {
        return a ^ b;
    }
};

// ==========================================
// 2. 辅助工具：矩阵转换与保存
// ==========================================

// 将 Eigen::MatrixXi (稠密) 转换为 SparseMatrix<int>，并随机赋予 GF 系数
SparseMatrix<int> denseToSparse(const MatrixXi &dense, SimpleGF &gf) {
    int rows = dense.rows();
    int cols = dense.cols();
    SparseMatrix<int> sparse(rows, cols);
    vector<Triplet<int>> t;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (dense(i, j) == 1) {
                // TIT 生成的是二进制矩阵 (0/1)
                // 这里我们给它赋予一个 GF(64) 的随机非零系数，以构造成非二进制码
                int val = gf.exp(rand() % (gf.getSize() - 1));
                t.emplace_back(i, j, val);
            }
        }
    }
    sparse.setFromTriplets(t.begin(), t.end());
    return sparse;
}

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

// ==========================================
// 3. 构造与计算逻辑
// ==========================================

// PEG 构造 H1
SparseMatrix<int> buildH_PEG(int M, int N, int dv, SimpleGF &gf) {
    cout << "  [PEG] Generating H1 (M=" << M << ", N=" << N << ")..." << endl;
    vector<int> degSeq(N, dv);
    BigGirth bg(M, N, degSeq.data(), "dummy.dat", 1, 10000, false);
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

// 计算生成矩阵 G1 (H1 的核)
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
    int rank = pivot_cols.size();
    SparseMatrix<int> G(n - rank, n);
    vector<Triplet<int>> tG;
    int f_idx = 0;
    for (int j = 0; j < n; ++j) {
        if (is_pivot[j])
            continue;
        tG.emplace_back(f_idx, j, 1);
        for (int r = 0; r < rank; ++r) {
            int val = H(r, j);
            if (val)
                tG.emplace_back(f_idx, pivot_cols[r], val);
        }
        f_idx++;
    }
    G.setFromTriplets(tG.begin(), tG.end());
    return G;
}

// 克罗内克积
SparseMatrix<int> kroneckerNB(const SparseMatrix<int> &A, const SparseMatrix<int> &B,
                              SimpleGF &gf) {
    SparseMatrix<int> C(A.rows() * B.rows(), A.cols() * B.cols());
    vector<Triplet<int>> t;
    t.reserve((long long)A.nonZeros() * B.nonZeros());
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
// 4. Main
// ==========================================
int main() {
    srand(time(0));
    SimpleGF gf(6); // GF(64)
    cout << "=== QLDPC: Hybrid (PEG + TIT) Construction ===\n";

    // --------------------------------------------------------
    // Step 1: 生成 H1 (PEG) 并计算 G1
    // -------------------------------------------------k-------
    // 参数: 30行 60列 (列重3)
    SparseMatrix<int> H1 = buildH_PEG(60, 120, 3, gf);
    SparseMatrix<int> G1 = computeGeneratorMatrix(H1, gf);

    // --------------------------------------------------------
    // Step 2: 使用 TIT 算法生成 H2 和 H3
    // --------------------------------------------------------
    cout << "  [TIT] Constructing H2 and H3 using binary_codes.cpp logic...\n";
    // 参数说明 (J, L, P, sigma, tau)
    // J=3, L=5, P=10 --> H2 大小 (3*10)x(5*10) = 30x50
    //                    H3 大小 (3*10)x(5*10) = 30x50
    // 你可以调整 P 来改变码长，比如 P=20 -> 60x100
    int J = 3, L = 5, P = 163;
    int sigma = 2, tau = 3; // TIT 论文中的移位参数

    // 调用库函数生成一对正交矩阵 (Hc, Hd)
    pair<MatrixXi, MatrixXi> tit_matrices = constructBinaryHMatrices(J, L, P, sigma, tau);

    // 转换并赋予随机 GF 系数
    // H2 = Hc
    SparseMatrix<int> H2 = denseToSparse(tit_matrices.first, gf);
    // H3 = Hd (注意：TIT 保证 Hc * Hd^T = 0，所以 H3 可以直接作为 H2 的正交对偶)
    SparseMatrix<int> H3 = denseToSparse(tit_matrices.second, gf);

    cout << "  [TIT] Generated H2: " << H2.rows() << "x" << H2.cols() << endl;
    cout << "  [TIT] Generated H3: " << H3.rows() << "x" << H3.cols() << endl;

    // --------------------------------------------------------
    // Step 3: 准备单位阵 (用于对齐)
    // --------------------------------------------------------
    // H1 (30x60) -> n1 = 60
    // H2 (30x50) -> n2 = 50 (如果 P=10, L=5)
    int n1 = H1.cols();
    int n2 = H2.cols();

    SparseMatrix<int> In2(n2, n2);
    for (int i = 0; i < n2; ++i)
        In2.insert(i, i) = 1;
    SparseMatrix<int> In1(n1, n1);
    for (int i = 0; i < n1; ++i)
        In1.insert(i, i) = 1;

    // --------------------------------------------------------
    // Step 4: 构造 HX (上下拼接)
    // --------------------------------------------------------
    // Part A: H1 (x) In2 -> (30x60) x (50x50) = 1500 x 3000
    // Part B: In1 (x) H2 -> (60x60) x (30x50) = 1800 x 3000
    // 列数 3000 = n1 * n2，完美对齐！

    cout << "  [Build] HX (Stacked)..." << endl;
    SparseMatrix<int> HX_top = kroneckerNB(H1, In2, gf);
    SparseMatrix<int> HX_bot = kroneckerNB(In1, H2, gf);

    SparseMatrix<int> HX(HX_top.rows() + HX_bot.rows(), HX_top.cols());
    vector<Triplet<int>> tX;
    for (int k = 0; k < HX_top.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(HX_top, k); it; ++it)
            tX.emplace_back(it.row(), it.col(), it.value());

    int row_off = HX_top.rows();
    for (int k = 0; k < HX_bot.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(HX_bot, k); it; ++it)
            tX.emplace_back(it.row() + row_off, it.col(), it.value());
    HX.setFromTriplets(tX.begin(), tX.end());

    // --------------------------------------------------------
    // Step 5: 构造 HZ (G1 x H3)
    // --------------------------------------------------------
    // G1 (30x60) x H3 (30x50) -> 900 x 3000
    // 列数 3000，与 HX 完美一致！

    cout << "  [Build] HZ (G1 x H3)..." << endl;
    SparseMatrix<int> HZ = kroneckerNB(G1, H3, gf);

    cout << "  [Result] HX: " << HX.rows() << "x" << HX.cols() << endl;
    cout << "  [Result] HZ: " << HZ.rows() << "x" << HZ.cols() << endl;

    // --------------------------------------------------------
    // Step 6: 保存
    // --------------------------------------------------------
    saveToEssaiFormat(HX, "HX_Safe_GF64.txt", gf);
    saveToEssaiFormat(HZ, "HZ_Safe_GF64.txt", gf);

    return 0;
}

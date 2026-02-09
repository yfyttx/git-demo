/**
 * main.cpp
 * * QLDPC Constructor for Essai Simulator (Non-Binary)
 * * Output format: N M GF (Header) + Index/Value pairs
 */

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

// PEG 头文件
#include "../peg-unige/BigGirth.h"
#include "../peg-unige/Random.h"

using namespace Eigen;

// ==========================================
// 1. GF(16) 算术与对数表
// ==========================================
class SimpleGF16 {
    // GF(2^4) primitive polynomial: x^4 + x + 1
    int exp_table[32] = {1, 2, 4, 8, 3, 6,  12, 11, 5,  10, 7,  14, 15, 13, 9, 1,
                         2, 4, 8, 3, 6, 12, 11, 5,  10, 7,  14, 15, 13, 9,  1, 2};
    int log_table[16] = {-1, 0, 1, 4, 2, 8, 5, 10, 3, 14, 9, 7, 6, 13, 11, 12};

  public:
    int getPower() const {
        return 4;
    }
    int getSize() const {
        return 16;
    }

    int mul(int a, int b) const {
        if (a == 0 || b == 0)
            return 0;
        return exp_table[(log_table[a] + log_table[b]) % 15];
    }

    int exp(int k) const {
        while (k < 0)
            k += 15;
        return exp_table[k % 15];
    }

    // 获取元素的指数形式 (用于写入文件)
    // 假设 essai 需要的是幂次 p (其中 val = alpha^p)
    int log(int val) const {
        if (val <= 0 || val >= 16)
            return 0;
        return log_table[val];
    }
};

// ==========================================
// 2. 专用保存函数 (Essai 格式)
// ==========================================
void saveToEssaiFormat(const SparseMatrix<int> &H, const std::string &filename, SimpleGF16 &gf) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Cannot open " << filename << "\n";
        return;
    }

    int M = H.rows(); // 校验节点数
    int N = H.cols(); // 变量节点数
    int GF = 16;      // 我们固定使用 GF(16)

    // 1. Header: N M GF
    out << N << " " << M << " " << GF << "\n";

    // 计算度分布
    std::vector<int> dv(N, 0);
    std::vector<int> dc(M, 0);
    for (int k = 0; k < H.outerSize(); ++k) {
        for (SparseMatrix<int>::InnerIterator it(H, k); it; ++it) {
            dv[it.col()]++;
            dc[it.row()]++;
        }
    }

    // 2. Variable Degrees
    for (int i = 0; i < N; ++i)
        out << dv[i] << " ";
    out << "\n";

    // 3. Check Degrees
    for (int i = 0; i < M; ++i)
        out << dc[i] << " ";
    out << "\n";

    // 4. Matrix Data (Row-wise: Index Value Index Value ...)
    // 注意：Essai 需要 1-based Index
    // Value 通常是元素的幂次 (0..14)，或者整数值。
    // 根据 init.c 的逻辑 "code->matValue[m][k] = temp_int + 1"，推测文件里存的是幂次 p
    // 使得存储值为 p+1。我们这里写入幂次 (log value)。

    // 为了按行写入，我们需要转置或者构建行列表
    std::vector<std::vector<std::pair<int, int>>> rows(M);
    for (int k = 0; k < H.outerSize(); ++k) { // H is column-major by default in Eigen
        for (SparseMatrix<int>::InnerIterator it(H, k); it; ++it) {
            int r = it.row();
            int c = it.col();
            int val = it.value();
            int power = gf.log(val);           // 获取幂次 0..14
            rows[r].push_back({c + 1, power}); // 1-based column index
        }
    }

    for (int i = 0; i < M; ++i) {
        for (auto &p : rows[i]) {
            out << p.first << " " << p.second << "   ";
        }
        out << "\n";
    }

    out.close();
    std::cout << "  -> Saved to " << filename << " (Essai format: N M GF)\n";
}

// ==========================================
// 3. 构造逻辑 (同前，但更精简)
// ==========================================
SparseMatrix<int> buildH1_via_OfficialPEG(int M, int N, int dv, SimpleGF16 &gf) {
    std::cout << "  [PEG] Generating H1 (M=" << M << ", N=" << N << ")..." << std::endl;
    std::vector<int> degSeq(N, dv);
    BigGirth bg(M, N, degSeq.data(), "dummy.dat", 1, 10000, false);
    bg.loadH();

    SparseMatrix<int> H1(M, N);
    std::vector<Triplet<int>> triplets;
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            if (bg.H[i][j] == 1) {
                int val = gf.exp(rand() % 15);
                triplets.emplace_back(i, j, val);
            }
        }
    }
    H1.setFromTriplets(triplets.begin(), triplets.end());
    return H1;
}

SparseMatrix<int> kroneckerNB(const SparseMatrix<int> &A, const SparseMatrix<int> &B,
                              SimpleGF16 &gf) {
    SparseMatrix<int> C(A.rows() * B.rows(), A.cols() * B.cols());
    std::vector<Triplet<int>> t;
    long long est = (long long)A.nonZeros() * B.nonZeros();
    if (est > 50000000)
        std::cout << "  [Warn] Large Matrix: " << est << " entries\n";
    t.reserve(est);

    for (int k = 0; k < A.outerSize(); ++k) {
        for (SparseMatrix<int>::InnerIterator itA(A, k); itA; ++itA) {
            for (int l = 0; l < B.outerSize(); ++l) {
                for (SparseMatrix<int>::InnerIterator itB(B, l); itB; ++itB) {
                    int val = gf.mul(itA.value(), itB.value());
                    if (val != 0)
                        t.emplace_back(itA.row() * B.rows() + itB.row(),
                                       itA.col() * B.cols() + itB.col(), val);
                }
            }
        }
    }
    C.setFromTriplets(t.begin(), t.end());
    return C;
}

int main() {
    srand(time(0));
    SimpleGF16 gf;
    std::cout << "=== QLDPC Constructor: Essai Format ===\n";

    // 1. 生成 H1 (PEG) & H2 (Systematic)
    // 稍微调小一点参数以便快速测试
    SparseMatrix<int> H1 = buildH1_via_OfficialPEG(30, 60, 3, gf);

    // H2 (Systematic)
    int m2 = 5, n2 = 10;
    SparseMatrix<int> H2(m2, n2);
    std::vector<Triplet<int>> tH2;
    for (int c = 0; c < n2 - m2; ++c)
        for (int d = 0; d < 2; ++d)
            tH2.emplace_back(rand() % m2, c, gf.exp(rand() % 15)); // P
    for (int r = 0; r < m2; ++r)
        tH2.emplace_back(r, r + (n2 - m2), 1); // I
    H2.setFromTriplets(tH2.begin(), tH2.end());

    // 2. 构造 HX (Non-Binary)
    // HX = [ H1 (x) I ]
    //      [ I (x) H2 ]
    std::cout << "  [Build] Constructing HX_NB...\n";
    SparseMatrix<int> I1(H1.cols(), H1.cols());
    for (int i = 0; i < H1.cols(); ++i)
        I1.insert(i, i) = 1;
    SparseMatrix<int> I2(H2.cols(), H2.cols());
    for (int i = 0; i < H2.cols(); ++i)
        I2.insert(i, i) = 1;

    SparseMatrix<int> HX_top = kroneckerNB(H1, I2, gf);
    SparseMatrix<int> HX_bot = kroneckerNB(I1, H2, gf);

    SparseMatrix<int> HX(HX_top.rows() + HX_bot.rows(), HX_top.cols());
    std::vector<Triplet<int>> tX;
    for (int k = 0; k < HX_top.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(HX_top, k); it; ++it)
            tX.emplace_back(it.row(), it.col(), it.value());
    int off = HX_top.rows();
    for (int k = 0; k < HX_bot.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(HX_bot, k); it; ++it)
            tX.emplace_back(it.row() + off, it.col(), it.value());
    HX.setFromTriplets(tX.begin(), tX.end());

    // 3. 保存为 Essai 格式
    // 注意：这里保存的是 HX (Non-Binary)，而不是展开后的 Binary
    // 这样 Essai 就能正确识别 GF=16 并运行
    saveToEssaiFormat(HX, "HX_NB_Essai.txt", gf);

    return 0;
}

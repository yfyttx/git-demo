/**
 * main.cpp
 * * Revised for Generalized Shor Code (Quantum Product Code)
 * * Implements the construction:
 * * HX = [ H1 (x) I ]  (Vertical Stack)
 * * [ I (x) H2 ]
 * * HZ = G1 (x) G2
 * *
 * * This structure guarantees HX * HZ^T = 0.
 */

#include "gf2p/gf2p.h"
#include "ldpc_codes/binary_codes.h"
#include "ldpc_codes/nonbinary_codes.h"
#include "utils/matrix_utils.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

using namespace Eigen;

// ==========================================
// 辅助函数：非二进制到二进制的展开
// ==========================================
SparseMatrix<int> binaryExpand(const SparseMatrix<int> &H_nb, GF2p &gf) {
    int p = 0;
    int size = gf.getSize();
    while ((size >>= 1) > 0)
        p++;

    int M = H_nb.rows();
    int N = H_nb.cols();

    SparseMatrix<int> H_bin(M * p, N * p);
    std::vector<Triplet<int>> triplets;

    for (int k = 0; k < H_nb.outerSize(); ++k) {
        for (SparseMatrix<int>::InnerIterator it(H_nb, k); it; ++it) {
            int r = it.row();
            int c = it.col();
            int val = it.value();

            if (val != 0) {
                for (int col_bit = 0; col_bit < p; ++col_bit) {
                    int basis_elem = (1 << col_bit);
                    int product = gf.mul(basis_elem, val);
                    for (int row_bit = 0; row_bit < p; ++row_bit) {
                        if ((product >> row_bit) & 1) {
                            triplets.emplace_back(r * p + row_bit, c * p + col_bit, 1);
                        }
                    }
                }
            }
        }
    }
    H_bin.setFromTriplets(triplets.begin(), triplets.end());
    return H_bin;
}

// ==========================================
// 辅助函数：二进制张量积 (Kronecker Product)
// ==========================================
SparseMatrix<int> kroneckerBinary(const SparseMatrix<int> &A, const SparseMatrix<int> &B) {
    int m1 = A.rows(), n1 = A.cols();
    int m2 = B.rows(), n2 = B.cols();

    SparseMatrix<int> C(m1 * m2, n1 * n2);
    std::vector<Triplet<int>> triplets;

    // 预估非零元数量以防重新分配内存
    // 注意: 对于生成矩阵 G，这个数可能很大
    long long est_nnz = (long long)A.nonZeros() * B.nonZeros();
    // 简单限制一下 reserve 大小防止 bad_alloc (视内存情况调整)
    if (est_nnz < 500000000)
        triplets.reserve(est_nnz);

    for (int k = 0; k < A.outerSize(); ++k) {
        for (SparseMatrix<int>::InnerIterator itA(A, k); itA; ++itA) {
            int rowA = itA.row();
            int colA = itA.col();
            // A 是二进制矩阵，value 必然是 1

            for (int l = 0; l < B.outerSize(); ++l) {
                for (SparseMatrix<int>::InnerIterator itB(B, l); itB; ++itB) {
                    int rowB = itB.row();
                    int colB = itB.col();
                    triplets.emplace_back(rowA * m2 + rowB, colA * n2 + colB, 1);
                }
            }
        }
    }
    C.setFromTriplets(triplets.begin(), triplets.end());
    return C;
}

// ==========================================
// 核心函数：计算生成矩阵 G (从校验矩阵 H)
// 原理：求解 H * G^T = 0
// ==========================================
SparseMatrix<int> computeGeneratorMatrix(const SparseMatrix<int> &H_sparse) {
    // 转换为稠密矩阵进行高斯消元
    MatrixXi H = MatrixXi(H_sparse);
    int rows = H.rows();
    int cols = H.cols();

    std::vector<int> col_perm(cols);
    std::iota(col_perm.begin(), col_perm.end(), 0);

    int pivot_row = 0;
    std::vector<int> pivot_cols;

    // 1. 高斯消元 (RREF)
    for (int j = 0; j < cols && pivot_row < rows; ++j) {
        int sel = -1;
        for (int i = pivot_row; i < rows; ++i) {
            if (H(i, j) == 1) {
                sel = i;
                break;
            }
        }

        if (sel == -1)
            continue;

        // 交换行
        if (sel != pivot_row) {
            for (int k = 0; k < cols; ++k)
                std::swap(H(sel, k), H(pivot_row, k));
        }

        pivot_cols.push_back(j);

        // 消元
        for (int i = 0; i < rows; ++i) {
            if (i != pivot_row && H(i, j) == 1) {
                for (int k = j; k < cols; ++k)
                    H(i, k) ^= H(pivot_row, k);
            }
        }
        pivot_row++;
    }

    int rank = pivot_row;
    int k_dim = cols - rank; // 码的维数 k

    std::cout << "    [G-Calc] Matrix Rank: " << rank << ", Dimension k: " << k_dim << std::endl;

    if (k_dim == 0)
        return SparseMatrix<int>(0, cols);

    // 2. 构造基础解系 (G^T)
    std::vector<bool> is_pivot(cols, false);
    for (int p : pivot_cols)
        is_pivot[p] = true;

    std::vector<int> free_cols;
    for (int j = 0; j < cols; ++j) {
        if (!is_pivot[j])
            free_cols.push_back(j);
    }

    MatrixXi G_T(cols, k_dim);
    G_T.setZero();

    for (int i = 0; i < k_dim; ++i) {
        int free_col_idx = free_cols[i];

        // 自由变量设为 1
        G_T(free_col_idx, i) = 1;

        // 回代求主元变量
        // Pivot equations: x_pivot + sum(x_free) = 0  => x_pivot = sum(x_free)
        for (int r = 0; r < rank; ++r) {
            int p_col = pivot_cols[r];
            if (H(r, free_col_idx) == 1) {
                G_T(p_col, i) = 1;
            }
        }
    }

    // 3. 转置得到 G
    MatrixXi G = G_T.transpose();
    return G.sparseView();
}

int main() {
    try {
        srand(time(0));

        // ==========================================
        // 参数设置 (因为 G 稠密，建议先用小参数测试)
        // ==========================================
        int P = 5; // 较小的单位阵大小
        int J = 3; // 列重
        int L = 4; // 行重
        int sigma = 0;
        int tau = 1;
        int gf_power = 4; // GF(16)

        std::cout << "=== Generalized Shor Code Constructor ===\n";
        std::cout << "=== Structure: HX=[H1xI; IxH2], HZ=G1xG2 ===\n";

        // ==========================================
        // Step 1-4: 构造基础码 C1 (H1, G1)
        // ==========================================
        std::cout << "\n[Step 1] Constructing Base Code C1..." << std::endl;
        auto [Hc_bin, Hd_bin] = constructBinaryHMatrices(J, L, P, sigma, tau);

        // 初始化 GF
        std::vector<int> poly;
        if (gf_power == 4)
            poly = {1, 1, 0, 0, 1}; // x^4+x+1
        else
            poly = {1, 1, 0, 0, 0, 0, 1}; // x^6+x+1 (Default)
        GF2p gf(gf_power, poly);

        // 构造非二进制 H
        std::vector<std::pair<int, int>> Hc_nonzero = getNonZeroPositions(Hc_bin);
        MatrixXi linear_system = buildLinearSystem(Hc_nonzero, Hc_bin, Hd_bin, gf.getMod());
        std::vector<int> solution = gaussElimination(linear_system, gf.getMod());

        SparseMatrix<int> H1_NB(Hc_bin.rows(), Hc_bin.cols());
        std::vector<Triplet<int>> h1_triplets;
        for (size_t i = 0; i < Hc_nonzero.size(); ++i) {
            auto [row, col] = Hc_nonzero[i];
            int log_val = (!solution.empty()) ? solution[i] : (rand() % (gf.getMod() - 1) + 1);
            if (gf.exp(log_val) != 0)
                h1_triplets.emplace_back(row, col, gf.exp(log_val));
        }
        H1_NB.setFromTriplets(h1_triplets.begin(), h1_triplets.end());

        // H1: Binary Expansion
        SparseMatrix<int> H1 = binaryExpand(H1_NB, gf);
        std::cout << "  -> H1 Size: " << H1.rows() << "x" << H1.cols() << std::endl;

        // G1: Compute Generator Matrix
        std::cout << "  -> Computing G1..." << std::endl;
        SparseMatrix<int> G1 = computeGeneratorMatrix(H1);
        std::cout << "  -> G1 Size: " << G1.rows() << "x" << G1.cols() << std::endl;

        // ==========================================
        // Step 5: 准备 C2 (这里设 C2 = C1)
        // ==========================================
        SparseMatrix<int> H2 = H1;
        SparseMatrix<int> G2 = G1;

        // ==========================================
        // Step 6: 构造量子码矩阵 HX, HZ
        // ==========================================
        std::cout << "\n[Step 6] Constructing Quantum Code Matrices..." << std::endl;

        int n1 = H1.cols();
        int n2 = H2.cols();
        SparseMatrix<int> I_n1(n1, n1);
        I_n1.setIdentity();
        SparseMatrix<int> I_n2(n2, n2);
        I_n2.setIdentity();

        // 1. 构造 HX (Vertical Stacking)
        // HX = [ H1 (x) I_n2 ]
        //      [ I_n1 (x) H2 ]
        std::cout << "  -> Building HX (Vertical Stack)..." << std::endl;
        SparseMatrix<int> HX_top = kroneckerBinary(H1, I_n2);
        SparseMatrix<int> HX_bot = kroneckerBinary(I_n1, H2);

        // 拼接 Triplet
        std::vector<Triplet<int>> tX;
        tX.reserve(HX_top.nonZeros() + HX_bot.nonZeros());

        // Top Part
        for (int k = 0; k < HX_top.outerSize(); ++k)
            for (SparseMatrix<int>::InnerIterator it(HX_top, k); it; ++it)
                tX.emplace_back(it.row(), it.col(), it.value());

        // Bottom Part (Row index offset by HX_top.rows())
        int row_offset = HX_top.rows();
        for (int k = 0; k < HX_bot.outerSize(); ++k)
            for (SparseMatrix<int>::InnerIterator it(HX_bot, k); it; ++it)
                tX.emplace_back(it.row() + row_offset, it.col(), it.value());

        SparseMatrix<int> HX(HX_top.rows() + HX_bot.rows(), HX_top.cols());
        HX.setFromTriplets(tX.begin(), tX.end());

        // 2. 构造 HZ = G1 (x) G2
        std::cout << "  -> Building HZ = G1 (x) G2 (Warning: Dense)..." << std::endl;
        SparseMatrix<int> HZ = kroneckerBinary(G1, G2);

        std::cout << "  -> Final Sizes: HX(" << HX.rows() << "x" << HX.cols() << "), HZ("
                  << HZ.rows() << "x" << HZ.cols() << ")" << std::endl;

        // ==========================================
        // Step 7: 保存矩阵
        // ==========================================
        saveToKN(HX, "HX_binary.txt", gf);
        saveToKN(HZ, "HZ_binary.txt", gf);
        std::cout << "  -> Saved to HX_binary.txt and HZ_binary.txt" << std::endl;

        // ==========================================
        // Step 8: 验证正交性
        // ==========================================
        std::cout << "\n[Step 8] Verifying Orthogonality (HX * HZ^T)..." << std::endl;

        // Product = HX * HZ^T
        // Using sparse multiplication. Result should be 0 mod 2.
        SparseMatrix<int> Product = HX * HZ.transpose();

        bool isOrthogonal = true;
        for (int k = 0; k < Product.outerSize(); ++k) {
            for (SparseMatrix<int>::InnerIterator it(Product, k); it; ++it) {
                if (it.value() % 2 != 0) {
                    isOrthogonal = false;
                    break;
                }
            }
            if (!isOrthogonal)
                break;
        }

        if (isOrthogonal) {
            std::cout << "  RESULT: \033[1;32mPASSED\033[0m (Matrices are Orthogonal)\n";
        } else {
            std::cout << "  RESULT: \033[1;31mFAILED\033[0m (Product is not zero matrix)\n";
        }

    } catch (const std::exception &e) {
        std::cerr << "\n[ERROR] " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

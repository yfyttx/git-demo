/**
 * main.cpp
 * * Revised for Binary Hypergraph Product Construction
 * * Strategy: Expand H1 to Binary FIRST, then construct HGP.
 * * This guarantees HX * HZ^T = 0 mod 2.
 */

#include "gf2p/gf2p.h"
#include "ldpc_codes/binary_codes.h"
#include "ldpc_codes/nonbinary_codes.h"
#include "utils/matrix_utils.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <iostream>
#include <vector>

using namespace Eigen;

// ==========================================
// 辅助函数：非二进制到二进制的展开 (Binary Image Expansion)
// ==========================================
SparseMatrix<int> binaryExpand(const SparseMatrix<int> &H_nb, GF2p &gf) {
    // 自动计算 p (根据 gf 大小)
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
            int val = it.value(); // GF element

            if (val != 0) {
                // 生成 val 对应的 p*p 二进制矩阵块 (列向量基底展开)
                for (int col_bit = 0; col_bit < p; ++col_bit) {
                    int basis_elem = (1 << col_bit);
                    int product = gf.mul(basis_elem, val); // GF乘法

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
// 辅助函数：二进制张量积 (Binary Kronecker Product)
// ==========================================
SparseMatrix<int> kroneckerBinary(const SparseMatrix<int> &A, const SparseMatrix<int> &B) {
    int m1 = A.rows(), n1 = A.cols();
    int m2 = B.rows(), n2 = B.cols();

    SparseMatrix<int> C(m1 * m2, n1 * n2);
    std::vector<Triplet<int>> triplets;

    for (int k = 0; k < A.outerSize(); ++k) {
        for (SparseMatrix<int>::InnerIterator itA(A, k); itA; ++itA) {
            int rowA = itA.row();
            int colA = itA.col();
            // A中非零元素必然是1 (二进制矩阵)

            // 复制整个矩阵 B 到对应块
            for (int l = 0; l < B.outerSize(); ++l) {
                for (SparseMatrix<int>::InnerIterator itB(B, l); itB; ++itB) {
                    int rowB = itB.row();
                    int colB = itB.col();

                    // 位置映射: Row = rA*m2 + rB
                    triplets.emplace_back(rowA * m2 + rowB, colA * n2 + colB, 1);
                }
            }
        }
    }
    C.setFromTriplets(triplets.begin(), triplets.end());
    return C;
}

int main() {
    try {
        // ==========================================
        // 参数设置
        // ==========================================
        int P = 7;
        int J = 2;
        int L = 6;
        int sigma = 2;
        int tau = 3;
        int gf_power = 6; // GF(64)

        std::cout << "=== Non-Binary Seed -> Binary HGP Constructor ===\n";

        // ==========================================
        // Step 1: 构造二进制骨架
        // ==========================================
        std::cout << "[Step 1] Constructing QC-LDPC Skeleton..." << std::endl;
        auto [Hc_bin, Hd_bin] = constructBinaryHMatrices(J, L, P, sigma, tau);

        // ==========================================
        // Step 2: 初始化 GF
        // ==========================================
        std::vector<int> irreducible_poly = {1, 1, 0, 0, 0, 0, 1};
        GF2p gf(gf_power, irreducible_poly);

        // ==========================================
        // Step 3: 构造非二进制 H1
        // ==========================================
        std::cout << "[Step 3] Constructing Seed Matrix H1 (Non-Binary)..." << std::endl;
        std::vector<std::pair<int, int>> Hc_nonzero = getNonZeroPositions(Hc_bin);
        MatrixXi linear_system = buildLinearSystem(Hc_nonzero, Hc_bin, Hd_bin, gf.getMod());
        std::vector<int> solution = gaussElimination(linear_system, gf.getMod());

        SparseMatrix<int> H1_NB(Hc_bin.rows(), Hc_bin.cols());
        std::vector<Triplet<int>> h1_triplets;
        for (size_t i = 0; i < Hc_nonzero.size(); ++i) {
            auto [row, col] = Hc_nonzero[i];
            int log_val = (!solution.empty()) ? solution[i] : (rand() % (gf.getMod() - 1) + 1);
            int poly_val = gf.exp(log_val);
            if (poly_val != 0)
                h1_triplets.emplace_back(row, col, poly_val);
        }
        H1_NB.setFromTriplets(h1_triplets.begin(), h1_triplets.end());
        std::cout << "  -> H1_NB created: " << H1_NB.rows() << "x" << H1_NB.cols() << std::endl;

        // ==========================================
        // Step 4: 立即展开为二进制矩阵 (关键修改!)
        // ==========================================
        std::cout << "\n[Step 4] Expanding H1 to Binary FIRST..." << std::endl;
        SparseMatrix<int> H1 = binaryExpand(H1_NB, gf);
        SparseMatrix<int> H2 = H1; // 对称构造
        std::cout << "  -> H1 Binary Size: " << H1.rows() << "x" << H1.cols() << std::endl;

        // ==========================================
        // Step 5: 构造超图积 (使用二进制矩阵)
        // H_X = [H1 (x) I,  I (x) H2^T]
        // H_Z = [I (x) H2,  H1^T (x) I]
        // ==========================================
        std::cout << "\n[Step 5] Constructing Hypergraph Product (Binary Domain)..." << std::endl;

        int r1 = H1.rows(), n1 = H1.cols();
        int r2 = H2.rows(), n2 = H2.cols();

        SparseMatrix<int> I_n1(n1, n1);
        I_n1.setIdentity();
        SparseMatrix<int> I_r1(r1, r1);
        I_r1.setIdentity();
        SparseMatrix<int> I_n2(n2, n2);
        I_n2.setIdentity();
        SparseMatrix<int> I_r2(r2, r2);
        I_r2.setIdentity();

        SparseMatrix<int> H1_T = H1.transpose();
        SparseMatrix<int> H2_T = H2.transpose();

        std::cout << "  -> Computing Kronecker products..." << std::endl;
        // 注意：这里使用 kroneckerBinary，不再涉及 GF 运算
        SparseMatrix<int> A = kroneckerBinary(H1, I_n2);
        SparseMatrix<int> B = kroneckerBinary(I_r1, H2_T);
        SparseMatrix<int> C = kroneckerBinary(I_n1, H2);
        SparseMatrix<int> D = kroneckerBinary(H1_T, I_r2);

        // 拼接
        std::cout << "  -> Assembling HX and HZ..." << std::endl;
        SparseMatrix<int> HX(A.rows(), A.cols() + B.cols());
        SparseMatrix<int> HZ(C.rows(), C.cols() + D.cols());

        // 简单的 Triplet 拼接
        std::vector<Triplet<int>> tX, tZ;

        // HX = [A | B]
        for (int k = 0; k < A.outerSize(); ++k)
            for (SparseMatrix<int>::InnerIterator it(A, k); it; ++it)
                tX.emplace_back(it.row(), it.col(), it.value());
        for (int k = 0; k < B.outerSize(); ++k)
            for (SparseMatrix<int>::InnerIterator it(B, k); it; ++it)
                tX.emplace_back(it.row(), it.col() + A.cols(), it.value());

        // HZ = [C | D]
        for (int k = 0; k < C.outerSize(); ++k)
            for (SparseMatrix<int>::InnerIterator it(C, k); it; ++it)
                tZ.emplace_back(it.row(), it.col(), it.value());
        for (int k = 0; k < D.outerSize(); ++k)
            for (SparseMatrix<int>::InnerIterator it(D, k); it; ++it)
                tZ.emplace_back(it.row(), it.col() + C.cols(), it.value());

        HX.setFromTriplets(tX.begin(), tX.end());
        HZ.setFromTriplets(tZ.begin(), tZ.end());

        std::cout << "  -> Final Quantum Code Sizes: HX(" << HX.rows() << "x" << HX.cols()
                  << "), HZ(" << HZ.rows() << "x" << HZ.cols() << ")" << std::endl;

        // ==========================================
        // Step 6: 保存
        // ==========================================
        saveToKN(HX, "HX_binary.txt", gf);
        saveToKN(HZ, "HZ_binary.txt", gf);
        std::cout << "  -> Files saved." << std::endl;

        // ==========================================
        // Step 7: 验证正交性
        // ==========================================
        std::cout << "\n[Step 7] Verifying Commutativity (HX * HZ^T = 0)..." << std::endl;
        SparseMatrix<int> Product = HX * HZ.transpose();

        bool clean = true;
        for (int k = 0; k < Product.outerSize(); ++k) {
            for (SparseMatrix<int>::InnerIterator it(Product, k); it; ++it) {
                if (it.value() % 2 != 0) { // 模2检查
                    clean = false;
                    break;
                }
            }
            if (!clean)
                break;
        }

        if (clean)
            std::cout << "  RESULT: \033[1;32mPASSED\033[0m (Valid CSS Code)\n";
        else
            std::cout << "  RESULT: \033[1;31mFAILED\033[0m (Check logic)\n";

    } catch (const std::exception &e) {
        std::cerr << "\n[ERROR] Exception: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

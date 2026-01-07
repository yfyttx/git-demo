/**
 * main.cpp
 * * 基于论文 "Quantum Error Correction Beyond the Bounded Distance Decoding Limit"
 * 构造 Non-Binary QC-LDPC CSS 码 (Example 5: N=42, M=14)
 * * [修改说明]
 * 使用 saveToKN 替代 saveToALIST 以适配 NB_LDPC_FB_public 仿真器的特殊格式要求
 */

#include "gf2p/gf2p.h"
#include "ldpc_codes/binary_codes.h"
#include "ldpc_codes/nonbinary_codes.h"
#include "utils/matrix_utils.h" // 包含 saveToKN 的声明
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>

using namespace Eigen;

int main() {
    try {
        // ==========================================
        // 参数设置 (对应论文 Example 5)
        // ==========================================
        // N = P * L = 7 * 6 = 42
        // M = P * J = 7 * 2 = 14
        int P = 7;
        int J = 2;
        int L = 6;
        int sigma = 2;
        int tau = 3;
        int gf_power = 6; // GF(2^6) = GF(64)

        std::cout << "=== Non-Binary Quantum LDPC Code Constructor ===\n";
        std::cout << "Target: N=" << (P * L) << ", M=" << (P * J) << ", GF(2^" << gf_power
                  << ")\n\n";

        // ==========================================
        // Step 1: 构造二进制骨架 (Binary Skeleton)
        // ==========================================
        std::cout << "[Step 1] Constructing Binary QC Matrices..." << std::endl;
        auto [Hc_bin, Hd_bin] = constructBinaryHMatrices(J, L, P, sigma, tau);
        std::cout << "  -> Skeleton size: " << Hc_bin.rows() << "x" << Hc_bin.cols() << "\n";

        // ==========================================
        // Step 2: 初始化有限域 GF(64)
        // ==========================================
        std::cout << "[Step 2] Initializing GF(64)..." << std::endl;
        // 使用本原多项式 p(x) = 1 + x + x^6 (对应的二进制为 1000011 -> array: 1,1,0,0,0,0,1)
        std::vector<int> irreducible_poly = {1, 1, 0, 0, 0, 0, 1};
        GF2p gf(gf_power, irreducible_poly);

        // ==========================================
        // Step 3: 求解 H_gamma (系数赋值)
        // ==========================================
        std::cout << "[Step 3] Solving H_gamma coefficients..." << std::endl;

        std::vector<std::pair<int, int>> Hc_nonzero = getNonZeroPositions(Hc_bin);
        // 建立线性方程组以消除短环 (Cycle constraints)
        MatrixXi linear_system = buildLinearSystem(Hc_nonzero, Hc_bin, Hd_bin, gf.getMod());
        // 高斯消元求解系数 (在 log 域)
        std::vector<int> log_gamma_solution = gaussElimination(linear_system, gf.getMod());

        // 存储 H_gamma
        // 1. Sparse Matrix 用于保存文件
        SparseMatrix<int> H_gamma_sparse(Hc_bin.rows(), Hc_bin.cols());
        std::vector<Triplet<int>> gamma_triplets;
        gamma_triplets.reserve(Hc_nonzero.size());

        // 2. Dense Matrix 用于后续求解 H_delta 和验证
        MatrixXi H_gamma_vals_dense = MatrixXi::Zero(Hc_bin.rows(), Hc_bin.cols());

        for (size_t i = 0; i < Hc_nonzero.size(); ++i) {
            auto [row, col] = Hc_nonzero[i];

            // 如果高斯消元返回空(无约束)，给予随机值，否则使用解出的值
            int log_val = (!log_gamma_solution.empty()) ? log_gamma_solution[i]
                                                        : (rand() % (gf.getMod() - 1) + 1);

            int poly_val = gf.exp(log_val); // 转换为多项式值 (0..63)

            if (poly_val != 0) {
                gamma_triplets.emplace_back(row, col, poly_val);
                H_gamma_vals_dense(row, col) = poly_val;
            }
        }
        H_gamma_sparse.setFromTriplets(gamma_triplets.begin(), gamma_triplets.end());

        // [关键修改] 使用 saveToKN 保存兼容仿真器的格式
        // 注意：这里需要传入 gf 对象来进行 多项式值 -> 幂次值 的转换
        saveToKN(H_gamma_sparse, "H_gamma.txt", gf);
        std::cout << "  -> H_gamma saved (KN Format)." << std::endl;

        // ==========================================
        // Step 4: 求解 H_delta (正交性约束)
        // ==========================================
        std::cout << "[Step 4] Solving H_delta coefficients (Orthogonality)..." << std::endl;

        SparseMatrix<int> H_delta_sparse(Hd_bin.rows(), Hd_bin.cols());
        std::vector<Triplet<int>> delta_triplets;
        delta_triplets.reserve(Hd_bin.rows() * L * 2);

        // 创建 H_delta 的密集矩阵用于 Step 5 的验证
        MatrixXi H_delta_dense = MatrixXi::Zero(Hd_bin.rows(), Hd_bin.cols());

        for (int i = 0; i < Hd_bin.rows(); ++i) {
            // 进度条
            if (i % 10 == 0)
                std::cout << "\r  -> Processing row " << i << "/" << Hd_bin.rows() << std::flush;

            // 对每一行求解局部方程组
            std::vector<int> delta_row_vals =
                solveDeltaRow(i, Hc_bin, Hd_bin, H_gamma_vals_dense, gf);

            for (int j = 0; j < Hd_bin.cols(); ++j) {
                int val = delta_row_vals[j];
                if (val != 0) {
                    delta_triplets.emplace_back(i, j, val);
                    H_delta_dense(i, j) = val; // 填充密集矩阵用于验证
                }
            }
        }
        std::cout << "\n  -> H_delta solved." << std::endl;

        H_delta_sparse.setFromTriplets(delta_triplets.begin(), delta_triplets.end());

        // [关键修改] 使用 saveToKN 保存兼容仿真器的格式
        saveToKN(H_delta_sparse, "H_delta.txt", gf);
        std::cout << "  -> H_delta saved (KN Format)." << std::endl;

        // ==========================================
        // Step 5: 验证正交性 (Verification)
        // ==========================================
        std::cout << "[Step 5] Verifying H_gamma * H_delta^T == 0 ..." << std::endl;

        bool is_orthogonal = verifyNonbinaryOrthogonality(H_gamma_vals_dense, H_delta_dense, gf);

        if (is_orthogonal) {
            std::cout << "\n============================================\n";
            std::cout << "  RESULT: \033[1;32mPASSED\033[0m\n";
            std::cout << "  The matrices satisfy the CSS constraint.\n";
            std::cout << "  Files generated in KN format for NB_LDPC_FB_public simulator.\n";
            std::cout << "============================================\n";
        } else {
            std::cout << "\n============================================\n";
            std::cout << "  RESULT: \033[1;31mFAILED\033[0m\n";
            std::cout << "  The matrices are NOT orthogonal.\n";
            std::cout << "  Please check the Gaussian Elimination logic.\n";
            std::cout << "============================================\n";
        }

    } catch (const std::exception &e) {
        std::cerr << "\n[ERROR] Exception caught: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

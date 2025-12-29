#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse> // 必须引入
#include <fstream> 
#include "gf2p/gf2p.h"
#include "ldpc_codes/binary_codes.h"
#include "ldpc_codes/nonbinary_codes.h"
#include "utils/matrix_utils.h"

using namespace Eigen;

int main() {
    try {
        // --- 参数设置 ---
        // 注意：如果是 1e6 级别，P, L 需要相应调整
        int P = 7; 
        int J = 2;
        int L = 6;
        int sigma = 2;
        int tau = 3;

        std::cout << "=== Non-Binary Quantum LDPC Code Constructor (Sparse Mode) ===\n";
        
        // 1. 二进制骨架 (Step 1)
        auto [Hc_bin, Hd_bin] = constructBinaryHMatrices(J, L, P, sigma, tau);
        std::cout << "[Step 1] Skeleton constructed. Size: " << Hc_bin.rows() << "x" << Hc_bin.cols() << "\n";

        // 2. GF 初始化 (Step 2)
        int p = 6;
        std::vector<int> irreducible_poly = {1, 1, 0, 0, 0, 0, 1}; 
        GF2p gf(p, irreducible_poly);
        
        // 3. 构造 H_gamma (Step 3) - 使用 Triplet 构建稀疏矩阵
        std::cout << "[Step 3] Solving H_gamma...\n";
        std::vector<std::pair<int, int>> Hc_nonzero = getNonZeroPositions(Hc_bin);
        MatrixXi linear_system = buildLinearSystem(Hc_nonzero, Hc_bin, Hd_bin, gf.getMod());
        std::vector<int> log_gamma_solution = gaussElimination(linear_system, gf.getMod());

        // 使用 SparseMatrix 存储 H_gamma_vals
        // 预估非零元数量 = Hc_nonzero.size()
        SparseMatrix<int> H_gamma_sparse(Hc_bin.rows(), Hc_bin.cols());
        std::vector<Triplet<int>> gamma_triplets;
        gamma_triplets.reserve(Hc_nonzero.size());

        // 临时密集矩阵用于 solveDeltaRow 查询 (如果 RAM 允许)
        // 警告：对于 1e6，这步仍然危险。理想情况应重写 solveDeltaRow 接受 Sparse。
        // 但为了兼容现有 dense solver，暂时保留。如果 OOM，需重写 solver。
        MatrixXi H_gamma_vals_dense = MatrixXi::Zero(Hc_bin.rows(), Hc_bin.cols());

        for (size_t i = 0; i < Hc_nonzero.size(); ++i) {
            auto [row, col] = Hc_nonzero[i];
            int log_val = log_gamma_solution[i];
            int poly_val = gf.exp(log_val); // 转换为多项式值 (0..63)
            
            // 0 是空元素，非0是有效值
            if (poly_val != 0) {
                gamma_triplets.push_back(Triplet<int>(row, col, poly_val));
                H_gamma_vals_dense(row, col) = poly_val; 
            }
        }
        H_gamma_sparse.setFromTriplets(gamma_triplets.begin(), gamma_triplets.end());
        
        // 直接保存 H_gamma.alist (不生成 txt)
        saveToALIST(H_gamma_sparse, "H_gamma.alist");

        // 4. 构造 H_delta (Step 4)
        std::cout << "[Step 4] Solving H_delta...\n";
        SparseMatrix<int> H_delta_sparse(Hd_bin.rows(), Hd_bin.cols());
        std::vector<Triplet<int>> delta_triplets;
        // 预估每行非零元，大致预留内存
        delta_triplets.reserve(Hd_bin.rows() * L * 2); 

        for (int i = 0; i < Hd_bin.rows(); ++i) {
            if (i % 100 == 0) std::cout << "\rProcessing row " << i << "/" << Hd_bin.rows() << std::flush;
            
            // 这里传入 H_gamma_vals_dense。
            // 对于极大矩阵，这里是瓶颈。如果崩溃，需要将 solveDeltaRow 改为接受 SparseMatrix。
            std::vector<int> delta_row_vals = solveDeltaRow(i, Hc_bin, Hd_bin, H_gamma_vals_dense, gf);
            
            for (int j = 0; j < Hd_bin.cols(); ++j) {
                int val = delta_row_vals[j];
                if (val != 0) {
                    delta_triplets.push_back(Triplet<int>(i, j, val));
                }
            }
        }
        std::cout << "\n";
        H_delta_sparse.setFromTriplets(delta_triplets.begin(), delta_triplets.end());

        // 直接保存 H_delta.alist
        saveToALIST(H_delta_sparse, "H_delta.alist");

        std::cout << "Success! ALIST files generated directly.\n";

    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

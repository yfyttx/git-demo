#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <fstream> 
#include "gf2p/gf2p.h"
#include "ldpc_codes/binary_codes.h"
#include "ldpc_codes/nonbinary_codes.h"
#include "utils/matrix_utils.h"

using namespace Eigen;

// (saveMatrixToFile 函数保持不变)
void saveMatrixToFile(const Eigen::MatrixXi& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                file << matrix(i, j) << " "; 
            }
            file << "\n"; 
        }
        file.close();
        std::cout << "\nSuccessfully saved LOG matrix to " << filename << std::endl;
    } else {
        std::cerr << "\nError: Could not open file " << filename << " for writing." << std::endl;
    }
}

int main() {
    try {
        
        int P = 7;
        int J = 2;
        int L = 6;
        int sigma = 2;
        int tau = 3;

        std::cout << "Quantum LDPC Code Constructor\n";
        std::cout << "Parameters: P=" << P << ", J=" << J << ", L=" << L 
                  << ", sigma=" << sigma << ", tau=" << tau << std::endl;
        std::cout << "-------------------------------------------\n";

        
        auto [Hc_bin, Hd_bin] = constructBinaryHMatrices(J, L, P, sigma, tau);
        std::cout << "Step 1: Binary skeleton matrices H_C and H_D constructed.\n";

        
        int p = 6;
        // 确保这是 GF(2^6) 的一个有效本原多项式
        // x^6 + x + 1 (系数为 [1, 1, 0, 0, 0, 0, 1] )
        std::vector<int> irreducible_poly = {1, 1, 0, 0, 0, 0, 1}; 
        GF2p gf(p, irreducible_poly);
        std::cout << "Step 2: Galois Field GF(2^" << p << ") initialized.\n";
        
       
        std::vector<std::pair<int, int>> Hc_nonzero = getNonZeroPositions(Hc_bin);

        
        MatrixXi linear_system = buildLinearSystem(Hc_nonzero, Hc_bin, Hd_bin, gf.getMod());
        std::cout << "Step 3a: Linear system for H_gamma constructed.\n";

       
        std::vector<int> log_gamma_solution = gaussElimination(linear_system, gf.getMod());
        
        // --- 核心修改从这里开始 ---

        // H_gamma 存储 *多项式值* (用于验证)
        MatrixXi H_gamma_vals = MatrixXi::Zero(Hc_bin.rows(), Hc_bin.cols());
        // H_gamma_logs 存储 *指数值* (用于保存到文件)
        MatrixXi H_gamma_logs = MatrixXi::Zero(Hc_bin.rows(), Hc_bin.cols());

        for (size_t i = 0; i < Hc_nonzero.size(); ++i) {
            auto [row, col] = Hc_nonzero[i];
            int log_val = log_gamma_solution[i];
            
            // 存储多项式值 (0-63)
            H_gamma_vals(row, col) = gf.exp(log_val); 
            
            // 存储指数值 (0-62, 对应 alpha^0 到 alpha^62)
            // MATLAB 工具箱约定：0 代表零元素，1 代表 alpha^0, 2 代表 alpha^1 ...
            // 但论文图2约定：'0' 代表 alpha^0。
            // 让我们遵循 "log_val" 原样存储，并假设 0 仍然是 0 元素。
            // 不，最安全的约定是：0=零元素，log+1 = alpha^log
            // ...
            // 让我们做一个更简单的假设：
            // H_gamma.txt 中，0 = 零元素, 并且 (log_val) = alpha^(log_val)
            // 不对，论文图2说 log(alpha^0) = 0。所以 '0' 存储为 0。
            // 让我们就存储 log_gamma_solution[i]
            
            H_gamma_logs(row, col) = log_val;
        }
        std::cout << "Step 3b: Non-binary matrix H_gamma (values and logs) constructed.\n";
        
        
        // H_delta_vals 存储 *多项式值* (用于验证)
        MatrixXi H_delta_vals = MatrixXi::Zero(Hd_bin.rows(), Hd_bin.cols());
        // H_delta_logs 存储 *指数值* (用于保存到文件)
        MatrixXi H_delta_logs = MatrixXi::Zero(Hd_bin.rows(), Hd_bin.cols());

        for (int i = 0; i < Hd_bin.rows(); ++i) {
            // solveDeltaRow 应该返回多项式值，因为它在内部进行 gf 运算
            std::vector<int> delta_row_vals_vec = solveDeltaRow(i, Hc_bin, Hd_bin, H_gamma_vals, gf);
            
            for (int j = 0; j < Hd_bin.cols(); ++j) {
                int val = delta_row_vals_vec[j];
                H_delta_vals(i, j) = val;

                // 将多项式值转换回指数 (log) 以便保存
                if (val == 0) {
                    H_delta_logs(i, j) = 0; // 假设 0 = 零元素
                } else {
                    H_delta_logs(i, j) = gf.log(val); // 存储指数
                }
            }
        }
        std::cout << "Step 4: Non-binary matrix H_delta (values and logs) constructed.\n";

        
        // 验证 *必须* 使用 *多项式值* 矩阵
        bool is_orthogonal = verifyNonbinaryOrthogonality(H_gamma_vals, H_delta_vals, gf);
        std::cout << "\nFinal Check: Orthogonality (H_gamma * H_delta^T == 0): " 
                  << (is_orthogonal ? "SUCCESS" : "FAILURE") << "\n";

        // (printAsHexLog 函数现在可能具有误导性，但我们不再需要它)
        
        // --- 核心修改在这里 ---
        // 保存 *指数* 矩阵 (H_gamma_logs) 到文件
        // (请确保路径正确)
        saveMatrixToFile(H_gamma_logs, "C:\\Users\\yfyttx\\Documents\\code\\quantum computing\\NB_LDPC\\NB-LDPC-toolbox-master\\H_gamma.txt");
        
        // 保存 *指数* 矩阵 (H_delta_logs) 到文件
        saveMatrixToFile(H_delta_logs, "C:\\Users\\yfyttx\\Documents\\code\\quantum computing\\NB_LDPC\\NB-LDPC-toolbox-master\\H_delta.txt");
        
    } catch (const std::exception& e) {
        std::cerr << "\nAn error occurred: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
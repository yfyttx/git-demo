#include "utils/matrix_utils.h"
#include "gf2p/gf2p.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

using namespace Eigen;

// --- 辅助函数 ---
std::vector<std::pair<int, int>> getNonZeroPositions(const MatrixXi &mat) {
    std::vector<std::pair<int, int>> positions;
    for (int i = 0; i < mat.rows(); ++i) {
        for (int j = 0; j < mat.cols(); ++j) {
            if (mat(i, j) != 0) {
                positions.emplace_back(i, j);
            }
        }
    }
    return positions;
}

bool verifyNonbinaryOrthogonality(const MatrixXi &H_gamma, const MatrixXi &H_delta,
                                  const GF2p &gf) {

    if (H_gamma.cols() != H_delta.cols() || H_gamma.rows() != H_delta.rows()) {
        return false;
    }

    MatrixXi H_delta_T = H_delta.transpose();

    for (int i = 0; i < H_gamma.rows(); ++i) {
        for (int j = 0; j < H_delta_T.cols(); ++j) {
            int sum = 0;
            for (int k = 0; k < H_gamma.cols(); ++k) {
                int term = gf.mul(H_gamma(i, k), H_delta_T(k, j));
                sum = gf.add(sum, term);
            }
            if (sum != 0) {
                return false;
            }
        }
    }
    return true;
}

// ... print 函数省略，如果需要可以保留 ...

// --- [核心修改] saveToKN 实现 ---
// 这种格式专门用于 NB_LDPC_FB_public 仿真器 (init.c with KN_matrix defined)
void saveToKN(const Eigen::SparseMatrix<int> &mat, const std::string &filename, const GF2p &gf) {
    int M = mat.rows();
    int N = mat.cols();
    int q = 64; // GF(64)

    // 1. 统计度数
    std::vector<int> dv(N, 0);
    std::vector<int> dc(M, 0);

    // 遍历稀疏矩阵非零元
    for (int k = 0; k < mat.outerSize(); ++k) {
        for (Eigen::SparseMatrix<int>::InnerIterator it(mat, k); it; ++it) {
            if (it.value() != 0) {
                dv[it.col()]++;
                dc[it.row()]++;
            }
        }
    }

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return;
    }

    std::cout << "Writing KN-Format to " << filename << " (N=" << N << ", M=" << M << ")..."
              << std::endl;

    // 2. 写入头部 (N M GF)
    file << N << " " << M << " " << q << "\n";

    // 3. 写入度数列表
    for (int i = 0; i < N; ++i)
        file << dv[i] << (i == N - 1 ? "" : " ");
    file << "\n";
    for (int i = 0; i < M; ++i)
        file << dc[i] << (i == M - 1 ? "" : " ");
    file << "\n";

    // 4. 写入行连接与数值 (Row Connections & Values)
    // 为了按行写入，我们需要转置或者行遍历。SparseMatrix 是列优先，转置后行遍历最高效。
    Eigen::SparseMatrix<int> mat_T = mat.transpose(); // mat_T 的列就是原矩阵的行

    for (int i = 0; i < M; ++i) {
        // 对于第 i 行 (原矩阵)
        std::vector<std::pair<int, int>> row_data;
        for (Eigen::SparseMatrix<int>::InnerIterator it(mat_T, i); it; ++it) {
            if (it.value() != 0) {
                int col_idx = it.row() + 1; // 1-based index
                int poly_val = it.value();
                int log_val = gf.log(poly_val); // 转换为幂次 (0..62)
                row_data.push_back({col_idx, log_val});
            }
        }

        // 写入文件
        for (size_t k = 0; k < row_data.size(); ++k) {
            file << row_data[k].first << " " << row_data[k].second;
            if (k < row_data.size() - 1)
                file << " ";
        }
        file << "\n";
    }

    file.close();
    std::cout << "Done (KN format)." << std::endl;
}

// 占位函数，为了保持兼容性
void saveToALIST(const Eigen::SparseMatrix<int> &mat, const std::string &filename) {
    std::cout << "Warning: saveToALIST is deprecated. Use saveToKN." << std::endl;
}

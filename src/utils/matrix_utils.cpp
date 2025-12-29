#include "utils/matrix_utils.h"
#include "gf2p/gf2p.h"
#include <iostream>
#include <fstream>      // [必须添加] 用于写文件
#include <vector>
#include <algorithm>    // [必须添加] 用于 std::max
#include <Eigen/Sparse> // [必须添加]
#include <Eigen/Dense>

using namespace Eigen;

// --- 原有函数保持不变 ---

std::vector<std::pair<int, int>> getNonZeroPositions(const MatrixXi& mat) {
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

bool verifyNonbinaryOrthogonality(
    const MatrixXi& H_gamma,
    const MatrixXi& H_delta,
    const GF2p& gf) {
    
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

void printNonbinaryMatrix(
    const MatrixXi& mat,
    const std::string& name,
    const GF2p& gf,
    int P) {
    
    std::cout << "\n=== " << name << " (" << mat.rows() << "x" << mat.cols() 
              << ", GF(2^" << gf.getSize() << ")) ===" << std::endl;
    
    for (int i = 0; i < mat.rows(); ++i) {
        for (int j = 0; j < mat.cols(); ++j) {
            std::cout << gf.toString(mat(i, j)) << " ";
            if ((j + 1) % P == 0 && j < mat.cols() -1) std::cout << "| ";
        }
        std::cout << std::endl;
        if ((i + 1) % P == 0 && i < mat.rows() -1) {
            for(int k=0; k < mat.cols() + (mat.cols()/P -1)*2; ++k) std::cout << "-";
            std::cout << std::endl;
        }
    }
}

void printAsHexLog(
    const Eigen::MatrixXi& mat,
    const std::string& name,
    const GF2p& gf,
    int P) {

    std::cout << "\n=== " << name << " (Hex Log Representation) ===" << std::endl;

    for (int i = 0; i < mat.rows(); ++i) {
        for (int j = 0; j < mat.cols(); ++j) {
            int element = mat(i, j);

            if (element == 0) {
                std::cout << " "; 
            } else {
                int log_val = gf.log(element);
                std::stringstream ss;
                ss << std::hex << log_val;
                std::cout << ss.str();
            }
            
            std::cout << " ";

            if ((j + 1) % P == 0 && j < mat.cols() - 1) {
                std::cout << "| ";
            }
        }
        std::cout << std::endl;

        if ((i + 1) % P == 0 && i < mat.rows() - 1) {
            for(int k=0; k < (mat.cols() * 2) + (mat.cols() / P - 1) * 2 -1; ++k) {
                std::cout << "-";
            }
            std::cout << std::endl;
        }
    }
}

// --- [必须添加] 新增的 saveToALIST 实现 ---

void saveToALIST(const Eigen::SparseMatrix<int>& mat, const std::string& filename) {
    int M = mat.rows();
    int N = mat.cols();
    int q = 64; // GF(64)

    // 1. 统计度数
    std::vector<int> dv(N, 0);
    std::vector<int> dc(M, 0);
    int max_dv = 0;
    int max_dc = 0;

    // 遍历稀疏矩阵非零元
    for (int k = 0; k < mat.outerSize(); ++k) {
        for (Eigen::SparseMatrix<int>::InnerIterator it(mat, k); it; ++it) {
            if (it.value() != 0) {
                dv[it.col()]++;
                dc[it.row()]++;
            }
        }
    }

    for (int d : dv) max_dv = std::max(max_dv, d);
    for (int d : dc) max_dc = std::max(max_dc, d);

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return;
    }

    std::cout << "Writing ALIST to " << filename << " (N=" << N << ", M=" << M << ")..." << std::endl;

    // 2. 写入头部
    file << N << " " << M << " " << q << "\n";
    file << max_dv << " " << max_dc << "\n";

    // 3. 写入度数列表
    for (int i = 0; i < N; ++i) file << dv[i] << (i == N - 1 ? "" : " ");
    file << "\n";
    for (int i = 0; i < M; ++i) file << dc[i] << (i == M - 1 ? "" : " ");
    file << "\n";

    // 4. 写入列连接 (Column Connections: RowIndex)
    // SparseMatrix 默认是 Column Major，遍历非常快
    for (int j = 0; j < N; ++j) {
        std::vector<int> col_rows; 
        for (Eigen::SparseMatrix<int>::InnerIterator it(mat, j); it; ++it) {
             if (it.value() != 0) {
                 col_rows.push_back(it.row() + 1); // ALIST 使用 1-based 索引
             }
        }
        for (int r : col_rows) file << r << " ";
        // 补0以对齐 max_dv
        for (size_t k = col_rows.size(); k < (size_t)max_dv; ++k) file << "0 ";
        file << "\n";
    }

    // 5. 写入行连接 (Row Connections: ColIndex)
    // 使用转置矩阵来加速行遍历
    Eigen::SparseMatrix<int> mat_T = mat.transpose();
    for (int i = 0; i < M; ++i) {
        std::vector<int> row_cols;
        for (Eigen::SparseMatrix<int>::InnerIterator it(mat_T, i); it; ++it) {
            if (it.value() != 0) {
                // 转置后，it.row() 是原矩阵的 col
                row_cols.push_back(it.row() + 1); 
            }
        }
        for (int c : row_cols) file << c << " ";
        for (size_t k = row_cols.size(); k < (size_t)max_dc; ++k) file << "0 ";
        file << "\n";
    }

    // 6. 写入列值 (Column Values) - NB-LDPC 扩展
    for (int j = 0; j < N; ++j) {
        std::vector<int> col_vals;
        for (Eigen::SparseMatrix<int>::InnerIterator it(mat, j); it; ++it) {
             if (it.value() != 0) col_vals.push_back(it.value());
        }
        for (int v : col_vals) file << v << " ";
        for (size_t k = col_vals.size(); k < (size_t)max_dv; ++k) file << "0 ";
        file << "\n";
    }
    
    // 7. 写入行值 (Row Values) - NB-LDPC 扩展
    for (int i = 0; i < M; ++i) {
        std::vector<int> row_vals;
        for (Eigen::SparseMatrix<int>::InnerIterator it(mat_T, i); it; ++it) {
            if (it.value() != 0) row_vals.push_back(it.value());
        }
        for (int v : row_vals) file << v << " ";
        for (size_t k = row_vals.size(); k < (size_t)max_dc; ++k) file << "0 ";
        file << "\n";
    }

    file.close();
    std::cout << "Done." << std::endl;
}

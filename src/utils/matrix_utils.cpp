#include "utils/matrix_utils.h"
#include "gf2p/gf2p.h"
#include <iostream>
using namespace Eigen;

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
                // 对于零元素，打印一个点来表示空白
                std::cout << " "; 
            } else {
                // 对于非零元素：
                // 1. 计算其在 GF(2^p) 中的对数
                int log_val = gf.log(element);

                // 2. 将对数值转换为十六进制字符串
                std::stringstream ss;
                ss << std::hex << log_val;
                std::cout << ss.str();
            }
            
            // 为了对齐，在每个元素后加一个空格
            std::cout << " ";

            // 在 PxP 块之间添加垂直分隔符
            if ((j + 1) % P == 0 && j < mat.cols() - 1) {
                std::cout << "| ";
            }
        }
        std::cout << std::endl;

        // 在 J/2 组行之后添加水平分隔符
        if ((i + 1) % P == 0 && i < mat.rows() - 1) {
            // 计算分隔线的长度以实现对齐
            for(int k=0; k < (mat.cols() * 2) + (mat.cols() / P - 1) * 2 -1; ++k) {
                std::cout << "-";
            }
            std::cout << std::endl;
        }
    }
}
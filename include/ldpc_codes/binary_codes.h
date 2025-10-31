#ifndef BINARY_CODES_H
#define BINARY_CODES_H

#include <Eigen/Dense>
#include <utility>

// 声明构造二进制Hc和Hd矩阵的函数
std::pair<Eigen::MatrixXi, Eigen::MatrixXi> constructBinaryHMatrices(int J, int L, int P, int sigma, int tau);

#endif // BINARY_CODES_H
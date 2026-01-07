#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include <Eigen/Dense>
#include <Eigen/Sparse> // [必须添加]
#include <string>
#include <utility>
#include <vector>

class GF2p; // 前向声明

std::vector<std::pair<int, int>> getNonZeroPositions(const Eigen::MatrixXi &mat);

bool verifyNonbinaryOrthogonality(const Eigen::MatrixXi &H_gamma, const Eigen::MatrixXi &H_delta,
                                  const GF2p &gf);

void printNonbinaryMatrix(const Eigen::MatrixXi &mat, const std::string &name, const GF2p &gf,
                          int P);

void printAsHexLog(const Eigen::MatrixXi &mat, const std::string &name, const GF2p &gf, int P);

// [必须添加] saveToALIST 的声明
void saveToALIST(const Eigen::SparseMatrix<int> &mat, const std::string &filename);
void saveToKN(const Eigen::SparseMatrix<int> &mat, const std::string &filename, const GF2p &gf);

// saveToALIST 的声明 (保留以兼容旧代码，虽然后面可能不用)
#endif // MATRIX_UTILS_H

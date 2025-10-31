#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <utility>

class GF2p; // 前向声明

std::vector<std::pair<int, int>> getNonZeroPositions(const Eigen::MatrixXi& mat);

bool verifyNonbinaryOrthogonality(
    const Eigen::MatrixXi& H_gamma,
    const Eigen::MatrixXi& H_delta,
    const GF2p& gf);

void printNonbinaryMatrix(
    const Eigen::MatrixXi& mat,
    const std::string& name,
    const GF2p& gf,
    int P);

void printAsHexLog(
    const Eigen::MatrixXi& mat,
    const std::string& name,
    const GF2p& gf,
    int P);
#endif // MATRIX_UTILS_H
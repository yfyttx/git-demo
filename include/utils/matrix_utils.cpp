#include <Eigen/Sparse>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

// 直接导出为非二进制 ALIST 格式 (兼容 NB-LDPC toolbox)
// 格式: N M GF
//       max_dv max_dc
//       dv_list
//       dc_list
//       Col Connections (Pos Val Pos Val ...)
//       Row Connections (Pos Val Pos Val ...)
void saveToALIST(const Eigen::SparseMatrix<int> &mat, const std::string &filename) {
    int M = mat.rows();
    int N = mat.cols();
    int q = 64; // GF(64)

    // 1. 统计度数 (Degrees)
    std::vector<int> dv(N, 0);
    std::vector<int> dc(M, 0);
    int max_dv = 0;
    int max_dc = 0;

    // 遍历非零元素统计度数
    for (int k = 0; k < mat.outerSize(); ++k) {
        for (Eigen::SparseMatrix<int>::InnerIterator it(mat, k); it; ++it) {
            // it.row(), it.col(), it.value()
            if (it.value() != 0) { // 0 是空元素
                dv[it.col()]++;
                dc[it.row()]++;
            }
        }
    }

    for (int d : dv)
        max_dv = std::max(max_dv, d);
    for (int d : dc)
        max_dc = std::max(max_dc, d);

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return;
    }

    std::cout << "Writing ALIST to " << filename << " (N=" << N << ", M=" << M << ")..."
              << std::endl;

    // 2. 写入头部
    file << N << " " << M << " " << q << "\n";
    file << max_dv << " " << max_dc << "\n";

    // 3. 写入度数列表
    for (int i = 0; i < N; ++i)
        file << dv[i] << (i == N - 1 ? "" : " ");
    file << "\n";
    for (int i = 0; i < M; ++i)
        file << dc[i] << (i == M - 1 ? "" : " ");
    file << "\n";

    // 4. 写入列连接 (Column Connections: RowIndex Value)
    // SparseMatrix 默认是 Column Major，这步非常快
    for (int j = 0; j < N; ++j) {
        // 查找第 j 列的所有非零元素
        // 注意：ALIST 使用 1-based indexing
        std::vector<std::pair<int, int>> col_entries;

        // Eigen SparseMatrix 迭代器
        // 必须确保矩阵已压缩 (makeCompressed)，通常构造后就是压缩的
        // 为了安全，我们可以使用 coeffRef 吗？不，太慢。
        // 如果 mat 是 ColMajor (默认)，直接迭代：
        for (Eigen::SparseMatrix<int>::InnerIterator it(mat, j); it; ++it) {
            if (it.value() != 0) {
                col_entries.push_back({it.row() + 1, it.value()});
            }
        }

        for (const auto &p : col_entries) {
            file << p.first << " " << p.second << " ";
        }
        // Padding (如果需要对齐，虽然 fscanf 通常不需要，但保持格式整洁)
        // for (size_t k = col_entries.size(); k < max_dv; ++k) file << "0 0 ";
        file << "\n";
    }

    // 5. 写入行连接 (Row Connections: ColIndex Value)
    // 需要行遍历。将矩阵转置或者构建行列表。
    // 为了效率，构建临时的行主序矩阵
    Eigen::SparseMatrix<int, Eigen::RowMajor> mat_row_major = mat;

    for (int i = 0; i < M; ++i) {
        for (Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(mat_row_major, i); it;
             ++it) {
            if (it.value() != 0) {
                file << (it.col() + 1) << " " << it.value() << " ";
            }
        }
        file << "\n";
    }

    file.close();
    std::cout << "Done." << std::endl;
}

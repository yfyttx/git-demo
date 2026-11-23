/**
 * Construct an orthogonal nonbinary TIT matrix pair and its binary expansion.
 * The four outputs are consumed together by simulate_qec.
 */

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>
#include <ctime>
#include "gf2p/gf2p.h"
#include<queue>
// 引入 PEG 和 TIT 构造头文件
#include "../peg-unige/BigGirth.h"
#include "../peg-unige/Random.h"
#include "ldpc_codes/binary_codes.h"
#include "ldpc_codes/nonbinary_codes.h" // TIT 构造函数在这里
using namespace Eigen;
using namespace std;
int calculate_girth(const SparseMatrix<int>& H) {
    int rows = H.rows();
    int cols = H.cols();
    if (H.nonZeros() == 0) return 0;

    int min_girth = 999;
    int total_nodes = rows + cols;

    // 构建邻接表
    std::vector<std::vector<int>> adj(total_nodes);
    for (int k = 0; k < H.outerSize();++k) { // 兼容 ColMajor 遍历
        for (SparseMatrix<int>::InnerIterator it(H, k); it; ++it) {
            int r = it.row();
            int c = it.col();
            adj[c].push_back(r + cols);
            adj[r + cols].push_back(c);
        }
    }

    // 对每个节点进行 BFS
    for (int start_node = 0; start_node < total_nodes; ++start_node) {
        if (adj[start_node].empty()) continue;

        std::vector<int> dist(total_nodes, -1);
        std::vector<int> parent(total_nodes, -1);
        std::queue<int> q;

        dist[start_node] = 0;
        q.push(start_node);

        while (!q.empty()) {
            int u = q.front();
            q.pop();

            for (int v : adj[u]) {
                if (dist[v] == -1) {
                    dist[v] = dist[u] + 1;
                    parent[v] = u;
                    q.push(v);
                } else if (v != parent[u]) {
                    int cycle_len = dist[u] + dist[v] + 1;
                    if (cycle_len < min_girth) min_girth = cycle_len;
                    if (min_girth == 4) return 4; // 已达最小可能围长，提前退出
                }
            }
        }
    }
    return (min_girth == 999) ? 0 : min_girth;
}
// ==========================================
// 2. 辅助工具：矩阵转换与保存
// ==========================================

// [核心新增]：通过行列随机缩放，生成严格正交的非二元 TIT 矩阵对
pair<SparseMatrix<int>, SparseMatrix<int>>
generateOrthogonalNB_TIT(const MatrixXi &Hc, const MatrixXi &Hd, GF2p &gf,
                         std::mt19937 &rng) {
    return constructKasaiNonbinaryPair(Hc, Hd, gf, rng);
}

void saveToEssaiFormat(const SparseMatrix<int> &H, const std::string &filename, GF2p &gf) {
    ofstream out(filename);
    if (!out) throw std::runtime_error("Cannot open matrix output: " + filename);
    int M = H.rows(), N = H.cols();
    out << N << " " << M << " " << gf.getSize() << "\n";
    vector<int> dv(N, 0), dc(M, 0);
    for (int k = 0; k < H.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(H, k); it; ++it) {
            dv[it.col()]++;
            dc[it.row()]++;
        }
    for (int i = 0; i < N; ++i)
        out << (i ? " " : "") << dv[i];
    out << "\n";
    for (int i = 0; i < M; ++i)
        out << (i ? " " : "") << dc[i];
    out << "\n";
    vector<vector<pair<int, int>>> rows(M);
    for (int k = 0; k < H.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(H, k); it; ++it)
            rows[it.row()].push_back({it.col() + 1, gf.log(it.value())});
    for (int i = 0; i < M; ++i) {
        for (size_t j = 0; j < rows[i].size(); ++j) {
            if (j) out << " ";
            out << rows[i][j].first << " " << rows[i][j].second;
        }
        out << "\n";
    }
    out.close();
    if (!out) throw std::runtime_error("Failed to write matrix output: " + filename);
    cout << "  -> Saved " << filename << " (" << N << "x" << M << ")\n";
}

// ==========================================
// 3. 构造与计算逻辑
// ==========================================

// PEG 构造 H1
SparseMatrix<int> buildH_PEG(int M, int N, int dv, GF2p &gf) {
    cout << "  [PEG] Generating H1 (M=" << M << ", N=" << N << ")..." << endl;
    vector<int> degSeq(N, dv);
    BigGirth bg(M, N, degSeq.data(), "dummy.dat", 1, 10000, false);
    bg.loadH();
    SparseMatrix<int> H(M, N);
    vector<Triplet<int>> t;
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            if (bg.H[i][j])
                t.emplace_back(i, j, gf.exp(rand() % (gf.getSize() - 1)));
    H.setFromTriplets(t.begin(), t.end());
    return H;
}

// 计算生成矩阵 G1 (H1 的核空间)
SparseMatrix<int> computeGeneratorMatrix(const SparseMatrix<int> &H_sparse, GF2p &gf) {
    int m = H_sparse.rows(), n = H_sparse.cols();
    // 转换为稠密矩阵处理，消元过程会产生大量填充
    MatrixXi H = MatrixXi::Zero(m, n);
    for (int k = 0; k < H_sparse.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(H_sparse, k); it; ++it)
            H(it.row(), it.col()) = it.value();

    int pivot_row = 0;
    vector<int> pivot_cols;
    vector<bool> is_pivot(n, false);
    vector<bool> col_used(n, false); // 记录哪些列已被选为主元

    // --- 核心优化：基于稀疏度的动态主元选择 ---
    while (pivot_row < m) {
        int best_col = -1;
        int min_degree = n + 1;

        // 寻找剩余列中非零元素最少（度数最小）且非零的列
        for (int j = 0; j < n; ++j) {
            if (col_used[j]) continue;

            int current_degree = 0;
            int first_nonzero_row = -1;
            for (int r = pivot_row; r < m; ++r) {
                if (H(r, j) != 0) {
                    current_degree++;
                    if (first_nonzero_row == -1) first_nonzero_row = r;
                }
            }

            if (current_degree > 0 && current_degree < min_degree) {
                min_degree = current_degree;
                best_col = j;
            }
        }

        // 如果找不到任何非零元素，说明剩下的行全为0，提前结束
        if (best_col == -1) break;

        // 找到该列中第一个非零元素所在的行，并交换到当前 pivot_row
        int sel_row = -1;
        for (int r = pivot_row; r < m; ++r) {
            if (H(r, best_col) != 0) {
                sel_row = r;
                break;
            }
        }

        if (sel_row != pivot_row)
            H.row(pivot_row).swap(H.row(sel_row));

        // 归一化主元行
        int iv = gf.inv(H(pivot_row, best_col));
        for (int k = 0; k < n; ++k) {
            if (H(pivot_row, k) != 0)
                H(pivot_row, k) = gf.mul(H(pivot_row, k), iv);
        }

        // 消去其他行在该列的元素（简化阶梯形）
        for (int r = 0; r < m; ++r) {
            if (r != pivot_row && H(r, best_col) != 0) {
                int factor = H(r, best_col);
                for (int k = 0; k < n; ++k) {
                    if (H(pivot_row, k) != 0 || H(r, k) != 0) {
                        H(r, k) = gf.sub(H(r, k), gf.mul(factor, H(pivot_row, k)));
                    }
                }
            }
        }

        pivot_cols.push_back(best_col);
        is_pivot[best_col] = true;
        col_used[best_col] = true;
        pivot_row++;
    }

    // --- 构造生成矩阵 G ---
    int rank = pivot_cols.size();
    SparseMatrix<int> G(n - rank, n);
    vector<Triplet<int>> tG;
    int g_row_idx = 0;

    // 根据 H * G^T = 0 构造基向量
    for (int j = 0; j < n; ++j) {
        if (is_pivot[j]) continue;

        // 自由变量位置设为 1
        tG.emplace_back(g_row_idx, j, 1);

        // 对应主元位置填入 H 矩阵消元后的系数
        for (int i = 0; i < rank; ++i) {
            int val = H(i, j);
            if (val != 0)
                tG.emplace_back(g_row_idx, pivot_cols[i], val);
        }
        g_row_idx++;
    }

    G.setFromTriplets(tG.begin(), tG.end());
    return G;
}

// 克罗内克积
SparseMatrix<int> kroneckerNB(const SparseMatrix<int> &A, const SparseMatrix<int> &B,
                              GF2p &gf) {
    SparseMatrix<int> C(A.rows() * B.rows(), A.cols() * B.cols());
    vector<Triplet<int>> t;
    t.reserve((long long)A.nonZeros() * B.nonZeros());
    for (int k = 0; k < A.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator itA(A, k); itA; ++itA)
            for (int l = 0; l < B.outerSize(); ++l)
                for (SparseMatrix<int>::InnerIterator itB(B, l); itB; ++itB) {
                    int val = gf.mul(itA.value(), itB.value());
                    if (val)
                        t.emplace_back(itA.row() * B.rows() + itB.row(),
                                       itA.col() * B.cols() + itB.col(), val);
                }
    C.setFromTriplets(t.begin(), t.end());
    return C;
}

SparseMatrix<int> makeIdentityNB(int n) {
    SparseMatrix<int> I(n, n);
    vector<Triplet<int>> entries;
    entries.reserve(n);
    for (int i = 0; i < n; ++i) entries.emplace_back(i, i, 1);
    I.setFromTriplets(entries.begin(), entries.end());
    return I;
}

SparseMatrix<int> makeSelectorWithTrailingZero(int k) {
    // [I_k | 0]: the final inner-code coordinate is absent from this block.
    SparseMatrix<int> S(k, k + 1);
    vector<Triplet<int>> entries;
    entries.reserve(k);
    for (int i = 0; i < k; ++i) entries.emplace_back(i, i, 1);
    S.setFromTriplets(entries.begin(), entries.end());
    return S;
}

SparseMatrix<int> makeAllOnesRow(int n) {
    SparseMatrix<int> H(1, n);
    vector<Triplet<int>> entries;
    entries.reserve(n);
    for (int j = 0; j < n; ++j) entries.emplace_back(0, j, 1);
    H.setFromTriplets(entries.begin(), entries.end());
    return H;
}

SparseMatrix<int> verticalStackNB(const SparseMatrix<int> &top,
                                  const SparseMatrix<int> &bottom) {
    if (top.cols() != bottom.cols())
        throw invalid_argument("Cannot stack matrices with different column counts");
    SparseMatrix<int> result(top.rows() + bottom.rows(), top.cols());
    vector<Triplet<int>> entries;
    entries.reserve(top.nonZeros() + bottom.nonZeros());
    for (int row = 0; row < top.outerSize(); ++row)
        for (SparseMatrix<int>::InnerIterator it(top, row); it; ++it)
            entries.emplace_back(it.row(), it.col(), it.value());
    for (int row = 0; row < bottom.outerSize(); ++row)
        for (SparseMatrix<int>::InnerIterator it(bottom, row); it; ++it)
            entries.emplace_back(top.rows() + it.row(), it.col(), it.value());
    result.setFromTriplets(entries.begin(), entries.end());
    return result;
}

int countZeroColumns(const SparseMatrix<int> &matrix) {
    vector<bool> has_nonzero(matrix.cols(), false);
    for (int outer = 0; outer < matrix.outerSize(); ++outer)
        for (SparseMatrix<int>::InnerIterator it(matrix, outer); it; ++it)
            has_nonzero[it.col()] = true;
    return static_cast<int>(count(has_nonzero.begin(), has_nonzero.end(), false));
}

// ==========================================
// 4. [重点新增]：GF域正交性严谨验证函数
// ==========================================
// 验证 HX * HZ^T = 0 在 GF(64) 下是否成立
int verifyOrthogonality(const SparseMatrix<int> &HX, const SparseMatrix<int> &HZ,
                        const GF2p &gf) {
    if (HX.cols() != HZ.cols())
        return -1; // 列数不匹配，无法进行内积

    int Mx = HX.rows();
    int Mz = HZ.rows();

    // 为了极速计算稀疏矩阵行的内积，转换为行主序 (RowMajor)
    SparseMatrix<int, RowMajor> HX_r = HX;
    SparseMatrix<int, RowMajor> HZ_r = HZ;

    int error_count = 0;

    // 遍历 HX 的每一行
    for (int i = 0; i < Mx; ++i) {
        // 遍历 HZ 的每一行 (相当于 HZ^T 的列)
        for (int j = 0; j < Mz; ++j) {
            int dot_product = 0;

            // 使用双指针法快速计算稀疏向量的点乘
            SparseMatrix<int, RowMajor>::InnerIterator itX(HX_r, i);
            SparseMatrix<int, RowMajor>::InnerIterator itZ(HZ_r, j);

            while (itX && itZ) {
                if (itX.col() == itZ.col()) {
                    int term = gf.mul(itX.value(), itZ.value());
                    dot_product = gf.add(dot_product, term); // GF加法
                    ++itX;
                    ++itZ;
                } else if (itX.col() < itZ.col()) {
                    ++itX;
                } else {
                    ++itZ;
                }
            }

            if (dot_product != 0) {
                error_count++;
            }
        }
    }
    return error_count;
}
int verifyBinaryOrthogonality(const SparseMatrix<int> &HX, const SparseMatrix<int> &HZ) {
    if (HX.cols() != HZ.cols()) return -1;
    SparseMatrix<int, RowMajor> HX_r = HX;
    SparseMatrix<int, RowMajor> HZ_r = HZ;
    int error_count = 0;
    for (int i = 0; i < HX_r.rows(); ++i) {
        for (int j = 0; j < HZ_r.rows(); ++j) {
            int parity = 0;
            SparseMatrix<int, RowMajor>::InnerIterator itX(HX_r, i);
            SparseMatrix<int, RowMajor>::InnerIterator itZ(HZ_r, j);
            while (itX && itZ) {
                if (itX.col() == itZ.col()) {
                    parity ^= 1;
                    ++itX;
                    ++itZ;
                } else if (itX.col() < itZ.col()) {
                    ++itX;
                } else {
                    ++itZ;
                }
            }
            error_count += parity;
        }
    }
    return error_count;
}
// GF(2^e) multiplication matrix in the field's polynomial basis.
MatrixXi getCompanionMatrix(int val, GF2p &gf) {
    const int e = gf.getPower();
    MatrixXi A = MatrixXi::Zero(e, e);
    for (int input_bit = 0; input_bit < e; ++input_bit) {
        const int product = gf.mul(val, 1 << input_bit);
        for (int output_bit = 0; output_bit < e; ++output_bit) {
            A(output_bit, input_bit) = (product >> output_bit) & 1;
        }
    }
    return A;
}

// ==========================================
// 2. 展开函数：将 GF(q) 矩阵转换为二进制矩阵
// ==========================================
SparseMatrix<int> convertToBinary(const SparseMatrix<int> &H_nb, GF2p &gf, bool is_transposed) {
    const int e = gf.getPower();
    int M_nb = H_nb.rows();
    int N_nb = H_nb.cols();

    SparseMatrix<int> H_bin(M_nb * e, N_nb * e);
    vector<Triplet<int>> triplets;

    for (int k = 0; k < H_nb.outerSize(); ++k) {
        for (SparseMatrix<int>::InnerIterator it(H_nb, k); it; ++it) {
            MatrixXi comp = getCompanionMatrix(it.value(), gf);
            // 论文提到 HZ 构造时可能需要转置处理 [cite: 248]
            if (is_transposed) comp = comp.transpose().eval();

            for (int i = 0; i < e; ++i) {
                for (int j = 0; j < e; ++j) {
                    if (comp(i, j) == 1) {
                        triplets.emplace_back(it.row() * e + i, it.col() * e + j, 1);
                    }
                }
            }
        }
    }
    H_bin.setFromTriplets(triplets.begin(), triplets.end());
    return H_bin;
}

// ==========================================
// 3. 专门保存物理二进制矩阵的函数
// ==========================================
void saveToBinaryFormat(const SparseMatrix<int> &H, const std::string &filename) {
    ofstream out(filename);
    if (!out) throw std::runtime_error("Cannot open matrix output: " + filename);
    int M = H.rows(), N = H.cols();
    // 物理二进制矩阵 Q=2
    out << N << " " << M << " 2\n";
    vector<int> dv(N, 0), dc(M, 0);
    for (int k = 0; k < H.outerSize(); ++k) {
        for (SparseMatrix<int>::InnerIterator it(H, k); it; ++it) {
            ++dv[it.col()];
            ++dc[it.row()];
        }
    }
    for (int i = 0; i < N; ++i) out << (i ? " " : "") << dv[i];
    out << '\n';
    for (int i = 0; i < M; ++i) out << (i ? " " : "") << dc[i];
    out << '\n';

    // 转换为 RowMajor 方便保存
    SparseMatrix<int, RowMajor> Hr = H;
    for (int i = 0; i < M; ++i) {
        bool first = true;
        for (SparseMatrix<int, RowMajor>::InnerIterator it(Hr, i); it; ++it) {
            // 保存列索引(从1开始)和值(1)
            if (!first) out << " ";
            out << it.col() + 1 << " 1";
            first = false;
        }
        out << "\n";
    }
    out.close();
    if (!out) throw std::runtime_error("Failed to write matrix output: " + filename);
    cout << "  -> Saved Binary Physical Matrix: " << filename << endl;
}
// ==========================================
// 5. Main
// ==========================================
int main(int argc, char *argv[]) {
    if (argc > 2) {
        cerr << "Usage: construct [seed]\n";
        return 1;
    }
    const unsigned int seed = argc == 2 ? static_cast<unsigned int>(stoul(argv[1]))
                                        : std::random_device{}();
    std::mt19937 rng(seed);
    srand(seed);
    cout << "seed=" << seed << '\n';
    GF2p gf(6, {1,0,0,0,0,1,1});
    cout << "=== GF(64) TIT Matrix Construction ===\n";

    // --------------------------------------------------------
    // Step 2: Generate the orthogonal TIT parity-check pair.
    // --------------------------------------------------------
    cout << "  [TIT] Generating the orthogonal HX/HZ pair...\n";
    const int J = 2;
    int L = 6, P = 31;
    int sigma = 5, tau = 2;

    // 1. 获取基础二进制 TIT 矩阵 (骨架)
    pair<MatrixXi, MatrixXi> tit_bin = constructBinaryHMatrices(J, L, P, sigma, tau);

    auto [hx_tit, hz_tit] = generateOrthogonalNB_TIT(tit_bin.first, tit_bin.second, gf, rng);

    cout << "  [Check] TIT HX non-zeros: " << hx_tit.nonZeros() << endl;
    cout << "  [Check] TIT HZ non-zeros: " << hz_tit.nonZeros() << endl;
    const int tit_orthogonality_errors = verifyOrthogonality(hx_tit, hz_tit, gf);
    cout << "  [Check] TIT Orthogonality (HX*HZ^T) error count: "
         << tit_orthogonality_errors << endl;
    if (tit_orthogonality_errors != 0) {
        cerr << "Refusing to save a non-orthogonal TIT matrix pair.\n";
        return 1;
    }

    // The current decoder consumes this directly generated orthogonal TIT pair.
    // Keep the logical matrices paired with their binary expansions below.
    cout << "  [Result] TIT HX: " << hx_tit.rows() << "x" << hx_tit.cols() << endl;
    cout << "  [Result] TIT HZ: " << hz_tit.rows() << "x" << hz_tit.cols() << endl;

    cout << "  [TIT HX Girth]: " << calculate_girth(hx_tit) << endl;
    cout << "  [TIT HZ Girth]: " << calculate_girth(hz_tit) << endl;
    // --------------------------------------------------------
    // Step 5: 核心验证 - GF(64) 正交性检查
    // --------------------------------------------------------
    cout << "\n=== Verifying Quantum Orthogonality ===" << endl;
    cout << "  Checking if HX * HZ^T = 0 in GF(64)..." << endl;

    int non_orthogonal_elements = verifyOrthogonality(hx_tit, hz_tit, gf);

    if (non_orthogonal_elements == 0) {
        cout << "  [SUCCESS] -> HX and HZ are perfectly orthogonal in GF(64)!" << endl;
    } else {
        cout << "  [FAILED]  -> Found " << non_orthogonal_elements
             << " non-zero elements in HX * HZ^T." << endl;
        cerr << "Refusing to save matrices that violate the quantum orthogonality check.\n";
        return 1;
    }
    cout << "=======================================\n" << endl;

    // --------------------------------------------------------
    // Step 6: 保存
    // --------------------------------------------------------

    // 在 main 函数末尾，保存之前调用：
    cout << "  [Convert] Expanding to Binary Physical Matrices (Companion Matrix)..." << endl;

    // H_Gamma 展开为 HX, H_Delta 展开为 HZ [cite: 125, 126]
    SparseMatrix<int> hx_tit_physical = convertToBinary(hx_tit, gf, false);
    SparseMatrix<int> hz_tit_physical = convertToBinary(hz_tit, gf, true); // Z侧通常取转置

    cout << "  [Result] Physical TIT HX size: " << hx_tit_physical.rows() << "x" << hx_tit_physical.cols() << endl;
    cout << "  [Result] Physical TIT HZ size: " << hz_tit_physical.rows() << "x" << hz_tit_physical.cols() << endl;
    const int binary_orthogonality_errors = verifyBinaryOrthogonality(hx_tit_physical, hz_tit_physical);
    cout << "  [Check] Physical binary orthogonality error count: "
         << binary_orthogonality_errors << endl;
    if (binary_orthogonality_errors != 0) {
        cerr << "Refusing to save non-commuting physical stabilizers.\n";
        return 1;
    }

    // 保存为二进制格式 (Q 设为 2)
    std::filesystem::create_directories("build");
    saveToBinaryFormat(hx_tit_physical, "build/HX_Physical_Binary.txt");
    saveToBinaryFormat(hz_tit_physical, "build/HZ_Physical_Binary.txt");
    saveToEssaiFormat(hx_tit, "build/HX_Safe_GF64.txt", gf);
    saveToEssaiFormat(hz_tit, "build/HZ_Safe_GF64.txt", gf);

    // --------------------------------------------------------
    // Optional comparison construction following the video/Eq. (23):
    // Original cascade assignment requested for this experiment:
    // H_X = hz_tit_outer_for_hx kron inner_generator, and
    // H_Z = [hx_tit_outer_for_hz kron selector; I_n kron inner_parity_check].
    // Keep these outputs separate
    // from the direct TIT baseline above.
    // --------------------------------------------------------
    cout << "\n=== Building cascade comparison matrices ===" << endl;
    const int inner_length = 6;
    const int inner_dimension = inner_length - 1;
    SparseMatrix<int> hz_tit_outer_for_hx = hz_tit;
    SparseMatrix<int> hx_tit_outer_for_hz = hx_tit;
    SparseMatrix<int> inner_parity_check = makeAllOnesRow(inner_length);
    SparseMatrix<int> inner_generator = computeGeneratorMatrix(inner_parity_check, gf);
    SparseMatrix<int> selector = makeSelectorWithTrailingZero(inner_dimension);
    SparseMatrix<int> hz_cascade_upper = kroneckerNB(hx_tit_outer_for_hz, selector, gf);
    SparseMatrix<int> hz_cascade_lower = kroneckerNB(
        makeIdentityNB(hx_tit_outer_for_hz.cols()), inner_parity_check, gf);
    SparseMatrix<int> hx_cascade = kroneckerNB(hz_tit_outer_for_hx, inner_generator, gf);
    SparseMatrix<int> hz_cascade = verticalStackNB(hz_cascade_upper, hz_cascade_lower);
    SparseMatrix<int> hz_video_upper = kroneckerNB(
        hx_tit_outer_for_hz, makeIdentityNB(inner_length), gf);
    SparseMatrix<int> hz_video = verticalStackNB(hz_video_upper, hz_cascade_lower);

    cout << "  [Cascade] inner parity check: " << inner_parity_check.rows() << "x"
         << inner_parity_check.cols() << ", inner generator: "
         << inner_generator.rows() << "x" << inner_generator.cols() << endl;
    cout << "  [Cascade] HX: " << hx_cascade.rows() << "x" << hx_cascade.cols()
         << ", HZ: " << hz_cascade.rows() << "x" << hz_cascade.cols() << endl;
    cout << "  [Theorem 23] HZ upper-block zero columns: "
         << countZeroColumns(hz_cascade_upper)
         << ", complete HZ zero columns: " << countZeroColumns(hz_cascade) << endl;
    cout << "  [Video identity-block variant] HZ: " << hz_video.rows() << "x"
         << hz_video.cols() << ", GF orthogonality errors: "
         << verifyOrthogonality(hx_cascade, hz_video, gf) << endl;
    if (inner_generator.rows() != inner_dimension || inner_generator.cols() != inner_length ||
        verifyOrthogonality(inner_parity_check, inner_generator, gf) != 0) {
        cerr << "Inner code generator validation failed.\n";
        return 1;
    }
    const int cascade_orthogonality_errors = verifyOrthogonality(hx_cascade, hz_cascade, gf);
    cout << "  [Check] Cascade GF orthogonality errors: "
         << cascade_orthogonality_errors << endl;
    if (cascade_orthogonality_errors != 0) {
        cerr << "Refusing to save non-orthogonal cascade comparison matrices.\n";
        return 1;
    }

    SparseMatrix<int> hx_cascade_physical = convertToBinary(hx_cascade, gf, false);
    SparseMatrix<int> hz_cascade_physical = convertToBinary(hz_cascade, gf, true);
    const int cascade_binary_errors =
        verifyBinaryOrthogonality(hx_cascade_physical, hz_cascade_physical);
    cout << "  [Check] Cascade binary orthogonality errors: " << cascade_binary_errors << endl;
    if (cascade_binary_errors != 0) {
        cerr << "Refusing to save non-commuting binary cascade matrices.\n";
        return 1;
    }
    saveToBinaryFormat(hx_cascade_physical, "build/HX_Cascade_Physical_Binary.txt");
    saveToBinaryFormat(hz_cascade_physical, "build/HZ_Cascade_Physical_Binary.txt");
    saveToEssaiFormat(hz_tit_outer_for_hx, "build/HZ_TIT_Outer_For_HX_GF64.txt", gf);
    saveToEssaiFormat(hx_tit_outer_for_hz, "build/HX_TIT_Outer_For_HZ_GF64.txt", gf);
    saveToEssaiFormat(hx_cascade, "build/HX_Cascade_GF64.txt", gf);
    saveToEssaiFormat(hz_cascade, "build/HZ_Cascade_GF64.txt", gf);

    SparseMatrix<int> hz_video_physical = convertToBinary(hz_video, gf, true);
    const int video_binary_errors =
        verifyBinaryOrthogonality(hx_cascade_physical, hz_video_physical);
    cout << "  [Check] Video identity-block binary orthogonality errors: "
         << video_binary_errors << endl;
    if (video_binary_errors != 0) {
        cerr << "Refusing to save non-commuting video-variant matrices.\n";
        return 1;
    }
    saveToBinaryFormat(hx_cascade_physical, "build/HX_Video_Physical_Binary.txt");
    saveToBinaryFormat(hz_video_physical, "build/HZ_Video_Physical_Binary.txt");
    saveToEssaiFormat(hx_cascade, "build/HX_Video_GF64.txt", gf);
    saveToEssaiFormat(hz_video, "build/HZ_Video_GF64.txt", gf);
    std::cout << "Matrix generation complete." << std::endl;
    return 0;
}

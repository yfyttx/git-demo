/**
 * main.cpp
 * * QLDPC Constructor: Hybrid (PEG + TIT QC-LDPC)
 * *
 * * Part 1: H1 (PEG) -> G1 (Kernel)
 * * Part 2: H2, H3 (TIT QC-LDPC Pair, where H2*H3^T = 0)
 * *
 * * HX = [ H1 (x) In2 ]
 * * [ In1 (x) H2 ]
 * *
 * * HZ = [ G1 (x) H3 ]
 */

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
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
generateOrthogonalNB_TIT(const MatrixXi &Hc, const MatrixXi &Hd, GF2p &gf) {
    int m2 = Hc.rows(), n2 = Hc.cols();
    int m3 = Hd.rows(), n3 = Hd.cols();

    // --- Step 1: 求解 H2 (H_gamma) ---
    // 提取 Hc 的非零位置用于建立线性方程组
    vector<pair<int, int>> Hc_nonzero_pos;
    for (int i = 0; i < m2; ++i)
        for (int j = 0; j < n2; ++j)
            if (Hc(i, j) == 1) Hc_nonzero_pos.emplace_back(i, j);

    // 求解对数域线性方程组 (mod q-1)
    // 注意：这里需要你之前写的 buildLinearSystem 和 gaussElimination
    MatrixXi system = buildLinearSystem(Hc_nonzero_pos, Hc, Hd, gf.getSize() - 1);
    vector<int> log_gamma = gaussElimination(system, gf.getSize() - 1);

    // 构造非二进制 H2 并转为稠密矩阵方便后续求解 H3
    SparseMatrix<int> H2(m2, n2);
    MatrixXi H2_dense = MatrixXi::Zero(m2, n2);
    vector<Triplet<int>> t2;
    for (size_t i = 0; i < Hc_nonzero_pos.size(); ++i) {
        int val = gf.exp(log_gamma[i]);
        int r = Hc_nonzero_pos[i].first;
        int c = Hc_nonzero_pos[i].second;
        t2.emplace_back(r, c, val);
        H2_dense(r, c) = val;
    }
    H2.setFromTriplets(t2.begin(), t2.end());

    // --- Step 2: 求解 H3 (H_delta) ---
    // 这一步确保每一行都与 H2 正交，且非零！
    SparseMatrix<int> H3(m3, n3);
    vector<Triplet<int>> t3;
    for (int i = 0; i < m3; ++i) {
        // 调用你新写的 solveDeltaRow
        vector<int> delta_row = solveDeltaRow(i, Hc, Hd, H2_dense, gf);
        for (int j = 0; j < n3; ++j) {
            if (delta_row[j] != 0) {
                t3.emplace_back(i, j, delta_row[j]);
            }
        }
    }
    H3.setFromTriplets(t3.begin(), t3.end());

    // 检查 H3 是否还是空的
    if (t3.empty()) {
        std::cerr << "CRITICAL ERROR: H3 has no non-zero elements! Check solveDeltaRow." << std::endl;
    }

    return {H2, H3};
}

void saveToEssaiFormat(const SparseMatrix<int> &H, const std::string &filename, GF2p &gf) {
    ofstream out(filename);
    if (!out)
        return;
    int M = H.rows(), N = H.cols();
    out << N << " " << M << " " << gf.getSize() << "\n";
    vector<int> dv(N, 0), dc(M, 0);
    for (int k = 0; k < H.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(H, k); it; ++it) {
            dv[it.col()]++;
            dc[it.row()]++;
        }
    for (int i = 0; i < N; ++i)
        out << dv[i] << " ";
    out << "\n";
    for (int i = 0; i < M; ++i)
        out << dc[i] << " ";
    out << "\n";
    vector<vector<pair<int, int>>> rows(M);
    for (int k = 0; k < H.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(H, k); it; ++it)
            rows[it.row()].push_back({it.col() + 1, gf.log(it.value())});
    for (int i = 0; i < M; ++i) {
        for (auto &p : rows[i])
            out << p.first << " " << p.second << "    ";
        out << "\n";
    }
    out.close();
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
// 获取 GF(64) 元素的友矩阵映射 (e=6)

    // 展开后的尺寸为 (M*e) x (N*e) 
// ==========================================
// 1. 友矩阵映射：GF(64) -> 6x6 Binary Matrix
// ==========================================
MatrixXi getCompanionMatrix(int val, GF2p &gf) {
    // 修复：手动指定 e 为 6 (针对 GF64)
    int e = 6; 
    MatrixXi A = MatrixXi::Zero(e, e);
    
    if (val == 0) return A; // 0 映射为全零矩阵 [cite: 126, 243]

    // 构造基础友矩阵 A(alpha)
    // 根据论文公式 (7)，基于本原多项式构造 [cite: 239, 240]
    // GF64 常用多项式: 1 + x + x^6
    MatrixXi A_alpha = MatrixXi::Zero(e, e);
    for (int j = 0; j < e - 1; ++j) {
        A_alpha(j + 1, j) = 1; // 次对角线设为 1
    }
    // 最后一列对应多项式系数 [cite: 239]
    A_alpha(0, 5) = 1; // a0 = 1
    A_alpha(1, 5) = 1; // a1 = 1
    A_alpha(2, 5) = 0;
    A_alpha(3, 5) = 0;
    A_alpha(4, 5) = 0;
    A_alpha(5, 5) = 0; // a5 = 0

    // 计算 A(alpha^l)，其中 l = log(val) [cite: 245]
    int l = gf.log(val);
    A = MatrixXi::Identity(e, e);
    for (int i = 0; i < l; ++i) {
        // 矩阵乘法并在 GF(2) 下取模
        A = (A * A_alpha).unaryExpr([](int x) { return x % 2; });
    }
    
    return A;
}

// ==========================================
// 2. 展开函数：将 GF(q) 矩阵转换为二进制矩阵 
// ==========================================
SparseMatrix<int> convertToBinary(const SparseMatrix<int> &H_nb, GF2p &gf, bool is_transposed) {
    int e = 6; // 硬编码针对 GF64
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
    if (!out) return;
    int M = H.rows(), N = H.cols();
    // 物理二进制矩阵 Q=2 
    out << N << " " << M << " 2\n";
    // 占位符：dv/dc 信息
    out << "\n\n"; 
    
    // 转换为 RowMajor 方便保存
    SparseMatrix<int, RowMajor> Hr = H;
    for (int i = 0; i < M; ++i) {
        for (SparseMatrix<int, RowMajor>::InnerIterator it(Hr, i); it; ++it) {
            // 保存列索引(从1开始)和值(1)
            out << it.col() + 1 << " 1    "; 
        }
        out << "\n";
    }
    out.close();
    cout << "  -> Saved Binary Physical Matrix: " << filename << endl;
}
// ==========================================
// 5. Main
// ==========================================
int main() {
    srand(time(0));
    GF2p gf(6, {1,0,0,0,0,1,1}); 
    cout << "=== QLDPC: Hybrid (PEG + TIT) Construction ===\n";

    // --------------------------------------------------------
    // Step 1: 生成 H1 (PEG) 并计算 G1
    // --------------------------------------------------------
    SparseMatrix<int> H1 = buildH_PEG(10, 20, 2, gf);
    cout << "  [Kernel] Computing G1 such that H1 * G1^T = 0 in GF(64)..." << endl;
    SparseMatrix<int> G1 = computeGeneratorMatrix(H1, gf);

    // --------------------------------------------------------
    // Step 2: 严格按照 Kasai 论文逻辑生成正交的 H2 和 H3
    // --------------------------------------------------------
    cout << "  [TIT] Generating H2 and H3 with Linear Constraint Solving...\n";
    const int J = 2;
    int L = 6, P = 31;
    int sigma = 5, tau = 2;
    
    // 1. 获取基础二进制 TIT 矩阵 (骨架)
    pair<MatrixXi, MatrixXi> tit_bin = constructBinaryHMatrices(J, L, P, sigma, tau);

    auto [H2, H3] = generateOrthogonalNB_TIT(tit_bin.first, tit_bin.second, gf);

    cout << "  [Check] H2 non-zeros: " << H2.nonZeros() << endl;
    cout << "  [Check] H3 non-zeros: " << H3.nonZeros() << endl;
    // 预验证：此时 error count 必须为 0
    cout << "  [Check] NB-TIT Orthogonality (H2*H3^T) error count: " 
         << verifyOrthogonality(H2, H3, gf) << endl;

    // --------------------------------------------------------
    // Step 3: 准备单位阵
    // --------------------------------------------------------
    int n1 = H1.cols();
    int n2 = H2.cols();

    SparseMatrix<int> In2(n2, n2);
    for (int i = 0; i < n2; ++i)
        In2.insert(i, i) = 1;
    SparseMatrix<int> In1(n1, n1);
    for (int i = 0; i < n1; ++i)
        In1.insert(i, i) = 1;

    // --------------------------------------------------------
    // Step 4: 构造 HX 和 HZ
    // --------------------------------------------------------
    cout << "  [Build] HX (Stacked)..." << endl;
    SparseMatrix<int> HX_top = kroneckerNB(H1, In2, gf);
    SparseMatrix<int> HX_bot = kroneckerNB(In1, H2, gf);

    SparseMatrix<int> HX(HX_top.rows() + HX_bot.rows(), HX_top.cols());
    vector<Triplet<int>> tX;
    for (int k = 0; k < HX_top.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(HX_top, k); it; ++it)
            tX.emplace_back(it.row(), it.col(), it.value());

    int row_off = HX_top.rows();
    for (int k = 0; k < HX_bot.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(HX_bot, k); it; ++it)
            tX.emplace_back(it.row() + row_off, it.col(), it.value());
    HX.setFromTriplets(tX.begin(), tX.end());


    
    cout << "  [Build] HZ (G1 x H3)..." << endl;
    SparseMatrix<int> HZ = kroneckerNB(G1, H3, gf);


    HX = H2;
    HZ = H3;
    cout << "  [Result] HX: " << HX.rows() << "x" << HX.cols() << endl;
    cout << "  [Result] HZ: " << HZ.rows() << "x" << HZ.cols() << endl;

    cout << "  [H1 Girth]: " << calculate_girth(H1) << endl;
    cout << "  [G1 Girth]: " << calculate_girth(G1) << endl;
    cout << "  [H2 Girth]: " << calculate_girth(H2) << endl;
    cout << "  [H3 Girth]: " << calculate_girth(H3) << endl;
    cout << "  [Hx Girth]: " << calculate_girth(HX) << endl;
    cout << "  [Hz Girth]: " << calculate_girth(HZ) << endl;
    // --------------------------------------------------------
    // Step 5: 核心验证 - GF(64) 正交性检查
    // --------------------------------------------------------
    cout << "\n=== Verifying Quantum Orthogonality ===" << endl;
    cout << "  Checking if HX * HZ^T = 0 in GF(64)..." << endl;
    
    int non_orthogonal_elements = verifyOrthogonality(HX, HZ, gf);

    if (non_orthogonal_elements == 0) {
        cout << "  [SUCCESS] -> HX and HZ are perfectly orthogonal in GF(64)!" << endl;
    } else {
        cout << "  [FAILED]  -> Found " << non_orthogonal_elements
             << " non-zero elements in HX * HZ^T." << endl;
    }
    cout << "=======================================\n" << endl;

    // --------------------------------------------------------
    // Step 6: 保存
    // --------------------------------------------------------

    // 在 main 函数末尾，保存之前调用：
    cout << "  [Convert] Expanding to Binary Physical Matrices (Companion Matrix)..." << endl;

    // H_Gamma 展开为 HX, H_Delta 展开为 HZ [cite: 125, 126]
    SparseMatrix<int> HX_physical = convertToBinary(HX, gf, false); 
    SparseMatrix<int> HZ_physical = convertToBinary(HZ, gf, true); // Z侧通常取转置 

    cout << "  [Result] Physical HX size: " << HX_physical.rows() << "x" << HX_physical.cols() << endl;
    cout << "  [Result] Physical HZ size: " << HZ_physical.rows() << "x" << HZ_physical.cols() << endl;

    // 保存为二进制格式 (Q 设为 2)
    saveToBinaryFormat(HX_physical, "HX_Physical_Binary.txt");
    saveToBinaryFormat(HZ_physical, "HZ_Physical_Binary.txt");
    saveToEssaiFormat(HX, "HX_Safe_GF64.txt", gf);
    saveToEssaiFormat(HZ, "HZ_Safe_GF64.txt", gf);
    std::cout << "CHECK: VERSION 2.0 - TESTING PARAMETERS" << std::endl;
    return 0;
}

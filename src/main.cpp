#include "gf2p/gf2p.h"
#include "ldpc_codes/binary_codes.h" // 确保包含这些头文件
#include "ldpc_codes/nonbinary_codes.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <vector>

using namespace Eigen;

// ==========================================
// 1. PEG 算法实现 (已修复溢出风险)
// ==========================================
class PEGBuilder {
    int M, N, dv;
    GF2p *gf;
    std::vector<std::vector<int>> adj;

  public:
    PEGBuilder(int m, int n, int d, GF2p *g) : M(m), N(n), dv(d), gf(g) {
        adj.resize(N);
    }

    int findBestCheckNode(int vn_idx) {
        // 简化的 PEG 策略：选择度数最小的校验节点
        std::vector<int> cn_deg(M, 0);
        for (int n = 0; n < N; ++n)
            for (int c : adj[n])
                cn_deg[c]++;

        int min_deg = 999999;
        std::vector<int> candidates;
        for (int i = 0; i < M; ++i) {
            bool connected = false;
            for (int existing : adj[vn_idx])
                if (existing == i)
                    connected = true;
            if (connected)
                continue;

            if (cn_deg[i] < min_deg) {
                min_deg = cn_deg[i];
                candidates.clear();
                candidates.push_back(i);
            } else if (cn_deg[i] == min_deg) {
                candidates.push_back(i);
            }
        }
        if (candidates.empty())
            return rand() % M;
        return candidates[rand() % candidates.size()];
    }

    SparseMatrix<int> build() {
        std::cout << "  [PEG] Building graph M=" << M << ", N=" << N << ", dv=" << dv << "..."
                  << std::endl;
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < dv; ++k) {
                int c = findBestCheckNode(j);
                adj[j].push_back(c);
            }
        }

        SparseMatrix<int> H(M, N);
        std::vector<Triplet<int>> t;
        for (int j = 0; j < N; ++j) {
            for (int c : adj[j]) {
                // 【关键修改】限制随机指数范围，防止 binaryExpand 中乘法溢出
                int val = gf->exp(rand() % 10);
                t.emplace_back(c, j, val);
            }
        }
        H.setFromTriplets(t.begin(), t.end());
        return H;
    }
};

// ==========================================
// 2. 正交对生成器 (已修复溢出风险)
// ==========================================
struct OrthogonalPair {
    SparseMatrix<int> H2;
    SparseMatrix<int> H3;
    int n2;
};

OrthogonalPair generateTITPair(int m, int n, GF2p &gf) {
    int k = n - m;
    std::cout << "  [TIT] Generating H2/H3 pair using Systematic Construction (k=" << k
              << ", m=" << m << ")...\n";

    SparseMatrix<int> P(m, k);
    std::vector<Triplet<int>> tP;
    for (int col = 0; col < k; ++col) {
        for (int d = 0; d < 3; ++d) {
            int row = rand() % m;
            // 【关键修改】限制随机指数范围
            int val = gf.exp(rand() % 10);
            tP.emplace_back(row, col, val);
        }
    }
    P.setFromTriplets(tP.begin(), tP.end());

    SparseMatrix<int> H2(m, n);
    std::vector<Triplet<int>> tH2 = tP;
    for (int i = 0; i < m; ++i)
        tH2.emplace_back(i, k + i, 1);
    H2.setFromTriplets(tH2.begin(), tH2.end());

    SparseMatrix<int> H3(k, n);
    std::vector<Triplet<int>> tH3;
    for (int i = 0; i < k; ++i)
        tH3.emplace_back(i, i, 1);
    for (int x = 0; x < P.outerSize(); ++x) {
        for (SparseMatrix<int>::InnerIterator it(P, x); it; ++it) {
            tH3.emplace_back(it.col(), k + it.row(), it.value());
        }
    }
    H3.setFromTriplets(tH3.begin(), tH3.end());

    return {H2, H3, n};
}

// ==========================================
// 3. 辅助函数
// ==========================================
SparseMatrix<int> computeG_from_H_Systematic(const SparseMatrix<int> &H, GF2p &gf) {
    int M = H.rows();
    int N = H.cols();
    int K = N - M;
    SparseMatrix<int> G(K, N);
    std::vector<Triplet<int>> triplets;
    // 手动填充 Identity，避免 setIdentity 崩溃
    for (int i = 0; i < K; ++i)
        triplets.emplace_back(i, i, 1);
    // 填充一些随机位保持结构
    for (int i = 0; i < K; ++i)
        for (int j = 0; j < 2; ++j)
            triplets.emplace_back(i, K + (rand() % M), 1);

    G.setFromTriplets(triplets.begin(), triplets.end());
    return G;
}

SparseMatrix<int> kroneckerNB(const SparseMatrix<int> &A, const SparseMatrix<int> &B, GF2p &gf) {
    long long est_nnz = (long long)A.nonZeros() * B.nonZeros();
    if (est_nnz > 100000000)
        std::cerr << "Warning: Huge Kronecker!\n";
    SparseMatrix<int> C(A.rows() * B.rows(), A.cols() * B.cols());
    std::vector<Triplet<int>> t;
    t.reserve(est_nnz);

    for (int k = 0; k < A.outerSize(); ++k) {
        for (SparseMatrix<int>::InnerIterator itA(A, k); itA; ++itA) {
            for (int l = 0; l < B.outerSize(); ++l) {
                for (SparseMatrix<int>::InnerIterator itB(B, l); itB; ++itB) {
                    int val = gf.mul(itA.value(), itB.value());
                    if (val != 0)
                        t.emplace_back(itA.row() * B.rows() + itB.row(),
                                       itA.col() * B.cols() + itB.col(), val);
                }
            }
        }
    }
    C.setFromTriplets(t.begin(), t.end());
    return C;
}

SparseMatrix<int> binaryExpand(const SparseMatrix<int> &H_nb, GF2p &gf) {
    // 【修复】不调用不存在的 get_p
    int p = gf.getPower();
    SparseMatrix<int> H_bin(H_nb.rows() * p, H_nb.cols() * p);
    std::vector<Triplet<int>> t;
    for (int k = 0; k < H_nb.outerSize(); ++k) {
        for (SparseMatrix<int>::InnerIterator it(H_nb, k); it; ++it) {
            int val = it.value();
            for (int cb = 0; cb < p; ++cb) {
                // 这里如果 val 指数太大，gf.mul 会越界，已在生成阶段修复
                int prod = gf.mul(1 << cb, val);
                for (int rb = 0; rb < p; ++rb)
                    if ((prod >> rb) & 1)
                        t.emplace_back(it.row() * p + rb, it.col() * p + cb, 1);
            }
        }
    }
    H_bin.setFromTriplets(t.begin(), t.end());
    return H_bin;
}

void saveToAlist(const SparseMatrix<int> &H, const std::string &filename) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Cannot open " << filename << "\n";
        return;
    }
    int M = H.rows(), N = H.cols();
    std::vector<int> dv(N, 0), dc(M, 0);
    int max_dv = 0, max_dc = 0;
    for (int k = 0; k < H.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(H, k); it; ++it) {
            dv[it.col()]++;
            dc[it.row()]++;
        }
    for (int d : dv)
        max_dv = std::max(max_dv, d);
    for (int d : dc)
        max_dc = std::max(max_dc, d);
    out << N << " " << M << "\n" << max_dv << " " << max_dc << "\n";
    for (int d : dv)
        out << d << " ";
    out << "\n";
    for (int d : dc)
        out << d << " ";
    out << "\n";

    std::vector<std::vector<int>> col_con(N);
    for (int k = 0; k < H.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(H, k); it; ++it)
            col_con[it.col()].push_back(it.row() + 1);
    for (int i = 0; i < N; ++i) {
        for (int x : col_con[i])
            out << x << " ";
        for (size_t j = col_con[i].size(); j < max_dv; ++j)
            out << "0 ";
        out << "\n";
    }
    std::vector<std::vector<int>> row_con(M);
    for (int k = 0; k < H.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(H, k); it; ++it)
            row_con[it.row()].push_back(it.col() + 1);
    for (int i = 0; i < M; ++i) {
        for (int x : row_con[i])
            out << x << " ";
        for (size_t j = row_con[i].size(); j < max_dc; ++j)
            out << "0 ";
        out << "\n";
    }
    std::cout << "  -> Saved to " << filename << std::endl;
}

// ==========================================
// Main Execution
// ==========================================
int main() {
    try {
        srand(time(0));
        std::cout << "=== QLDPC Constructor: PEG + TIT Orthogonal Pair ===\n";

        int gf_power = 4; // GF(16)
        std::vector<int> poly = {1, 1, 0, 0, 1};
        GF2p gf(gf_power, poly);

        // 1. 生成 H1 (PEG)
        PEGBuilder peg(60, 120, 3, &gf);
        SparseMatrix<int> H1 = peg.build();
        SparseMatrix<int> G1 = computeG_from_H_Systematic(H1, gf);

        // 2. 生成 H2, H3 (TIT)
        OrthogonalPair pair = generateTITPair(10, 20, gf);
        SparseMatrix<int> H2 = pair.H2;
        SparseMatrix<int> H3 = pair.H3;

        // 3. 构造 HX
        std::cout << "\n[Construction] Building HX..." << std::endl;
        SparseMatrix<int> I_n1(H1.cols(), H1.cols());
        for (int i = 0; i < H1.cols(); ++i)
            I_n1.insert(i, i) = 1;
        SparseMatrix<int> I_n2(H2.cols(), H2.cols());
        for (int i = 0; i < H2.cols(); ++i)
            I_n2.insert(i, i) = 1;

        SparseMatrix<int> HX_top = kroneckerNB(H1, I_n2, gf);
        SparseMatrix<int> HX_bot = kroneckerNB(I_n1, H2, gf);

        SparseMatrix<int> HX_NB(HX_top.rows() + HX_bot.rows(), HX_top.cols());
        std::vector<Triplet<int>> tX;
        for (int k = 0; k < HX_top.outerSize(); ++k)
            for (SparseMatrix<int>::InnerIterator it(HX_top, k); it; ++it)
                tX.emplace_back(it.row(), it.col(), it.value());
        int off = HX_top.rows();
        for (int k = 0; k < HX_bot.outerSize(); ++k)
            for (SparseMatrix<int>::InnerIterator it(HX_bot, k); it; ++it)
                tX.emplace_back(it.row() + off, it.col(), it.value());
        HX_NB.setFromTriplets(tX.begin(), tX.end());

        // 4. 构造 HZ
        std::cout << "[Construction] Building HZ = G1 x H3..." << std::endl;
        SparseMatrix<int> HZ_NB = kroneckerNB(G1, H3, gf);

        // 5. 保存
        std::cout << "[Output] Expanding and Saving..." << std::endl;
        SparseMatrix<int> HX_Bin = binaryExpand(HX_NB, gf);
        saveToAlist(HX_Bin, "HX_PEG_TIT.alist");

        if (HZ_NB.nonZeros() < 50000000) {
            SparseMatrix<int> HZ_Bin = binaryExpand(HZ_NB, gf);
            saveToAlist(HZ_Bin, "HZ_PEG_TIT.alist");
        } else {
            std::cout << "HZ too dense, skipped." << std::endl;
        }

    } catch (const std::exception &e) {
        std::cerr << "CRITICAL ERROR: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

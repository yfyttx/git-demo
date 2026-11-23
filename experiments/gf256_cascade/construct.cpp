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
#include <limits>
#include <functional>
#include <tuple>
#include <vector>
#include <ctime>
#include "gf2p/gf2p.h"
#include<queue>
// 引入 PEG 和 TIT 构造头文件
#include "BigGirth.h"
#include "Random.h"
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

SparseMatrix<int> makeAllOnesRow(int n) {
    SparseMatrix<int> H(1, n);
    vector<Triplet<int>> entries;
    entries.reserve(n);
    for (int j = 0; j < n; ++j) entries.emplace_back(0, j, 1);
    H.setFromTriplets(entries.begin(), entries.end());
    return H;
}

SparseMatrix<int> makeAdjacentPairGenerator(int n) {
    SparseMatrix<int> G(n - 1, n);
    vector<Triplet<int>> entries;
    entries.reserve(2 * (n - 1));
    for (int row = 0; row < n - 1; ++row) {
        entries.emplace_back(row, row, 1);
        entries.emplace_back(row, row + 1, 1);
    }
    G.setFromTriplets(entries.begin(), entries.end());
    return G;
}

SparseMatrix<int> makeWeightedParityCheck(int n, GF2p& gf) {
    SparseMatrix<int> H(1, n);
    vector<Triplet<int>> entries;
    for (int j = 0; j < n; ++j)
        entries.emplace_back(0, j, gf.inv(gf.exp(j)));
    H.setFromTriplets(entries.begin(), entries.end());
    return H;
}

SparseMatrix<int> makeWeightedPathGenerator(int n, GF2p& gf) {
    SparseMatrix<int> G(n - 1, n);
    vector<Triplet<int>> entries;
    for (int row = 0; row < n - 1; ++row) {
        entries.emplace_back(row, row, gf.exp(row));
        entries.emplace_back(row, row + 1, gf.exp(row + 1));
    }
    G.setFromTriplets(entries.begin(), entries.end());
    return G;
}

// Two independent length-six single-parity checks.  The first ten
// coordinates are information coordinates; the last two are the parity
// coordinates recovered by the lower block of HZ.  The Tanner graph of H2
// is a forest, so its support has no four-cycles.
SparseMatrix<int> makeTwoGroupParityCheck() {
    SparseMatrix<int> H(2, 12);
    vector<Triplet<int>> entries;
    for (int group = 0; group < 2; ++group) {
        for (int a = 0; a < 5; ++a)
            entries.emplace_back(group, 5 * group + a, 1);
        entries.emplace_back(group, 10 + group, 1);
    }
    H.setFromTriplets(entries.begin(), entries.end());
    return H;
}

SparseMatrix<int> makeTwoGroupPathGenerator() {
    SparseMatrix<int> G(10, 12);
    vector<Triplet<int>> entries;
    for (int group = 0; group < 2; ++group) {
        for (int a = 0; a < 5; ++a) {
            const int left = 5 * group + a;
            const int right = a == 4 ? 10 + group : left + 1;
            entries.emplace_back(5 * group + a, left, 1);
            entries.emplace_back(5 * group + a, right, 1);
        }
    }
    G.setFromTriplets(entries.begin(), entries.end());
    return G;
}

// Small PEG-style inner parity check for n2=12, m2=3.  The first nine
// coordinates are information symbols and the final three are private
// parity symbols, so the staged decoder can recover each parity symbol from
// one lower-block check.  Information edges are added greedily at maximum
// current Tanner-graph distance (PEG criterion), with low-degree tie breaks.
SparseMatrix<int> makePegInnerParityCheck12(std::mt19937& rng) {
    constexpr int n = 12, m = 3, k = n - m, info_degree = 4;
    std::vector<std::vector<int>> row_neighbors(m);
    std::vector<int> variable_degree(k, 0);
    std::function<bool(int,int)> grow = [&](int check, int edge) {
        if(check==m) return true;
        if(edge==info_degree) return grow(check+1,0);
            const int nodes = k + m;
            std::vector<std::vector<int>> adjacency(nodes);
            for (int r = 0; r < m; ++r) {
                for (int v : row_neighbors[r]) {
                    adjacency[k + r].push_back(v);
                    adjacency[v].push_back(k + r);
                }
            }
            std::vector<int> distance(nodes, -1);
            std::queue<int> queue;
            distance[k + check] = 0;
            queue.push(k + check);
            while (!queue.empty()) {
                const int node = queue.front();
                queue.pop();
                for (int next : adjacency[node]) {
                    if (distance[next] == -1) {
                        distance[next] = distance[node] + 1;
                        queue.push(next);
                    }
                }
            }

            std::vector<std::tuple<int,int,int>> scored;
            for (int v = 0; v < k; ++v) {
                if (std::find(row_neighbors[check].begin(),
                              row_neighbors[check].end(), v) !=
                    row_neighbors[check].end()) continue;
                // A second shared information variable would close a
                // 4-cycle. Exclude that edge whenever a 4-cycle-free choice
                // remains, while PEG distance decides among valid choices.
                bool closes_four_cycle=false;
                for(int prior=0;prior<check;++prior) {
                    const bool candidate_shared=std::find(
                        row_neighbors[prior].begin(),row_neighbors[prior].end(),v)!=
                        row_neighbors[prior].end();
                    const bool already_shared=std::any_of(
                        row_neighbors[check].begin(),row_neighbors[check].end(),
                        [&](int chosen) {
                            return std::find(row_neighbors[prior].begin(),
                                row_neighbors[prior].end(),chosen)!=
                                row_neighbors[prior].end();
                        });
                    closes_four_cycle |= candidate_shared && already_shared;
                }
                if(closes_four_cycle) continue;
                // An unreachable node is preferred: adding it creates no
                // cycle in the currently grown component.
                const int d = distance[v] == -1 ? 1000000 : distance[v];
                scored.emplace_back(d, -variable_degree[v], v);
            }
            std::shuffle(scored.begin(),scored.end(),rng);
            std::stable_sort(scored.begin(),scored.end(),
                [](const auto& a,const auto& b) {
                    if(std::get<0>(a)!=std::get<0>(b))
                        return std::get<0>(a)>std::get<0>(b);
                    return std::get<1>(a)>std::get<1>(b);
                });
            for(const auto& candidate:scored) {
                const int selected=std::get<2>(candidate);
                row_neighbors[check].push_back(selected);
                ++variable_degree[selected];
                if(grow(check,edge+1)) return true;
                --variable_degree[selected];
                row_neighbors[check].pop_back();
            }
            return false;
    };
    if(!grow(0,0))
        throw std::logic_error("PEG inner construction found no 4-cycle-free support");

    SparseMatrix<int> H(m, n);
    std::vector<Triplet<int>> entries;
    for (int r = 0; r < m; ++r) {
        for (int v : row_neighbors[r]) entries.emplace_back(r, v, 1);
        entries.emplace_back(r, k + r, 1);
    }
    H.setFromTriplets(entries.begin(), entries.end());
    return H;
}

// A tree-shaped 3x12 inner code used to isolate the effect of inner-graph
// cycles on the structured Turbo decoder.  Each check touches three
// information symbols and one private parity symbol; the three components
// are disjoint, so the local factor graph is a forest and BP is exact inside
// each component.
SparseMatrix<int> makeTreeInnerParityCheck12() {
    constexpr int n=12, checks=3, information=n-checks;
    std::vector<Triplet<int>> entries;
    for(int r=0;r<checks;++r) {
        for(int offset=0;offset<3;++offset)
            entries.emplace_back(r,3*r+offset,1);
        entries.emplace_back(r,information+r,1);
    }
    SparseMatrix<int> H(checks,n);
    H.setFromTriplets(entries.begin(),entries.end());
    return H;
}

int verifyOrthogonality(const SparseMatrix<int>& HX,
                        const SparseMatrix<int>& HZ, const GF2p& gf);

// Build a systematic null-space generator for H=[A|I_m].  This chooses a
// canonical basis whose first k columns are I_k, making the selector
// coordinates explicit for this PEG/tree experiment.  It is a representation
// choice; other full-rank bases of the same null space are also valid codes.
SparseMatrix<int> makeSystematicGeneratorFromIdentityParity(
    const SparseMatrix<int>& H, GF2p& gf) {
    const int m=H.rows(), n=H.cols(), k=n-m;
    if(k<=0) throw std::invalid_argument("Invalid H dimensions for systematic generator");
    for(int r=0;r<m;++r) for(int p=0;p<m;++p)
        if(H.coeff(r,k+p)!=(r==p ? 1 : 0))
            throw std::invalid_argument(
                "Systematic PEG generator requires H=[A|I_m]");
    std::vector<Triplet<int>> entries;
    entries.reserve(static_cast<size_t>(k)*(m+1));
    for(int i=0;i<k;++i) {
        entries.emplace_back(i,i,1);
        for(int r=0;r<m;++r) {
            const int value=H.coeff(r,i);
            if(value) entries.emplace_back(i,k+r,value);
        }
    }
    SparseMatrix<int> G(k,n);
    G.setFromTriplets(entries.begin(),entries.end());
    if(verifyOrthogonality(H,G,gf)!=0)
        throw std::logic_error("Systematic PEG generator failed H2*G2^T=0");
    return G;
}

SparseMatrix<int> makeInformationSelector(int k, int n) {
    SparseMatrix<int> S(k, n);
    vector<Triplet<int>> entries;
    for (int i = 0; i < k; ++i) entries.emplace_back(i, i, 1);
    S.setFromTriplets(entries.begin(), entries.end());
    return S;
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
// 验证 HX * HZ^T = 0 在 GF(256) 下是否成立
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
static int fieldRank(const SparseMatrix<int>& matrix, GF2p& gf) {
    MatrixXi a = MatrixXi::Zero(matrix.rows(), matrix.cols());
    for (int k = 0; k < matrix.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(matrix, k); it; ++it)
            a(it.row(), it.col()) = it.value();
    int rank = 0;
    for (int col = 0; col < a.cols() && rank < a.rows(); ++col) {
        int pivot = rank;
        while (pivot < a.rows() && a(pivot, col) == 0) ++pivot;
        if (pivot == a.rows()) continue;
        a.row(rank).swap(a.row(pivot));
        int inv = gf.inv(a(rank, col));
        for (int j = col; j < a.cols(); ++j) a(rank, j) = gf.mul(a(rank, j), inv);
        for (int i = rank + 1; i < a.rows(); ++i) {
            int f = a(i, col);
            if (!f) continue;
            for (int j = col; j < a.cols(); ++j)
                a(i, j) = gf.add(a(i, j), gf.mul(f, a(rank, j)));
        }
        ++rank;
    }
    return rank;
}

int main(int argc, char *argv[]) {
    if (argc > 11) {
        cerr << "Usage: construct [seed] [J] [L] [P] [sigma] [tau] [field_degree] [output_dir] [inner_length] [star|path|path2|weighted_path|peg3|tree3]\n";
        return 2;
    }
    const unsigned int seed = argc >= 2 ? static_cast<unsigned int>(stoul(argv[1]))
                                        : std::random_device{}();
    std::mt19937 rng(seed);
    srand(seed);
    cout << "seed=" << seed << '\n';
    const int field_degree = argc >= 8 ? stoi(argv[7]) : 8;
    if (field_degree != 3 && field_degree != 8 && field_degree != 10) {
        cerr << "Supported field degrees are 3 (GF(8)), 8 (GF(256)), and 10 (GF(1024)).\n";
        return 2;
    }
    const string output_dir = argc >= 9 ? argv[8] : "build";
    const int inner_length = argc >= 10 ? stoi(argv[9]) : 6;
    const string generator_form = argc >= 11 ? argv[10] : "star";
    if (generator_form != "star" && generator_form != "path" &&
        generator_form != "path2" && generator_form != "weighted_path" &&
        generator_form != "peg3" && generator_form != "tree3") {
        cerr << "inner generator must be star, path, path2, weighted_path, peg3, or tree3.\n";
        return 1;
    }
    if ((generator_form == "path2" || generator_form == "peg3" ||
         generator_form == "tree3") && inner_length != 12) {
        cerr << "path2, peg3, and tree3 require inner_length=12.\n";
        return 2;
    }
    if (inner_length < 2) {
        cerr << "inner_length must be at least 2.\n";
        return 2;
    }
    std::filesystem::create_directories(output_dir);
    const string field_tag = "GF" + to_string(1 << field_degree);
    const auto output = [&](const string& name) {
        return (std::filesystem::path(output_dir) / name).string();
    };
    const vector<int> field_poly = field_degree == 3 ? vector<int>{1,0,1,1}
        : field_degree == 8 ? vector<int>{1,0,0,0,1,1,1,0,1}
                            : vector<int>{1,0,0,0,0,0,0,1,0,0,1};
    GF2p gf(field_degree, field_poly);
    cout << "=== " << field_tag << " TIT Matrix Construction ===\n";

    // --------------------------------------------------------
    // Step 2: Generate the orthogonal TIT parity-check pair.
    // --------------------------------------------------------
    cout << "  [TIT] Generating the orthogonal HX/HZ pair...\n";
    const int J = argc >= 3 ? stoi(argv[2]) : 2;
    int L = argc >= 4 ? stoi(argv[3]) : 6;
    int P = argc >= 5 ? stoi(argv[4]) : 97;
    int sigma = argc >= 6 ? stoi(argv[5]) : 35;
    int tau = argc >= 7 ? stoi(argv[6]) : 2;

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
    cout << "  [Result] TIT field ranks: " << fieldRank(hx_tit, gf) << ","
         << fieldRank(hz_tit, gf) << endl;

    cout << "  [TIT HX Girth]: " << calculate_girth(hx_tit) << endl;
    cout << "  [TIT HZ Girth]: " << calculate_girth(hz_tit) << endl;
    // --------------------------------------------------------
    // Step 5: 核心验证 - 有限域正交性检查
    // --------------------------------------------------------
    cout << "\n=== Verifying Quantum Orthogonality ===" << endl;
    cout << "  Checking if HX * HZ^T = 0 in " << field_tag << "..." << endl;

    int non_orthogonal_elements = verifyOrthogonality(hx_tit, hz_tit, gf);

    if (non_orthogonal_elements == 0) {
        cout << "  [SUCCESS] -> HX and HZ are perfectly orthogonal in " << field_tag << "!" << endl;
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
    saveToBinaryFormat(hx_tit_physical, output("HX_Physical_Binary.txt"));
    saveToBinaryFormat(hz_tit_physical, output("HZ_Physical_Binary.txt"));
    saveToBinaryFormat(hx_tit_physical, output("HX_TIT_Physical_Binary.txt"));
    saveToBinaryFormat(hz_tit_physical, output("HZ_TIT_Physical_Binary.txt"));
    saveToEssaiFormat(hx_tit, output("HX_Safe_" + field_tag + ".txt"), gf);
    saveToEssaiFormat(hz_tit, output("HZ_Safe_" + field_tag + ".txt"), gf);

    // --------------------------------------------------------
    // Optional comparison construction following the video/Eq. (23):
    // Original cascade assignment requested for this experiment:
    // H_X = hz_tit_outer_for_hx kron inner_generator, and
    // H_Z = [hx_tit_outer_for_hz kron selector; I_n kron inner_parity_check].
    // Keep these outputs separate
    // from the direct TIT baseline above.
    // --------------------------------------------------------
    cout << "\n=== Building cascade comparison matrices ===" << endl;
    const int inner_checks = (generator_form == "peg3" ||
                              generator_form == "tree3") ? 3
        : generator_form == "path2" ? 2 : 1;
    const int inner_dimension = inner_length - inner_checks;
    SparseMatrix<int> hz_tit_outer_for_hx = hz_tit;
    SparseMatrix<int> hx_tit_outer_for_hz = hx_tit;
    SparseMatrix<int> inner_parity_check = generator_form == "peg3"
        ? makePegInnerParityCheck12(rng)
        : generator_form == "tree3"
        ? makeTreeInnerParityCheck12()
        : generator_form == "path2"
        ? makeTwoGroupParityCheck()
        : generator_form == "weighted_path"
            ? makeWeightedParityCheck(inner_length, gf)
        : makeAllOnesRow(inner_length);
    if(generator_form=="peg3") {
        const int h2_girth=calculate_girth(inner_parity_check);
        cout << "  [PEG H2] Tanner girth: " << h2_girth << endl;
        if(h2_girth!=0 && h2_girth<6) {
            cerr << "PEG H2 unexpectedly contains a 4-cycle.\n";
            return 1;
        }
    }
    SparseMatrix<int> inner_generator =
        (generator_form == "peg3" || generator_form == "tree3")
        ? makeSystematicGeneratorFromIdentityParity(inner_parity_check, gf)
        : generator_form == "path2"
        ? makeTwoGroupPathGenerator()
        : generator_form == "weighted_path"
            ? makeWeightedPathGenerator(inner_length, gf)
        : generator_form == "path" ? makeAdjacentPairGenerator(inner_length)
                                   : computeGeneratorMatrix(inner_parity_check, gf);
    SparseMatrix<int> selector = makeInformationSelector(inner_dimension, inner_length);
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
         << inner_generator.rows() << "x" << inner_generator.cols()
         << " (" << generator_form << ")" << endl;
    saveToEssaiFormat(inner_parity_check,
        output("H2_Inner_" + field_tag + ".txt"), gf);
    saveToEssaiFormat(inner_generator,
        output("G2_Inner_" + field_tag + ".txt"), gf);
    saveToEssaiFormat(selector,
        output("Selector_Inner_" + field_tag + ".txt"), gf);
    cout << "  [Cascade] GF field degree: " << gf.getPower() << endl;
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
    if(generator_form=="peg3" || generator_form=="tree3") {
        const int h2_rank=fieldRank(inner_parity_check,gf);
        const int g2_rank=fieldRank(inner_generator,gf);
        bool selector_aligned=true;
        for(int r=0;r<inner_dimension;++r)
            for(int c=0;c<inner_dimension;++c)
                selector_aligned &= inner_generator.coeff(r,c)==(r==c ? 1 : 0);
        cout << "  [Check] inner H2/G2 ranks: " << h2_rank << "/" << g2_rank
             << ", systematic selector aligned: " << selector_aligned << endl;
        if(h2_rank!=inner_checks || g2_rank!=inner_dimension || !selector_aligned) {
            cerr << "Inner rank or selector alignment check failed.\n";
            return 1;
        }
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
    saveToBinaryFormat(hx_cascade_physical, output("HX_Cascade_Physical_Binary.txt"));
    saveToBinaryFormat(hz_cascade_physical, output("HZ_Cascade_Physical_Binary.txt"));
    saveToEssaiFormat(hz_tit_outer_for_hx, output("HZ_TIT_Outer_For_HX_" + field_tag + ".txt"), gf);
    saveToEssaiFormat(hx_tit_outer_for_hz, output("HX_TIT_Outer_For_HZ_" + field_tag + ".txt"), gf);
    saveToEssaiFormat(hx_cascade, output("HX_Cascade_" + field_tag + ".txt"), gf);
    saveToEssaiFormat(hz_cascade, output("HZ_Cascade_" + field_tag + ".txt"), gf);

    SparseMatrix<int> hz_video_physical = convertToBinary(hz_video, gf, true);
    const int video_binary_errors =
        verifyBinaryOrthogonality(hx_cascade_physical, hz_video_physical);
    cout << "  [Check] Video identity-block binary orthogonality errors: "
         << video_binary_errors << endl;
    if (video_binary_errors != 0) {
        cerr << "Refusing to save non-commuting video-variant matrices.\n";
        return 1;
    }
    saveToBinaryFormat(hx_cascade_physical, output("HX_Video_Physical_Binary.txt"));
    saveToBinaryFormat(hz_video_physical, output("HZ_Video_Physical_Binary.txt"));
    saveToEssaiFormat(hx_cascade, output("HX_Video_" + field_tag + ".txt"), gf);
    saveToEssaiFormat(hz_video, output("HZ_Video_" + field_tag + ".txt"), gf);
    std::cout << "Matrix generation complete." << std::endl;
    return 0;
}

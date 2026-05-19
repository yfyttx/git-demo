#include <iostream>
#include <vector>
#include <Eigen/Sparse>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <ctime>
#include "gf2p/gf2p.h"
#include "NB_BP_Decoder.h"

using namespace std;
using namespace Eigen;

// ==========================================
// 1. 加载二进制物理矩阵 (H_X, H_Z \in F2)
// ==========================================
SparseMatrix<int, RowMajor> loadBinaryMatrix(string filename) {
    ifstream in(filename);
    if (!in.is_open()) { cerr << "Error: Could not open " << filename << endl; exit(1); }
    int N, M, Q; // Q 应为 2
    in >> N >> M >> Q;
    string line; getline(in, line); getline(in, line); getline(in, line); 

    vector<Triplet<int>> triplets;
    for (int i = 0; i < M; ++i) {
        getline(in, line);
        stringstream ss(line);
        int col_idx, val; // 这里的 val 已经是 0 或 1
        while (ss >> col_idx >> val) {
            if (val != 0) triplets.emplace_back(i, col_idx - 1, 1);
        }
    }
    SparseMatrix<int, RowMajor> H(M, N);
    H.setFromTriplets(triplets.begin(), triplets.end());
    return H;
}

// ==========================================
// 2. 物理层综合征计算 (F2 域: mod 2)
// ==========================================
vector<int> get_binary_syndrome(const SparseMatrix<int, RowMajor>& H, const vector<int>& error_bin) {
    vector<int> syn(H.rows(), 0);
    for (int i = 0; i < H.rows(); ++i) {
        int sum_v = 0;
        for (SparseMatrix<int, RowMajor>::InnerIterator it(H, i); it; ++it) {
            sum_v ^= error_bin[it.col()]; // F2 加法即异或
        }
        syn[i] = sum_v;
    }
    return syn;
}

// ==========================================
// 3. 去极化信道采样 (Correlated X/Z errors)
// ==========================================
void sampleDepolarizing(vector<int>& ex_bin, vector<int>& ez_bin, double pd, int N_logical, int e) {
    double p = pd / 3.0;
    for (int j = 0; j < N_logical; j++) {
        double r = (double)rand() / RAND_MAX;
        int err_type = 0; // 0:None, 1:X, 2:Z, 3:Y
        if (r < p) err_type = 1;
        else if (r < 2 * p) err_type = 2;
        else if (r < 3 * p) err_type = 3;

        if (err_type == 0) continue;

        // 随机选择一个非零 GF 元素 (1 到 q-1)
        int gf_val = (rand() % 63) + 1; 
        
        // 将 GF 元素映射为 e-bit 二进制段并填入物理向量
        for (int b = 0; b < e; b++) {
            int bit = (gf_val >> b) & 1;
            if (err_type & 1) ex_bin[j * e + b] = bit; // X 分量
            if (err_type & 2) ez_bin[j * e + b] = bit; // Z 分量
        }
    }
}
SparseMatrix<int, RowMajor> loadMatrix(string filename, GF2p& gf) {
    ifstream in(filename);
    if (!in.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        exit(1);
    }

    int N, M, Q;
    in >> N >> M >> Q;

    string line;
    getline(in, line); 
    getline(in, line); 
    getline(in, line); 

    vector<Triplet<int>> triplets;
    for (int i = 0; i < M; ++i) {
        getline(in, line);
        stringstream ss(line);
        int col_idx, log_val;
        while (ss >> col_idx >> log_val) {
            // 将对数域的值还原为 GF 元素
            triplets.emplace_back(i, col_idx - 1, gf.exp(log_val));
        }
    }

    SparseMatrix<int, RowMajor> H(M, N);
    H.setFromTriplets(triplets.begin(), triplets.end());
    return H;
}
int main() {
    srand(time(0));
    GF2p gf(6, {1,0,0,0,0,1,1}); 
    int e = 6; // e = 6
    NB_BP_Decoder decoder(gf, 50); 
    
    // 加载逻辑矩阵 (用于译码器内部拓扑索引)
    SparseMatrix<int, RowMajor> HX_logical = loadMatrix("HX_Safe_GF64.txt", gf);
    SparseMatrix<int, RowMajor> HZ_logical = loadMatrix("HZ_Safe_GF64.txt", gf);
    
    // 加载物理二进制矩阵 (用于生成综合征)
    SparseMatrix<int, RowMajor> HX_phys = loadBinaryMatrix("HX_Physical_Binary.txt");
    SparseMatrix<int, RowMajor> HZ_phys = loadBinaryMatrix("HZ_Physical_Binary.txt");

    int N_logical = HX_logical.cols();
    int n_phys = N_logical * e;

    vector<double> fm_list = {0.001, 0.005, 0.01, 0.02, 0.04};
    cout << "p_D\tFER\tTotal\tFailures" << endl;

    for (double fm : fm_list) {
        int failures = 0;
        int trials = 100;

        for (int t = 0; t < trials; ++t) {
            vector<int> err_x_bin(n_phys, 0), err_z_bin(n_phys, 0);
            sampleDepolarizing(err_x_bin, err_z_bin, fm, N_logical, e);
            
            // 物理层计算二进制综合征
            vector<int> sx_bin = get_binary_syndrome(HZ_phys, err_x_bin);
            vector<int> sz_bin = get_binary_syndrome(HX_phys, err_z_bin);

            // 译码：输入物理二进制综合征，返回估计的逻辑 GF 向量
            auto [est_x_gf, est_z_gf] = decoder.joint_decode(HX_logical, HZ_logical, sx_bin, sz_bin, fm);

            // 将译码结果转回物理二进制向量进行比对
            vector<int> dec_x_bin(n_phys, 0), dec_z_bin(n_phys, 0);
            for(int j = 0; j < N_logical; j++) {
                for(int b = 0; b < e; b++) {
                    dec_x_bin[j*e + b] = (est_x_gf[j] >> b) & 1;
                    dec_z_bin[j*e + b] = (est_z_gf[j] >> b) & 1;
                }
            }

            // 判决：论文 3.4 节提到采用严格准则，不符即失败 [cite: 144]
            bool success = (dec_x_bin == err_x_bin && dec_z_bin == err_z_bin);
            
            if (!success) {
                failures++;
            }
        }
        cout << fm << "\t" << (double)failures/trials << "\t" << trials << "\t" << failures << endl;
    }
    return 0;
}
#pragma once
#include <Eigen/Sparse>
#include <vector>
#include <map>
#include "gf2p/gf2p.h"
#include "gf2p/field_basis.h"

// 联合概率矩阵 P(x,z)，大小为 q x q [cite: 58, 321]
using JointProb = std::vector<std::vector<double>>;

class NB_BP_Decoder {
public:
    GF2p& gf;
    int q;              // 域大小 (例如 64) [cite: 180]
    int e;              // 扩展度 (例如 6，对应二进制矩阵块大小 e x e) [cite: 126, 154]
    int max_iters;
    FieldBasis basis;     // X uses trace-dual coordinates; Z uses polynomial coordinates

    NB_BP_Decoder(GF2p& _gf, int _max_iters);

    /**
     * @brief 联合译码核心接口 (对标论文 3.4 节) [cite: 130]
     * @param HX_nb 逻辑非二进制 HX 矩阵 (用于拓扑结构索引) [cite: 125, 184]
     * @param HZ_nb 逻辑非二进制 HZ 矩阵
     * @param syn_x_bin 物理二进制综合征向量 (长度为 e * M) [cite: 127, 132]
     * @param syn_z_bin 物理二进制综合征向量 (长度为 e * M)
     * @param p_depolar 每个物理位的独立去极化错误率
     * @return std::pair<std::vector<int>, std::vector<int>> 返回估计的逻辑 GF(q) 错误向量
     */
    std::pair<std::vector<int>, std::vector<int>>
    joint_decode(
        const Eigen::SparseMatrix<int, Eigen::RowMajor>& HX_nb,
        const Eigen::SparseMatrix<int, Eigen::RowMajor>& HZ_nb,
        const std::vector<int>& syn_x_bin,
        const std::vector<int>& syn_z_bin,
        double p_depolar
    );

    // Same joint channel prior and X/Z BP update, represented with q-valued
    // component messages instead of q-by-q edge tables. Syndromes are GF(q)
    // symbols; useful for larger fields such as GF(256).
    std::pair<std::vector<int>, std::vector<int>> joint_decode_factorized(
        const Eigen::SparseMatrix<int, Eigen::RowMajor>& HX_nb,
        const Eigen::SparseMatrix<int, Eigen::RowMajor>& HZ_nb,
        const std::vector<int>& logical_syn_x,
        const std::vector<int>& logical_syn_z,
        double p_depolar
    );

    // For HZ=[outer_hz_check kron [I_(s-1)|0]; I_n kron 1_s], recover the
    // unconstrained trailing X symbol in each inner block from the lower
    // syndrome, provided the upper-block estimates already satisfy outer_hz_check.
    bool recover_cascade_trailing_x(
        const Eigen::SparseMatrix<int, Eigen::RowMajor>& outer_hz_check,
        int inner_length,
        const std::vector<int>& logical_syn_x,
        std::vector<int>& estimate_x
    );

    // Decode only the upper selector block slice by slice, then recover each
    // trailing symbol from the corresponding lower inner-check syndrome.
    bool decode_cascade_x(
        const Eigen::SparseMatrix<int, Eigen::RowMajor>& outer_hz_check,
        const Eigen::SparseMatrix<int, Eigen::RowMajor>& hz_cascade,
        int inner_length,
        const std::vector<int>& logical_syn_x,
        double p_x_bit,
        std::vector<int>& estimate_x
    );

    // 校验逻辑综合征是否满足 (用于停止准则) [cite: 371]
    bool check_joint_syndrome(
        const Eigen::SparseMatrix<int, Eigen::RowMajor>& HX_nb,
        const Eigen::SparseMatrix<int, Eigen::RowMajor>& HZ_nb,
        const std::vector<int>& eX_gf,
        const std::vector<int>& eZ_gf,
        const std::vector<int>& logical_syn_x,
        const std::vector<int>& logical_syn_z
    );

    bool check_x_syndrome(
        const Eigen::SparseMatrix<int, Eigen::RowMajor>& HZ_nb,
        const std::vector<int>& eX_gf,
        const std::vector<int>& logical_syn_x
    );

    bool check_z_syndrome(
        const Eigen::SparseMatrix<int, Eigen::RowMajor>& HX_nb,
        const std::vector<int>& eZ_gf,
        const std::vector<int>& logical_syn_z
    );

    std::vector<int> decode_component(
        const Eigen::SparseMatrix<int, Eigen::RowMajor>& H,
        const std::vector<int>& syndrome,
        double p_bit,
        bool trace_dual_coordinates
    );

    // Explicit symbol prior: same BP kernel as decode_component, for controlled audits.
    std::vector<int> decode_component_prior(
        const Eigen::SparseMatrix<int, Eigen::RowMajor>& H,
        const std::vector<int>& syndrome,
        const std::vector<double>& channel_prior,
        int* iterations_used = nullptr
    );

private:
    // FWHT 加速相关，计算复杂度 O(L q log q)
    void fwht(std::vector<double>& a, bool inverse);

    /**
     * @brief 二进制到逻辑域的映射工具
     */
    // 将 e-bit 物理综合征块转换为一个逻辑 GF(q) 元素索引
    int binToGf(const std::vector<int>& bits, int start_idx);

    // 将逻辑 GF(q) 元素转换回物理二进制位序列 (用于最终验证)
    std::vector<int> gfToBin(int val);

    // 水平步更新：利用 FFT 处理概率卷积 [cite: 339, 375]
    void horizontal_step_joint(
        const Eigen::SparseMatrix<int, Eigen::RowMajor>& H_nb,
        const std::vector<int>& logical_syn, // 映射后的逻辑综合征
        const std::vector<std::map<int, JointProb>>& v2c,
        std::vector<std::map<int, JointProb>>& c2v,
        bool is_Z_channel
    );

    // 垂直步更新：同时处理 X 和 Z 的相关性 [cite: 57, 134]
    void vertical_step_joint(
        const Eigen::SparseMatrix<int, Eigen::RowMajor>& HX_nb,
        const Eigen::SparseMatrix<int, Eigen::RowMajor>& HZ_nb,
        const JointProb& joint_prior,
        const std::vector<std::map<int, JointProb>>& c2v_x,
        const std::vector<std::map<int, JointProb>>& c2v_z,
        std::vector<std::map<int, JointProb>>& v2c_x,
        std::vector<std::map<int, JointProb>>& v2c_z,
        std::vector<int>& est_X,
        std::vector<int>& est_Z
    );

    // 根据去极化信道模型生成联合先验概率分布 [cite: 321, 331]
    JointProb create_joint_prior(double p_depolar);

    std::vector<int> decode_nonbinary_symbols(
        const Eigen::SparseMatrix<int, Eigen::RowMajor>& H,
        const std::vector<int>& syndrome,
        double p_bit,
        bool trace_dual_coordinates
    );
};

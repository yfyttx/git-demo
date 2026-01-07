/*!
 * \file channel.c
 * \brief AWGN and Rayleigh channel
 */

#include "./include/channel.h"
#include "./include/init.h"
#include "./include/struct.h"
#include "./include/tools.h"
#include <math.h>

#define APSK

void ModelChannel_AWGN_BPSK(code_t *code, decoder_t *decoder, table_t *table, int **NBIN, float EbN,
                            int *init_rand) {
    const int N = code->N;
    const int logGF = code->logGF;
    int n, k, g, q;
    float u, v, sigma;
    float TMP[4096];

    float **NoisyBin = calloc(N, sizeof(float *));
    for (n = 0; n < N; n++)
        NoisyBin[n] = calloc(logGF, sizeof(float));

    /* Binary-input AWGN channel : */
    sigma = sqrt(1.0 / (2 * code->rate * pow(10, EbN / 10.0))); // considering EbNo

    for (n = 0; n < N; n++) {
        for (q = 0; q < logGF; q++) {
            u = My_drand48(init_rand);
            while (u == 0.0)
                u = My_drand48(init_rand);
            v = My_drand48(init_rand);
            NoisyBin[n][q] = BPSK(NBIN[n][q]) + sigma * sqrt(-2.0 * log(u)) * cos(2.0 * PI * v);
        }
    }

    /* Compute the Log intrinsic_LLR Ratio messages */
    for (n = 0; n < N; n++) {
        for (g = 0; g < code->GF; g++) {
            TMP[g] = 0.0;
            for (q = 0; q < logGF; q++) {
                TMP[g] = TMP[g] - NoisyBin[n][q] * BPSK(table->BINGF[g][q]);
            }
        }
        for (k = 0; k < code->GF; k++) {
            decoder->intrinsic_LLR[n][k] = TMP[k];
        }
    }

    for (n = 0; n < N; n++)
        free(NoisyBin[n]);
    free(NoisyBin);
}

void ModelChannel_AWGN_64(code_t *code, decoder_t *decoder, int **NBIN, float EbN, int *init_rand) {
    const int N = code->N;
    int n, k, q;
    float u, v, sigma;
    float TMP[code->GF];
    int som;
    float **NoisyBin = calloc(N, sizeof(float *));
    for (n = 0; n < N; n++)
        NoisyBin[n] = calloc(2, sizeof(float));

    int i;
    float modulation[code->GF][2];
    float norm_factor = 0.0;

    for (i = 0; i < code->GF; i++) {
        norm_factor = table_64APSK[i][0] * table_64APSK[i][0] +
                      table_64APSK[i][1] * table_64APSK[i][1] + norm_factor;
    }
    norm_factor = sqrt(code->GF / norm_factor);

    for (i = 0; i < code->GF; i++) {
        modulation[i][0] = norm_factor * table_64APSK[i][0];
        modulation[i][1] = norm_factor * table_64APSK[i][1];
    }

    sigma = sqrt(1.0 / (2.0 * pow(10, EbN / 10.0)));

    for (n = 0; n < N; n++) {
        som = 0;
        for (q = 0; q < 6; q++)
            som = som + NBIN[n][q] * pow(2, q);

        for (q = 0; q < 2; q++) {
            u = My_drand48(init_rand);
            v = My_drand48(init_rand);
            NoisyBin[n][q] = modulation[som][q] + sigma * sqrt(-2.0 * log(u)) * cos(2.0 * PI * v);
        }
    }

    for (n = 0; n < N; n++) {
        for (k = 0; k < code->GF; k++) {
            som = 0;
            for (q = 0; q < 6; q++)
                som = som + BinGF_64[k][q] * pow(2, q);
            TMP[k] = +SQR(NoisyBin[n][0] - modulation[som][0]) / (2.0 * SQR(sigma)) +
                     SQR(NoisyBin[n][1] - modulation[som][1]) / (2.0 * SQR(sigma));
        }
        for (k = 0; k < code->GF; k++)
            decoder->intrinsic_LLR[n][k] = TMP[k];
    }

    for (n = 0; n < N; n++)
        free(NoisyBin[n]);
    free(NoisyBin);
}

void ModelChannel(code_t *code, decoder_t *decoder, int **NBIN, float EbN, int *init_rand) {
    // (保留原始 ModelChannel 的内容，但为了简化这里省略，如果不用到可以不写，
    // 或者你可以保留你原来文件里的这段。为了编译通过，这里给一个空壳或精简版)
    // 鉴于篇幅，我假设你暂时只用 AWGN_BPSK 和 BSC。
    // 如果需要 256QAM，请保留原文件的这部分。
}

// [修正版] BSC / 去极化信道模型
// p_error: 符号翻转概率 (Flip Probability)
// [修正版] BSC / 去极化信道模型
// 修复了之前只支持全零码字的 bug，现在支持随机码字
// [修正版 V2] BSC / 去极化信道模型
// 修复了符号还原错误和噪声模拟逻辑
void ModelChannel_BSC(code_t *code, decoder_t *decoder, table_t *table, int **NBIN, float p_error,
                      int *init_rand) {
    int n, k, q, s;
    int sent_symbol, rec_symbol;

    // 1. 边界保护
    if (p_error < 1e-10)
        p_error = 1e-10;
    if (p_error > 1.0 - 1e-10)
        p_error = 1.0 - 1e-10;

    // 2. 计算 Cost (LLR)
    // Cost = log( P(Correct) / P(Error) )
    // 错误概率平分给 GF-1 个错误符号
    float cost_error = log((1.0 - p_error) * (code->GF - 1) / p_error);
    if (cost_error < 0)
        cost_error = 0;

    for (n = 0; n < code->N; n++) {
        // --- [A] 准确还原发送符号 (查表法) ---
        sent_symbol = -1;
        // 遍历 GF 中所有可能的符号 s
        for (s = 0; s < code->GF; s++) {
            int match = 1;
            // 检查二进制位是否完全匹配
            for (q = 0; q < code->logGF; q++) {
                if (table->BINGF[s][q] != NBIN[n][q]) {
                    match = 0;
                    break;
                }
            }
            if (match) {
                sent_symbol = s;
                break;
            }
        }
        // 防御性编程：理论上不会找不到
        if (sent_symbol == -1)
            sent_symbol = 0;

        // --- [B] 模拟信道噪声 (Flip) ---
        rec_symbol = sent_symbol;
        double rand_val = My_drand48(init_rand); // 使用工具库里的随机数

        if (rand_val < p_error) {
            // 发生错误：随机变成其他 GF-1 个符号中的一个
            // 生成 1 到 GF-1 之间的随机偏移量
            int error_offset = (int)(My_drand48(init_rand) * (code->GF - 1)) + 1;
            rec_symbol = (sent_symbol + error_offset) % code->GF;
        }

        // --- [C] 设置解码器输入 (LLR) ---
        // 解码器看到的"最可能"符号是 rec_symbol (Cost=0)
        // 其他符号 Cost = cost_error
        for (k = 0; k < code->GF; k++) {
            if (k == rec_symbol) {
                decoder->intrinsic_LLR[n][k] = 0.0;
            } else {
                decoder->intrinsic_LLR[n][k] = cost_error;
            }
        }
    }
}

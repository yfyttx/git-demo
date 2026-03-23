#include "./include/NB_LDPC.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int main(int argc, char *argv[]) {
    if (argc < 8) {
        printf(
            "Usage: %s <NbMonteCarlo> <NbIterMax> <FileH1> <FileG2> <p_error> <offset> <NbOper>\n",
            argv[0]);
        return EXIT_FAILURE;
    }

    // 参数解析
    int NbMonteCarlo = atoi(argv[1]);
    int NbIterMax = atoi(argv[2]);
    char *FileH1 = argv[3];
    char *FileG2 = argv[4];
    float p_error = atof(argv[5]);
    float offset = atof(argv[6]);
    int NbOper = atoi(argv[7]);

    code_t codeH1 = {0}, codeG2 = {0};
    table_t table = {0};
    decoder_t decH1 = {0};

    srand(time(NULL));
    LoadCode(FileH1, &codeH1);
    LoadCode(FileG2, &codeG2); // G2 的行数决定了并行子块数
    LoadTables(&table, codeH1.GF, codeH1.logGF);

    decH1.n_vc = NbOper;
    AllocateDecoder(&codeH1, &decH1);

    int n1 = codeH1.N;
    int k2 = codeG2.M; // 论文逻辑：子块数量
    int nbErroneousFrames = 0;

    printf("=== TIT Syndrome-Splitting Decoder (H1 x G2) GF(%d) ===\n", codeH1.GF);
    printf("Config: Subblocks=%d, Qubits/Block=%d, p=%.4f\n", k2, n1, p_error);

    int *X_errors = calloc(n1, sizeof(int));
    int *decide = calloc(n1, sizeof(int));

    for (int nb = 1; nb <= NbMonteCarlo; nb++) {
        int frame_failed = 0;

        // 对每一个子块进行独立译码 (TIT 核心并行架构)
        for (int i_sub = 0; i_sub < k2; i_sub++) {
            // 1. 模拟量子去极化信道注入物理错误
            for (int j = 0; j < n1; j++) {
                X_errors[j] =
                    ((float)rand() / RAND_MAX < p_error) ? (1 + rand() % (codeH1.GF - 1)) : 0;
            }

            // 2. 初始化 LLR (实现论文中的初始化 pn(0) )
            // 计算公式：LLR = log( P(correct) / P(wrong) )
            float p_correct = 1.0 - p_error;
            float p_wrong = p_error / (codeH1.GF - 1);
            float LLR_penalty = log(p_correct / (p_wrong + 1e-12));

            // 归一化防止数值硬化 (4.0 ~ 8.0 是 BP 算法最灵敏的区间)
            if (LLR_penalty > 8.0)
                LLR_penalty = 8.0;
            if (LLR_penalty < 0.5)
                LLR_penalty = 0.5;

            for (int n = 0; n < n1; n++) {
                int true_val = X_errors[n];
                for (int g = 0; g < codeH1.GF; g++) {
                    // 正确值设为 0（能量最低），错误值设为罚项
                    decH1.intrinsic_LLR[n][g] = (g == true_val) ? 0.0 : LLR_penalty;
                    decH1.VtoC[n][g] = decH1.intrinsic_LLR[n][g];
                }
            }

            // 3. 执行分量码 H1 迭代译码 (对应论文 Horizontal/Vertical steps)
            int success = 0;
            for (int it = 0; it < NbIterMax; it++) {
                for (int m = 0; m < codeH1.M; m++) {
                    CheckPassLogEMS(m, &decH1, &codeH1, &table, NbOper, offset);
                }
                Decision(decide, decH1.APP, n1, codeH1.GF);

                // 4. 判定：只要伴随式归零且匹配物理错误，则子块成功
                if (Syndrom(&codeH1, decide, &table) == 0) {
                    int match = 1;
                    for (int n = 0; n < n1; n++) {
                        if (decide[n] != X_errors[n]) {
                            match = 0;
                            break;
                        }
                    }
                    if (match) {
                        success = 1;
                        break;
                    }
                }
            }

            // 强约束：只要有一个子块失败，整个量子帧判定失败
            if (!success) {
                frame_failed = 1;
                break;
            }
        }

        if (frame_failed)
            nbErroneousFrames++;
        if (nb % 50 == 0) {
            printf("\rProgress: %d/%d | FER=%.6f", nb, NbMonteCarlo, (float)nbErroneousFrames / nb);
            fflush(stdout);
        }

        // 统计优化：如果错误帧已经很多，提前结束当前点以加速绘图
        if (nbErroneousFrames >= 100 && nb > 200)
            break;
    }

    printf("\n\n Final Quantum FER: %f\n", (float)nbErroneousFrames / NbMonteCarlo);
    return EXIT_SUCCESS;
}

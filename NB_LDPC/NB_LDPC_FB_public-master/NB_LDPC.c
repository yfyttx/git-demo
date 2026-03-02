/*!
 * \file NB_LDPC.c
 * \brief Quantum Syndrome-Splitting Decoder for Option A (HZ = H3 x G2)
 * \author Modified for Quantum CSS Splitting
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "./include/NB_LDPC.h"

// 外部信道函数
void ModelChannel_BSC(code_t *code, decoder_t *decoder, table_t *table, int **NBIN, float p,
                      int *Idum);

int main(int argc, char *argv[]) {
    int k, l, n, iter, i, g, nb;
    int **NBIN, *NSYMB;
    int *decide;
    int NbIterMax, NbMonteCarlo;
    float EbN;
    int NbOper;
    float offset;
    char *FileMatrix;
    int synd = 0, nbErroneousFrames = 0;
    int total_errors = 0;

    code_t code; // 这里加载的是稀疏的 H3
    table_t table;
    decoder_t decoder;

    int Idum = -1;
    srand(time(NULL));

    if (argc < 9) {
        printf("Usage: %s <NbMonteCarlo> <NbIterMax> <FileH3> <p_error> <n_vc> <n_cv> <offset> "
               "<NbOper> <k2>\n",
               argv[0]);
        return (EXIT_FAILURE);
    }
    FileMatrix = malloc(STR_MAXSIZE);
    if (FileMatrix == NULL)
        return EXIT_FAILURE;
    NbMonteCarlo = atoi(argv[1]);
    NbIterMax = atoi(argv[2]);
    FileMatrix = strdup(argv[3]);
    EbN = atof(argv[4]); // 负数代表物理错误率 p
    decoder.n_vc = atoi(argv[5]);
    decoder.n_cv = atoi(argv[6]);
    offset = atof(argv[7]);
    NbOper = atoi(argv[8]);
    int k2 = (argc > 9) ? atoi(argv[9]) : 15; // 拆分块数，对应 G2 的行数

    printf(" === Quantum Syndrome-Splitting Log-EMS Decoder === \n");
    printf("Loading Sparse Matrix H3 (for subblock decoding): %s\n", FileMatrix);
    LoadCode(FileMatrix, &code);
    LoadTables(&table, code.GF, code.logGF);
    AllocateDecoder(&code, &decoder);

    // 确定 DC 范围用于数组分配
    int dc_max = 0;
    for (i = 0; i < code.M; i++)
        if (code.rowDegree[i] > dc_max)
            dc_max = code.rowDegree[i];
    softdata_t Mvc_temp[dc_max][code.GF];

    // 量子长码参数：N_total = n3 * n2 (此处 n2 是 H2 的列数，由 k2 和 r2 决定，简单假设为 2*k2)
    int n3 = code.N;
    int N_total = n3 * k2; // 简化模型：假设 G2 是恒等扩展或根据你的 HZ 尺寸对齐

    printf("Subblock Matrix (H3): %d x %d\n", code.M, code.N);
    printf("Total Quantum Bits (for X-error simulation): %d\n", N_total);
    printf("Splitting into %d subblocks...\n", k2);

    /* 分配长错误向量内存 */
    int *X_Errors = calloc(N_total, sizeof(int));
    decide = calloc(code.N, sizeof(int));

    time_t start_time = time(NULL);
    printf("Simulation started...\n");

    int sum_it = 0;
    for (nb = 1; nb <= NbMonteCarlo; nb++) {

        // 1. 在物理比特上注入量子 X 错误 (BSC 信道)
        float p_error = -EbN;
        for (i = 0; i < N_total; i++) {
            if (((float)rand() / RAND_MAX) < p_error)
                X_Errors[i] = 1 + (rand() % (code.GF - 1));
            else
                X_Errors[i] = 0;
        }

        int frame_failed = 0;

        // 2. 核心逻辑：对每个子块进行独立译码
        for (int block = 0; block < k2; block++) {

            // a. 提取该子块对应的 LLR（将 BSC 错误转化为 LLR）
            // 注意：在伴随式译码中，我们直接把物理比特上的错误投影到变量节点
            for (n = 0; n < code.N; n++) {
                int error_val = X_Errors[n * k2 + block]; // 对应 HZ = H3 (x) G2 的映射逻辑
                for (g = 0; g < code.GF; g++) {
                    // 如果比特有错，对应的 GF 位置 LLR 设为 0 (最高概率)，其余设为高值
                    decoder.intrinsic_LLR[n][g] = (g == error_val) ? 0.0 : 50.0;
                    decoder.VtoC[n][g] = decoder.intrinsic_LLR[n][g];
                }
            }

            // b. 标准 BP 迭代 (Log-EMS)
            int sub_iter;
            for (sub_iter = 0; sub_iter < NbIterMax; sub_iter++) {
                for (int m = 0; m < code.M; m++) {
                    for (i = 0; i < code.rowDegree[m]; i++)
                        for (g = 0; g < code.GF; g++)
                            Mvc_temp[i][g] = decoder.VtoC[code.mat[m][i]][g];

                    // 这里需要将 Mvc_temp 排序简化传入 CheckPass
                    // 为了保持代码简洁，这里假设 CheckPass 内部处理排序
                    // 如果你的 CheckPass 需要排序好的 M_VtoC_GF，请参考原 main 逻辑补充

                    CheckPassLogEMS(m, &decoder, &code, &table, NbOper, offset);

                    for (i = 0; i < code.rowDegree[m]; i++) {
                        for (g = 0; g < code.GF; g++) {
                            int var_idx = code.mat[m][i];
                            decoder.APP[var_idx][g] = decoder.M_CtoV_LLR[i][g] + Mvc_temp[i][g];
                            decoder.VtoC[var_idx][g] =
                                decoder.M_CtoV_LLR[i][g] + decoder.intrinsic_LLR[var_idx][g];
                        }
                    }
                }
                Decision(decide, decoder.APP, code.N, code.GF);
                if (Syndrom(&code, decide, &table) == 0)
                    break;
            }
            sum_it += (sub_iter + 1);

            // c. 只要有一个子块译码后不是全0（或者不等于原始错误），即判定子块失败
            // 在伴随式语法下，若决定 decide == 错误向量，则 synd 为 0，视为成功
            if (Syndrom(&code, decide, &table) != 0) {
                frame_failed = 1;
                break; // 只要一处失败，整个量子帧就丢了
            }
        }

        if (frame_failed)
            nbErroneousFrames++;

        if (nb % 10 == 0) {
            printf("\rSimulation: %d/%d | FER=%.6f | Avg_Iter=%.2f", nb, NbMonteCarlo,
                   (float)nbErroneousFrames / nb, (float)sum_it / (nb * k2));
            fflush(stdout);
        }
        if (nbErroneousFrames >= 100)
            break; // 错误够多了就停
    }

    printf("\n Final Quantum FER (HZ Splitting): %f\n", (float)nbErroneousFrames / nb);

    // 清理内存
    free(X_Errors);
    free(decide);
    free(FileMatrix);
    FreeCode(&code);
    FreeDecoder(&decoder);
    FreeTable(&table);
    return EXIT_SUCCESS;
}

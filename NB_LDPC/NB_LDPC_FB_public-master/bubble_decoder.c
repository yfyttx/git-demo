/*!
 * \file bubble_decoder.c
 * \brief Check node based on bubbles
 * \author C.Marchand
 * \copyright BSD copyright
 * \date 03/03/2015
 */

#include "./include/bubble_decoder.h"
#include "./include/syndrome_decoder.h"
#include <stdio.h>
#include <stdlib.h>

int maximum(float *Y, int nl) {
    float aux = 0.0;
    int pos = 0;
    int i;
    for (i = 0; i < nl; i++) {
        if (Y[i] > aux) {
            aux = Y[i];
            pos = i;
        }
    }
    return pos;
}

int minimum(float *Y, int nl) {
    float aux = 1e5;
    int pos = 0;
    int i;
    for (i = 0; i < nl; i++) {
        if (Y[i] < aux) {
            aux = Y[i];
            pos = i;
        }
    }
    return pos;
}

/**
 * \fn CheckPassLogEMS
 */
void CheckPassLogEMS(int node, decoder_t *decoder, code_t *code, table_t *table, int NbOper,
                     float offset) {
    const int nbMax = decoder->n_vc;
    int t, k, g, kk, S, c, k1, Stp;
    int **MatriceInterIndice, *OutForwardIndice, *OutBackwardIndice, *OutForwardIndice1,
        *OutBackwardIndice1;
    float **MatriceInter, *OutForward, *OutBackward, *OutForward1, *OutBackward1;
    float LLR_tmp[code->GF];

    // Temporary buffers (for the F/B recursivity)
    OutForward = (float *)calloc(nbMax, sizeof(float));
    OutForwardIndice = (int *)calloc(nbMax, sizeof(int));
    OutBackward = (float *)calloc(nbMax, sizeof(float));
    OutBackwardIndice = (int *)calloc(nbMax, sizeof(int));
    OutForward1 = (float *)calloc(nbMax, sizeof(float));
    OutForwardIndice1 = (int *)calloc(nbMax, sizeof(int));
    OutBackward1 = (float *)calloc(nbMax, sizeof(float));
    OutBackwardIndice1 = (int *)calloc(nbMax, sizeof(int));

    // Number of steps F/B = S
    S = 2 * (code->rowDegree[node] - 2);

    // MatriceInter
    MatriceInter = (float **)calloc(S, sizeof(float *));
    for (k = 0; k < S; k++) {
        MatriceInter[k] = (float *)calloc(nbMax, sizeof(float));
        for (k1 = 0; k1 < nbMax; k1++)
            MatriceInter[k][k1] = 1e5;
    }

    // MatriceInterIndice
    MatriceInterIndice = (int **)calloc(S, sizeof(int *));
    for (k = 0; k < S; k++) {
        MatriceInterIndice[k] = (int *)calloc(nbMax, sizeof(int));
        for (k1 = 0; k1 < nbMax; k1++)
            MatriceInterIndice[k][k1] = -1;
    }

    // Initialization
    for (k = 0; k < nbMax; k++) {
        OutForward[k] = 1e5;
        OutBackward[k] = 1e5;
        OutForwardIndice[k] = -1;
        OutBackwardIndice[k] = -1;
    }

    // Rotation VtoC
    for (t = 0; t < code->rowDegree[node]; t++) {
        for (k = 0; k < nbMax; k++) {
            if (decoder->M_VtoC_GF[t][k] != -1) {
                c = decoder->M_VtoC_GF[t][k];
                c = table->MULDEC[c][code->matValue[node][t]];
                decoder->M_VtoC_GF[t][k] = c;
            }
            // else printf("M_VtoC_GF[%d][%d]=%d\n",t,k,decoder->M_VtoC_GF[t][k]); // Commented out
            // noise
        }
    }

    // Initialisation algorithm
    for (k = 0; k < nbMax; k++) {
        OutForward[k] = decoder->M_VtoC_LLR[0][k];
        OutForwardIndice[k] = decoder->M_VtoC_GF[0][k];
        OutBackward[k] = decoder->M_VtoC_LLR[code->rowDegree[node] - 1][k];
        OutBackwardIndice[k] = decoder->M_VtoC_GF[code->rowDegree[node] - 1][k];
    }

    // Start of the recusivity
    for (kk = 1; kk < (code->rowDegree[node] - 1); kk++) {
        for (k = 0; k < nbMax; k++) {
            // forward step
            OutForward1[k] = decoder->M_VtoC_LLR[kk][k];
            OutForwardIndice1[k] = decoder->M_VtoC_GF[kk][k];
            // backward step
            OutBackward1[k] = decoder->M_VtoC_LLR[code->rowDegree[node] - kk - 1][k];
            OutBackwardIndice1[k] = decoder->M_VtoC_GF[code->rowDegree[node] - kk - 1][k];
        }

        for (k = 0; k < nbMax; k++) {
            // Storage of the intermediate vectors
            MatriceInter[kk - 1][k] = OutForward[k];
            MatriceInterIndice[kk - 1][k] = OutForwardIndice[k];
            MatriceInter[2 * (code->rowDegree[node] - 2) - kk][k] = OutBackward[k];
            MatriceInterIndice[2 * (code->rowDegree[node] - 2) - kk][k] = OutBackwardIndice[k];
        }

        if (kk < (code->rowDegree[node] - 1))
            ElementaryStep(OutForward, OutForward1, OutForwardIndice, OutForwardIndice1, OutForward,
                           OutForwardIndice, table->ADDGF, code->GF, nbMax, NbOper);

        if (kk < (code->rowDegree[node] - 1))
            ElementaryStep(OutBackward, OutBackward1, OutBackwardIndice, OutBackwardIndice1,
                           OutBackward, OutBackwardIndice, table->ADDGF, code->GF, nbMax, NbOper);
    }

    // Update of vectors CtoV (first and last)
    for (k = 0; k < nbMax; k++) {
        decoder->M_CtoV_LLR[code->rowDegree[node] - 1][k] = OutForward[k];
        decoder->M_CtoV_GF[code->rowDegree[node] - 1][k] = OutForwardIndice[k];
        decoder->M_CtoV_LLR[0][k] = OutBackward[k];
        decoder->M_CtoV_GF[0][k] = OutBackwardIndice[k];
    }

    // Update of the others vectors CtoV
    for (k = 0; k < (code->rowDegree[node] - 2); k++) {
        ElementaryStep(MatriceInter[k], MatriceInter[(code->rowDegree[node] - 2) + k],
                       MatriceInterIndice[k], MatriceInterIndice[(code->rowDegree[node] - 2) + k],
                       OutForward, OutForwardIndice, table->ADDGF, code->GF, nbMax, NbOper);
        for (g = 0; g < nbMax; g++) {
            decoder->M_CtoV_LLR[k + 1][g] = OutForward[g];
            decoder->M_CtoV_GF[k + 1][g] = OutForwardIndice[g];
        }
    }

    // Rotation CtoV
    for (t = 0; t < code->rowDegree[node]; t++) {
        for (k = 0; k < nbMax; k++) {
            if (decoder->M_CtoV_GF[t][k] == -1) {
                Stp = k;
                break;
            } else
                Stp = nbMax;
        }

        for (k = 0; k < Stp; k++) {
            decoder->M_CtoV_GF[t][k] =
                table->DIVDEC[decoder->M_CtoV_GF[t][k]][code->matValue[node][t]];
        }

        for (k = 0; k < code->GF; k++)
            LLR_tmp[k] = decoder->M_CtoV_LLR[t][Stp - 1] + offset;
        for (k = 0; k < Stp; k++)
            LLR_tmp[decoder->M_CtoV_GF[t][k]] = decoder->M_CtoV_LLR[t][k];
        for (k = 0; k < code->GF; k++) {
            decoder->M_CtoV_GF[t][k] = k;
            decoder->M_CtoV_LLR[t][k] = LLR_tmp[k];
        }
    }

    for (k = 0; k < S; k++)
        free(MatriceInterIndice[k]);
    free(MatriceInterIndice);
    for (k = 0; k < S; k++)
        free(MatriceInter[k]);
    free(MatriceInter);
    free(OutForward);
    free(OutForwardIndice);
    free(OutBackward);
    free(OutBackwardIndice);
    free(OutForward1);
    free(OutForwardIndice1);
    free(OutBackward1);
    free(OutBackwardIndice1);
}

/**
 * \fn CheckPassLogEMS_dc3
 */
void CheckPassLogEMS_dc3(int node, decoder_t *decoder, code_t *code, table_t *table, int NbOper,
                         float offset) {
    const int n_vc = decoder->n_vc;
    const int n_cv = decoder->n_cv;
    int t, c, k, Stp;
    float LLR_tmp[code->GF];

    for (t = 0; t < code->rowDegree[node]; t++) {
        for (k = 0; k < n_vc; k++) {
            if (decoder->M_VtoC_GF[t][k] != -1) {
                c = decoder->M_VtoC_GF[t][k];
                c = table->MULDEC[c][code->matValue[node][t]];
                decoder->M_VtoC_GF[t][k] = c;
            }
        }
    }

    ElementaryStep_nm(decoder->M_VtoC_LLR[0], decoder->M_VtoC_LLR[1], decoder->M_VtoC_GF[0],
                      decoder->M_VtoC_GF[1], decoder->M_CtoV_LLR[2], decoder->M_CtoV_GF[2],
                      table->ADDGF, code->GF, n_vc, n_vc, n_cv, NbOper); // forward
    ElementaryStep_nm(decoder->M_VtoC_LLR[1], decoder->M_VtoC_LLR[2], decoder->M_VtoC_GF[1],
                      decoder->M_VtoC_GF[2], decoder->M_CtoV_LLR[0], decoder->M_CtoV_GF[0],
                      table->ADDGF, code->GF, n_vc, n_vc, n_cv, NbOper); // backward
    ElementaryStep_nm(decoder->M_VtoC_LLR[0], decoder->M_VtoC_LLR[2], decoder->M_VtoC_GF[0],
                      decoder->M_VtoC_GF[2], decoder->M_CtoV_LLR[1], decoder->M_CtoV_GF[1],
                      table->ADDGF, code->GF, n_vc, n_vc, n_cv, NbOper); // merge

    for (t = 0; t < code->rowDegree[node]; t++) {
        for (k = 0; k < n_cv; k++) {
            if (decoder->M_CtoV_GF[t][k] == -1) {
                Stp = k;
                break;
            } else
                Stp = n_cv;
        }

        for (k = 0; k < Stp; k++) {
            decoder->M_CtoV_GF[t][k] =
                table->DIVDEC[decoder->M_CtoV_GF[t][k]][code->matValue[node][t]];
        }

        for (k = 0; k < code->GF; k++)
            LLR_tmp[k] = decoder->M_CtoV_LLR[t][Stp - 1] + offset;
        for (k = 0; k < Stp; k++)
            LLR_tmp[decoder->M_CtoV_GF[t][k]] = decoder->M_CtoV_LLR[t][k];
        for (k = 0; k < code->GF; k++) {
            decoder->M_CtoV_GF[t][k] = k;
            decoder->M_CtoV_LLR[t][k] = LLR_tmp[k];
        }
    }
}

/**
 * \fn ElementaryStep
 * \brief [Modified] Removed printf noise for -1 inputs
 */
int ElementaryStep(float *Input1, float *Input2, int *IndiceInput1, int *IndiceInput2,
                   float *Output, int *IndiceOut, int **ADDGF, int GF, int nbMax, int nbOper) {
    int nmU = nbMax;
    int nmV = nbMax;
    int nmS = nbMax;
    float loc_Output[nbMax];
    int loc_IndiceOut[nbMax];
    int i, j, s, ss, Indice_aux = 0;
    int nb_bubble = 8;
    int pos;
    float tab_aux[nbMax][nbMax];
    float tab_comp[3][nb_bubble];
    int GFvalues_already_in_Out[GF];

    // init
    for (i = 0; i < GF; i++)
        GFvalues_already_in_Out[i] = -1;

    for (i = 0; i < nbMax; i++) {
        loc_Output[i] = 1e5;
        loc_IndiceOut[i] = -1;
    }
    for (i = 0; i < nmU; i++) {
        for (j = 0; j < nmV; j++)
            tab_aux[i][j] = 0;
    }

    // pre-compute the addition of bubble check
    // horizontal
    for (j = 0; j < nb_bubble >> 1; j++) {
        for (i = 0; i < nmU; i++)
            tab_aux[i][j] = Input1[i] + Input2[j];
    }
    // vertical
    for (i = 0; i < nb_bubble >> 1; i++) {
        for (j = nb_bubble >> 1; j < nmV; j++)
            tab_aux[i][j] = Input1[i] + Input2[j];
    }

    // tab_comp init
    for (j = 0; j < nb_bubble / 2; j++) {
        tab_comp[0][j] = tab_aux[j][0];
        tab_comp[1][j] = j;
        tab_comp[2][j] = 0;
    }
    for (j = 0; j < nb_bubble / 2; j++) {
        tab_comp[0][j + nb_bubble / 2] = tab_aux[nb_bubble / 2][j];
        tab_comp[1][j + nb_bubble / 2] = nb_bubble / 2;
        tab_comp[2][j + nb_bubble / 2] = j;
    }

    // filling Out and IndiceOut
    s = 0;
    for (ss = 0; ss < nbOper; ss++) {
        pos = minimum(tab_comp[0], nb_bubble);

        // [修复核心]
        // 遇到 -1 表示已经没有更多有效的组合了。
        // 原代码在这里 printf("out of bounds") 其实是误报，这里我们直接安静地 break 即可。
        if ((IndiceInput1[(int)(tab_comp[1][pos])] == -1) ||
            (IndiceInput2[(int)(tab_comp[2][pos])] == -1)) {
            // break 表示停止填充，剩下的 output 保持为 -1 (Init状态)，这是正常的
            break;
        }

        Indice_aux = IndiceInput1[(int)(tab_comp[1][pos])] ^ IndiceInput2[(int)(tab_comp[2][pos])];

        if (GFvalues_already_in_Out[Indice_aux] == -1) {
            loc_Output[s] = tab_comp[0][pos];
            loc_IndiceOut[s] = Indice_aux;
            GFvalues_already_in_Out[Indice_aux] = 1;
            s++;
        }

        if (s == nmS)
            break;

        // control limits of tab_aux
        if ((tab_comp[1][pos] >= nmU - 1) || (tab_comp[2][pos] >= nmV - 1)) {
            break;
        }

        // update tab_comp with next value from tab_aux
        if (pos > nb_bubble / 2 - 1) {
            tab_comp[1][pos] = tab_comp[1][pos] + 1;
        } else {
            tab_comp[2][pos] = tab_comp[2][pos] + 1;
        }
        tab_comp[0][pos] = tab_aux[(int)tab_comp[1][pos]][(int)tab_comp[2][pos]];
    }

    for (i = 0; i < nbMax; i++) {
        Output[i] = loc_Output[i];
        IndiceOut[i] = loc_IndiceOut[i];
    }

    return (0);
}

/**
 * \fn ElementaryStep_nm
 * \brief [Modified] Removed printf noise for -1 inputs
 */
int ElementaryStep_nm(float *Input1, float *Input2, int *IndiceInput1, int *IndiceInput2,
                      float *Output, int *IndiceOut, int **ADDGF, int GF, int nmU, int nmV, int nmS,
                      int nbOper) {
    float loc_Output[nmS];
    int loc_IndiceOut[nmS];
    int i, j, s, ss, Indice_aux = 0;
    int nb_bubble = 8;
    int pos;
    float tab_aux[nmU + 1][nmV + 1];
    float tab_comp[3][nb_bubble];
    int GFvalues_already_Out[GF];

    // init
    for (i = 0; i < GF; i++)
        GFvalues_already_Out[i] = -1;

    for (i = 0; i < nmS; i++) {
        loc_Output[i] = 1e5;
        loc_IndiceOut[i] = -1;
    }
    for (i = 0; i < nmU + 1; i++) {
        for (j = 0; j < nmV + 1; j++)
            tab_aux[i][j] = 0;
    }

    // pre-compute
    for (j = 0; j < nb_bubble >> 1; j++) {
        for (i = 0; i < nmU; i++)
            tab_aux[i][j] = Input1[i] + Input2[j];
        tab_aux[i][j] = 100;
    }
    for (i = 0; i < nb_bubble >> 1; i++) {
        for (j = nb_bubble >> 1; j < nmV; j++)
            tab_aux[i][j] = Input1[i] + Input2[j];
        tab_aux[i][j] = 100;
    }

    // tab_comp init
    for (j = 0; j < nb_bubble / 2; j++) {
        tab_comp[0][j] = tab_aux[j][0];
        tab_comp[1][j] = j;
        tab_comp[2][j] = 0;
    }
    for (j = 0; j < nb_bubble / 2; j++) {
        tab_comp[0][j + nb_bubble / 2] = tab_aux[nb_bubble / 2][j];
        tab_comp[1][j + nb_bubble / 2] = nb_bubble / 2;
        tab_comp[2][j + nb_bubble / 2] = j;
    }

    // filling Out
    s = 0;
    for (ss = 0; ss < nbOper; ss++) {
        pos = minimum(tab_comp[0], nb_bubble);

        // [修复核心]
        // 同样，这里遇到 -1 直接 break，不再报错。
        if ((IndiceInput1[(int)(tab_comp[1][pos])] == -1) ||
            (IndiceInput2[(int)(tab_comp[2][pos])] == -1)) {
            break;
        }

        Indice_aux = IndiceInput1[(int)(tab_comp[1][pos])] ^ IndiceInput2[(int)(tab_comp[2][pos])];

        if (GFvalues_already_Out[Indice_aux] == -1) {
            loc_Output[s] = tab_comp[0][pos];
            loc_IndiceOut[s] = Indice_aux;
            GFvalues_already_Out[Indice_aux] = 1;
            s++;
        }

        if (s == nmS)
            break;

        if (pos > nb_bubble / 2 - 1) {
            tab_comp[1][pos] = tab_comp[1][pos] + 1;
        } else {
            tab_comp[2][pos] = tab_comp[2][pos] + 1;
        }
        tab_comp[0][pos] = tab_aux[(int)tab_comp[1][pos]][(int)tab_comp[2][pos]];
    }

    for (i = 0; i < nmS; i++) {
        Output[i] = loc_Output[i];
        IndiceOut[i] = loc_IndiceOut[i];
    }

    return (0);
}

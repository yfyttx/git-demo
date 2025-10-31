#include "ldpc_codes/binary_codes.h"

using namespace Eigen;

// --- 辅助函数：计算 a 在模 m 下的乘法逆元 ---
// 使用扩展欧几里得算法，只对 m 是素数时保证正确
long long power(long long base, long long exp, long long mod) {
    long long res = 1;
    base %= mod;
    while (exp > 0) {
        if (exp % 2 == 1) res = (res * base) % mod;
        base = (base * base) % mod;
        exp /= 2;
    }
    return res;
}

long long modInverse(long long n, long long mod) {
    return power(n, mod - 2, mod);
}


// --- 能够正确处理负指数的 powMod 函数 ---
int powMod(int base, int exp, int modulus) {
    if (exp == 0) return 1;
    if (exp > 0) {
        return (int)power(base, exp, modulus);
    }
    // exp < 0
    long long pos_exp = -exp;
    long long base_pow = power(base, pos_exp, modulus);
    return (int)modInverse(base_pow, modulus);
}


// --- 辅助函数：创建循环置换矩阵 ---
MatrixXi createCircularShift(int P, int k) {
    MatrixXi mat = MatrixXi::Zero(P, P);
    k = (k % P + P) % P; // 确保 k 在 [0, P-1]
    
    for (int i = 0; i < P; ++i) {
        mat(i, (i + k) % P) = 1;
    }
    return mat;
}


// --- 主函数：构造二进制矩阵 ---
std::pair<MatrixXi, MatrixXi> constructBinaryHMatrices(int J, int L, int P, int sigma, int tau) {
    int rows = J * P;
    int cols = L * P;
    MatrixXi Hc = MatrixXi::Zero(rows, cols);
    MatrixXi Hd = MatrixXi::Zero(rows, cols);

    for (int j = 0; j < J; ++j) {
        for (int l = 0; l < L; ++l) {
            int c_shift, d_shift;

            if (l < L / 2) {
                c_shift = powMod(sigma, -j + l, P);
            } else {
                long long temp = powMod(sigma, -j + l, P);
                c_shift = (int)(( (long long)tau * temp) % P);
            }

            if (l < L / 2) {
                long long temp = powMod(sigma, j - l, P);
                long long val = ( (long long)tau * temp) % P;
                d_shift = (P - val) % P;
            } else {
                long long temp = powMod(sigma, j - l, P);
                d_shift = (P - temp) % P;
            }

            Hc.block(j * P, l * P, P, P) = createCircularShift(P, c_shift);
            Hd.block(j * P, l * P, P, P) = createCircularShift(P, d_shift);
        }
    }

    return {Hc, Hd};
}
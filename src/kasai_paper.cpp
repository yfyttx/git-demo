// Reproduce the construction and separate component decoder in Kasai et al.,
// IEEE Transactions on Information Theory 58(2), 1223-1230 (2012), Sec. II-IV.
#include <Eigen/Sparse>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

#include "NB_BP_Decoder.h"
#include "ldpc_codes/binary_codes.h"
#include "ldpc_codes/nonbinary_codes.h"

using Matrix = Eigen::SparseMatrix<int, Eigen::RowMajor>;

static std::vector<int> polynomial(int degree) {
    // Irreducible polynomials; the GF class selects a primitive root itself.
    int bits = 0;
    switch (degree) {
        case 6: bits = 0x43; break;   // x^6+x+1
        case 8: bits = 0x11d; break;  // x^8+x^4+x^3+x^2+1
        case 9: bits = 0x211; break;  // x^9+x^4+1
        case 10: bits = 0x409; break; // x^10+x^3+1
        default: throw std::invalid_argument("Supported field degrees: 6, 8, 9, 10");
    }
    std::vector<int> coefficients(degree + 1);
    for (int i = 0; i <= degree; ++i)
        coefficients[i] = (bits >> (degree - i)) & 1;
    return coefficients;
}

static std::vector<int> syndrome(const Matrix& check, const std::vector<int>& error,
                                 const GF2p& gf) {
    std::vector<int> result(check.rows(), 0);
    for (int row = 0; row < check.rows(); ++row)
        for (Matrix::InnerIterator edge(check, row); edge; ++edge)
            result[row] = gf.add(result[row], gf.mul(edge.value(), error[edge.col()]));
    return result;
}

static Matrix binary_expansion(const Matrix& source, int degree, const GF2p& gf,
                               bool transpose_blocks) {
    std::vector<Eigen::Triplet<int>> entries;
    for (int row = 0; row < source.rows(); ++row) {
        for (Matrix::InnerIterator edge(source, row); edge; ++edge) {
            for (int input = 0; input < degree; ++input) {
                const int product = gf.mul(edge.value(), 1 << input);
                for (int output = 0; output < degree; ++output) {
                    if ((product >> output) & 1)
                        entries.emplace_back(row * degree + (transpose_blocks ? input : output),
                                             edge.col() * degree + (transpose_blocks ? output : input), 1);
                }
            }
        }
    }
    Matrix expanded(source.rows() * degree, source.cols() * degree);
    expanded.setFromTriplets(entries.begin(), entries.end());
    return expanded;
}

static void save_market(const Matrix& matrix, const std::string& filename) {
    std::ofstream output(filename);
    if (!output) throw std::runtime_error("Cannot write " + filename);
    output << "%%MatrixMarket matrix coordinate integer general\n";
    output << matrix.rows() << ' ' << matrix.cols() << ' ' << matrix.nonZeros() << '\n';
    for (int row = 0; row < matrix.rows(); ++row)
        for (Matrix::InnerIterator entry(matrix, row); entry; ++entry)
            output << row + 1 << ' ' << entry.col() + 1 << ' ' << entry.value() << '\n';
}

static void verify_pair(const Matrix& gamma, const Matrix& delta, const GF2p& gf) {
    if (gamma.cols() != delta.cols()) throw std::runtime_error("Different matrix widths");
    Matrix delta_rows = delta;
    for (int r = 0; r < gamma.rows(); ++r) {
        for (int s = 0; s < delta.rows(); ++s) {
            int product = 0;
            for (Matrix::InnerIterator a(gamma, r); a; ++a) {
                const int b = delta_rows.coeff(s, a.col());
                if (b) product = gf.add(product, gf.mul(a.value(), b));
            }
            if (product) throw std::runtime_error("Gamma Delta^T is nonzero");
        }
    }
}

static int field_rank(const Matrix& matrix, const GF2p& gf) {
    Eigen::MatrixXi work = Eigen::MatrixXi::Zero(matrix.rows(), matrix.cols());
    for (int row = 0; row < matrix.rows(); ++row)
        for (Matrix::InnerIterator edge(matrix, row); edge; ++edge)
            work(row, edge.col()) = edge.value();
    int rank = 0;
    for (int col = 0; col < matrix.cols() && rank < matrix.rows(); ++col) {
        int pivot = rank;
        while (pivot < matrix.rows() && work(pivot, col) == 0) ++pivot;
        if (pivot == matrix.rows()) continue;
        work.row(rank).swap(work.row(pivot));
        const int inverse = gf.inv(work(rank, col));
        for (int j = col; j < matrix.cols(); ++j)
            work(rank, j) = gf.mul(work(rank, j), inverse);
        for (int i = rank + 1; i < matrix.rows(); ++i) {
            const int factor = work(i, col);
            if (!factor) continue;
            for (int j = col; j < matrix.cols(); ++j)
                work(i, j) = gf.add(work(i, j), gf.mul(factor, work(rank, j)));
        }
        ++rank;
    }
    return rank;
}

int main(int argc, char** argv) {
    try {
        if (argc > 10) throw std::invalid_argument(
            "Usage: kasai_paper [trials] [seed] [f_m] [degree] [P] [sigma] [tau] [L] [output-prefix]");
        const int trials = argc > 1 ? std::stoi(argv[1]) : 20;
        const unsigned seed = argc > 2 ? std::stoul(argv[2]) : 20261018;
        const double fm = argc > 3 ? std::stod(argv[3]) : 0.02;
        const int degree = argc > 4 ? std::stoi(argv[4]) : 8;
        const int P = argc > 5 ? std::stoi(argv[5]) : 31;
        const int sigma = argc > 6 ? std::stoi(argv[6]) : 5;
        const int tau = argc > 7 ? std::stoi(argv[7]) : 2;
        const int L = argc > 8 ? std::stoi(argv[8]) : 6;
        constexpr int J = 2;
        if (trials <= 0 || fm < 0 || fm > 2.0 / 3.0)
            throw std::invalid_argument("Invalid trial count or marginal flip probability");
        GF2p gf(degree, polynomial(degree));
        std::mt19937 construction_rng(seed);
        auto supports = constructBinaryHMatrices(J, L, P, sigma, tau);
        auto generated = constructKasaiNonbinaryPair(
            supports.first, supports.second, gf, construction_rng);
        Matrix gamma = generated.first, delta = generated.second;
        verify_pair(gamma, delta, gf);
        const Matrix hc = binary_expansion(gamma, degree, gf, false);
        const Matrix hd = binary_expansion(delta, degree, gf, true);
        const Matrix product = hc * hd.transpose();
        for (int row = 0; row < product.rows(); ++row)
            for (Matrix::InnerIterator entry(product, row); entry; ++entry)
                if (entry.value() & 1)
                    throw std::runtime_error("Binary H_C H_D^T is nonzero over GF(2)");
        if (argc > 9) {
            const std::string prefix = argv[9];
            save_market(gamma, prefix + "_Gamma.mtx");
            save_market(delta, prefix + "_Delta.mtx");
            save_market(hc, prefix + "_HC.mtx");
            save_market(hd, prefix + "_HD.mtx");
        }
        const int gamma_rank = field_rank(gamma, gf);
        const int delta_rank = field_rank(delta, gf);

        NB_BP_Decoder decoder(gf, 50);
        std::mt19937 channel_rng(seed ^ 0x9e3779b9u);
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        const double pd = 1.5 * fm;
        const int N = gamma.cols();
        const FieldBasis basis(gf);
        int gamma_block_errors = 0, delta_block_errors = 0;
        int gamma_syndrome_failures = 0, delta_syndrome_failures = 0;
        int physical_errors = 0;
        for (int trial = 0; trial < trials; ++trial) {
            std::vector<int> z(N, 0), x(N, 0);
            for (int symbol = 0; symbol < N; ++symbol) {
                int dual_x_bits = 0;
                for (int bit = 0; bit < degree; ++bit) {
                    const double draw = uniform(channel_rng);
                    const bool x_bit = draw < pd / 3 || (draw >= 2 * pd / 3 && draw < pd);
                    const bool z_bit = draw >= pd / 3 && draw < pd;
                    dual_x_bits |= static_cast<int>(x_bit) << bit;
                    z[symbol] |= static_cast<int>(z_bit ||
                        (draw >= 2 * pd / 3 && draw < pd)) << bit;
                    physical_errors += draw < pd;
                }
                x[symbol] = basis.fromDual(dual_x_bits);
            }
            // H_C=A(Gamma), H_D=A(Delta)^T. Their syndrome symbols are
            // Gamma*z (polynomial basis) and Delta*x (trace-dual basis).
            const auto sz = syndrome(gamma, z, gf);
            const auto sx = syndrome(delta, x, gf);
            const auto zhat = decoder.decode_component(gamma, sz, fm, false);
            const auto xhat = decoder.decode_component(delta, sx, fm, true);
            gamma_syndrome_failures += syndrome(gamma, zhat, gf) != sz;
            delta_syndrome_failures += syndrome(delta, xhat, gf) != sx;
            gamma_block_errors += zhat != z;
            delta_block_errors += xhat != x;
            if (trials >= 1000 && (trial + 1) % 500 == 0)
                std::cerr << "completed=" << trial + 1 << '/' << trials << '\n';
        }
        std::cout << std::setprecision(8)
                  << "Kasai2012 J=" << J << " L=" << L << " P=" << P
                  << " sigma=" << sigma << " tau=" << tau
                  << " q=" << gf.getSize() << " N_qubits=" << N * degree
                  << " binary_checks=" << hc.rows() << 'x' << hc.cols()
                  << " nominal_RQ=" << 1.0 - 2.0 * J / L
                  << " actual_RQ=" << 1.0 - double(gamma_rank + delta_rank) / N
                  << " ranks=" << gamma_rank << ',' << delta_rank
                  << " seed=" << seed << " max_iters=50\n"
                  << "f_m=" << fm << " p_D=" << pd << " trials=" << trials
                  << " mean_physical_errors=" << double(physical_errors) / trials << '\n'
                  << "Gamma_block_error=" << double(gamma_block_errors) / trials
                  << " Delta_block_error=" << double(delta_block_errors) / trials
                  << " Gamma_block_errors=" << gamma_block_errors
                  << " Delta_block_errors=" << delta_block_errors
                  << " Gamma_syndrome_fail=" << double(gamma_syndrome_failures) / trials
                  << " Delta_syndrome_fail=" << double(delta_syndrome_failures) / trials
                  << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Kasai reproduction error: " << error.what() << '\n';
        return 1;
    }
}

#include <Eigen/Sparse>

#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "NB_BP_Decoder.h"
#include "gf2p/gf2p.h"

using namespace Eigen;
using namespace std;

SparseMatrix<int, RowMajor> loadMatrix(const string& filename, GF2p& gf) {
    ifstream in(filename);
    if (!in) {
        throw runtime_error("Could not open matrix file: " + filename);
    }

    int cols = 0, rows = 0, field_size = 0;
    if (!(in >> cols >> rows >> field_size) || cols <= 0 || rows <= 0 || field_size != gf.getSize()) {
        throw runtime_error("Invalid GF matrix header in " + filename);
    }

    string line;
    getline(in, line);  // finish header line
    getline(in, line);  // variable degrees
    getline(in, line);  // check degrees

    vector<Triplet<int>> triplets;
    for (int row = 0; row < rows; ++row) {
        if (!getline(in, line)) {
            throw runtime_error("Unexpected end of matrix file: " + filename);
        }
        istringstream row_stream(line);
        int col_index = 0, log_value = 0;
        while (row_stream >> col_index >> log_value) {
            if (col_index < 1 || col_index > cols || log_value < 0 || log_value >= gf.getMod()) {
                throw runtime_error("Invalid entry in matrix file: " + filename);
            }
            triplets.emplace_back(row, col_index - 1, gf.exp(log_value));
        }
        if (!row_stream.eof()) {
            throw runtime_error("Malformed row in matrix file: " + filename);
        }
    }

    SparseMatrix<int, RowMajor> matrix(rows, cols);
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    return matrix;
}

SparseMatrix<int, RowMajor> loadBinaryMatrix(const string& filename) {
    ifstream in(filename);
    if (!in) {
        throw runtime_error("Could not open matrix file: " + filename);
    }

    int cols = 0, rows = 0, field_size = 0;
    if (!(in >> cols >> rows >> field_size) || cols <= 0 || rows <= 0 || field_size != 2) {
        throw runtime_error("Invalid binary matrix header in " + filename);
    }

    string line;
    getline(in, line);  // finish header line
    getline(in, line);  // optional variable degrees
    getline(in, line);  // optional check degrees

    vector<Triplet<int>> triplets;
    for (int row = 0; row < rows; ++row) {
        if (!getline(in, line)) {
            throw runtime_error("Unexpected end of matrix file: " + filename);
        }
        istringstream row_stream(line);
        int col_index = 0, value = 0;
        while (row_stream >> col_index >> value) {
            if (col_index < 1 || col_index > cols || (value != 0 && value != 1)) {
                throw runtime_error("Invalid entry in binary matrix file: " + filename);
            }
            if (value == 1) {
                triplets.emplace_back(row, col_index - 1, 1);
            }
        }
        if (!row_stream.eof()) {
            throw runtime_error("Malformed row in binary matrix file: " + filename);
        }
    }

    SparseMatrix<int, RowMajor> matrix(rows, cols);
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    return matrix;
}

vector<int> getBinarySyndrome(const SparseMatrix<int, RowMajor>& matrix,
                              const vector<int>& error_bits) {
    if (matrix.cols() != static_cast<int>(error_bits.size())) {
        throw invalid_argument("Binary error length does not match parity-check columns");
    }

    vector<int> syndrome(matrix.rows(), 0);
    for (int row = 0; row < matrix.rows(); ++row) {
        for (SparseMatrix<int, RowMajor>::InnerIterator it(matrix, row); it; ++it) {
            syndrome[row] ^= error_bits[it.col()];
        }
    }
    return syndrome;
}

void sampleDepolarizing(vector<int>& error_x, vector<int>& error_z, double p) {
    if (p < 0.0 || p > 1.0 || error_x.size() != error_z.size()) {
        throw invalid_argument("Invalid depolarizing channel parameters or vector dimensions");
    }

    const double component_probability = p / 3.0;
    for (size_t bit = 0; bit < error_x.size(); ++bit) {
        const double sample = static_cast<double>(rand()) / (RAND_MAX + 1.0);
        if (sample < component_probability) error_x[bit] = 1;
        else if (sample < 2.0 * component_probability) error_z[bit] = 1;
        else if (sample < p) error_x[bit] = error_z[bit] = 1;
    }
}

vector<int> binaryToFieldSymbols(const vector<int>& bits, int field_degree,
                                 const FieldBasis& basis, bool trace_dual) {
    if (field_degree <= 0 || bits.size() % static_cast<size_t>(field_degree) != 0) {
        throw invalid_argument("Binary error vector is not a whole number of GF symbols");
    }
    vector<int> symbols(bits.size() / field_degree, 0);
    for (size_t symbol = 0; symbol < symbols.size(); ++symbol) {
        for (int bit = 0; bit < field_degree; ++bit) {
            symbols[symbol] |= (bits[symbol * field_degree + bit] & 1) << bit;
        }
        if (trace_dual) symbols[symbol] = basis.fromDual(symbols[symbol]);
    }
    return symbols;
}

void validateBinaryExpansion(const SparseMatrix<int, RowMajor>& logical,
                             const SparseMatrix<int, RowMajor>& physical,
                             const GF2p& gf, const FieldBasis& basis,
                             bool transposed_blocks) {
    const int degree = gf.getPower();
    if (physical.rows() != logical.rows() * degree ||
        physical.cols() != logical.cols() * degree) {
        throw runtime_error("Binary and GF matrices have incompatible dimensions");
    }
    for (int row = 0; row < logical.rows(); ++row) {
        for (int output_bit = 0; output_bit < degree; ++output_bit) {
            vector<int> expected(physical.cols(), 0);
            for (SparseMatrix<int, RowMajor>::InnerIterator entry(logical, row); entry; ++entry) {
                for (int input_bit = 0; input_bit < degree; ++input_bit) {
                    const int input = 1 << input_bit;
                    const int output = transposed_blocks
                        ? basis.toDual(gf.mul(entry.value(), basis.fromDual(input)))
                        : gf.mul(entry.value(), input);
                    expected[entry.col() * degree + input_bit] = (output >> output_bit) & 1;
                }
            }
            for (SparseMatrix<int, RowMajor>::InnerIterator entry(physical,
                    row * degree + output_bit); entry; ++entry) {
                if (expected[entry.col()] != 1 || entry.value() != 1) {
                    throw runtime_error("Binary matrix does not match its GF matrix");
                }
                expected[entry.col()] = 0;
            }
            for (int bit : expected) {
                if (bit != 0) throw runtime_error("Binary matrix does not match its GF matrix");
            }
        }
    }
}

bool isInStabilizerRowSpace(const SparseMatrix<int, RowMajor>& checks,
                            const vector<int>& residual, const GF2p& gf) {
    if (checks.cols() != static_cast<int>(residual.size())) {
        throw invalid_argument("Residual error length does not match stabilizer width");
    }

    // Solve lambda^T * checks = residual over GF(q).
    MatrixXi augmented = MatrixXi::Zero(checks.cols(), checks.rows() + 1);
    for (int row = 0; row < checks.rows(); ++row) {
        for (SparseMatrix<int, RowMajor>::InnerIterator it(checks, row); it; ++it) {
            augmented(it.col(), row) = it.value();
        }
    }
    for (int col = 0; col < checks.cols(); ++col) {
        augmented(col, checks.rows()) = residual[col];
    }

    int pivot_row = 0;
    for (int variable = 0; variable < checks.rows() && pivot_row < checks.cols(); ++variable) {
        int selected = pivot_row;
        while (selected < checks.cols() && augmented(selected, variable) == 0) ++selected;
        if (selected == checks.cols()) continue;
        augmented.row(pivot_row).swap(augmented.row(selected));

        const int inverse = gf.inv(augmented(pivot_row, variable));
        for (int col = variable; col <= checks.rows(); ++col) {
            augmented(pivot_row, col) = gf.mul(augmented(pivot_row, col), inverse);
        }
        for (int row = 0; row < checks.cols(); ++row) {
            if (row == pivot_row || augmented(row, variable) == 0) continue;
            const int factor = augmented(row, variable);
            for (int col = variable; col <= checks.rows(); ++col) {
                augmented(row, col) = gf.add(
                    augmented(row, col), gf.mul(factor, augmented(pivot_row, col)));
            }
        }
        ++pivot_row;
    }

    for (int row = pivot_row; row < checks.cols(); ++row) {
        if (augmented(row, checks.rows()) != 0) return false;
    }
    return true;
}

int main(int argc, char* argv[]) {
    try {
        if (argc > 5) {
            throw invalid_argument(
                "Usage: simulate_qec [trials-per-error-rate] [seed] [error-rate] "
                "[direct|cascade|video]");
        }
        const int trials = argc >= 2 ? stoi(argv[1]) : 100;
        if (trials <= 0) {
            throw invalid_argument("The trial count must be positive");
        }

        const unsigned int seed = argc >= 3 ? static_cast<unsigned int>(stoul(argv[2]))
                                            : static_cast<unsigned int>(time(nullptr));
        const string construction = argc >= 5 ? argv[4] : "direct";
        if (construction != "direct" && construction != "cascade" && construction != "video") {
            throw invalid_argument("Construction must be 'direct', 'cascade', or 'video'");
        }
        srand(seed);
        cerr << "seed=" << seed << " construction=" << construction << '\n';
        GF2p gf(6, {1, 0, 0, 0, 0, 1, 1});
        const int field_degree = 6;
        NB_BP_Decoder decoder(gf, 50);
        FieldBasis basis(gf);

        string direct_dir = "build/";
        if (construction == "direct") {
            const vector<string> files = {"HX_Safe_GF64.txt", "HZ_Safe_GF64.txt",
                                          "HX_Physical_Binary.txt", "HZ_Physical_Binary.txt"};
            for (const string& file : files) {
                if (!filesystem::is_regular_file(direct_dir + file)) {
                    direct_dir = "matrices/gf64_tit/";
                    break;
                }
            }
        }
        const string hx_prefix = construction == "direct" ? direct_dir + "HX" :
            (construction == "cascade" ? "build/HX_Cascade" : "build/HX_Video");
        const string hz_prefix = construction == "direct" ? direct_dir + "HZ" :
            (construction == "cascade" ? "build/HZ_Cascade" : "build/HZ_Video");
        const string suffix = construction == "direct" ? "_Safe_GF64.txt" : "_GF64.txt";
        const auto hx_logical = loadMatrix(hx_prefix + suffix, gf);
        const auto hz_logical = loadMatrix(hz_prefix + suffix, gf);
        SparseMatrix<int, RowMajor> cascade_outer_hz_check;
        SparseMatrix<int, RowMajor> cascade_outer_hx_check;
        if (construction == "cascade") {
            cascade_outer_hz_check = loadMatrix("build/HX_TIT_Outer_For_HZ_GF64.txt", gf);
            cascade_outer_hx_check = loadMatrix("build/HZ_TIT_Outer_For_HX_GF64.txt", gf);
        }
        const auto hx_physical = loadBinaryMatrix(hx_prefix + "_Physical_Binary.txt");
        const auto hz_physical = loadBinaryMatrix(hz_prefix + "_Physical_Binary.txt");

        if (hx_logical.cols() != hz_logical.cols()) {
            throw runtime_error("Logical HX and HZ must have the same number of columns");
        }
        if (hx_physical.cols() != hx_logical.cols() * field_degree ||
            hz_physical.cols() != hz_logical.cols() * field_degree ||
            hx_physical.rows() != hx_logical.rows() * field_degree ||
            hz_physical.rows() != hz_logical.rows() * field_degree) {
            throw runtime_error("Logical and binary matrices have incompatible dimensions. "
                                "Regenerate all four matrices with the same constructor run.");
        }
        validateBinaryExpansion(hx_logical, hx_physical, gf, basis, false);
        validateBinaryExpansion(hz_logical, hz_physical, gf, basis, true);

        const int symbol_count = hx_logical.cols();
        const int physical_bits = symbol_count * field_degree;
        const vector<double> error_rates = argc >= 4 ? vector<double>{stod(argv[3])}
                                                   : vector<double>{0.001, 0.005, 0.01, 0.02, 0.04};
        cout << "p_D\tLogicalFailureRate\tTotal\tFailures\tSyndromeFailures"
                "\tXSyndromeFailures\tZSyndromeFailures"
                "\tMeanPhysicalErrors\tMeanResidualSymbols\tUndetectedLogicalFailures"
                "\tExactRecoveryFrames\tCascadeTrailingRecoveries"
                "\tCascadeStagedAttempts\tCascadeStagedValid"
                "\tZSoftFallbackAttempts\tZSoftFallbackSuccesses\n";

        for (double p : error_rates) {
            int failures = 0;
            int syndrome_failures = 0;
            int x_syndrome_failures = 0;
            int z_syndrome_failures = 0;
            int undetected_logical_failures = 0;
            int exact_recovery_frames = 0;
            int cascade_trailing_recoveries = 0;
            int cascade_staged_attempts = 0;
            int cascade_staged_successes = 0;
            int z_soft_fallback_attempts = 0;
            int z_soft_fallback_successes = 0;
            int physical_errors = 0;
            int residual_symbols = 0;
            for (int trial = 0; trial < trials; ++trial) {
                vector<int> error_x(physical_bits, 0), error_z(physical_bits, 0);
                sampleDepolarizing(error_x, error_z, p);
                for (int bit = 0; bit < physical_bits; ++bit) {
                    physical_errors += error_x[bit] || error_z[bit];
                }

                const vector<int> syndrome_x = getBinarySyndrome(hz_physical, error_x);
                const vector<int> syndrome_z = getBinarySyndrome(hx_physical, error_z);
                const vector<int> logical_syndrome_x =
                    binaryToFieldSymbols(syndrome_x, field_degree, basis, true);
                const vector<int> logical_syndrome_z =
                    binaryToFieldSymbols(syndrome_z, field_degree, basis, false);
                vector<int> corrected_x;
                vector<int> corrected_z;
                if (construction == "cascade") {
                    // Evaluate the requested staged decoder on H_Z: decode
                    // only the five upper-selector slices, then recover the
                    // sixth symbol per inner block from the lower H_2 check.
                    ++cascade_staged_attempts;
                    vector<int> staged_x;
                    if (decoder.decode_cascade_x(cascade_outer_hz_check, hz_logical, 6,
                            logical_syndrome_x, p, staged_x)) {
                        ++cascade_staged_successes;
                        ++cascade_trailing_recoveries;
                        corrected_x = std::move(staged_x);
                    } else {
                        corrected_x.assign(symbol_count, 0);
                    }
                    // Keep the other component on its ordinary H_X decoder;
                    // this isolates the effect of the requested H_Z strategy.
                    corrected_z = decoder.decode_component(
                        hx_logical, logical_syndrome_z, 2.0 * p / 3.0, false);
                } else {
                    const auto [estimate_x, estimate_z] = decoder.joint_decode(
                        hx_logical, hz_logical, syndrome_x, syndrome_z, p);
                    corrected_x = estimate_x;
                    corrected_z = estimate_z;
                }
                if (!decoder.check_z_syndrome(hx_logical, corrected_z, logical_syndrome_z)) {
                    ++z_soft_fallback_attempts;
                    const vector<int> soft_z = decoder.decode_component(
                        hx_logical, logical_syndrome_z, 2.0 * p / 3.0, false);
                    if (decoder.check_z_syndrome(hx_logical, soft_z, logical_syndrome_z)) {
                        corrected_z = soft_z;
                        ++z_soft_fallback_successes;
                    }
                }
                const bool x_syndrome_ok = decoder.check_x_syndrome(
                    hz_logical, corrected_x, logical_syndrome_x);
                const bool z_syndrome_ok = decoder.check_z_syndrome(
                    hx_logical, corrected_z, logical_syndrome_z);
                const bool syndrome_ok = x_syndrome_ok && z_syndrome_ok;
                if (!x_syndrome_ok) ++x_syndrome_failures;
                if (!z_syndrome_ok) ++z_syndrome_failures;
                if (!syndrome_ok) {
                    ++syndrome_failures;
                }

                const vector<int> actual_x = binaryToFieldSymbols(error_x, field_degree, basis, true);
                const vector<int> actual_z = binaryToFieldSymbols(error_z, field_degree, basis, false);
                vector<int> residual_x(symbol_count), residual_z(symbol_count);
                bool exact_recovery = true;
                for (int symbol = 0; symbol < symbol_count; ++symbol) {
                    residual_x[symbol] = gf.add(actual_x[symbol], corrected_x[symbol]);
                    residual_z[symbol] = gf.add(actual_z[symbol], corrected_z[symbol]);
                    residual_symbols += residual_x[symbol] != 0 || residual_z[symbol] != 0;
                    exact_recovery &= residual_x[symbol] == 0 && residual_z[symbol] == 0;
                }
                exact_recovery_frames += exact_recovery;

                const bool x_is_stabilizer = isInStabilizerRowSpace(hx_logical, residual_x, gf);
                const bool z_is_stabilizer = isInStabilizerRowSpace(hz_logical, residual_z, gf);
                if (!x_is_stabilizer || !z_is_stabilizer) {
                    ++failures;
                    if (syndrome_ok) ++undetected_logical_failures;
                }
            }
            cout << p << '\t' << static_cast<double>(failures) / trials << '\t'
                 << trials << '\t' << failures << '\t' << syndrome_failures << '\t'
                 << x_syndrome_failures << '\t' << z_syndrome_failures << '\t'
                 << static_cast<double>(physical_errors) / trials << '\t'
                 << static_cast<double>(residual_symbols) / trials << '\t'
                 << undetected_logical_failures << '\t' << exact_recovery_frames << '\t'
                 << cascade_trailing_recoveries << '\t' << cascade_staged_attempts << '\t'
                 << cascade_staged_successes << '\t' << z_soft_fallback_attempts << '\t'
                 << z_soft_fallback_successes << '\n';
        }
    } catch (const exception& error) {
        cerr << "Simulation error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}

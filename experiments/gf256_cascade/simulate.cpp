#include <Eigen/Sparse>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_set>
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

struct PauliChannel {
    double p_i = 1.0;
    double p_x = 0.0;
    double p_y = 0.0;
    double p_z = 0.0;

    static PauliChannel zBiased(double p, double bias) {
        if (!std::isfinite(p) || !std::isfinite(bias) ||
            p <= 0.0 || p >= 1.0 || bias <= 0.0)
            throw invalid_argument("Require 0 < p < 1 and finite Z bias > 0");
        const double low = p / (bias + 2.0);
        return {1.0 - p, low, low, bias * low};
    }

    double marginalX() const { return p_x + p_y; }
    double marginalZ() const { return p_z + p_y; }
    double zGivenX(bool x_bit) const {
        return x_bit ? p_y / (p_x + p_y) : p_z / (p_i + p_z);
    }
};

PauliChannel pauli_channel;

void samplePauli(vector<int>& error_x, vector<int>& error_z) {
    if (error_x.size() != error_z.size())
        throw invalid_argument("Pauli channel vector dimensions differ");
    for (size_t bit = 0; bit < error_x.size(); ++bit) {
        const double sample = static_cast<double>(rand()) / (RAND_MAX + 1.0);
        if (sample < pauli_channel.p_x) error_x[bit] = 1;
        else if (sample < pauli_channel.p_x + pauli_channel.p_z) error_z[bit] = 1;
        else if (sample < pauli_channel.p_x + pauli_channel.p_z + pauli_channel.p_y)
            error_x[bit] = error_z[bit] = 1;
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

// Row echelon basis for repeated GF(q) row-space membership queries.
class FieldRowSpace {
public:
    FieldRowSpace(const SparseMatrix<int, RowMajor>& matrix, GF2p& field,
                  const SparseMatrix<int, RowMajor>* inner_h2=nullptr,
                  const SparseMatrix<int, RowMajor>* inner_g2=nullptr)
        : width_(matrix.cols()), gf_(field) {
        if(inner_h2) inner_h2_=*inner_h2;
        if(inner_g2) inner_g2_=*inner_g2;
        vector<vector<int>> dense(matrix.rows(), vector<int>(width_, 0));
        for (int r=0; r<matrix.rows(); ++r)
            for (SparseMatrix<int, RowMajor>::InnerIterator it(matrix,r); it; ++it)
                dense[r][it.col()]=it.value();
        int rank=0;
        for (int col=0; col<width_ && rank<matrix.rows(); ++col) {
            int pivot=rank;
            while (pivot<matrix.rows() && dense[pivot][col]==0) ++pivot;
            if (pivot==matrix.rows()) continue;
            std::swap(dense[rank],dense[pivot]);
            const int inverse=gf_.inv(dense[rank][col]);
            for (int j=col; j<width_; ++j)
                dense[rank][j]=gf_.mul(dense[rank][j],inverse);
            for (int r=rank+1; r<matrix.rows(); ++r) {
                const int factor=dense[r][col];
                if (!factor) continue;
                for (int j=col; j<width_; ++j)
                    dense[r][j]=gf_.add(dense[r][j],gf_.mul(factor,dense[rank][j]));
            }
            pivots_.push_back(col);
            basis_.push_back(std::move(dense[rank]));
            ++rank;
        }
    }

    vector<int> reduce(vector<int> value) const {
        if (value.size()!=static_cast<size_t>(width_))
            throw invalid_argument("Row-space vector has incorrect length");
        for (size_t r=0; r<basis_.size(); ++r) {
            const int factor=value[pivots_[r]];
            if (!factor) continue;
            for (int j=pivots_[r]; j<width_; ++j)
                value[j]=gf_.add(value[j],gf_.mul(factor,basis_[r][j]));
        }
        return value;
    }

    bool contains(vector<int> value) const {
        value=reduce(std::move(value));
        return std::all_of(value.begin(),value.end(),[](int x){return x==0;});
    }

    const SparseMatrix<int, RowMajor>* innerH2() const {
        return inner_h2_ ? &*inner_h2_ : nullptr;
    }

    const SparseMatrix<int, RowMajor>* innerG2() const {
        return inner_g2_ ? &*inner_g2_ : nullptr;
    }

private:
    int width_;
    GF2p& gf_;
    vector<int> pivots_;
    vector<vector<int>> basis_;
    std::optional<SparseMatrix<int, RowMajor>> inner_h2_;
    std::optional<SparseMatrix<int, RowMajor>> inner_g2_;
};

bool isCascadeXStabilizer(const vector<int>& residual,
                          const FieldRowSpace& outer_space, GF2p& gf,
                          int n2, int k2, const vector<int>& scales) {
    const int blocks=static_cast<int>(residual.size())/n2;
    if (n2 == 12 && k2 == 10) {
        for (int group=0; group<2; ++group) {
            const int tail=10+group;
            for (int j=0; j<blocks; ++j) {
                int sum=residual[j*n2+tail];
                for (int a=0; a<5; ++a)
                    sum=gf.add(sum,residual[j*n2+5*group+a]);
                if (sum != 0) return false;
            }
            for (int a=0; a<5; ++a) {
                vector<int> coordinate(blocks);
                for (int j=0; j<blocks; ++j)
                    coordinate[j]=residual[j*n2+5*group+a];
                if (!outer_space.contains(std::move(coordinate))) return false;
            }
        }
        return true;
    }
    // A single-check inner code may store a path basis rather than [I|A].
    // Its basis-independent membership check uses H2-derived scales below.
    if(const auto* g2=outer_space.innerG2(); g2 && k2 < n2-1) {
        if(g2->rows()==k2 && g2->cols()==n2) {
            bool systematic=true;
            for(int a=0;a<k2;++a) for(int b=0;b<k2;++b)
                systematic &= g2->coeff(a,b)==(a==b ? 1 : 0);
            if(!systematic)
                throw invalid_argument(
                    "Generic cascade X-stabilizer check requires G2=[I|A]");
            for(int block=0;block<blocks;++block)
                for(int p=k2;p<n2;++p) {
                    int expected=0;
                    for(int a=0;a<k2;++a)
                        expected=gf.add(expected,gf.mul(
                            residual[block*n2+a],g2->coeff(a,p)));
                    if(expected!=residual[block*n2+p]) return false;
                }
            for(int a=0;a<k2;++a) {
                vector<int> coordinate(blocks);
                for(int block=0;block<blocks;++block)
                    coordinate[block]=residual[block*n2+a];
                if(!outer_space.contains(std::move(coordinate))) return false;
            }
            return true;
        }
    }
    for (int j=0; j<blocks; ++j) {
        int sum=0;
        for (int a=0; a<k2; ++a)
            sum=gf.add(sum,gf.mul(residual[j*n2+a+1],gf.inv(scales[a+1])));
        if (sum!=gf.mul(residual[j*n2],gf.inv(scales[0]))) return false;
    }
    for (int a=0; a<k2; ++a) {
        vector<int> coordinate(blocks);
        for (int j=0; j<blocks; ++j)
            coordinate[j]=gf.mul(residual[j*n2+a+1],gf.inv(scales[a+1]));
        if (!outer_space.contains(std::move(coordinate))) return false;
    }
    return true;
}

bool isCascadeZStabilizer(const vector<int>& residual,
                          const FieldRowSpace& outer_space, GF2p& gf,
                          int n2, int k2, const vector<int>& scales) {
    const int blocks=static_cast<int>(residual.size())/n2;
    if (n2 == 12 && k2 == 10) {
        for (int group=0; group<2; ++group)
            for (int a=0; a<5; ++a) {
                vector<int> coordinate(blocks);
                for (int j=0; j<blocks; ++j)
                    coordinate[j]=gf.add(residual[j*n2+5*group+a],
                                         residual[j*n2+10+group]);
                if (!outer_space.contains(std::move(coordinate))) return false;
            }
        return true;
    }
    if(const auto* h2=outer_space.innerH2()) {
        const int checks=h2->rows();
        if(h2->cols()==n2 && checks==n2-k2) {
            bool systematic_parity=true;
            for(int a=0;a<checks;++a) for(int b=0;b<checks;++b)
                systematic_parity &= h2->coeff(a,k2+b)==(a==b ? 1 : 0);
            if(!systematic_parity)
                throw invalid_argument(
                    "Generic cascade Z-stabilizer check requires H2=[A|I]");
            for(int a=0;a<k2;++a) {
                vector<int> coordinate(blocks);
                for(int block=0;block<blocks;++block) {
                    int value=residual[block*n2+a];
                    for(int t=0;t<checks;++t)
                        value=gf.add(value,gf.mul(
                            residual[block*n2+k2+t],h2->coeff(t,a)));
                    coordinate[block]=value;
                }
                if(!outer_space.contains(std::move(coordinate))) return false;
            }
            return true;
        }
    }
    for (int a=0; a<k2; ++a) {
        vector<int> coordinate(blocks);
        for (int j=0; j<blocks; ++j)
            coordinate[j]=gf.add(gf.mul(residual[j*n2+a],scales[a]),
                                 gf.mul(residual[j*n2+k2],scales[k2]));
        if (!outer_space.contains(std::move(coordinate))) return false;
    }
    return true;
}

// Canonical quotient-space label for the cascade Z-stabilizer coset. Equal
// labels are equivalent iff their difference is a Z stabilizer. This is for
// candidate diagnostics/selection; it never uses the sampled error.
vector<int> cascadeZCosetSignature(const vector<int>& word,
                                   const FieldRowSpace& outer_space,
                                   GF2p& gf, int n2, int k2,
                                   const vector<int>& scales) {
    if (word.size()%static_cast<size_t>(n2)!=0)
        throw invalid_argument("Cascade word has invalid inner length");
    const int blocks=static_cast<int>(word.size())/n2;
    vector<int> signature;
    if (n2==12 && k2==10) {
        signature.reserve(10*blocks);
        for (int group=0; group<2; ++group)
            for (int a=0; a<5; ++a) {
                vector<int> coordinate(blocks);
                for (int j=0; j<blocks; ++j)
                    coordinate[j]=gf.add(word[j*n2+5*group+a],
                                         word[j*n2+10+group]);
                coordinate=outer_space.reduce(std::move(coordinate));
                signature.insert(signature.end(),coordinate.begin(),coordinate.end());
        }
        return signature;
    }
    if(const auto* h2=outer_space.innerH2()) {
        const int checks=h2->rows();
        if(h2->cols()==n2 && checks==n2-k2) {
            for(int a=0;a<checks;++a) for(int b=0;b<checks;++b)
                if(h2->coeff(a,k2+b)!=(a==b ? 1 : 0))
                    throw invalid_argument(
                        "Generic cascade Z-coset key requires H2=[A|I]");
            signature.reserve(k2*blocks);
            for(int a=0;a<k2;++a) {
                vector<int> coordinate(blocks);
                for(int block=0;block<blocks;++block) {
                    int value=word[block*n2+a];
                    for(int t=0;t<checks;++t)
                        value=gf.add(value,gf.mul(
                            word[block*n2+k2+t],h2->coeff(t,a)));
                    coordinate[block]=value;
                }
                coordinate=outer_space.reduce(std::move(coordinate));
                signature.insert(signature.end(),coordinate.begin(),coordinate.end());
            }
            return signature;
        }
    }
    signature.reserve(k2*blocks);
    for (int a=0; a<k2; ++a) {
        vector<int> coordinate(blocks);
        for (int j=0; j<blocks; ++j)
            coordinate[j]=gf.add(gf.mul(word[j*n2+a],scales[a]),
                                 gf.mul(word[j*n2+k2],scales[k2]));
        coordinate=outer_space.reduce(std::move(coordinate));
        signature.insert(signature.end(),coordinate.begin(),coordinate.end());
    }
    return signature;
}

string cascadeZCosetKey(const vector<int>& word,
                       const FieldRowSpace& outer_space, GF2p& gf,
                       int n2, int k2, const vector<int>& scales) {
    const auto signature=cascadeZCosetSignature(
        word,outer_space,gf,n2,k2,scales);
    string key;
    key.reserve(signature.size()*2);
    for (int value:signature) {
        key.push_back(static_cast<char>(value&0xff));
        key.push_back(static_cast<char>((value>>8)&0xff));
    }
    return key;
}

using SoftVec = std::vector<double>;
using SoftMatrix = std::vector<std::vector<SoftVec>>;

void fwhtLocal(SoftVec& a, bool inverse) {
    const int n = static_cast<int>(a.size());
    for (int len = 1; len < n; len <<= 1) {
        for (int i = 0; i < n; i += 2 * len) {
            for (int j = 0; j < len; ++j) {
                const double u = a[i + j], v = a[i + len + j];
                a[i + j] = u + v;
                a[i + len + j] = u - v;
            }
        }
    }
    if (inverse) for (double& x : a) x /= n;
}

void normalizeLocal(SoftVec& v) {
    constexpr double floor = 1e-300;
    double sum = 0.0;
    for (double& x : v) { x = std::max(floor, x); sum += x; }
    for (double& x : v) x /= sum;
}

SoftVec xorConvolution(const SoftVec& a, const SoftVec& b) {
    SoftVec fa = a, fb = b;
    fwhtLocal(fa, false); fwhtLocal(fb, false);
    for (size_t i = 0; i < fa.size(); ++i) fa[i] *= fb[i];
    fwhtLocal(fa, true);
    for (double& x : fa) if (x < 0.0) x = 0.0;
    normalizeLocal(fa);
    return fa;
}

struct OuterSoftResult {
    SoftMatrix extrinsic;
    std::vector<int> hard;
    std::vector<SoftVec> posterior;
    bool valid;
};

struct RedundantOuterCheck {
    int first_row;
    int second_row;
    int second_scale;
};

// Add low-weight GF(q)-linear combinations of pairs of outer checks that
// overlap in a column. Choosing the coefficient ratio cancels that shared
// symbol, so each generated row is a valid redundant check and usually stays
// substantially sparser than an arbitrary row combination.
std::vector<RedundantOuterCheck> buildPairwiseRedundantOuterChecks(
    const SparseMatrix<int,RowMajor>& H,
    GF2p& gf,
    size_t max_additional_rows) {
    struct Candidate {
        RedundantOuterCheck transform;
        int weight;
        int first;
        int second;
    };
    const int rows=H.rows();
    std::vector<std::vector<std::pair<int,int>>> row_entries(rows);
    std::vector<std::vector<int>> column_rows(H.cols());
    for(int r=0;r<rows;++r)
        for(SparseMatrix<int,RowMajor>::InnerIterator it(H,r);it;++it) {
            row_entries[r].emplace_back(it.col(),it.value());
            column_rows[it.col()].push_back(r);
        }
    std::vector<Candidate> candidates;
    std::unordered_set<std::string> seen;
    auto encode=[](const std::map<int,int>& entries) {
        std::string key;
        key.reserve(entries.size()*6);
        for(const auto& [col,value]:entries) {
            for(int shift=0;shift<4;++shift)
                key.push_back(static_cast<char>((static_cast<uint32_t>(col)>>(shift*8))&255));
            key.push_back(static_cast<char>(value&255));
            key.push_back(static_cast<char>((value>>8)&255));
        }
        return key;
    };
    for(int col=0;col<H.cols();++col) {
        const auto& incident=column_rows[col];
        for(size_t a=0;a<incident.size();++a)
            for(size_t b=a+1;b<incident.size();++b) {
                const int first=incident[a], second=incident[b];
                int first_value=0, second_value=0;
                for(const auto& [c,value]:row_entries[first]) if(c==col) first_value=value;
                for(const auto& [c,value]:row_entries[second]) if(c==col) second_value=value;
                if(first_value==0 || second_value==0) continue;
                const int scale=gf.mul(first_value,gf.inv(second_value));
                std::map<int,int> combined;
                for(const auto& [c,value]:row_entries[first]) combined[c]=value;
                for(const auto& [c,value]:row_entries[second]) {
                    const int scaled=gf.mul(scale,value);
                    const auto found=combined.find(c);
                    if(found==combined.end()) combined.emplace(c,scaled);
                    else {
                        found->second=gf.add(found->second,scaled);
                        if(found->second==0) combined.erase(found);
                    }
                }
                if(combined.empty()) continue;
                const std::string key=encode(combined);
                if(!seen.insert(key).second) continue;
                candidates.push_back({{first,second,scale},
                    static_cast<int>(combined.size()),first,second});
            }
    }
    std::stable_sort(candidates.begin(),candidates.end(),
        [](const Candidate& a,const Candidate& b) {
            if(a.weight!=b.weight) return a.weight<b.weight;
            if(a.first!=b.first) return a.first<b.first;
            return a.second<b.second;
        });
    if(candidates.size()>max_additional_rows) {
        std::vector<int> row_usage(static_cast<size_t>(rows),0);
        std::vector<char> selected(candidates.size(),0);
        std::vector<Candidate> balanced;
        balanced.reserve(max_additional_rows);
        for(size_t step=0;step<max_additional_rows;++step) {
            size_t best=candidates.size();
            int best_usage=std::numeric_limits<int>::max();
            for(size_t index=0;index<candidates.size();++index) if(!selected[index]) {
                const auto& candidate=candidates[index];
                const int usage=row_usage[candidate.first]+row_usage[candidate.second];
                if(usage<best_usage) { best=index; best_usage=usage; }
            }
            if(best==candidates.size()) break;
            selected[best]=1;
            balanced.push_back(candidates[best]);
            ++row_usage[candidates[best].first];
            ++row_usage[candidates[best].second];
        }
        candidates=std::move(balanced);
    }
    std::vector<RedundantOuterCheck> result;
    result.reserve(candidates.size());
    for(const auto& candidate:candidates) result.push_back(candidate.transform);
    return result;
}

struct ProjectionDiagnostics {
    int bp_valid = 0;
    int raw_valid = 0;
    int final_valid = 0;
    int outer_refreshed = 0;
    int repaired = 0;
    int repair_failed = 0;
    int hard_disagreements = 0;
    int joint_candidates_tested = 0;
    int joint_updates = 0;
    int inactivation_attempts = 0;
    int inactivation_successes = 0;
    int inactivated_checks = 0;
    std::vector<std::vector<std::vector<int>>> outer_candidate_list;
    std::vector<std::vector<std::vector<int>>> neighbor_candidate_list;
    std::vector<std::vector<std::vector<int>>> pair_candidate_list;
    int pair_audit_left = -1;
    int pair_audit_candidates = 0;
    double joint_log_gain = 0.0;
};

OuterSoftResult runOuterSoftBP(
    NB_BP_Decoder& decoder,
    const SparseMatrix<int, RowMajor>& H,
    const std::vector<int>& syndrome,
    const std::vector<SoftVec>& priors,
    int iterations,
    GF2p& gf,
    bool layered_schedule=false) {
    const int q = gf.getSize(), rows = H.rows(), vars = H.cols();
    std::vector<std::vector<int>> row_cols(rows);
    std::vector<std::vector<std::pair<int,int>>> var_edges(vars);
    SoftMatrix v2c(rows), c2v(rows);
    for (int r = 0; r < rows; ++r) {
        for (SparseMatrix<int, RowMajor>::InnerIterator it(H, r); it; ++it) {
            int e = static_cast<int>(row_cols[r].size());
            row_cols[r].push_back(it.col()); var_edges[it.col()].emplace_back(r,e);
            v2c[r].push_back(priors[it.col()]);
            c2v[r].push_back(SoftVec(q, 1.0/q));
        }
    }
    std::vector<int> hard(vars, 0);
    std::vector<SoftVec> posterior(vars, SoftVec(q, 1.0/q));
    bool valid = false;
    for (int iter = 0; iter < iterations; ++iter) {
        auto update_check=[&](int r) {
            const int d = row_cols[r].size();
            if (d == 0) return;
            std::vector<SoftVec> transformed(d, SoftVec(q));
            for (int e = 0; e < d; ++e) {
                SoftVec scaled(q, 0.0);
                const int h = H.coeff(r, row_cols[r][e]);
                for (int x = 0; x < q; ++x) scaled[gf.mul(h,x)] += v2c[r][e][x];
                fwhtLocal(scaled, false); transformed[e] = std::move(scaled);
            }
            std::vector<SoftVec> prefix(d+1, SoftVec(q,1.0)), suffix(d+1, SoftVec(q,1.0));
            for (int e=0;e<d;++e) for(int x=0;x<q;++x) prefix[e+1][x]=prefix[e][x]*transformed[e][x];
            for (int e=d-1;e>=0;--e) for(int x=0;x<q;++x) suffix[e][x]=suffix[e+1][x]*transformed[e][x];
            for (int e=0;e<d;++e) {
                SoftVec other(q);
                for(int x=0;x<q;++x) other[x]=prefix[e][x]*suffix[e+1][x];
                fwhtLocal(other,true);
                const int h=H.coeff(r,row_cols[r][e]);
                SoftVec out(q);
                for(int x=0;x<q;++x) out[x]=std::max(0.0,other[gf.sub(syndrome[r],gf.mul(h,x))]);
                normalizeLocal(out); c2v[r][e]=std::move(out);
            }
        };
        auto update_variable_messages=[&](int v) {
            if (var_edges[v].empty()) return;
            for (const auto& [tr,te]:var_edges[v]) {
                SoftVec out=priors[v];
                for (const auto& [r,e]:var_edges[v]) if(r!=tr || e!=te) {
                    for(int x=0;x<q;++x) out[x]*=c2v[r][e][x];
                    normalizeLocal(out);
                }
                v2c[tr][te]=std::move(out);
            }
        };
        if (layered_schedule) {
            // Layered nonbinary BP: after each check update, immediately
            // refresh outgoing variable messages for its neighbors. Later
            // checks in the sweep therefore consume the newest information.
            for (int r=0;r<rows;++r) {
                update_check(r);
                for (int v:row_cols[r]) update_variable_messages(v);
            }
        } else {
            for (int r=0;r<rows;++r) update_check(r);
            for (int v=0;v<vars;++v) update_variable_messages(v);
        }
        for (int v=0;v<vars;++v) {
            SoftVec post=priors[v];
            for (const auto& [r,e]:var_edges[v]) {
                for(int x=0;x<q;++x) post[x]*=c2v[r][e][x];
                normalizeLocal(post);
            }
            posterior[v]=post;
            hard[v]=std::distance(post.begin(),std::max_element(post.begin(),post.end()));
        }
        valid=decoder.check_z_syndrome(H,hard,syndrome);
        if(valid) break;
    }
    SoftMatrix ext(vars, std::vector<SoftVec>());
    for(int v=0;v<vars;++v) {
        SoftVec x(q,1.0);
        for(const auto& [r,e]:var_edges[v]) {
            for(int a=0;a<q;++a) x[a]*=c2v[r][e][a];
            normalizeLocal(x);
        }
        ext[v].push_back(std::move(x));
    }
    return {std::move(ext),std::move(hard),std::move(posterior),valid};
}

struct InactivationResult {
    std::vector<int> estimate;
    int attempts = 0;
    int successful_check = -1;
    bool syndrome_ok = false;
};

// GF(q) adaptation of stabilizer-inactivation for one outer projection.
// `dual_checks` must contain low-weight rows in ker(H), so inactivating their
// support preserves a small exact syndrome-recovery system.
InactivationResult decodeOuterWithStabilizerInactivation(
    NB_BP_Decoder& decoder,
    const SparseMatrix<int,RowMajor>& H,
    const SparseMatrix<int,RowMajor>& dual_checks,
    const std::vector<int>& syndrome,
    const std::vector<SoftVec>& priors,
    const std::vector<SoftVec>& posterior,
    int iterations,
    int max_inactivations,
    GF2p& gf) {
    const int rows=H.rows(), vars=H.cols(), q=gf.getSize();
    if (dual_checks.cols()!=vars || syndrome.size()!=static_cast<size_t>(rows) ||
        priors.size()!=static_cast<size_t>(vars) ||
        posterior.size()!=static_cast<size_t>(vars) || iterations<=0 ||
        max_inactivations<0)
        throw std::invalid_argument("Invalid stabilizer-inactivation input");

    struct RankedCheck { double reliability; int row; };
    std::vector<RankedCheck> order;
    order.reserve(dual_checks.rows());
    for (int check=0;check<dual_checks.rows();++check) {
        double reliability=0.0;
        bool nonempty=false;
        for (SparseMatrix<int,RowMajor>::InnerIterator it(dual_checks,check);it;++it) {
            nonempty=true;
            double first=0.0,second=0.0;
            for (double probability:posterior[it.col()]) {
                if(probability>first) { second=first; first=probability; }
                else if(probability>second) second=probability;
            }
            reliability+=std::log(std::max(1e-300,static_cast<double>(first))) -
                         std::log(std::max(1e-300,static_cast<double>(second)));
        }
        if(nonempty) order.push_back({reliability,check});
    }
    std::stable_sort(order.begin(),order.end(),[](const auto& a,const auto& b) {
        if(a.reliability!=b.reliability) return a.reliability<b.reliability;
        return a.row<b.row;
    });

    InactivationResult result;
    result.estimate.assign(vars,0);
    const int trials=std::min(max_inactivations,static_cast<int>(order.size()));
    for(int attempt=0;attempt<trials;++attempt) {
        const int dual_row=order[attempt].row;
        std::vector<int> support;
        std::vector<char> inactivated(static_cast<size_t>(vars),0);
        for(SparseMatrix<int,RowMajor>::InnerIterator it(dual_checks,dual_row);it;++it) {
            support.push_back(it.col());
            inactivated[it.col()]=1;
        }
        if(support.empty()) continue;

        // The selected dual check must commute with every row of H. Verify
        // this invariant at runtime before using it to remove variables.
        bool orthogonal=true;
        for(int r=0;r<rows && orthogonal;++r) {
            int dot=0;
            for(SparseMatrix<int,RowMajor>::InnerIterator it(H,r);it;++it)
                if(inactivated[it.col()])
                    dot=gf.add(dot,gf.mul(it.value(),dual_checks.coeff(dual_row,it.col())));
            orthogonal=(dot==0);
        }
        if(!orthogonal)
            throw std::runtime_error("Selected inactivation support is not in ker(H)");

        std::vector<int> outside_columns;
        std::vector<int> old_to_new(static_cast<size_t>(vars),-1);
        for(int v=0;v<vars;++v) if(!inactivated[v]) {
            old_to_new[v]=static_cast<int>(outside_columns.size());
            outside_columns.push_back(v);
        }
        std::vector<int> kept_rows;
        for(int r=0;r<rows;++r) {
            bool touches=false;
            for(SparseMatrix<int,RowMajor>::InnerIterator it(H,r);it;++it)
                if(inactivated[it.col()]) { touches=true; break; }
            if(!touches) kept_rows.push_back(r);
        }
        if(kept_rows.empty()) continue;

        SparseMatrix<int,RowMajor> reduced_h(
            static_cast<int>(kept_rows.size()),static_cast<int>(outside_columns.size()));
        std::vector<Triplet<int>> reduced_entries;
        std::vector<int> reduced_syndrome(kept_rows.size());
        std::vector<SoftVec> reduced_priors(outside_columns.size());
        for(size_t j=0;j<outside_columns.size();++j)
            reduced_priors[j]=priors[outside_columns[j]];
        for(size_t r=0;r<kept_rows.size();++r) {
            const int old_row=kept_rows[r];
            reduced_syndrome[r]=syndrome[old_row];
            for(SparseMatrix<int,RowMajor>::InnerIterator it(H,old_row);it;++it)
                reduced_entries.emplace_back(static_cast<int>(r),old_to_new[it.col()],it.value());
        }
        reduced_h.setFromTriplets(reduced_entries.begin(),reduced_entries.end());
        reduced_h.makeCompressed();
        const auto outside=runOuterSoftBP(
            decoder,reduced_h,reduced_syndrome,reduced_priors,iterations,gf);
        ++result.attempts;
        if(!outside.valid) continue;

        // Solve every check touching the inactivated support, then choose
        // among the affine solutions using only decoder posterior values.
        std::vector<std::vector<int>> local_matrix;
        std::vector<int> local_rhs;
        local_matrix.reserve(static_cast<size_t>(rows));
        local_rhs.reserve(static_cast<size_t>(rows));
        std::vector<int> support_position(static_cast<size_t>(vars),-1);
        for(size_t j=0;j<support.size();++j) support_position[support[j]]=static_cast<int>(j);
        std::vector<int> outside_estimate(static_cast<size_t>(vars),0);
        for(size_t j=0;j<outside_columns.size();++j)
            outside_estimate[outside_columns[j]]=outside.hard[j];
        bool inconsistent=false;
        for(int r=0;r<rows;++r) {
            std::vector<int> equation(support.size(),0);
            bool touches=false;
            int rhs=syndrome[r];
            for(SparseMatrix<int,RowMajor>::InnerIterator it(H,r);it;++it) {
                const int position=support_position[it.col()];
                if(position>=0) {
                    equation[position]=it.value();
                    touches=true;
                } else {
                    rhs=gf.add(rhs,gf.mul(it.value(),outside_estimate[it.col()]));
                }
            }
            if(touches) {
                local_matrix.push_back(std::move(equation));
                local_rhs.push_back(rhs);
            } else if(rhs!=0) {
                inconsistent=true;
                break;
            }
        }
        if(inconsistent || local_matrix.empty()) continue;

        int rank=0;
        std::vector<int> pivot_columns;
        for(int col=0;col<static_cast<int>(support.size()) && rank<static_cast<int>(local_matrix.size());++col) {
            int pivot=rank;
            while(pivot<static_cast<int>(local_matrix.size()) && local_matrix[pivot][col]==0) ++pivot;
            if(pivot==static_cast<int>(local_matrix.size())) continue;
            std::swap(local_matrix[rank],local_matrix[pivot]);
            std::swap(local_rhs[rank],local_rhs[pivot]);
            const int inverse=gf.inv(local_matrix[rank][col]);
            for(size_t j=0;j<support.size();++j)
                local_matrix[rank][j]=gf.mul(local_matrix[rank][j],inverse);
            local_rhs[rank]=gf.mul(local_rhs[rank],inverse);
            for(int r=0;r<static_cast<int>(local_matrix.size());++r) if(r!=rank) {
                const int factor=local_matrix[r][col];
                if(!factor) continue;
                for(size_t j=0;j<support.size();++j)
                    local_matrix[r][j]=gf.add(
                        local_matrix[r][j],gf.mul(factor,local_matrix[rank][j]));
                local_rhs[r]=gf.add(local_rhs[r],gf.mul(factor,local_rhs[rank]));
            }
            pivot_columns.push_back(col);
            ++rank;
        }
        for(int r=rank;r<static_cast<int>(local_matrix.size());++r) {
            const bool zero=std::all_of(local_matrix[r].begin(),local_matrix[r].end(),
                                        [](int x){return x==0;});
            if(zero && local_rhs[r]!=0) inconsistent=true;
        }
        if(inconsistent) continue;
        std::vector<char> is_pivot(support.size(),0);
        std::vector<int> pivot_row_for_column(support.size(),-1);
        for(int r=0;r<rank;++r) {
            is_pivot[pivot_columns[r]]=1;
            pivot_row_for_column[pivot_columns[r]]=r;
        }
        std::vector<int> free_columns;
        for(int col=0;col<static_cast<int>(support.size());++col)
            if(!is_pivot[col]) free_columns.push_back(col);
        size_t solution_count=1;
        for(size_t i=0;i<free_columns.size();++i) {
            if(solution_count>65536/static_cast<size_t>(q)) { solution_count=65537; break; }
            solution_count*=static_cast<size_t>(q);
        }
        if(solution_count>65536) continue;

        double best_score=-std::numeric_limits<double>::infinity();
        std::vector<int> best_candidate;
        for(size_t code=0;code<solution_count;++code) {
            std::vector<int> values(support.size(),0);
            size_t remainder=code;
            for(int col:free_columns) {
                values[col]=static_cast<int>(remainder%static_cast<size_t>(q));
                remainder/=static_cast<size_t>(q);
            }
            for(int col:pivot_columns) {
                const int row=pivot_row_for_column[col];
                int value=local_rhs[row];
                for(int free_col:free_columns)
                    value=gf.add(value,gf.mul(local_matrix[row][free_col],values[free_col]));
                values[col]=value;
            }
            std::vector<int> candidate=outside_estimate;
            double score=0.0;
            for(size_t j=0;j<support.size();++j) {
                candidate[support[j]]=values[j];
                score+=std::log(std::max(1e-300,posterior[support[j]][values[j]]));
            }
            if(score>best_score && decoder.check_z_syndrome(H,candidate,syndrome)) {
                best_score=score;
                best_candidate=std::move(candidate);
            }
        }
        if(!best_candidate.empty()) {
            result.estimate=std::move(best_candidate);
            result.successful_check=dual_row;
            result.syndrome_ok=true;
            return result;
        }
    }
    return result;
}

// Structured soft X decoding on the H_Z cascade graph. Each selector
// coordinate is an outer LDPC graph; each physical block is an inner H2 graph.
// Only component-to-component extrinsic messages are exchanged, so the channel
// prior is counted once at each variable.
OuterSoftResult runCascadeXSoftBP(
    NB_BP_Decoder& decoder,
    const SparseMatrix<int,RowMajor>& outer_h,
    const SparseMatrix<int,RowMajor>& hz_cascade,
    const std::vector<int>& syndrome,
    const std::vector<SoftVec>& channel_priors,
    int turbo_rounds,
    int outer_iterations,
    int inner_iterations,
    GF2p& gf) {
    const int n1=outer_h.cols(), m1=outer_h.rows();
    const int n2=hz_cascade.cols()/std::max(1,n1), q=gf.getSize();
    if(n1<=m1 || n2<=1 || hz_cascade.cols()!=n1*n2 || turbo_rounds<=0 ||
       outer_iterations<=0 || inner_iterations<=0 ||
       channel_priors.size()!=static_cast<size_t>(n1*n2))
        throw std::invalid_argument("Invalid structured soft-X dimensions or iterations");
    const int numerator=hz_cascade.rows()-m1*n2;
    const int denominator=n1-m1;
    if(numerator<=0 || numerator%denominator!=0)
        throw std::invalid_argument("Cannot infer H2 dimensions for structured soft-X BP");
    const int inner_checks=numerator/denominator;
    const int k=n2-inner_checks;
    if(inner_checks<1 || k<1 || hz_cascade.rows()!=m1*k+n1*inner_checks ||
       syndrome.size()!=static_cast<size_t>(hz_cascade.rows()))
        throw std::invalid_argument("Unsupported H_Z cascade dimensions for structured soft-X BP");

    SparseMatrix<int,RowMajor> h2(inner_checks,n2);
    for(int r=0;r<inner_checks;++r)
        for(SparseMatrix<int,RowMajor>::InnerIterator it(hz_cascade,m1*k+r);it;++it) {
            if(it.col()>=n2)
                throw std::invalid_argument("H2 row couples distinct physical blocks");
            h2.insert(r,it.col())=it.value();
        }
    h2.makeCompressed();
    for(int block=0;block<n1;++block)
        for(int r=0;r<inner_checks;++r) {
            const int full_row=m1*k+block*inner_checks+r;
            for(SparseMatrix<int,RowMajor>::InnerIterator it(hz_cascade,full_row);it;++it) {
                if(it.col()/n2!=block || h2.coeff(r,it.col()%n2)!=it.value())
                    throw std::invalid_argument("H2 blocks are not repeated block-diagonal copies");
            }
            for(int col=0;col<n2;++col)
                if(hz_cascade.coeff(full_row,block*n2+col)!=h2.coeff(r,col))
                    throw std::invalid_argument("H2 block coefficients differ across blocks");
        }
    for(int r=0;r<m1;++r) for(int a=0;a<k;++a) {
        const int upper_row=r*k+a;
        for(int block=0;block<n1;++block)
            for(int col=0;col<n2;++col) {
                const int expected=col==a ? outer_h.coeff(r,block) : 0;
                if(hz_cascade.coeff(upper_row,block*n2+col)!=expected)
                    throw std::invalid_argument("H_Z upper block is not outer_h tensor selector");
            }
    }

    std::vector<SoftVec> outer_to_inner(n1*k,SoftVec(q,1.0/q));
    std::vector<SoftVec> inner_to_variable(n1*n2,SoftVec(q,1.0/q));
    std::vector<OuterSoftResult> outer_results(k);
    std::vector<int> hard(n1*n2,0);
    std::vector<SoftVec> posterior(n1*n2,SoftVec(q,1.0/q));
    bool valid=false;
    for(int round=0;round<turbo_rounds;++round) {
        for(int a=0;a<k;++a) {
            std::vector<SoftVec> priors(n1);
            std::vector<int> coordinate_syndrome(m1);
            for(int r=0;r<m1;++r) coordinate_syndrome[r]=syndrome[r*k+a];
            for(int block=0;block<n1;++block) {
                const int v=block*n2+a;
                priors[block]=channel_priors[v];
                for(int x=0;x<q;++x) priors[block][x]*=inner_to_variable[block*n2+a][x];
                normalizeLocal(priors[block]);
            }
            outer_results[a]=runOuterSoftBP(
                decoder,outer_h,coordinate_syndrome,priors,outer_iterations,gf);
            for(int block=0;block<n1;++block)
                outer_to_inner[block*k+a]=outer_results[a].extrinsic[block][0];
        }
        for(int block=0;block<n1;++block) {
            std::vector<SoftVec> priors(n2);
            std::vector<int> inner_syndrome(inner_checks);
            for(int r=0;r<inner_checks;++r)
                inner_syndrome[r]=syndrome[m1*k+block*inner_checks+r];
            for(int a=0;a<n2;++a) {
                const int v=block*n2+a;
                priors[a]=channel_priors[v];
                if(a<k) for(int x=0;x<q;++x)
                    priors[a][x]*=outer_to_inner[block*k+a][x];
                normalizeLocal(priors[a]);
            }
            const auto inner_result=runOuterSoftBP(
                decoder,h2,inner_syndrome,priors,inner_iterations,gf);
            for(int a=0;a<n2;++a)
                inner_to_variable[block*n2+a]=inner_result.extrinsic[a][0];
        }
        for(int block=0;block<n1;++block) for(int a=0;a<n2;++a) {
            const int v=block*n2+a;
            posterior[v]=channel_priors[v];
            if(a<k) for(int x=0;x<q;++x)
                posterior[v][x]*=outer_to_inner[block*k+a][x];
            for(int x=0;x<q;++x)
                posterior[v][x]*=inner_to_variable[block*n2+a][x];
            normalizeLocal(posterior[v]);
            hard[v]=std::distance(posterior[v].begin(),
                std::max_element(posterior[v].begin(),posterior[v].end()));
        }
        valid=decoder.check_x_syndrome(hz_cascade,hard,syndrome);
        if(valid) break;
    }
    SoftMatrix extrinsic(n1*n2,std::vector<SoftVec>(1,SoftVec(q,1.0/q)));
    for(int v=0;v<n1*n2;++v) extrinsic[v][0]=posterior[v];
    return {std::move(extrinsic),std::move(hard),std::move(posterior),valid};
}

// Push a posterior over GF(q) X-error symbols through the physical Pauli
// channel, retaining correlations among the bits of each field symbol.
// The bitwise channel tensor product is applied in O(e*2^e), not O(q^2).
std::vector<SoftVec> softXPosteriorToZPrior(
    const std::vector<SoftVec>& x_posterior,
    const FieldBasis& basis,
    const GF2p& gf,
    const PauliChannel& channel) {
    const int q=gf.getSize(), degree=gf.getPower();
    std::vector<SoftVec> z_prior(x_posterior.size(),SoftVec(q,0.0));
    const double z_if_x0=channel.zGivenX(false);
    const double z_if_x1=channel.zGivenX(true);
    for(size_t variable=0;variable<x_posterior.size();++variable) {
        if(x_posterior[variable].size()!=static_cast<size_t>(q))
            throw std::invalid_argument("Invalid X-symbol posterior size");
        SoftVec bits(q,0.0);
        for(int x=0;x<q;++x) {
            const double probability=x_posterior[variable][x];
            if(!std::isfinite(probability)||probability<0.0)
                throw std::invalid_argument("Invalid X-symbol posterior value");
            bits[basis.toDual(x)]+=probability;
        }
        normalizeLocal(bits);
        for(int bit=0;bit<degree;++bit) {
            const int mask=1<<bit;
            for(int pattern=0;pattern<q;++pattern) {
                if(pattern&mask) continue;
                const double x0=bits[pattern],x1=bits[pattern|mask];
                bits[pattern]=(1.0-z_if_x0)*x0+(1.0-z_if_x1)*x1;
                bits[pattern|mask]=z_if_x0*x0+z_if_x1*x1;
            }
        }
        normalizeLocal(bits);
        z_prior[variable]=std::move(bits);
    }
    return z_prior;
}

struct SyndromeRepairResult {
    std::vector<int> estimate;
    bool valid = false;
    int depth = 0;
    int unsatisfied = 0;
    double prior_cost = 0.0;
};

// Diagnostic-only bounded syndrome-guided search. It uses only H, the target
// syndrome, and the supplied channel priors; the transmitted error is never
// an input. This is not enabled in the normal simulation path.
SyndromeRepairResult repairSyndromeBeam(
    const SparseMatrix<int, RowMajor>& H,
    const std::vector<int>& syndrome,
    const std::vector<int>& initial,
    const std::vector<SoftVec>& priors,
    GF2p& gf,
    int beam_width=48,
    int max_depth=8,
    const std::vector<std::pair<int,int>>* physical_bit_constraints=nullptr,
    const std::vector<int>* inner_scales=nullptr,
    int inner_length=0) {
    struct Node {
        std::vector<int> word;
        std::vector<int> residual;
        int unsatisfied = 0;
        double cost = 0.0;
        int depth = 0;
    };
    const int rows=H.rows(), vars=H.cols();
    if (syndrome.size()!=static_cast<size_t>(rows) ||
        initial.size()!=static_cast<size_t>(vars) ||
        priors.size()!=static_cast<size_t>(vars))
        throw std::invalid_argument("Invalid syndrome-repair input dimensions");
    if (physical_bit_constraints &&
        (physical_bit_constraints->size()!=static_cast<size_t>(vars) ||
         !inner_scales || inner_length<=0 ||
         inner_scales->size()!=static_cast<size_t>(inner_length)))
        throw std::invalid_argument("Invalid constrained syndrome-repair inputs");
    std::vector<std::vector<std::pair<int,int>>> var_checks(vars);
    for(int r=0;r<rows;++r)
        for(SparseMatrix<int,RowMajor>::InnerIterator it(H,r);it;++it)
            var_checks[it.col()].emplace_back(r,it.value());
    auto evaluate=[&](const std::vector<int>& word) {
        std::vector<int> residual(rows);
        int unsat=0;
        for(int r=0;r<rows;++r) {
            int value=0;
            for(SparseMatrix<int,RowMajor>::InnerIterator it(H,r);it;++it)
                value=gf.add(value,gf.mul(it.value(),word[it.col()]));
            residual[r]=gf.sub(syndrome[r],value);
            unsat+=residual[r]!=0;
        }
        return std::pair<std::vector<int>,int>{std::move(residual),unsat};
    };
    auto prior_cost=[&](const std::vector<int>& word) {
        double cost=0.0;
        for(int v=0;v<vars;++v)
            cost-=std::log(std::max(1e-300,priors[v][word[v]]));
        return cost;
    };
    auto [initial_residual,initial_unsatisfied]=evaluate(initial);
    Node root{initial,std::move(initial_residual),initial_unsatisfied,
              prior_cost(initial),0};
    SyndromeRepairResult best{initial,initial_unsatisfied==0,0,
                              initial_unsatisfied,0.0};
    if(best.valid) return best;
    std::vector<Node> beam{std::move(root)};
    std::unordered_set<std::string> visited;
    auto keyOf=[](const std::vector<int>& word) {
        std::string key; key.resize(word.size()*2);
        for(size_t i=0;i<word.size();++i) {
            key[2*i]=static_cast<char>(word[i]&255);
            key[2*i+1]=static_cast<char>((word[i]>>8)&255);
        }
        return key;
    };
    visited.insert(keyOf(initial));
    const double base_cost=beam.front().cost;
    for(int depth=1;depth<=max_depth && !beam.empty();++depth) {
        std::vector<Node> next;
        for(const Node& node:beam) {
            int target=-1, target_degree=std::numeric_limits<int>::max();
            for(int r=0;r<rows;++r) if(node.residual[r]!=0) {
                int degree=0;
                for(SparseMatrix<int,RowMajor>::InnerIterator it(H,r);it;++it) ++degree;
                if(degree<target_degree) { target=r; target_degree=degree; }
            }
            if(target<0) {
                best={node.word,true,node.depth,node.unsatisfied,node.cost-base_cost};
                return best;
            }
            for(SparseMatrix<int,RowMajor>::InnerIterator it(H,target);it;++it) {
                const int v=it.col(), h=it.value();
                const int delta=gf.mul(node.residual[target],gf.inv(h));
                if(delta==0) continue;
                Node child=node;
                child.word[v]=gf.add(child.word[v],delta);
                if (physical_bit_constraints) {
                    const auto [mask,value]=(*physical_bit_constraints)[v];
                    if(mask) {
                        const int physical=gf.mul(child.word[v],
                            gf.inv((*inner_scales)[v%inner_length]));
                        if(((physical^value)&mask)!=0) continue;
                    }
                }
                for(const auto& [r,coefficient]:var_checks[v])
                    child.residual[r]=gf.sub(child.residual[r],
                                              gf.mul(coefficient,delta));
                child.unsatisfied=0;
                for(int value:child.residual) child.unsatisfied+=value!=0;
                child.depth=depth;
                child.cost=prior_cost(child.word);
                std::string key=keyOf(child.word);
                if(!visited.insert(std::move(key)).second) continue;
                if(child.unsatisfied==0)
                    return {std::move(child.word),true,depth,0,
                            child.cost-base_cost};
                next.push_back(std::move(child));
            }
        }
        std::sort(next.begin(),next.end(),[&](const Node& a,const Node& b) {
            // Allow temporary syndrome growth when a low-prior-cost change
            // may unlock a multi-check correction later in the search.
            const double score_a=10.0*a.unsatisfied+(a.cost-base_cost);
            const double score_b=10.0*b.unsatisfied+(b.cost-base_cost);
            if(score_a!=score_b) return score_a<score_b;
            return a.cost<b.cost;
        });
        if(static_cast<int>(next.size())>beam_width) next.resize(beam_width);
        beam=std::move(next);
        if(!beam.empty() && beam.front().unsatisfied<best.unsatisfied)
            best={beam.front().word,false,depth,beam.front().unsatisfied,
                  beam.front().cost-base_cost};
    }
    return best;
}

// Localized order-1 OSD around the currently unsatisfied checks. It forms a
// small GF(256) system from those rows, selects a reliability-ordered basis,
// then tests single-symbol deviations on the least reliable free variables.
SyndromeRepairResult localizedOsd(
    const SparseMatrix<int, RowMajor>& H,
    const std::vector<int>& syndrome,
    const std::vector<int>& initial,
    const std::vector<std::vector<double>>& posterior,
    GF2p& gf,
    int free_variables_to_test=12) {
    const int rows=H.rows(), vars=H.cols(), q=gf.getSize();
    std::vector<int> residual(rows), bad_rows;
    std::vector<std::vector<int>> row_vars(rows);
    std::vector<std::vector<std::pair<int,int>>> var_edges(vars);
    for(int r=0;r<rows;++r) {
        int value=0;
        for(SparseMatrix<int,RowMajor>::InnerIterator it(H,r);it;++it) {
            value=gf.add(value,gf.mul(it.value(),initial[it.col()]));
            row_vars[r].push_back(it.col());
            var_edges[it.col()].emplace_back(r,it.value());
        }
        residual[r]=gf.sub(syndrome[r],value);
        if(residual[r]) bad_rows.push_back(r);
    }
    SyndromeRepairResult best{initial,bad_rows.empty(),0,
                              static_cast<int>(bad_rows.size()),0.0};
    if(bad_rows.empty()) return best;
    std::vector<int> variables;
    std::vector<char> in_support(vars,0);
    for(int r:bad_rows) for(int v:row_vars[r]) if(!in_support[v]) {
        in_support[v]=1; variables.push_back(v);
    }
    auto confidence=[&](int v) {
        double first=-1.0, second=-1.0;
        for(double p:posterior[v]) {
            if(p>first) { second=first; first=p; }
            else if(p>second) second=p;
        }
        return std::log(std::max(1e-300,first))-
               std::log(std::max(1e-300,second));
    };
    std::sort(variables.begin(),variables.end(),[&](int a,int b) {
        return confidence(a)<confidence(b);
    });
    const int m=static_cast<int>(bad_rows.size()), n=static_cast<int>(variables.size());
    std::vector<std::vector<int>> a(m,std::vector<int>(n));
    std::vector<int> b(m);
    for(int i=0;i<m;++i) {
        b[i]=residual[bad_rows[i]];
        for(int j=0;j<n;++j) a[i][j]=H.coeff(bad_rows[i],variables[j]);
    }
    std::vector<int> pivot_cols;
    int rank=0;
    for(int col=0;col<n && rank<m;++col) {
        int pivot=rank;
        while(pivot<m && a[pivot][col]==0) ++pivot;
        if(pivot==m) continue;
        std::swap(a[rank],a[pivot]); std::swap(b[rank],b[pivot]);
        const int inverse=gf.inv(a[rank][col]);
        for(int j=col;j<n;++j) a[rank][j]=gf.mul(a[rank][j],inverse);
        b[rank]=gf.mul(b[rank],inverse);
        for(int i=0;i<m;++i) if(i!=rank && a[i][col]) {
            const int factor=a[i][col];
            for(int j=col;j<n;++j)
                a[i][j]=gf.sub(a[i][j],gf.mul(factor,a[rank][j]));
            b[i]=gf.sub(b[i],gf.mul(factor,b[rank]));
        }
        pivot_cols.push_back(col); ++rank;
    }
    if(rank<m) return best;
    std::vector<char> is_pivot(n,0);
    for(int c:pivot_cols) is_pivot[c]=1;
    std::vector<int> free_cols;
    for(int c=0;c<n;++c) if(!is_pivot[c]) free_cols.push_back(c);
    std::sort(free_cols.begin(),free_cols.end(),[&](int a_col,int b_col) {
        return confidence(variables[a_col])<confidence(variables[b_col]);
    });
    if(static_cast<int>(free_cols.size())>free_variables_to_test)
        free_cols.resize(free_variables_to_test);
    double best_score=-std::numeric_limits<double>::infinity();
    int tested=0;
    auto consider=[&](int free_col,int free_delta) {
        std::vector<int> delta(n,0);
        if(free_col>=0) delta[free_col]=free_delta;
        for(int i=0;i<m;++i) {
            int value=b[i];
            if(free_col>=0)
                value=gf.sub(value,gf.mul(a[i][free_col],free_delta));
            delta[pivot_cols[i]]=value;
        }
        std::vector<int> candidate=initial;
        for(int j=0;j<n;++j) if(delta[j])
            candidate[variables[j]]=gf.add(candidate[variables[j]],delta[j]);
        ++tested;
        int unsat=0;
        for(int r=0;r<rows;++r) {
            int value=0;
            for(SparseMatrix<int,RowMajor>::InnerIterator it(H,r);it;++it)
                value=gf.add(value,gf.mul(it.value(),candidate[it.col()]));
            unsat+=value!=syndrome[r];
        }
        double score=0.0;
        for(int v=0;v<vars;++v)
            score+=std::log(std::max(1e-300,posterior[v][candidate[v]]));
        if(unsat<best.unsatisfied || (unsat==best.unsatisfied && score>best_score)) {
            best={std::move(candidate),unsat==0,tested,unsat,-score};
            best_score=score;
        }
    };
    consider(-1,0); // OSD-0 candidate.
    for(int col:free_cols) for(int value=1;value<q;++value)
        consider(col,value); // OSD-1: one free symbol deviates from its hard value.
    return best;
}

struct FullOsdResult {
    std::vector<int> estimate;
    bool syndrome_ok = false;
    int rank = 0;
    int candidates_tested = 0;
    double log_posterior = -std::numeric_limits<double>::infinity();
};

// Full-matrix GF(q) OSD. In the legacy modes reliable symbols form the pivot basis; BP hard
// decisions are assigned to the free symbols, then one-symbol deviations are
// enumerated over the weakest free symbols. The complete syndrome is verified
// for each retained candidate before it can be returned.
FullOsdResult fullMatrixOsd1(
    const SparseMatrix<int, RowMajor>& H,
    const std::vector<int>& syndrome,
    const std::vector<int>& bp_hard,
    const std::vector<SoftVec>& posterior,
    GF2p& gf,
    int weakest_free_symbols=64,
    int osd_order=1,
    int alternatives_per_symbol=8,
    bool unreliable_pivots=false,
    std::vector<std::vector<int>>* candidate_list=nullptr,
    uint32_t basis_tie_seed=0,
    int basis_randomization_window=0,
    const std::vector<int>* audit_truth=nullptr,
    int* audit_free_deviations=nullptr,
    int* audit_selected_free_deviations=nullptr) {
    const int rows=H.rows(), vars=H.cols(), q=gf.getSize();
    if(syndrome.size()!=static_cast<size_t>(rows) ||
       bp_hard.size()!=static_cast<size_t>(vars) ||
       posterior.size()!=static_cast<size_t>(vars) || osd_order<0 || osd_order>2 ||
       alternatives_per_symbol<1)
        throw std::invalid_argument("Invalid full OSD input dimensions");
    std::vector<double> reliability(vars);
    std::vector<int> order(vars), inverse_order(vars);
    std::iota(order.begin(),order.end(),0);
    std::mt19937 tie_rng(basis_tie_seed);
    for(int v=0;v<vars;++v) {
        double first=-1.0, second=-1.0;
        for(double probability:posterior[v]) {
            if(probability>first) { second=first; first=probability; }
            else if(probability>second) second=probability;
        }
        reliability[v]=std::log(std::max(1e-300,first))-
                       std::log(std::max(1e-300,second));
    }
    // For projection repair, solve less reliable independent symbols and
    // preserve the hard decisions on the more reliable free symbols.
    const auto pivot_order=[&](int a,int b) {
        return unreliable_pivots ? reliability[a]<reliability[b]
                                 : reliability[a]>reliability[b];
    };
    if (basis_tie_seed==0) {
        std::sort(order.begin(),order.end(),pivot_order);
    } else if (basis_randomization_window>1) {
        std::sort(order.begin(),order.end(),pivot_order);
        for (int begin=0; begin<vars; begin+=basis_randomization_window) {
            const int end=std::min(vars,begin+basis_randomization_window);
            std::shuffle(order.begin()+begin,order.begin()+end,tie_rng);
        }
    } else {
        std::shuffle(order.begin(),order.end(),tie_rng);
        std::stable_sort(order.begin(),order.end(),pivot_order);
    }
    for(int pos=0;pos<vars;++pos) inverse_order[order[pos]]=pos;
    std::vector<int> a(static_cast<size_t>(rows)*vars,0), rhs=syndrome;
    for(int r=0;r<rows;++r)
        for(SparseMatrix<int,RowMajor>::InnerIterator it(H,r);it;++it)
            a[static_cast<size_t>(r)*vars+inverse_order[it.col()]]=it.value();
    std::vector<int> multiplication(static_cast<size_t>(q)*q);
    for(int x=0;x<q;++x) for(int y=0;y<q;++y)
        multiplication[static_cast<size_t>(x)*q+y]=gf.mul(x,y);
    std::vector<int> pivot_positions;
    int rank=0;
    for(int col=0;col<vars && rank<rows;++col) {
        int pivot=rank;
        while(pivot<rows && a[static_cast<size_t>(pivot)*vars+col]==0) ++pivot;
        if(pivot==rows) continue;
        if(pivot!=rank) {
            for(int j=0;j<vars;++j)
                std::swap(a[static_cast<size_t>(rank)*vars+j],
                          a[static_cast<size_t>(pivot)*vars+j]);
            std::swap(rhs[rank],rhs[pivot]);
        }
        const int inverse=gf.inv(a[static_cast<size_t>(rank)*vars+col]);
        for(int j=col;j<vars;++j)
            a[static_cast<size_t>(rank)*vars+j]=multiplication[
                static_cast<size_t>(a[static_cast<size_t>(rank)*vars+j])*q+inverse];
        rhs[rank]=multiplication[static_cast<size_t>(rhs[rank])*q+inverse];
        for(int r=0;r<rows;++r) if(r!=rank) {
            const int factor=a[static_cast<size_t>(r)*vars+col];
            if(!factor) continue;
            const size_t row_offset=static_cast<size_t>(r)*vars;
            const size_t pivot_offset=static_cast<size_t>(rank)*vars;
            for(int j=col;j<vars;++j) {
                const int product=multiplication[
                    static_cast<size_t>(factor)*q+a[pivot_offset+j]];
                a[row_offset+j]^=product; // GF(2^8) addition is XOR.
            }
            rhs[r]^=multiplication[static_cast<size_t>(factor)*q+rhs[rank]];
        }
        pivot_positions.push_back(col);
        ++rank;
    }
    for(int r=rank;r<rows;++r) if(rhs[r]!=0)
        return {bp_hard,false,rank,0,-std::numeric_limits<double>::infinity()};
    std::vector<char> is_pivot(vars,0);
    for(int col:pivot_positions) is_pivot[col]=1;
    std::vector<int> free_positions;
    for(int pos=0;pos<vars;++pos) if(!is_pivot[pos]) free_positions.push_back(pos);
    if (audit_truth) {
        if (audit_truth->size()!=static_cast<size_t>(vars))
            throw std::invalid_argument("OSD audit truth has incorrect length");
        int deviations=0;
        for (int pos:free_positions)
            deviations+=(*audit_truth)[order[pos]]!=bp_hard[order[pos]];
        if (audit_free_deviations) *audit_free_deviations=deviations;
    }
    const auto free_order=[&](int a_pos,int b_pos) {
        return reliability[order[a_pos]]<reliability[order[b_pos]];
    };
    if (basis_tie_seed==0)
        std::sort(free_positions.begin(),free_positions.end(),free_order);
    else std::stable_sort(free_positions.begin(),free_positions.end(),free_order);
    if(static_cast<int>(free_positions.size())>weakest_free_symbols)
        free_positions.resize(weakest_free_symbols);
    if (audit_truth && audit_selected_free_deviations) {
        int deviations=0;
        for (int pos:free_positions)
            deviations+=(*audit_truth)[order[pos]]!=bp_hard[order[pos]];
        *audit_selected_free_deviations=deviations;
    }
    std::vector<int> candidate=bp_hard;
    // Set every free variable to its BP hard value, then solve all pivot values.
    for(int row=0;row<rank;++row) {
        int value=rhs[row];
        const size_t row_offset=static_cast<size_t>(row)*vars;
        for(int pos=0;pos<vars;++pos) if(!is_pivot[pos])
            value^=multiplication[static_cast<size_t>(a[row_offset+pos])*q+
                                  bp_hard[order[pos]]];
        candidate[order[pivot_positions[row]]]=value;
    }
    auto scoreWord=[&](const std::vector<int>& word) {
        double score=0.0;
        for(int v=0;v<vars;++v)
            score+=std::log(std::max(1e-300,posterior[v][word[v]]));
        return score;
    };
    const std::vector<int> osd0=candidate;
    FullOsdResult best{osd0,false,rank,0,scoreWord(osd0)};
    auto valid=[&](const std::vector<int>& word) {
        for(int r=0;r<rows;++r) {
            int value=0;
            for(SparseMatrix<int,RowMajor>::InnerIterator it(H,r);it;++it)
                value=gf.add(value,gf.mul(it.value(),word[it.col()]));
            if(value!=syndrome[r]) return false;
        }
        return true;
    };
    best.syndrome_ok=valid(best.estimate);
    best.candidates_tested=1;
    if (candidate_list) {
        candidate_list->clear();
        if (best.syndrome_ok) candidate_list->push_back(best.estimate);
    }
    if (osd_order >= 1) for(int free_pos:free_positions) {
        const int free_var=order[free_pos], old_value=osd0[free_var];
        std::vector<int> alternatives;
        if (unreliable_pivots) {
            // Bounded candidate lists keep GF(1024) joint updates practical.
            alternatives.resize(q);
            std::iota(alternatives.begin(),alternatives.end(),0);
            const auto better=[&](int x,int y) {
                if(posterior[free_var][x]!=posterior[free_var][y])
                    return posterior[free_var][x]>posterior[free_var][y];
                return x<y;
            };
            const int keep=std::min(q,alternatives_per_symbol+1);
            std::partial_sort(alternatives.begin(),alternatives.begin()+keep,
                              alternatives.end(),better);
            alternatives.resize(keep);
        } else {
            alternatives.resize(q);
            std::iota(alternatives.begin(),alternatives.end(),0);
        }
        int considered=0;
        for(int new_value:alternatives) if(new_value!=old_value) {
            if(unreliable_pivots && considered++>=alternatives_per_symbol) break;
            const int delta=gf.sub(new_value,old_value);
            std::vector<int> trial=osd0;
            trial[free_var]=new_value;
            for(int row=0;row<rank;++row) {
                const int coefficient=a[static_cast<size_t>(row)*vars+free_pos];
                trial[order[pivot_positions[row]]]=gf.add(
                    osd0[order[pivot_positions[row]]],
                    gf.mul(coefficient,delta));
            }
            ++best.candidates_tested;
            if(!valid(trial)) continue;
            if (candidate_list) candidate_list->push_back(trial);
            const double score=scoreWord(trial);
            if(!best.syndrome_ok || score>best.log_posterior) {
                best.estimate=std::move(trial);
                best.syndrome_ok=true;
                best.log_posterior=score;
            }
        }
    }
    if (osd_order >= 2) {
        // OSD-2 over a bounded reliability-ordered subset. For each selected
        // pair, enumerate only the highest-posterior alternatives at each
        // free symbol, then update pivot symbols by linearity of H.
        const int pair_symbols=std::min(weakest_free_symbols,
                                       static_cast<int>(free_positions.size()));
        std::vector<std::vector<int>> alternatives(pair_symbols);
        for (int i=0; i<pair_symbols; ++i) {
            const int pos=free_positions[i], var=order[pos], old_value=osd0[var];
            std::vector<int> values(q);
            std::iota(values.begin(),values.end(),0);
            std::sort(values.begin(),values.end(),[&](int a_value,int b_value) {
                return posterior[var][a_value] > posterior[var][b_value];
            });
            for (int value:values) if (value!=old_value) {
                alternatives[i].push_back(value);
                if (static_cast<int>(alternatives[i].size())>=alternatives_per_symbol) break;
            }
        }
        for (int i=0; i<pair_symbols; ++i) {
            const int pos_i=free_positions[i], var_i=order[pos_i];
            const int old_i=osd0[var_i];
            for (int j=i+1; j<pair_symbols; ++j) {
                const int pos_j=free_positions[j], var_j=order[pos_j];
                const int old_j=osd0[var_j];
                for (int value_i:alternatives[i]) for (int value_j:alternatives[j]) {
                    const int delta_i=gf.sub(value_i,old_i);
                    const int delta_j=gf.sub(value_j,old_j);
                    std::vector<int> trial=osd0;
                    trial[var_i]=value_i;
                    trial[var_j]=value_j;
                    if (unreliable_pivots) {
                        // Preserve the exact syndrome after changing two
                        // free variables. The legacy branch is unchanged.
                        for(int row=0;row<rank;++row) {
                            const int coefficient_i=a[static_cast<size_t>(row)*vars+pos_i];
                            const int coefficient_j=a[static_cast<size_t>(row)*vars+pos_j];
                            trial[order[pivot_positions[row]]]=gf.add(
                                osd0[order[pivot_positions[row]]],
                                gf.add(gf.mul(coefficient_i,delta_i),
                                       gf.mul(coefficient_j,delta_j)));
                        }
                    }
                    for (int row=0; row<rank; ++row) {
                        const size_t row_offset=static_cast<size_t>(row)*vars;
                        const int pivot_delta=
                            gf.mul(a[row_offset+pos_i],delta_i) ^
                            gf.mul(a[row_offset+pos_j],delta_j);
                        trial[order[pivot_positions[row]]]=gf.add(
                            osd0[order[pivot_positions[row]]],pivot_delta);
                    }
                    ++best.candidates_tested;
                    if (!valid(trial)) continue;
                    if (candidate_list) candidate_list->push_back(trial);
                    const double score=scoreWord(trial);
                    if (!best.syndrome_ok || score>best.log_posterior) {
                        best.estimate=std::move(trial);
                        best.syndrome_ok=true;
                        best.log_posterior=score;
                    }
                }
            }
        }
    }
    return best;
}

struct GibbsCosetResult {
    struct LogicalDirection {
        int variable = -1;
        std::vector<std::pair<int,int>> pivot_updates;
    };
    std::vector<int> estimate;
    std::vector<std::string> sampled_class_keys;
    std::vector<LogicalDirection> logical_directions;
    std::vector<int> free_variables;
    std::vector<std::vector<std::pair<int,int>>> stabilizer_echelon;
    bool syndrome_ok = false;
    int samples = 0;
    int distinct_classes = 0;
    int winning_class_samples = 0;
    long long updates = 0;
    long long changed_updates = 0;
    long long exchange_attempts = 0;
    long long exchange_accepts = 0;
    int kernel_dimension = 0;
    int stabilizer_rank = 0;
    int logical_dimension = 0;
    long long logical_updates = 0;
    std::vector<double> replica_change_rates;
    double score_lag1_correlation = 0.0;
    double class_persistence = 0.0;
    double representative_log_posterior = -std::numeric_limits<double>::infinity();
};

// Gibbs sampling on the affine syndrome space H z = s. A single free
// polynomial-basis bit is updated at a time; pivot symbols are changed by
// the corresponding nullspace direction, so every transition stays on the
// exact syndrome coset. Samples are grouped by the exact cascade logical
// quotient signature, approximating posterior mass per logical class.
GibbsCosetResult sampleGfSyndromeCosets(
    const SparseMatrix<int, RowMajor>& H,
    const std::vector<int>& syndrome,
    const std::vector<int>& hard,
    const std::vector<SoftVec>& priors,
    const FieldRowSpace& z_stabilizer_space,
    int inner_length,
    int inner_dimension,
    const std::vector<int>& inner_scales,
    GF2p& gf,
    int sweeps,
    int burn_in_sweeps,
    int sample_every,
    uint32_t seed,
    int block_updates_per_sweep=0,
    int logical_updates_per_sweep=0,
    const std::vector<double>& inverse_temperatures={1.0},
    const SparseMatrix<int, RowMajor>* stabilizer_rows=nullptr,
    int stabilizer_updates_per_sweep=0) {
    const int rows=H.rows(), vars=H.cols(), q=gf.getSize();
    if (syndrome.size()!=static_cast<size_t>(rows) ||
        hard.size()!=static_cast<size_t>(vars) ||
        priors.size()!=static_cast<size_t>(vars) || sweeps<=0 ||
        burn_in_sweeps<0 || burn_in_sweeps>=sweeps || sample_every<=0 ||
        block_updates_per_sweep<0 ||
        logical_updates_per_sweep<0 ||
        stabilizer_updates_per_sweep<0 ||
        q<2 || (q&(q-1))!=0 || inverse_temperatures.empty() ||
        inverse_temperatures.front()!=1.0)
        throw std::invalid_argument("Invalid GF syndrome Gibbs input");
    for (size_t i=0;i<inverse_temperatures.size();++i)
        if (inverse_temperatures[i]<=0.0 || inverse_temperatures[i]>1.0 ||
            (i && inverse_temperatures[i]>=inverse_temperatures[i-1]))
            throw std::invalid_argument("Inverse temperatures must descend from 1 to positive values");
    if(stabilizer_updates_per_sweep>0 &&
       (!stabilizer_rows || stabilizer_rows->cols()!=vars))
        throw std::invalid_argument("Stabilizer row matrix is missing or has wrong width");
    if(logical_updates_per_sweep>0 &&
       (!stabilizer_rows || stabilizer_rows->cols()!=vars))
        throw std::invalid_argument("Logical coordinate moves require the stabilizer matrix");

    std::vector<double> reliability(vars);
    std::vector<int> order(vars), inverse_order(vars);
    std::iota(order.begin(),order.end(),0);
    for(int v=0;v<vars;++v) {
        double first=-1.0, second=-1.0;
        for(double probability:priors[v]) {
            if(probability>first) { second=first; first=probability; }
            else if(probability>second) second=probability;
        }
        reliability[v]=std::log(std::max(first,1e-300))-
                       std::log(std::max(second,1e-300));
    }
    std::sort(order.begin(),order.end(),[&](int a,int b) {
        return reliability[a]>reliability[b];
    });
    for(int pos=0;pos<vars;++pos) inverse_order[order[pos]]=pos;

    std::vector<int> a(static_cast<size_t>(rows)*vars,0), rhs=syndrome;
    for(int r=0;r<rows;++r)
        for(SparseMatrix<int,RowMajor>::InnerIterator it(H,r);it;++it)
            a[static_cast<size_t>(r)*vars+inverse_order[it.col()]]=it.value();
    std::vector<int> multiplication(static_cast<size_t>(q)*q);
    for(int x=0;x<q;++x) for(int y=0;y<q;++y)
        multiplication[static_cast<size_t>(x)*q+y]=gf.mul(x,y);
    std::vector<int> pivot_positions;
    int rank=0;
    for(int col=0;col<vars && rank<rows;++col) {
        int pivot=rank;
        while(pivot<rows && a[static_cast<size_t>(pivot)*vars+col]==0) ++pivot;
        if(pivot==rows) continue;
        if(pivot!=rank) {
            for(int j=0;j<vars;++j)
                std::swap(a[static_cast<size_t>(rank)*vars+j],
                          a[static_cast<size_t>(pivot)*vars+j]);
            std::swap(rhs[rank],rhs[pivot]);
        }
        const int inverse=gf.inv(a[static_cast<size_t>(rank)*vars+col]);
        for(int j=col;j<vars;++j)
            a[static_cast<size_t>(rank)*vars+j]=multiplication[
                static_cast<size_t>(a[static_cast<size_t>(rank)*vars+j])*q+inverse];
        rhs[rank]=multiplication[static_cast<size_t>(rhs[rank])*q+inverse];
        for(int r=0;r<rows;++r) if(r!=rank) {
            const int factor=a[static_cast<size_t>(r)*vars+col];
            if(!factor) continue;
            const size_t row_offset=static_cast<size_t>(r)*vars;
            const size_t pivot_offset=static_cast<size_t>(rank)*vars;
            for(int j=col;j<vars;++j)
                a[row_offset+j]^=multiplication[
                    static_cast<size_t>(factor)*q+a[pivot_offset+j]];
            rhs[r]^=multiplication[static_cast<size_t>(factor)*q+rhs[rank]];
        }
        pivot_positions.push_back(col);
        ++rank;
    }
    for(int r=rank;r<rows;++r) if(rhs[r]!=0) return {};
    std::vector<char> is_pivot(vars,0);
    for(int pos:pivot_positions) is_pivot[pos]=1;
    std::vector<int> free_positions;
    for(int pos=0;pos<vars;++pos) if(!is_pivot[pos]) free_positions.push_back(pos);

    std::vector<int> initial_state(vars,0);
    for(int pos:free_positions) initial_state[order[pos]]=hard[order[pos]];
    for(int r=0;r<rank;++r) {
        int value=rhs[r];
        const size_t offset=static_cast<size_t>(r)*vars;
        for(int pos:free_positions)
            value^=multiplication[static_cast<size_t>(a[offset+pos])*q+
                                  initial_state[order[pos]]];
        initial_state[order[pivot_positions[r]]]=value;
    }
    auto syndrome_valid=[&](const std::vector<int>& word) {
        for(int r=0;r<rows;++r) {
            int value=0;
            for(SparseMatrix<int,RowMajor>::InnerIterator it(H,r);it;++it)
                value=gf.add(value,gf.mul(it.value(),word[it.col()]));
            if(value!=syndrome[r]) return false;
        }
        return true;
    };
    if(!syndrome_valid(initial_state)) return {};

    std::vector<std::vector<std::pair<int,int>>> pivot_updates(free_positions.size());
    for(size_t f=0;f<free_positions.size();++f) {
        const int pos=free_positions[f];
        for(int r=0;r<rank;++r) {
            const int coefficient=a[static_cast<size_t>(r)*vars+pos];
            if(coefficient)
                pivot_updates[f].emplace_back(order[pivot_positions[r]],coefficient);
        }
    }
    std::vector<std::vector<double>> log_prior(vars,std::vector<double>(q));
    for(int v=0;v<vars;++v) for(int value=0;value<q;++value)
        log_prior[v][value]=std::log(std::max(priors[v][value],1e-300));
    std::vector<std::vector<std::pair<int,int>>> stabilizers;
    if(stabilizer_rows) {
        stabilizers.resize(stabilizer_rows->rows());
        const std::string zero_key=cascadeZCosetKey(std::vector<int>(vars,0),
            z_stabilizer_space,gf,inner_length,inner_dimension,inner_scales);
        for(int r=0;r<stabilizer_rows->rows();++r) {
            for(SparseMatrix<int,RowMajor>::InnerIterator it(*stabilizer_rows,r);it;++it)
                stabilizers[r].emplace_back(it.col(),it.value());
            if(!stabilizers[r].empty()) {
                std::vector<int> row_word(vars,0);
                for(const auto& [variable,coefficient]:stabilizers[r])
                    row_word[variable]=coefficient;
                if(cascadeZCosetKey(row_word,z_stabilizer_space,gf,inner_length,
                    inner_dimension,inner_scales)!=zero_key)
                    throw std::runtime_error(
                        "A supplied Z-stabilizer row changes the cascade coset key");
                for(int check=0;check<rows;++check) {
                    int value=0;
                    for(SparseMatrix<int,RowMajor>::InnerIterator it(H,check);it;++it)
                        value=gf.add(value,gf.mul(it.value(),row_word[it.col()]));
                    if(value!=0)
                        throw std::runtime_error(
                            "A supplied Z-stabilizer row does not preserve the syndrome");
                }
            }
        }
        stabilizers.erase(std::remove_if(stabilizers.begin(),stabilizers.end(),
            [](const auto& row){return row.empty();}),stabilizers.end());
    }
    if(stabilizer_updates_per_sweep>0 && stabilizers.empty())
        throw std::invalid_argument("No nonzero stabilizer rows to sample");

    // In the free-coordinate parameterization of ker(H), the restriction of
    // row(H_Z) gives the stabilizer subspace directly. Sparse elimination
    // splits the free coordinates into stabilizer pivots and independent
    // logical coordinates; moves on the latter change logical class.
    std::vector<int> free_index_by_variable(vars,-1);
    for(size_t f=0;f<free_positions.size();++f)
        free_index_by_variable[order[free_positions[f]]]=static_cast<int>(f);
    std::vector<std::map<int,int>> stabilizer_echelon(free_positions.size());
    int stabilizer_rank=0;
    if(stabilizer_rows) {
        for(int r=0;r<stabilizer_rows->rows();++r) {
            std::map<int,int> row;
            for(SparseMatrix<int,RowMajor>::InnerIterator it(*stabilizer_rows,r);it;++it) {
                const int f=free_index_by_variable[it.col()];
                if(f>=0) row[f]=gf.add(row[f],it.value());
            }
            for(auto it=row.begin();it!=row.end();) {
                if(it->second==0) it=row.erase(it); else ++it;
            }
            while(!row.empty()) {
                const int pivot=row.begin()->first;
                if(stabilizer_echelon[pivot].empty()) {
                    const int inverse=gf.inv(row.begin()->second);
                    for(auto& [column,value]:row)
                        value=gf.mul(value,inverse);
                    stabilizer_echelon[pivot]=std::move(row);
                    ++stabilizer_rank;
                    break;
                }
                const int factor=row.begin()->second;
                for(const auto& [column,value]:stabilizer_echelon[pivot]) {
                    const int updated=gf.add(row[column],gf.mul(factor,value));
                    if(updated==0) row.erase(column);
                    else row[column]=updated;
                }
            }
        }
    }
    std::vector<size_t> logical_free_positions;
    logical_free_positions.reserve(free_positions.size()-stabilizer_rank);
    for(size_t f=0;f<free_positions.size();++f)
        if(stabilizer_echelon[f].empty()) logical_free_positions.push_back(f);
    if(stabilizer_rank>static_cast<int>(free_positions.size()))
        throw std::runtime_error("Stabilizer rank exceeds syndrome kernel dimension");
    const int logical_dimension=static_cast<int>(logical_free_positions.size());
    std::vector<GibbsCosetResult::LogicalDirection> logical_directions;
    logical_directions.reserve(logical_free_positions.size());
    for(size_t f:logical_free_positions) {
        GibbsCosetResult::LogicalDirection direction;
        direction.variable=order[free_positions[f]];
        direction.pivot_updates=pivot_updates[f];
        logical_directions.push_back(std::move(direction));
    }

    struct ClassRecord { int count=0; double best_score=-1e300; std::vector<int> word; };
    std::map<std::string,ClassRecord> classes;
    const size_t replica_count=inverse_temperatures.size();
    std::vector<std::vector<int>> states(replica_count,initial_state);
    std::vector<std::mt19937> rngs;
    rngs.reserve(replica_count);
    for(size_t replica=0;replica<replica_count;++replica)
        rngs.emplace_back(seed+static_cast<uint32_t>(replica*0x9e3779b9U));
    std::vector<long long> updates(replica_count,0), changed_updates(replica_count,0);
    long long exchange_attempts=0, exchange_accepts=0;
    std::vector<std::string> cold_sample_keys;
    std::vector<double> cold_sample_scores;
    std::vector<size_t> scan(free_positions.size());
    std::iota(scan.begin(),scan.end(),0);
    auto score=[&](const std::vector<int>& word) {
        double total=0.0;
        for(int v=0;v<vars;++v) total+=log_prior[v][word[v]];
        return total;
    };
    auto heatbath_free_symbol=[&](std::vector<int>& state,size_t f,double beta,
                                  std::mt19937& rng,long long& update_count,
                                  long long& changed_count) {
        const int pos=free_positions[f], variable=order[pos];
        const int old_value=state[variable];
        ++update_count;
        std::vector<double> log_weights(q);
        double maximum=-std::numeric_limits<double>::infinity();
        for(int value=0;value<q;++value) {
            const int delta=gf.add(old_value,value);
            double candidate_score=log_prior[variable][value];
            for(const auto& [pivot_var,coefficient]:pivot_updates[f]) {
                const int changed=gf.add(state[pivot_var],gf.mul(coefficient,delta));
                candidate_score+=log_prior[pivot_var][changed];
            }
            log_weights[value]=beta*candidate_score;
            maximum=std::max(maximum,log_weights[value]);
        }
        double total=0.0;
        for(double& value:log_weights) {
            value=std::exp(std::max(-700.0,value-maximum));
            total+=value;
        }
        std::uniform_real_distribution<double> uniform(0.0,1.0);
        double draw=uniform(rng)*total;
        int new_value=q-1;
        for(int value=0;value<q;++value) {
            draw-=log_weights[value];
            if(draw<=0.0) { new_value=value; break; }
        }
        const int delta=gf.add(old_value,new_value);
        if(delta!=0) {
            ++changed_count;
            state[variable]=new_value;
            for(const auto& [pivot_var,coefficient]:pivot_updates[f])
                state[pivot_var]=gf.add(state[pivot_var],gf.mul(coefficient,delta));
        }
    };
    auto update_replica=[&](std::vector<int>& state, double beta,
                            std::mt19937& rng, long long& update_count,
                            long long& changed_count) {
        std::uniform_real_distribution<double> uniform(0.0,1.0);
        std::vector<size_t> replica_scan=scan;
        std::shuffle(replica_scan.begin(),replica_scan.end(),rng);
        if (block_updates_per_sweep>0) {
            const int block_count=std::min<int>(block_updates_per_sweep,replica_scan.size());
            for(int index=0;index<block_count;++index)
                heatbath_free_symbol(state,replica_scan[index],beta,rng,
                                     update_count,changed_count);
        } else for(size_t f:replica_scan) {
            const int pos=free_positions[f], variable=order[pos];
            for(int bit=0;bit<static_cast<int>(std::log2(q));++bit) {
                ++update_count;
                const int delta=1<<bit;
                const int old_free=state[variable], new_free=gf.add(old_free,delta);
                double old_score=log_prior[variable][old_free];
                double new_score=log_prior[variable][new_free];
                for(const auto& [pivot_var,coefficient]:pivot_updates[f]) {
                    const int old_value=state[pivot_var];
                    const int new_value=gf.add(old_value,gf.mul(coefficient,delta));
                    old_score+=log_prior[pivot_var][old_value];
                    new_score+=log_prior[pivot_var][new_value];
                }
                const double probability_new=1.0/(1.0+std::exp(
                    std::clamp(beta*(old_score-new_score),-700.0,700.0)));
                if(uniform(rng)<probability_new) {
                    ++changed_count;
                    state[variable]=new_free;
                    for(const auto& [pivot_var,coefficient]:pivot_updates[f])
                        state[pivot_var]=gf.add(state[pivot_var],
                                                gf.mul(coefficient,delta));
                }
            }
        }
        if(stabilizer_updates_per_sweep>0) {
            std::uniform_int_distribution<size_t> choose_row(0,stabilizers.size()-1);
            for(int update=0;update<stabilizer_updates_per_sweep;++update) {
                const auto& row=stabilizers[choose_row(rng)];
                ++update_count;
                std::vector<double> log_weights(q);
                double maximum=-std::numeric_limits<double>::infinity();
                for(int alpha=0;alpha<q;++alpha) {
                    double candidate_score=0.0;
                    for(const auto& [variable,coefficient]:row) {
                        const int value=gf.add(state[variable],gf.mul(alpha,coefficient));
                        candidate_score+=log_prior[variable][value];
                    }
                    log_weights[alpha]=beta*candidate_score;
                    maximum=std::max(maximum,log_weights[alpha]);
                }
                double total=0.0;
                for(double& value:log_weights) {
                    value=std::exp(std::max(-700.0,value-maximum));
                    total+=value;
                }
                std::uniform_real_distribution<double> uniform(0.0,1.0);
                double draw=uniform(rng)*total;
                int alpha=q-1;
                for(int value=0;value<q;++value) {
                    draw-=log_weights[value];
                    if(draw<=0.0) { alpha=value; break; }
                }
                if(alpha!=0) {
                    ++changed_count;
                    for(const auto& [variable,coefficient]:row)
                        state[variable]=gf.add(state[variable],gf.mul(alpha,coefficient));
                }
            }
        }
        if(logical_updates_per_sweep>0) {
            if(logical_free_positions.empty())
                throw std::runtime_error("No logical coordinate moves are available");
            std::uniform_int_distribution<size_t> choose_logical(
                0,logical_free_positions.size()-1);
            for(int update=0;update<logical_updates_per_sweep;++update)
                heatbath_free_symbol(state,
                    logical_free_positions[choose_logical(rng)],beta,rng,
                    update_count,changed_count);
        }
    };
    for(int sweep=0;sweep<sweeps;++sweep) {
        for(size_t replica=0;replica<replica_count;++replica)
            update_replica(states[replica],inverse_temperatures[replica],
                rngs[replica],updates[replica],changed_updates[replica]);
        const size_t parity=static_cast<size_t>(sweep%2);
        for(size_t left=parity;left+1<replica_count;left+=2) {
            ++exchange_attempts;
            const size_t right=left+1;
            const double log_accept=(inverse_temperatures[left]-
                inverse_temperatures[right])*(score(states[right])-score(states[left]));
            std::uniform_real_distribution<double> uniform(0.0,1.0);
            if(std::log(std::max(uniform(rngs[left]),1e-300))<
               std::min(0.0,log_accept)) {
                std::swap(states[left],states[right]);
                ++exchange_accepts;
            }
        }
        if(sweep>=burn_in_sweeps && (sweep-burn_in_sweeps)%sample_every==0) {
            const auto& cold_state=states.front();
            const std::string key=cascadeZCosetKey(cold_state,z_stabilizer_space,gf,
                inner_length,inner_dimension,inner_scales);
            auto& record=classes[key];
            ++record.count;
            const double value=score(cold_state);
            cold_sample_keys.push_back(key);
            cold_sample_scores.push_back(value);
            if(value>record.best_score) {
                record.best_score=value;
                record.word=cold_state;
            }
        }
    }
    GibbsCosetResult result;
    result.kernel_dimension=static_cast<int>(free_positions.size());
    result.stabilizer_rank=stabilizer_rank;
    result.logical_dimension=logical_dimension;
    result.logical_directions=std::move(logical_directions);
    result.free_variables.reserve(free_positions.size());
    for(int pos:free_positions) result.free_variables.push_back(order[pos]);
    result.stabilizer_echelon.resize(stabilizer_echelon.size());
    for(size_t pivot=0;pivot<stabilizer_echelon.size();++pivot)
        for(const auto& [column,value]:stabilizer_echelon[pivot])
            result.stabilizer_echelon[pivot].emplace_back(column,value);
    result.logical_updates=static_cast<long long>(logical_updates_per_sweep)*
        sweeps*static_cast<long long>(replica_count);
    result.samples=sweeps-burn_in_sweeps==0 ? 0 :
        1+(sweeps-burn_in_sweeps-1)/sample_every;
    result.updates=std::accumulate(updates.begin(),updates.end(),0LL);
    result.changed_updates=std::accumulate(changed_updates.begin(),changed_updates.end(),0LL);
    result.exchange_attempts=exchange_attempts;
    result.exchange_accepts=exchange_accepts;
    if(cold_sample_keys.size()>1) {
        size_t same_class=0;
        for(size_t i=1;i<cold_sample_keys.size();++i)
            same_class+=cold_sample_keys[i]==cold_sample_keys[i-1];
        result.class_persistence=static_cast<double>(same_class)/
                                 (cold_sample_keys.size()-1);
        double mean=std::accumulate(cold_sample_scores.begin(),
            cold_sample_scores.end(),0.0)/cold_sample_scores.size();
        double numerator=0.0,denominator=0.0;
        for(double value:cold_sample_scores) denominator+=(value-mean)*(value-mean);
        for(size_t i=1;i<cold_sample_scores.size();++i)
            numerator+=(cold_sample_scores[i-1]-mean)*(cold_sample_scores[i]-mean);
        if(denominator>0.0)
            result.score_lag1_correlation=numerator/denominator;
    }
    for(size_t replica=0;replica<replica_count;++replica)
        result.replica_change_rates.push_back(updates[replica] ?
            static_cast<double>(changed_updates[replica])/updates[replica] : 0.0);
    result.distinct_classes=static_cast<int>(classes.size());
    for(const auto& [key,record]:classes) {
        result.sampled_class_keys.push_back(key);
        if(record.count>result.winning_class_samples ||
           (record.count==result.winning_class_samples &&
            record.best_score>result.representative_log_posterior)) {
            result.winning_class_samples=record.count;
            result.estimate=record.word;
            result.representative_log_posterior=record.best_score;
        }
    }
    result.syndrome_ok=syndrome_valid(result.estimate);
    return result;
}

// Turbo-style outer/inner BP for [1|I_k], adjacent-pair G2, and the
// coefficient-scaled adjacent-pair basis. The two-group path2 code is
// decoded as two independent six-symbol chains before this core routine.
std::vector<int> decodeHxCascadeTurbo(
    NB_BP_Decoder& decoder,
    const SparseMatrix<int, RowMajor>& outer_h,
    const SparseMatrix<int, RowMajor>& hx_cascade,
    const std::vector<int>& syndrome,
    double p_depolar,
    const std::vector<int>* x_estimate,
    const FieldBasis& basis,
    GF2p& gf,
    int turbo_rounds=10,
    int outer_iterations=10,
    double damping=0.5,
    bool stop_on_valid=true,
    std::vector<SoftVec>* symbol_posterior_out=nullptr,
    bool projection_repair=false,
    ProjectionDiagnostics* projection_diagnostics=nullptr,
    bool joint_candidate_selection=false,
    double* joint_score_out=nullptr,
    bool refresh_outer_at_end=false,
    bool joint_list_selection=false,
    bool independent_outer=false,
    const std::vector<std::pair<int,int>>* physical_bit_constraints=nullptr,
    bool cavity_symbol_posterior=false,
    const std::vector<SoftVec>* z_symbol_priors=nullptr,
    const std::vector<int>* duplicated_outer_checks=nullptr,
    bool layered_outer_schedule=false,
    const std::vector<RedundantOuterCheck>* redundant_outer_checks=nullptr,
    const SparseMatrix<int,RowMajor>* dual_outer_checks=nullptr,
    int max_inactivation_checks=0,
    ProjectionDiagnostics* inactivation_diagnostics=nullptr) {
    const int q=gf.getSize(), n1=outer_h.cols(), m1=outer_h.rows();
    if (n1 <= 0 || m1 <= 0 || hx_cascade.cols() % n1 != 0)
        throw std::invalid_argument("Invalid outer dimensions for structured decoder");
    const int n2=hx_cascade.cols()/n1;
    double constraint_penalty_exponent=25.0;
    if (const char* value=std::getenv(
            "GF256_CASCADE_CONSTRAINT_PENALTY_EXPONENT"))
        constraint_penalty_exponent=std::stod(value);
    if (!std::isfinite(constraint_penalty_exponent) ||
        constraint_penalty_exponent<=0.0)
        throw std::invalid_argument("Constraint penalty exponent must be positive");
    if (physical_bit_constraints && physical_bit_constraints->size() !=
            static_cast<size_t>(n1*n2))
        throw std::invalid_argument("Invalid physical constraint size");
    if (physical_bit_constraints && hx_cascade.rows() != m1*(n2-1))
        throw std::invalid_argument("Guided constraints currently require one inner check");
    // The path2 inner code consists of two independent six-symbol paths.
    // Splitting here is exact: each outer syndrome slice touches only one
    // group. Reuse the same turbo factor for each path, then reassemble the
    // physical symbols in the original order.
    if (n2 == 12 && hx_cascade.rows() == m1 * 10) {
        if (syndrome.size() != static_cast<size_t>(m1 * 10))
            throw std::invalid_argument("Invalid path2 syndrome dimensions");
        std::vector<int> combined(n1 * n2, 0);
        if (symbol_posterior_out)
            symbol_posterior_out->assign(n1 * n2, SoftVec(q));
        double combined_score = 0.0;
        for (int group = 0; group < 2; ++group) {
            const auto physical_col = [group](int local) {
                return local == 5 ? 10 + group : 5 * group + local;
            };
            std::vector<Eigen::Triplet<int>> entries;
            for (int row = 0; row < m1 * 5; ++row) {
                const int full_row = (row / 5) * 10 + group * 5 + row % 5;
                for (SparseMatrix<int,RowMajor>::InnerIterator it(hx_cascade,full_row);
                     it; ++it) {
                    const int block = it.col() / n2;
                    const int inner_col = it.col() % n2;
                    for (int local = 0; local < 6; ++local)
                        if (physical_col(local) == inner_col)
                            entries.emplace_back(row, block * 6 + local, it.value());
                }
            }
            SparseMatrix<int,RowMajor> sub_hx(m1 * 5, n1 * 6);
            sub_hx.setFromTriplets(entries.begin(), entries.end());
            if (sub_hx.nonZeros() * 2 !=
                hx_cascade.nonZeros())
                throw std::invalid_argument("Unsupported path2 HX support");
            std::vector<int> sub_syndrome(m1 * 5), sub_x;
            std::vector<std::pair<int,int>> sub_constraints;
            for (int row = 0; row < m1 * 5; ++row)
                sub_syndrome[row] = syndrome[(row / 5) * 10 + group * 5 + row % 5];
            if (physical_bit_constraints) {
                sub_constraints.resize(n1 * 6);
                for (int block = 0; block < n1; ++block)
                    for (int local = 0; local < 6; ++local)
                        sub_constraints[block * 6 + local] =
                            (*physical_bit_constraints)[block * n2 +
                                physical_col(local)];
            }
            if (x_estimate) {
                sub_x.resize(n1 * 6);
                for (int block = 0; block < n1; ++block)
                    for (int local = 0; local < 6; ++local)
                        sub_x[block * 6 + local] =
                            (*x_estimate)[block * n2 + physical_col(local)];
            }
            std::vector<SoftVec> sub_z_priors;
            if(z_symbol_priors) {
                sub_z_priors.resize(static_cast<size_t>(n1*6),SoftVec(q));
                for(int block=0;block<n1;++block)
                    for(int local=0;local<6;++local)
                        sub_z_priors[block*6+local]=
                            (*z_symbol_priors)[block*n2+physical_col(local)];
            }
            std::vector<SoftVec> sub_posterior;
            double sub_score = 0.0;
            auto decoded = decodeHxCascadeTurbo(
                decoder, outer_h, sub_hx, sub_syndrome, p_depolar,
                x_estimate ? &sub_x : nullptr, basis, gf, turbo_rounds,
                outer_iterations, damping, stop_on_valid,
                symbol_posterior_out ? &sub_posterior : nullptr,
                projection_repair, nullptr, joint_candidate_selection,
                joint_score_out ? &sub_score : nullptr, refresh_outer_at_end,
                joint_list_selection, independent_outer,
                physical_bit_constraints ? &sub_constraints : nullptr,
                cavity_symbol_posterior,
                z_symbol_priors ? &sub_z_priors : nullptr,
                duplicated_outer_checks,layered_outer_schedule,
                redundant_outer_checks,dual_outer_checks,
                max_inactivation_checks,inactivation_diagnostics);
            combined_score += sub_score;
            for (int block = 0; block < n1; ++block)
                for (int local = 0; local < 6; ++local) {
                    const int destination = block * n2 + physical_col(local);
                    combined[destination] = decoded[block * 6 + local];
                    if (symbol_posterior_out)
                        (*symbol_posterior_out)[destination] =
                            std::move(sub_posterior[block * 6 + local]);
                }
        }
        if (joint_score_out) *joint_score_out = combined_score;
        if (projection_diagnostics)
            projection_diagnostics->final_valid =
                decoder.check_z_syndrome(hx_cascade, combined, syndrome);
        return combined;
    }
    if(hx_cascade.rows()%m1!=0)
        throw std::invalid_argument("Cascade HX row count is not divisible by outer checks");
    const int k2=hx_cascade.rows()/m1;
    if(hx_cascade.cols()!=n1*n2 || k2<=0 || k2>=n2 ||
       syndrome.size()!=static_cast<size_t>(m1*k2))
        throw std::invalid_argument("Invalid dimensions for structured turbo decoder");
    // Refuse any unsupported G2 instead of applying the wrong local factor.
    const std::vector<int> unit_scales(n2,1);
    std::vector<int> weighted_scales(n2,1);
    for (int j=0;j<n2;++j) weighted_scales[j]=gf.exp(j);
    const auto matches_generator=[&](bool path,
                                      const std::vector<int>& candidate_scales) {
        for (int r=0; r<m1; ++r) for (int a=0; a<k2; ++a) {
            int expected_entries=0;
            for (SparseMatrix<int, RowMajor>::InnerIterator it(outer_h,r); it; ++it) {
                const int left=path ? a : 0;
                if (hx_cascade.coeff(r*k2+a,it.col()*n2+left) !=
                        gf.mul(it.value(),candidate_scales[left]) ||
                    hx_cascade.coeff(r*k2+a,it.col()*n2+a+1) !=
                        gf.mul(it.value(),candidate_scales[a+1]))
                    return false;
                expected_entries+=2;
            }
            int actual_entries=0;
            for (SparseMatrix<int, RowMajor>::InnerIterator it(hx_cascade,r*k2+a); it; ++it)
                ++actual_entries;
            if (actual_entries != expected_entries) return false;
        }
        return true;
    };
    const bool plain_path=matches_generator(true,unit_scales);
    const bool weighted_path=!plain_path &&
        matches_generator(true,weighted_scales);
    const bool path_generator=plain_path || weighted_path;
    const bool star_generator=!path_generator &&
        matches_generator(false,unit_scales);
    const bool generic_generator=!path_generator && !star_generator;
    if (generic_generator && (physical_bit_constraints || joint_candidate_selection ||
            joint_list_selection || independent_outer || projection_repair ||
            duplicated_outer_checks || redundant_outer_checks || dual_outer_checks ||
            max_inactivation_checks>0))
        throw std::invalid_argument(
            "Generic G2 Turbo supports the basic structured decoder options only");
    SparseMatrix<int,RowMajor> outer_bp_h;
    std::vector<int> duplicate_rows;
    if(duplicated_outer_checks || redundant_outer_checks) {
        if(duplicated_outer_checks) {
        duplicate_rows=*duplicated_outer_checks;
        for(int row:duplicate_rows)
            if(row<0 || row>=m1)
                throw std::invalid_argument("Augmented outer-check index is out of range");
        }
        if(redundant_outer_checks)
            for(const auto& check:*redundant_outer_checks)
                if(check.first_row<0 || check.first_row>=m1 ||
                   check.second_row<0 || check.second_row>=m1 ||
                   check.second_scale==0 || check.second_scale>=q)
                    throw std::invalid_argument("Invalid redundant outer-check transform");
        const size_t extra_rows=duplicate_rows.size()+
            (redundant_outer_checks ? redundant_outer_checks->size() : 0);
        std::vector<Eigen::Triplet<int>> augmented_entries;
        augmented_entries.reserve(outer_h.nonZeros()*(1+extra_rows/std::max(1,m1)*2));
        for(int r=0;r<m1;++r)
            for(SparseMatrix<int,RowMajor>::InnerIterator it(outer_h,r);it;++it)
                augmented_entries.emplace_back(r,it.col(),it.value());
        int destination=m1;
        for(int source:duplicate_rows) {
            for(SparseMatrix<int,RowMajor>::InnerIterator it(outer_h,source);it;++it)
                augmented_entries.emplace_back(destination,it.col(),it.value());
            ++destination;
        }
        if(redundant_outer_checks) for(const auto& check:*redundant_outer_checks) {
            std::map<int,int> combined;
            for(SparseMatrix<int,RowMajor>::InnerIterator it(outer_h,check.first_row);it;++it)
                combined[it.col()]=it.value();
            for(SparseMatrix<int,RowMajor>::InnerIterator it(outer_h,check.second_row);it;++it) {
                const int scaled=gf.mul(check.second_scale,it.value());
                const auto found=combined.find(it.col());
                if(found==combined.end()) combined.emplace(it.col(),scaled);
                else {
                    found->second=gf.add(found->second,scaled);
                    if(found->second==0) combined.erase(found);
                }
            }
            for(const auto& [col,value]:combined)
                augmented_entries.emplace_back(destination,col,value);
            ++destination;
        }
        outer_bp_h.resize(m1+static_cast<int>(extra_rows),n1);
        outer_bp_h.setFromTriplets(augmented_entries.begin(),augmented_entries.end());
        outer_bp_h.makeCompressed();
    } else {
        outer_bp_h=outer_h;
    }
    const auto& scales=weighted_path ? weighted_scales : unit_scales;
    std::vector<int> inverse_scales(n2);
    for (int j=0;j<n2;++j) inverse_scales[j]=gf.inv(scales[j]);
    if (p_depolar <= 0.0 || p_depolar >= 1.0 || turbo_rounds <= 0 ||
        outer_iterations <= 0 || damping <= 0.0 || damping > 1.0)
        throw std::invalid_argument("Invalid channel or BP parameters");
    if (joint_candidate_selection && !projection_repair)
        throw std::invalid_argument("Joint outer candidates require projection repair");
    if (x_estimate && x_estimate->size() != static_cast<size_t>(n1*n2))
        throw std::invalid_argument("Invalid X estimate for conditional Z prior");
    if (z_symbol_priors && z_symbol_priors->size()!=static_cast<size_t>(n1*n2))
        throw std::invalid_argument("Invalid soft Z-prior dimensions");
    SoftMatrix channel(n1*n2, std::vector<SoftVec>(1, SoftVec(q)));
    for (int v=0; v<n1*n2; ++v) {
        const int xbits=x_estimate ? basis.toDual((*x_estimate)[v]) : 0;
        for(int z=0; z<q; ++z) {
            const int physical_z=gf.mul(z,inverse_scales[v%n2]);
            double probability=0.0;
            if(z_symbol_priors) {
                if((*z_symbol_priors)[v].size()!=static_cast<size_t>(q))
                    throw std::invalid_argument("Invalid soft Z-prior alphabet size");
                probability=(*z_symbol_priors)[v][physical_z];
            } else {
                probability=1.0;
                for(int b=0; b<gf.getPower(); ++b) {
                    const double pz=x_estimate
                        ? pauli_channel.zGivenX((xbits >> b)&1) : pauli_channel.marginalZ();
                    probability*=((physical_z>>b)&1) ? pz : 1.0-pz;
                }
            }
            if (physical_bit_constraints) {
                const auto [mask,value]=(*physical_bit_constraints)[v];
                // Keep constrained symbols numerically supported through the
                // FWHT while strongly suppressing values that violate a branch.
                if ((physical_z & mask) != value)
                    probability*=std::exp(-constraint_penalty_exponent);
            }
            channel[v][0][z]=probability;
        }
        normalizeLocal(channel[v][0]);
    }
    SoftMatrix inner(k2,std::vector<SoftVec>(n1,SoftVec(q,1.0/q)));
    if (!generic_generator) {
        for(int block=0; block<n1; ++block)
            for(int a=0; a<k2; ++a)
                inner[a][block]=xorConvolution(
                    channel[block*n2+(path_generator ? a : 0)][0],
                    channel[block*n2+a+1][0]);
    }
    std::vector<OuterSoftResult> outer(k2);
    std::vector<std::vector<int>> outer_syn(k2,
        std::vector<int>(m1+duplicate_rows.size()+
            (redundant_outer_checks ? redundant_outer_checks->size() : 0)));
    for(int a=0;a<k2;++a) {
        for(int r=0;r<m1;++r) outer_syn[a][r]=syndrome[r*k2+a];
        for(size_t copy=0;copy<duplicate_rows.size();++copy)
            outer_syn[a][m1+copy]=syndrome[duplicate_rows[copy]*k2+a];
        if(redundant_outer_checks)
            for(size_t index=0;index<redundant_outer_checks->size();++index) {
                const auto& check=(*redundant_outer_checks)[index];
                outer_syn[a][m1+duplicate_rows.size()+index]=gf.add(
                    syndrome[check.first_row*k2+a],
                    gf.mul(check.second_scale,
                           syndrome[check.second_row*k2+a]));
            }
    }
    if (generic_generator) {
        // Recover G2 from H_X = H_outer kron G2 and validate the Kronecker
        // factorization.  The inner local factor graph represents
        // y + G2 e = 0; y is decoded by the outer code while e carries the
        // physical channel prior.  This supports sparse systematic G2 rows
        // beyond the historical star/path special cases.
        int reference_row=-1, reference_block=-1, reference_value=0;
        for(int r=0;r<m1 && reference_row<0;++r)
            for(SparseMatrix<int,RowMajor>::InnerIterator it(outer_h,r);it;++it)
                if(it.value()!=0) {
                    reference_row=r; reference_block=it.col();
                    reference_value=it.value(); break;
                }
        if(reference_row<0)
            throw std::invalid_argument("Cannot infer G2 from a zero outer matrix");
        const int inverse_reference=gf.inv(reference_value);
        std::vector<int> g2(static_cast<size_t>(k2*n2),0);
        for(int a=0;a<k2;++a) for(int c=0;c<n2;++c)
            g2[static_cast<size_t>(a*n2+c)]=gf.mul(inverse_reference,
                hx_cascade.coeff(reference_row*k2+a,reference_block*n2+c));
        for(int r=0;r<m1;++r) for(int block=0;block<n1;++block)
            for(int a=0;a<k2;++a) for(int c=0;c<n2;++c) {
                const int expected=gf.mul(outer_h.coeff(r,block),
                    g2[static_cast<size_t>(a*n2+c)]);
                if(hx_cascade.coeff(r*k2+a,block*n2+c)!=expected)
                    throw std::invalid_argument("HX is not the expected outer_h kron G2");
            }
        std::vector<Eigen::Triplet<int>> local_entries;
        for(int a=0;a<k2;++a) {
            for(int c=0;c<n2;++c) {
                const int value=g2[static_cast<size_t>(a*n2+c)];
                if(value) local_entries.emplace_back(a,c,value);
            }
            local_entries.emplace_back(a,n2+a,1);
        }
        SparseMatrix<int,RowMajor> local_constraints(k2,n2+k2);
        local_constraints.setFromTriplets(local_entries.begin(),local_entries.end());
        local_constraints.makeCompressed();
        const std::vector<int> zero_inner_syndrome(k2,0);
        int local_iterations=5;
        if(const char* value=std::getenv("GF256_CASCADE_GENERIC_INNER_ITERS"))
            local_iterations=std::stoi(value);
        if(local_iterations<1 || local_iterations>50)
            throw std::invalid_argument(
                "GF256_CASCADE_GENERIC_INNER_ITERS must be in [1,50]");
        double generic_damping=damping;
        if(const char* value=std::getenv("GF256_CASCADE_GENERIC_DAMPING"))
            generic_damping=std::stod(value);
        if(!std::isfinite(generic_damping) || generic_damping<=0.0 ||
           generic_damping>1.0)
            throw std::invalid_argument(
                "GF256_CASCADE_GENERIC_DAMPING must be in (0,1]");

        auto run_inner_block=[&](int block, const SoftMatrix* y_priors,
                                 std::vector<SoftVec>* physical_posterior) {
            std::vector<SoftVec> local_priors(static_cast<size_t>(n2+k2),
                                               SoftVec(q,1.0/q));
            for(int c=0;c<n2;++c)
                local_priors[c]=channel[block*n2+c][0];
            if(y_priors)
                for(int a=0;a<k2;++a)
                    local_priors[n2+a]=(*y_priors)[a][block];
            const auto local=runOuterSoftBP(decoder,local_constraints,
                zero_inner_syndrome,local_priors,local_iterations,gf);
            if(physical_posterior)
                for(int c=0;c<n2;++c)
                    (*physical_posterior)[c]=local.posterior[c];
            return local;
        };

        const bool generic_trace=
            std::getenv("GF256_CASCADE_GENERIC_TRACE")!=nullptr;

        // Initialize outer-variable priors using the inner factor and the
        // channel only, then exchange extrinsic distributions by Turbo rounds.
        for(int block=0;block<n1;++block) {
            const auto local=run_inner_block(block,nullptr,nullptr);
            for(int a=0;a<k2;++a)
                inner[a][block]=local.extrinsic[n2+a][0];
        }
        std::vector<OuterSoftResult> outer(k2);
        std::vector<int> outer_syndrome(m1);
        std::vector<std::vector<SoftVec>> previous_outer_extrinsic;
        for(int round=0;round<turbo_rounds;++round) {
            for(int a=0;a<k2;++a) {
                for(int r=0;r<m1;++r) outer_syndrome[r]=syndrome[r*k2+a];
                outer[a]=runOuterSoftBP(decoder,outer_bp_h,outer_syndrome,
                    inner[a],outer_iterations,gf,layered_outer_schedule);
            }
            std::vector<int> projection_unsatisfied;
            if(generic_trace) {
                projection_unsatisfied.assign(k2,0);
                for(int a=0;a<k2;++a) for(int r=0;r<m1;++r) {
                    int check_value=0;
                    for(SparseMatrix<int,RowMajor>::InnerIterator it(outer_h,r);
                        it;++it)
                        check_value=gf.add(check_value,gf.mul(it.value(),
                            outer[a].hard[it.col()]));
                    projection_unsatisfied[a]+=
                        check_value!=syndrome[r*k2+a];
                }
            }
            double outer_tv_max=0.0,outer_tv_sum=0.0;
            int outer_tv_count=0;
            if(generic_trace && previous_outer_extrinsic.size()==
                    static_cast<size_t>(k2)) {
                for(int a=0;a<k2;++a) for(int block=0;block<n1;++block) {
                    double tv=0.0;
                    for(int x=0;x<q;++x)
                        tv+=std::abs(outer[a].extrinsic[block][0][x]-
                                     previous_outer_extrinsic[a][block][x]);
                    tv*=0.5;
                    outer_tv_max=std::max(outer_tv_max,tv);
                    outer_tv_sum+=tv;
                    ++outer_tv_count;
                }
            }
            SoftMatrix y_priors(k2,std::vector<SoftVec>(n1));
            for(int a=0;a<k2;++a)
                for(int block=0;block<n1;++block)
                    y_priors[a][block]=outer[a].extrinsic[block][0];
            if(generic_trace) {
                previous_outer_extrinsic.assign(k2,std::vector<SoftVec>(n1));
                for(int a=0;a<k2;++a) for(int block=0;block<n1;++block)
                    previous_outer_extrinsic[a][block]=outer[a].extrinsic[block][0];
            }
            double inner_tv_max=0.0,inner_tv_sum=0.0,inner_entropy_sum=0.0;
            int inner_tv_count=0;
            for(int block=0;block<n1;++block) {
                const auto local=run_inner_block(block,&y_priors,nullptr);
                for(int a=0;a<k2;++a) {
                    // Match the path/star Turbo schedule: damping is the
                    // fraction of the newly proposed local extrinsic message.
                    SoftVec proposed=local.extrinsic[n2+a][0];
                    for(int x=0;x<q;++x)
                        proposed[x]=generic_damping*proposed[x]+
                            (1.0-generic_damping)*inner[a][block][x];
                    normalizeLocal(proposed);
                    if(generic_trace) {
                        double tv=0.0,entropy=0.0;
                        for(int x=0;x<q;++x) {
                            tv+=std::abs(proposed[x]-inner[a][block][x]);
                            if(proposed[x]>0.0)
                                entropy-=proposed[x]*std::log(proposed[x]);
                        }
                        tv*=0.5;
                        inner_tv_max=std::max(inner_tv_max,tv);
                        inner_tv_sum+=tv;
                        inner_entropy_sum+=entropy;
                        ++inner_tv_count;
                    }
                    inner[a][block]=std::move(proposed);
                }
            }
            if(generic_trace) {
                double outer_entropy_sum=0.0;
                for(int a=0;a<k2;++a) for(int block=0;block<n1;++block)
                    for(double probability:outer[a].extrinsic[block][0])
                        if(probability>0.0)
                            outer_entropy_sum-=probability*std::log(probability);
                std::cerr << "[generic-g2-trace] round=" << (round+1)
                    << " outer_tv_mean="
                    << (outer_tv_count ? outer_tv_sum/outer_tv_count : -1.0)
                    << " outer_tv_max=" << outer_tv_max
                    << " inner_tv_mean="
                    << (inner_tv_count ? inner_tv_sum/inner_tv_count : 0.0)
                    << " inner_tv_max=" << inner_tv_max
                    << " outer_entropy_mean="
                    << outer_entropy_sum/(k2*n1)
                    << " inner_entropy_mean="
                    << (inner_tv_count ? inner_entropy_sum/inner_tv_count : 0.0)
                    << " projection_unsatisfied=[";
                for(int a=0;a<k2;++a) {
                    if(a) std::cerr << ',';
                    std::cerr << projection_unsatisfied[a];
                }
                std::cerr << "]\n";
            }
        }
        // Refresh the outer projections after the final inner update, then
        // form symbol APPs from the complete local inner factor graph.
        for(int a=0;a<k2;++a) {
            for(int r=0;r<m1;++r) outer_syndrome[r]=syndrome[r*k2+a];
            outer[a]=runOuterSoftBP(decoder,outer_bp_h,outer_syndrome,
                inner[a],outer_iterations,gf,layered_outer_schedule);
        }
        if(projection_diagnostics && std::getenv(
                "GF256_CASCADE_GENERIC_CANDIDATE_AUDIT")) {
            int audit_free_symbols=32,audit_order=1,audit_alternatives=4;
            if(const char* value=std::getenv(
                    "GF256_CASCADE_GENERIC_CANDIDATE_FREE"))
                audit_free_symbols=std::stoi(value);
            if(const char* value=std::getenv(
                    "GF256_CASCADE_GENERIC_CANDIDATE_ORDER"))
                audit_order=std::stoi(value);
            if(const char* value=std::getenv(
                    "GF256_CASCADE_GENERIC_CANDIDATE_ALTERNATIVES"))
                audit_alternatives=std::stoi(value);
            if(audit_free_symbols<1 || audit_free_symbols>256 ||
               audit_order<1 || audit_order>2 ||
               audit_alternatives<1 || audit_alternatives>8)
                throw std::invalid_argument(
                    "Generic candidate audit requires FREE in [1,256], ORDER in [1,2], and ALTERNATIVES in [1,8]");
            projection_diagnostics->outer_candidate_list.clear();
            projection_diagnostics->outer_candidate_list.reserve(k2);
            for(int a=0;a<k2;++a) {
                std::vector<int> projection_syndrome(m1);
                for(int r=0;r<m1;++r)
                    projection_syndrome[r]=syndrome[r*k2+a];
                std::vector<std::vector<int>> candidates;
                (void)fullMatrixOsd1(outer_h,projection_syndrome,
                    outer[a].hard,outer[a].posterior,gf,
                    std::min(n1,audit_free_symbols),audit_order,
                    audit_alternatives,true,&candidates);
                projection_diagnostics->outer_candidate_list.push_back(
                    std::move(candidates));
            }
        }
        const char* exact_recovery_env=
            std::getenv("GF256_CASCADE_GENERIC_EXACT_RECOVERY");
        const bool use_exact_projected_recovery=exact_recovery_env &&
            std::string(exact_recovery_env)=="1";
        if(!use_exact_projected_recovery) {
            // Preserve the original Turbo APP hard decision as the default;
            // exact projected recovery is an opt-in experiment because
            // syndrome feasibility alone does not choose the logical class.
            std::vector<int> estimate(static_cast<size_t>(n1*n2),0);
            if(symbol_posterior_out)
                symbol_posterior_out->assign(static_cast<size_t>(n1*n2),SoftVec(q));
            for(int block=0;block<n1;++block) {
                SoftMatrix y_priors(k2,std::vector<SoftVec>(n1));
                for(int a=0;a<k2;++a)
                    for(int j=0;j<n1;++j)
                        y_priors[a][j]=outer[a].extrinsic[j][0];
                std::vector<SoftVec> physical_posterior(n2,SoftVec(q));
                (void)run_inner_block(block,&y_priors,&physical_posterior);
                for(int c=0;c<n2;++c) {
                    estimate[block*n2+c]=std::distance(
                        physical_posterior[c].begin(),
                        std::max_element(physical_posterior[c].begin(),
                                         physical_posterior[c].end()));
                    if(symbol_posterior_out)
                        (*symbol_posterior_out)[block*n2+c]=
                            std::move(physical_posterior[c]);
                }
            }
            if(generic_trace) {
                std::vector<int> hard_consistency_mismatches(
                    static_cast<size_t>(k2),0);
                for(int block=0;block<n1;++block) for(int a=0;a<k2;++a) {
                    int projected_symbol=0;
                    for(int c=0;c<n2;++c)
                        projected_symbol=gf.add(projected_symbol,gf.mul(
                            g2[static_cast<size_t>(a*n2+c)],
                            estimate[block*n2+c]));
                    hard_consistency_mismatches[a]+=
                        projected_symbol!=outer[a].hard[block];
                }
                std::cerr << "[generic-g2-hard-consistency] mismatches_by_projection=[";
                for(int a=0;a<k2;++a) {
                    if(a) std::cerr << ',';
                    std::cerr << hard_consistency_mismatches[a];
                }
                std::cerr << "]\n";
            }
            if(projection_diagnostics)
                projection_diagnostics->final_valid=
                    decoder.check_z_syndrome(hx_cascade,estimate,syndrome);
            return estimate;
        }
        // The Turbo APP hard decisions above need not satisfy the coupling
        // equations jointly: independently hard-decoding e and y can break
        // G2 e = y even when both soft decoders look locally confident.
        // First project each outer coordinate onto its exact syndrome coset,
        // then solve each inner block exactly while searching only its three
        // free symbols.  We select an invertible k2-column minor of G2; each
        // assignment to the remaining coordinates determines the pivots.
        std::vector<std::vector<int>> projected_y(
            static_cast<size_t>(k2), std::vector<int>(static_cast<size_t>(n1)));
        for(int a=0;a<k2;++a) {
            std::vector<int> outer_syndrome(static_cast<size_t>(m1));
            for(int r=0;r<m1;++r) outer_syndrome[r]=syndrome[r*k2+a];
            const auto projected=fullMatrixOsd1(outer_h,outer_syndrome,
                outer[a].hard,outer[a].posterior,gf,
                std::min(n1,64),0,1,false);
            if(!projected.syndrome_ok)
                throw std::runtime_error("Could not project generic G2 outer coordinate onto its syndrome");
            for(int block=0;block<n1;++block)
                projected_y[a][block]=projected.estimate[block];
        }
        if(n2-k2!=3)
            throw std::invalid_argument(
                "Exact generic-G2 hard recovery currently requires exactly three free symbols");
        // Find k2 independent columns of G2, then invert that square minor.
        // The remaining three symbol positions are enumerated as free values.
        std::vector<int> echelon=g2;
        std::vector<int> pivot_columns;
        int pivot_row=0;
        for(int c=0;c<n2 && pivot_row<k2;++c) {
            int selected=pivot_row;
            while(selected<k2 && echelon[static_cast<size_t>(selected*n2+c)]==0)
                ++selected;
            if(selected==k2) continue;
            for(int j=0;j<n2;++j)
                std::swap(echelon[static_cast<size_t>(pivot_row*n2+j)],
                          echelon[static_cast<size_t>(selected*n2+j)]);
            const int inverse=gf.inv(echelon[static_cast<size_t>(pivot_row*n2+c)]);
            for(int j=c;j<n2;++j)
                echelon[static_cast<size_t>(pivot_row*n2+j)]=gf.mul(
                    echelon[static_cast<size_t>(pivot_row*n2+j)],inverse);
            for(int r=0;r<k2;++r) if(r!=pivot_row) {
                const int factor=echelon[static_cast<size_t>(r*n2+c)];
                if(!factor) continue;
                for(int j=c;j<n2;++j)
                    echelon[static_cast<size_t>(r*n2+j)]^=gf.mul(factor,
                        echelon[static_cast<size_t>(pivot_row*n2+j)]);
            }
            pivot_columns.push_back(c);
            ++pivot_row;
        }
        if(pivot_row!=k2)
            throw std::invalid_argument("G2 does not have full row rank");
        std::vector<char> is_pivot_column(static_cast<size_t>(n2),0);
        for(int c:pivot_columns) is_pivot_column[c]=1;
        std::vector<int> free_columns;
        for(int c=0;c<n2;++c) if(!is_pivot_column[c]) free_columns.push_back(c);
        std::vector<int> inverse_minor(static_cast<size_t>(k2*k2),0);
        std::vector<int> augmented(static_cast<size_t>(k2*2*k2),0);
        for(int r=0;r<k2;++r) {
            for(int c=0;c<k2;++c)
                augmented[static_cast<size_t>(r*2*k2+c)]=
                    g2[static_cast<size_t>(r*n2+pivot_columns[c])];
            augmented[static_cast<size_t>(r*2*k2+k2+r)]=1;
        }
        for(int c=0;c<k2;++c) {
            int selected=c;
            while(selected<k2 && augmented[static_cast<size_t>(selected*2*k2+c)]==0)
                ++selected;
            if(selected==k2) throw std::logic_error("Selected G2 minor is singular");
            for(int j=0;j<2*k2;++j)
                std::swap(augmented[static_cast<size_t>(c*2*k2+j)],
                          augmented[static_cast<size_t>(selected*2*k2+j)]);
            const int inverse=gf.inv(augmented[static_cast<size_t>(c*2*k2+c)]);
            for(int j=0;j<2*k2;++j)
                augmented[static_cast<size_t>(c*2*k2+j)]=gf.mul(
                    augmented[static_cast<size_t>(c*2*k2+j)],inverse);
            for(int r=0;r<k2;++r) if(r!=c) {
                const int factor=augmented[static_cast<size_t>(r*2*k2+c)];
                if(!factor) continue;
                for(int j=0;j<2*k2;++j)
                    augmented[static_cast<size_t>(r*2*k2+j)]^=gf.mul(factor,
                        augmented[static_cast<size_t>(c*2*k2+j)]);
            }
        }
        for(int r=0;r<k2;++r) for(int c=0;c<k2;++c)
            inverse_minor[static_cast<size_t>(r*k2+c)]=
                augmented[static_cast<size_t>(r*2*k2+k2+c)];

        std::vector<int> estimate(static_cast<size_t>(n1*n2),0);
        if(symbol_posterior_out)
            symbol_posterior_out->assign(static_cast<size_t>(n1*n2),SoftVec(q));
        for(int block=0;block<n1;++block) {
            SoftMatrix y_priors(k2,std::vector<SoftVec>(n1));
            for(int a=0;a<k2;++a)
                for(int j=0;j<n1;++j)
                    y_priors[a][j]=outer[a].extrinsic[j][0];
            std::vector<SoftVec> physical_posterior(n2,SoftVec(q));
            (void)run_inner_block(block,&y_priors,&physical_posterior);
            std::array<std::vector<int>,3> parity_candidates;
            for(int t=0;t<3;++t) {
                const int c=free_columns[t];
                parity_candidates[t].resize(static_cast<size_t>(q));
                std::iota(parity_candidates[t].begin(),parity_candidates[t].end(),0);
                std::partial_sort(parity_candidates[t].begin(),
                    parity_candidates[t].begin()+std::min(4,q),
                    parity_candidates[t].end(),[&](int x,int y) {
                        return physical_posterior[c][x]>physical_posterior[c][y];
                    });
                parity_candidates[t].resize(static_cast<size_t>(std::min(4,q)));
            }
            std::vector<int> best_word(static_cast<size_t>(n2),0);
            double best_score=-std::numeric_limits<double>::infinity();
            for(int p0:parity_candidates[0])
                for(int p1:parity_candidates[1])
                    for(int p2:parity_candidates[2]) {
                        const std::array<int,3> parity{p0,p1,p2};
                        std::vector<int> candidate(static_cast<size_t>(n2),0);
                        for(int t=0;t<3;++t) candidate[free_columns[t]]=parity[t];
                        std::vector<int> rhs(static_cast<size_t>(k2));
                        for(int a=0;a<k2;++a) {
                            int value=projected_y[a][block];
                            for(int t=0;t<3;++t)
                                value=gf.add(value,gf.mul(
                                    g2[static_cast<size_t>(a*n2+free_columns[t])],parity[t]));
                            rhs[a]=value;
                        }
                        for(int i=0;i<k2;++i) {
                            int value=0;
                            for(int j=0;j<k2;++j)
                                value=gf.add(value,gf.mul(
                                    inverse_minor[static_cast<size_t>(i*k2+j)],rhs[j]));
                            candidate[pivot_columns[i]]=value;
                        }
                        double score=0.0;
                        for(int c=0;c<n2;++c)
                            score+=std::log(std::max(1e-300,
                                channel[block*n2+c][0][candidate[c]]));
                        if(score>best_score) {
                            best_score=score;
                            best_word=std::move(candidate);
                        }
                    }
            for(int c=0;c<n2;++c) {
                estimate[block*n2+c]=best_word[c];
                if(symbol_posterior_out)
                    (*symbol_posterior_out)[block*n2+c]=
                        std::move(physical_posterior[c]);
            }
        }
        if(projection_diagnostics)
            projection_diagnostics->final_valid=
                decoder.check_z_syndrome(hx_cascade,estimate,syndrome);
        return estimate;
    }
    bool inner_updated_since_outer=false;
    for(int round=0;round<turbo_rounds;++round) {
        for(int a=0;a<k2;++a)
            outer[a]=runOuterSoftBP(decoder,outer_bp_h,outer_syn[a],inner[a],outer_iterations,gf,
                                    layered_outer_schedule);
        inner_updated_since_outer=false;
        if (independent_outer) break;
        if (stop_on_valid && std::all_of(outer.begin(), outer.end(),
                [](const OuterSoftResult& result) { return result.valid; })) break;
        for(int block=0;block<n1;++block) {
            std::vector<SoftVec> ext(k2);
            for(int a=0;a<k2;++a) ext[a]=outer[a].extrinsic[block][0];
            if (path_generator) {
                // Left[a] and right[a+1] exclude edge a's outer message.
                std::vector<SoftVec> left(n2),right(n2);
                left[0]=channel[block*n2][0];
                for(int a=0;a<k2;++a) {
                    left[a+1]=xorConvolution(left[a],ext[a]);
                    for(int t=0;t<q;++t)
                        left[a+1][t]*=channel[block*n2+a+1][0][t];
                    normalizeLocal(left[a+1]);
                }
                right[k2]=channel[block*n2+k2][0];
                for(int a=k2-1;a>=0;--a) {
                    right[a]=xorConvolution(right[a+1],ext[a]);
                    for(int t=0;t<q;++t)
                        right[a][t]*=channel[block*n2+a][0][t];
                    normalizeLocal(right[a]);
                }
                for(int a=0;a<k2;++a) {
                    SoftVec proposed=xorConvolution(left[a],right[a+1]);
                    for(int x=0;x<q;++x)
                        proposed[x]=damping*proposed[x]+(1.0-damping)*inner[a][block][x];
                    normalizeLocal(proposed);
                    inner[a][block]=std::move(proposed);
                }
            } else {
                std::vector<SoftVec> corr(k2);
                for(int a=0;a<k2;++a)
                    corr[a]=xorConvolution(ext[a],channel[block*n2+a+1][0]);
                for(int a=0;a<k2;++a) {
                    SoftVec weight=channel[block*n2][0];
                    for(int b=0;b<k2;++b) if(b!=a)
                        for(int t=0;t<q;++t) weight[t]*=corr[b][t];
                    normalizeLocal(weight);
                    SoftVec proposed=xorConvolution(weight,channel[block*n2+a+1][0]);
                    for(int x=0;x<q;++x)
                        proposed[x]=damping*proposed[x]+(1.0-damping)*inner[a][block][x];
                    normalizeLocal(proposed);
                    inner[a][block]=std::move(proposed);
                }
            }
        }
        inner_updated_since_outer=true;
    }
    // The last Turbo update changes the priors. Refresh the outer messages
    // before forming a hard decision from them.
    if ((projection_repair || refresh_outer_at_end) && inner_updated_since_outer) {
        for(int a=0;a<k2;++a)
            outer[a]=runOuterSoftBP(decoder,outer_bp_h,outer_syn[a],inner[a],outer_iterations,gf,
                                    layered_outer_schedule);
        if (projection_diagnostics) projection_diagnostics->outer_refreshed=1;
    }
    std::vector<std::vector<int>> projected(k2,std::vector<int>(n1));
    SoftMatrix projected_posterior(k2,std::vector<SoftVec>(n1,SoftVec(q)));
    const bool path_candidate_audit=projection_diagnostics &&
        std::getenv("GF256_CASCADE_PATH_CANDIDATE_AUDIT");
    int audit_free_symbols=32,audit_order=1,audit_alternatives=4;
    bool audit_joint_source=false;
    if(path_candidate_audit) {
        if(const char* value=std::getenv("GF256_CASCADE_PATH_CANDIDATE_FREE"))
            audit_free_symbols=std::stoi(value);
        if(const char* value=std::getenv("GF256_CASCADE_PATH_CANDIDATE_ORDER"))
            audit_order=std::stoi(value);
        if(const char* value=std::getenv("GF256_CASCADE_PATH_CANDIDATE_ALTERNATIVES"))
            audit_alternatives=std::stoi(value);
        if(const char* value=std::getenv("GF256_CASCADE_PATH_CANDIDATE_SOURCE")) {
            const std::string source=value;
            if(source!="outer" && source!="joint")
                throw std::invalid_argument(
                    "PATH_CANDIDATE_SOURCE must be outer or joint");
            audit_joint_source=source=="joint";
        }
        if(audit_free_symbols<1 || audit_free_symbols>256 ||
           audit_order<1 || audit_order>2 ||
           audit_alternatives<1 || audit_alternatives>8)
            throw std::invalid_argument(
                "Path candidate audit requires FREE in [1,256], ORDER in [1,2], and ALTERNATIVES in [1,8]");
        projection_diagnostics->outer_candidate_list.clear();
        projection_diagnostics->outer_candidate_list.reserve(k2);
    }
    for(int a=0;a<k2;++a) {
        for(int j=0;j<n1;++j) {
            SoftVec post=inner[a][j];
            for(int x=0;x<q;++x) post[x]*=outer[a].extrinsic[j][0][x];
            normalizeLocal(post);
            projected_posterior[a][j]=post;
            projected[a][j]=std::distance(post.begin(),std::max_element(post.begin(),post.end()));
        }
        if (projection_diagnostics) {
            projection_diagnostics->bp_valid += outer[a].valid;
            projection_diagnostics->raw_valid +=
                decoder.check_z_syndrome(outer_h,projected[a],outer_syn[a]);
            for(int j=0;j<n1;++j)
                projection_diagnostics->hard_disagreements +=
                    projected[a][j]!=outer[a].hard[j];
        }
        if(path_candidate_audit) {
            // Offline coverage audit only. Both sources use decoder beliefs;
            // the sampled error is consulted later by the diagnostic logger.
            std::vector<std::vector<int>> candidates;
            (void)fullMatrixOsd1(outer_h,outer_syn[a],
                audit_joint_source ? projected[a] : outer[a].hard,
                audit_joint_source ? projected_posterior[a] : outer[a].posterior,
                gf,std::min(n1,audit_free_symbols),audit_order,
                audit_alternatives,true,&candidates);
            projection_diagnostics->outer_candidate_list.push_back(
                std::move(candidates));
        }
        if (projection_repair) {
            if (outer[a].valid) {
                projected[a]=outer[a].hard;
            } else {
                // OSD-0 solves the exact outer syndrome while keeping the
                // most reliable free-symbol decisions. It has no access to
                // the sampled error and cannot guarantee logical success.
                const auto repaired=fullMatrixOsd1(
                    outer_h,outer_syn[a],outer[a].hard,
                    projected_posterior[a],gf,0,0,1,true);
                if (repaired.syndrome_ok) {
                    projected[a]=repaired.estimate;
                    if (projection_diagnostics) ++projection_diagnostics->repaired;
                } else if (projection_diagnostics) {
                    ++projection_diagnostics->repair_failed;
                }
            }
        }
    }
    // Diagnostic-only list feedback: generate new centers using neighboring
    // projection candidates and the shared physical-symbol chain. No truth,
    // no modification of the online estimate, and no claim of exact BP:
    // list evidence is deliberately tempered to limit feedback overconfidence.
    if(path_candidate_audit && path_generator &&
       std::getenv("GF256_CASCADE_NEIGHBOR_CANDIDATE_AUDIT")) {
        auto lists=projection_diagnostics->outer_candidate_list;
        projection_diagnostics->neighbor_candidate_list=lists;
        for(int pass=0;pass<2;++pass) {
            SoftMatrix messages(k2,std::vector<SoftVec>(n1,SoftVec(q)));
            for(int a=0;a<k2;++a) {
                auto ranked=lists[a];
                auto score=[&](const std::vector<int>& word) {
                    double value=0;
                    for(int j=0;j<n1;++j)
                        value+=std::log(std::max(1e-300,projected_posterior[a][j][word[j]]));
                    return value;
                };
                std::stable_sort(ranked.begin(),ranked.end(),[&](const auto& x,const auto& y) {
                    return score(x)>score(y);
                });
                if(ranked.size()>4) ranked.resize(4);
                for(int j=0;j<n1;++j) {
                    messages[a][j]=outer[a].extrinsic[j][0];
                    if(!ranked.empty()) {
                        for(double& value:messages[a][j]) value*=0.5;
                        for(const auto& word:ranked)
                            messages[a][j][word[j]]+=0.5/ranked.size();
                    }
                    normalizeLocal(messages[a][j]);
                }
            }
            SoftMatrix priors(k2,std::vector<SoftVec>(n1,SoftVec(q)));
            for(int block=0;block<n1;++block) {
                std::vector<SoftVec> left(n2),right(n2);
                left[0]=channel[block*n2][0];
                for(int a=0;a<k2;++a) {
                    left[a+1]=xorConvolution(left[a],messages[a][block]);
                    for(int t=0;t<q;++t)
                        left[a+1][t]*=channel[block*n2+a+1][0][t];
                    normalizeLocal(left[a+1]);
                }
                right[k2]=channel[block*n2+k2][0];
                for(int a=k2-1;a>=0;--a) {
                    right[a]=xorConvolution(right[a+1],messages[a][block]);
                    for(int t=0;t<q;++t)
                        right[a][t]*=channel[block*n2+a][0][t];
                    normalizeLocal(right[a]);
                    priors[a][block]=xorConvolution(left[a],right[a+1]);
                }
            }
            for(int a=0;a<k2;++a) {
                const auto updated=runOuterSoftBP(decoder,outer_h,outer_syn[a],
                    priors[a],outer_iterations,gf,layered_outer_schedule);
                std::vector<std::vector<int>> generated;
                (void)fullMatrixOsd1(outer_h,outer_syn[a],updated.hard,
                    updated.posterior,gf,std::min(n1,audit_free_symbols),
                    audit_order,audit_alternatives,true,&generated);
                for(const auto& word:generated)
                    if(!decoder.check_z_syndrome(outer_h,word,outer_syn[a]))
                        throw std::logic_error("Neighbor candidate violates syndrome");
                auto& accumulated=projection_diagnostics->neighbor_candidate_list[a];
                accumulated.insert(accumulated.end(),generated.begin(),generated.end());
                std::sort(accumulated.begin(),accumulated.end());
                accumulated.erase(std::unique(accumulated.begin(),accumulated.end()),accumulated.end());
                lists[a]=std::move(generated);
            }
        }
    }
    // Offline physical-window OSD. A candidate is one word over the three
    // physical positions touched by two neighboring G2 rows. Both projected
    // outer syndromes are solved in the SAME linear system, so their words
    // always agree on the shared physical position. This is diagnostic only.
    static bool pair_audit_run=false;
    if(path_candidate_audit && path_generator && !pair_audit_run &&
       std::getenv("GF256_CASCADE_PAIR_CANDIDATE_AUDIT")) {
        pair_audit_run=true;
        projection_diagnostics->pair_candidate_list=
            projection_diagnostics->outer_candidate_list;
        int chosen=-1, largest=-1;
        for(int a=0;a+1<k2;++a) {
            int residual=0;
            for(int c=a;c<=a+1;++c)
                for(int r=0;r<m1;++r) {
                    int observed=0;
                    for(SparseMatrix<int,RowMajor>::InnerIterator it(outer_h,r);it;++it)
                        observed^=gf.mul(it.value(),outer[c].hard[it.col()]);
                    residual+=observed!=outer_syn[c][r];
                }
            if(residual>largest) { largest=residual; chosen=a; }
        }
        const bool all_pairs=std::getenv("GF256_CASCADE_PAIR_ALL_AUDIT")!=nullptr;
        const int first_pair=all_pairs ? 0 : chosen;
        const int last_pair=all_pairs ? k2-2 : chosen;
        if(first_pair>=0) for(int left_coordinate=first_pair;
                              left_coordinate<=last_pair;++left_coordinate) {
            projection_diagnostics->pair_audit_left=left_coordinate;
            std::vector<Eigen::Triplet<int>> edges;
            for(int r=0;r<m1;++r) for(int offset=0;offset<2;++offset) {
                const int original_row=r*k2+left_coordinate+offset;
                for(SparseMatrix<int,RowMajor>::InnerIterator it(hx_cascade,original_row);it;++it) {
                    const int local=it.col()%n2;
                    if(local<left_coordinate || local>left_coordinate+2)
                        throw std::logic_error("Pair audit requires adjacent-pair G2 support");
                    edges.emplace_back(r*2+offset,(it.col()/n2)*3+
                        local-left_coordinate,it.value());
                }
            }
            SparseMatrix<int,RowMajor> pair_h(2*m1,3*n1);
            pair_h.setFromTriplets(edges.begin(),edges.end());
            std::vector<int> pair_syndrome(2*m1),pair_hard(3*n1);
            std::vector<SoftVec> pair_posterior(3*n1,SoftVec(q));
            for(int r=0;r<m1;++r) for(int offset=0;offset<2;++offset)
                pair_syndrome[r*2+offset]=outer_syn[left_coordinate+offset][r];
            for(int block=0;block<n1;++block) {
                std::vector<SoftVec> left(n2),right(n2);
                left[0]=channel[block*n2][0];
                for(int a=0;a<k2;++a) {
                    left[a+1]=xorConvolution(left[a],outer[a].extrinsic[block][0]);
                    for(int x=0;x<q;++x)
                        left[a+1][x]*=channel[block*n2+a+1][0][x];
                    normalizeLocal(left[a+1]);
                }
                right[k2]=channel[block*n2+k2][0];
                for(int a=k2-1;a>=0;--a) {
                    right[a]=xorConvolution(right[a+1],outer[a].extrinsic[block][0]);
                    for(int x=0;x<q;++x)
                        right[a][x]*=channel[block*n2+a][0][x];
                    normalizeLocal(right[a]);
                }
                for(int offset=0;offset<3;++offset) {
                    const int physical=left_coordinate+offset;
                    auto& posterior=pair_posterior[block*3+offset];
                    for(int x=0;x<q;++x)
                        posterior[x]=left[physical][x]*right[physical][x]/
                            std::max(1e-300,channel[block*n2+physical][0][x]);
                    normalizeLocal(posterior);
                    pair_hard[block*3+offset]=std::distance(
                        posterior.begin(),std::max_element(posterior.begin(),posterior.end()));
                }
            }
            std::vector<std::vector<int>> physical_candidates;
            const auto pair_result=fullMatrixOsd1(pair_h,pair_syndrome,pair_hard,
                pair_posterior,gf,std::min(24,audit_free_symbols),1,
                std::min(4,audit_alternatives),true,&physical_candidates);
            if(!pair_result.syndrome_ok || physical_candidates.empty())
                throw std::logic_error("Pair OSD did not solve both outer syndromes");
            projection_diagnostics->pair_audit_candidates+=physical_candidates.size();
            for(const auto& physical:physical_candidates) for(int offset=0;offset<2;++offset) {
                std::vector<int> word(n1);
                for(int block=0;block<n1;++block)
                    word[block]=gf.add(physical[block*3+offset],
                                       physical[block*3+offset+1]);
                if(!decoder.check_z_syndrome(outer_h,word,
                                             outer_syn[left_coordinate+offset]))
                    throw std::logic_error("Pair candidate violates outer syndrome");
                projection_diagnostics->pair_candidate_list[left_coordinate+offset].push_back(
                    std::move(word));
            }
            for(int offset=0;offset<2;++offset) {
                auto& list=projection_diagnostics->pair_candidate_list[left_coordinate+offset];
                std::sort(list.begin(),list.end());
                list.erase(std::unique(list.begin(),list.end()),list.end());
            }
        }
        if(all_pairs) projection_diagnostics->pair_audit_left=-2;
    }
    if (dual_outer_checks && max_inactivation_checks>0) {
        if(dual_outer_checks->cols()!=outer_h.cols())
            throw std::invalid_argument("Dual outer check dimensions do not match");
        for(int a=0;a<k2;++a) {
            const std::vector<int> original_syndrome(
                outer_syn[a].begin(),outer_syn[a].begin()+m1);
            if(decoder.check_z_syndrome(outer_h,projected[a],original_syndrome)) continue;
            const auto repaired=decodeOuterWithStabilizerInactivation(
                decoder,outer_h,*dual_outer_checks,original_syndrome,
                inner[a],projected_posterior[a],outer_iterations,
                max_inactivation_checks,gf);
            if(inactivation_diagnostics) {
                inactivation_diagnostics->inactivation_attempts+=repaired.attempts;
                inactivation_diagnostics->inactivation_successes+=repaired.syndrome_ok;
                inactivation_diagnostics->inactivated_checks+=
                    repaired.successful_check>=0 ? 1 : 0;
            }
            if(repaired.syndrome_ok) projected[a]=repaired.estimate;
        }
    }
    if (joint_candidate_selection) {
        // Each set of k projections leaves one free physical symbol per block.
        auto joint_log_likelihood=[&](const std::vector<std::vector<int>>& candidate) {
            double total=0.0;
            for(int block=0;block<n1;++block) {
                double probability=0.0;
                for(int t=0;t<q;++t) {
                    double term=channel[block*n2][0][t];
                    int symbol=t;
                    for(int a=0;a<k2;++a) {
                        symbol=gf.add(path_generator ? symbol : t,candidate[a][block]);
                        term*=channel[block*n2+a+1][0][symbol];
                    }
                    probability+=term;
                }
                total+=std::log(std::max(1e-300,probability));
            }
            return total;
        };
        double score=joint_log_likelihood(projected);
        const double initial_score=score;
        for(int pass=0;pass<2;++pass) {
            bool changed=false;
            for(int a=0;a<k2;++a) {
                std::vector<SoftVec> conditional(n1,SoftVec(q));
                for(int block=0;block<n1;++block) {
                    if (path_generator) {
                        // All other fixed projections split the physical
                        // chain into two independent sides around edge a.
                        SoftVec left(q),right(q);
                        for(int t=0;t<q;++t) {
                            int symbol=t;
                            left[t]=channel[block*n2+a][0][symbol];
                            for(int b=a-1;b>=0;--b) {
                                symbol=gf.add(symbol,projected[b][block]);
                                left[t]*=channel[block*n2+b][0][symbol];
                            }
                            symbol=t;
                            right[t]=channel[block*n2+a+1][0][symbol];
                            for(int b=a+1;b<k2;++b) {
                                symbol=gf.add(symbol,projected[b][block]);
                                right[t]*=channel[block*n2+b+1][0][symbol];
                            }
                        }
                        normalizeLocal(left);
                        normalizeLocal(right);
                        conditional[block]=xorConvolution(left,right);
                        continue;
                    }
                    SoftVec common(q);
                    for(int t=0;t<q;++t) {
                        common[t]=channel[block*n2][0][t];
                        for(int b=0;b<k2;++b) if(b!=a)
                            common[t]*=channel[block*n2+b+1][0][
                                gf.add(t,projected[b][block])];
                    }
                    normalizeLocal(common);
                    conditional[block]=xorConvolution(
                        common,channel[block*n2+a+1][0]);
                }
                // Solve pivot symbols exactly. The optional joint-list mode
                // uses a bounded OSD-2 candidate list; other modes retain
                // their smaller existing searches.
                std::vector<std::vector<int>> candidate_list;
                const int list_free_symbols = joint_list_selection
                    ? 10 : (k2==2 ? 8 : 4);
                const int list_osd_order = joint_list_selection
                    ? 2 : (k2==2 ? 2 : 1);
                const int list_alternatives = joint_list_selection
                    ? 2 : (k2==2 ? 6 : 4);
                const auto candidate=fullMatrixOsd1(
                    outer_h,outer_syn[a],projected[a],conditional,gf,
                    list_free_symbols,list_osd_order,list_alternatives,true,
                    joint_list_selection ? &candidate_list : nullptr);
                if(!candidate.syndrome_ok) continue;
                if (!joint_list_selection)
                    candidate_list.push_back(candidate.estimate);
                if (joint_list_selection) {
                    candidate_list.push_back(candidate.estimate);
                    std::sort(candidate_list.begin(),candidate_list.end());
                    candidate_list.erase(
                        std::unique(candidate_list.begin(),candidate_list.end()),
                        candidate_list.end());
                    auto projection_log_score=[&](const std::vector<int>& word) {
                        double score=0.0;
                        for(int j=0;j<n1;++j)
                            score+=std::log(std::max(1e-300,conditional[j][word[j]]));
                        return score;
                    };
                    std::stable_sort(candidate_list.begin(),candidate_list.end(),
                        [&](const auto& lhs,const auto& rhs) {
                            return projection_log_score(lhs)>
                                   projection_log_score(rhs);
                        });
                    // OSD-2 can produce hundreds of outer candidates. Keep
                    // only a bounded high-posterior list for the joint score.
                    if(candidate_list.size()>16) candidate_list.resize(16);
                }
                if (projection_diagnostics)
                    projection_diagnostics->joint_candidates_tested +=
                        static_cast<int>(candidate_list.size());
                for (const auto& alternative:candidate_list) {
                    auto trial=projected;
                    trial[a]=alternative;
                    const double trial_score=joint_log_likelihood(trial);
                    if(trial_score>score+1e-9) {
                        projected=std::move(trial);
                        score=trial_score;
                        changed=true;
                        if(projection_diagnostics) ++projection_diagnostics->joint_updates;
                    }
                }
            }
            if(!changed) break;
        }
        if(projection_diagnostics)
            projection_diagnostics->joint_log_gain=score-initial_score;
        if(joint_score_out) *joint_score_out=score;
    }
    std::vector<int> estimate(n1*n2,0);
    for(int block=0;block<n1;++block) {
        int best_tail=0; double best_score=-std::numeric_limits<double>::infinity();
        for(int tail=0;tail<q;++tail) {
            double score=std::log(std::max(1e-300,channel[block*n2][0][tail]));
            int symbol=tail;
            for(int a=0;a<k2;++a) {
                symbol=gf.add(path_generator ? symbol : tail,projected[a][block]);
                score+=std::log(std::max(1e-300,channel[block*n2+a+1][0][symbol]));
            }
            if(score>best_score){best_score=score;best_tail=tail;}
        }
        estimate[block*n2]=best_tail;
        for(int a=0;a<k2;++a)
            estimate[block*n2+a+1]=gf.add(
                projected[a][block],
                path_generator ? estimate[block*n2+a] : best_tail);
    }
    std::vector<int> physical_estimate=estimate;
    for (int v=0;v<n1*n2;++v)
        physical_estimate[v]=gf.mul(estimate[v],inverse_scales[v%n2]);
    if (projection_diagnostics)
        projection_diagnostics->final_valid=
            decoder.check_z_syndrome(hx_cascade,physical_estimate,syndrome);
    if (cavity_symbol_posterior) {
        // Physical priors are already present in the chain. Reusing a full
        // projection posterior here would count those priors a second time.
        for(int a=0;a<k2;++a) for(int j=0;j<n1;++j)
            projected_posterior[a][j]=outer[a].extrinsic[j][0];
    }
    if (symbol_posterior_out) {
        symbol_posterior_out->assign(static_cast<size_t>(n1*n2), SoftVec(q, 0.0));
        for (int block=0; block<n1; ++block) {
            if (path_generator) {
                // Marginalize the soft projection distributions along the
                // chain to give each physical symbol a posterior for OSD.
                std::vector<SoftVec> left(n2),right(n2);
                left[0]=channel[block*n2][0];
                for(int a=0;a<k2;++a) {
                    left[a+1]=xorConvolution(left[a],projected_posterior[a][block]);
                    for(int t=0;t<q;++t)
                        left[a+1][t]*=channel[block*n2+a+1][0][t];
                    normalizeLocal(left[a+1]);
                }
                right[k2]=channel[block*n2+k2][0];
                for(int a=k2-1;a>=0;--a) {
                    right[a]=xorConvolution(right[a+1],projected_posterior[a][block]);
                    for(int t=0;t<q;++t)
                        right[a][t]*=channel[block*n2+a][0][t];
                    normalizeLocal(right[a]);
                }
                for(int a=0;a<n2;++a) {
                    SoftVec posterior(q);
                    // Avoid dividing by tiny channel entries by rebuilding
                    // one cavity message from the adjacent branch.
                    if (a==0) posterior=right[0];
                    else if (a==k2) posterior=left[k2];
                    else {
                        SoftVec right_cavity=xorConvolution(
                            right[a+1],projected_posterior[a][block]);
                        for(int t=0;t<q;++t)
                            posterior[t]=left[a][t]*right_cavity[t];
                        normalizeLocal(posterior);
                    }
                    (*symbol_posterior_out)[block*n2+a]=std::move(posterior);
                }
                continue;
            }
            std::vector<SoftVec> projected_channel(k2);
            for (int a=0; a<k2; ++a)
                projected_channel[a]=xorConvolution(
                    projected_posterior[a][block],channel[block*n2+a+1][0]);
            SoftVec tail_posterior(q,0.0);
            for (int tail=0; tail<q; ++tail) {
                double weight=channel[block*n2][0][tail];
                for (int a=0; a<k2; ++a) weight*=projected_channel[a][tail];
                tail_posterior[tail]=weight;
            }
            normalizeLocal(tail_posterior);
            (*symbol_posterior_out)[block*n2]=tail_posterior;
            for (int a=0; a<k2; ++a) {
                SoftVec symbol_posterior(q,0.0);
                for (int tail=0; tail<q; ++tail) {
                    double weight=channel[block*n2][0][tail];
                    for (int b=0; b<k2; ++b)
                        if (b!=a) weight*=projected_channel[b][tail];
                    for (int symbol=0; symbol<q; ++symbol) {
                        const int projection=gf.add(symbol,tail);
                        symbol_posterior[symbol]+=weight*
                            channel[block*n2+a+1][0][symbol]*
                            projected_posterior[a][block][projection];
                    }
                }
                normalizeLocal(symbol_posterior);
                (*symbol_posterior_out)[block*n2+a+1]=std::move(symbol_posterior);
            }
        }
    }
    if (weighted_path && symbol_posterior_out) {
        for (int v=0;v<n1*n2;++v) {
            SoftVec physical(q);
            for (int z=0;z<q;++z)
                physical[z]=(*symbol_posterior_out)[v][gf.mul(z,scales[v%n2])];
            (*symbol_posterior_out)[v]=std::move(physical);
        }
    }
    return physical_estimate;
}

int main(int argc, char* argv[]) {
    try {
        if (argc > 9) {
            throw invalid_argument(
                "Usage: simulate_gf256_cascade [trials] [seed] [p_D] "
                "[tit|cascade|structured|conditioned|conditioned_projection_independent|conditioned_robust|conditioned_robust_softx|conditioned_robust_structured_softx|conditioned_robust_augmented|conditioned_robust_layered|conditioned_robust_overcomplete|conditioned_robust_si|conditioned_fullhx_bp|joint_factorized|conditioned_robust_refresh|conditioned_robust_joint_projection|conditioned_robust_joint_projection_safe|conditioned_robust_joint_list|conditioned_robust_outer_osd|conditioned_robust_joint_osd|conditioned_robust_fullhz|conditioned_robust_osd|conditioned_robust_osd2_weak|conditioned_robust_long|conditioned_robust_tuned|conditioned_robust_adaptive|conditioned_robust_restart_belief|conditioned_robust_restart_belief_global|conditioned_robust_tuned_restart_belief_global|conditioned_robust_tuned_syndrome_flip|conditioned_robust_tuned_syndrome_flip_layered] [target-trial] [field_degree] [matrix_dir] [inner_length]");
        }
        const int trials = argc >= 2 ? stoi(argv[1]) : 100;
        if (trials <= 0) {
            throw invalid_argument("The trial count must be positive");
        }

        const unsigned int seed = argc >= 3 ? static_cast<unsigned int>(stoul(argv[2]))
                                            : static_cast<unsigned int>(time(nullptr));
        const string construction = argc >= 5 ? argv[4] : "direct";
        const int target_trial = argc >= 6 ? stoi(argv[5]) : -1;
        const bool target_diagnostics = target_trial >= 0 &&
            std::getenv("GF256_CASCADE_SKIP_TARGET_DIAG") == nullptr;
        const int field_degree = argc >= 7 ? stoi(argv[6]) : 8;
        const string matrix_dir = argc >= 8 ? argv[7] : "build";
        const int inner_length = argc >= 9 ? stoi(argv[8]) : 6;
        const char* z_bias_setting = std::getenv("GF256_Z_BIAS");
        const double z_bias = z_bias_setting ? stod(z_bias_setting) : 1.0;
        if ((field_degree != 3 && field_degree != 8 && field_degree != 10) ||
            inner_length < 2)
            throw invalid_argument("Supported field degrees: 3, 8, 10; inner_length >= 2");
        if (target_trial >= trials)
            throw invalid_argument("Target trial must be less than the frame count");
        if (construction != "direct" && construction != "tit" && construction != "cascade" && construction != "structured" && construction != "conditioned" && construction != "conditioned_projection_independent" && construction != "conditioned_robust" && construction != "conditioned_robust_softx" && construction != "conditioned_robust_structured_softx" && construction != "conditioned_robust_augmented" && construction != "conditioned_robust_layered" && construction != "conditioned_robust_overcomplete" && construction != "conditioned_robust_si" && construction != "conditioned_fullhx_bp" && construction != "joint_factorized" && construction != "conditioned_robust_refresh" && construction != "conditioned_robust_joint_projection" && construction != "conditioned_robust_joint_projection_safe" && construction != "conditioned_robust_joint_list" && construction != "conditioned_robust_outer_osd" && construction != "conditioned_robust_joint_osd" && construction != "conditioned_robust_fullhz" && construction != "conditioned_robust_osd" && construction != "conditioned_robust_osd2_weak" && construction != "conditioned_robust_long" && construction != "conditioned_robust_tuned" && construction != "conditioned_robust_adaptive" && construction != "conditioned_robust_guess" && construction != "conditioned_robust_bpgd" && construction != "conditioned_robust_restart_belief" && construction != "conditioned_robust_restart_belief_global" && construction != "conditioned_robust_tuned_restart_belief_global" && construction != "conditioned_robust_tuned_syndrome_flip" && construction != "conditioned_robust_tuned_syndrome_flip_layered" && construction != "video") {
            throw invalid_argument("Unknown decoder mode");
        }
        const bool conditioned_mode = construction == "conditioned" ||
                                      construction == "conditioned_projection_independent" ||
                                      construction == "conditioned_robust" ||
                                      construction == "conditioned_robust_softx" ||
                                      construction == "conditioned_robust_structured_softx" ||
                                      construction == "conditioned_robust_augmented" ||
                                      construction == "conditioned_robust_layered" ||
                                      construction == "conditioned_robust_overcomplete" ||
                                      construction == "conditioned_robust_si" ||
                                      construction == "conditioned_robust_joint_projection_safe" ||
                                      construction == "conditioned_fullhx_bp" ||
                                      construction == "conditioned_robust_refresh" ||
                                      construction == "conditioned_robust_joint_projection" ||
                                      construction == "conditioned_robust_joint_list" ||
                                      construction == "conditioned_robust_outer_osd" ||
                                      construction == "conditioned_robust_joint_osd" ||
                                      construction == "conditioned_robust_fullhz" ||
                                      construction == "conditioned_robust_osd" ||
                                      construction == "conditioned_robust_osd2_weak" ||
                                      construction == "conditioned_robust_long" ||
                                      construction == "conditioned_robust_tuned" ||
                                      construction == "conditioned_robust_tuned_syndrome_flip" ||
                                      construction == "conditioned_robust_tuned_syndrome_flip_layered" ||
                                      construction == "conditioned_robust_adaptive" ||
                                      construction == "conditioned_robust_guess" ||
                                      construction == "conditioned_robust_bpgd" ||
                                      construction == "conditioned_robust_restart_belief" ||
                                      construction == "conditioned_robust_restart_belief_global" ||
                                      construction == "conditioned_robust_tuned_restart_belief_global";
        const bool restart_belief_global_mode =
            construction == "conditioned_robust_restart_belief_global" ||
            construction == "conditioned_robust_tuned_restart_belief_global";
        const bool syndrome_flip_mode =
            construction == "conditioned_robust_tuned_syndrome_flip" ||
            construction == "conditioned_robust_tuned_syndrome_flip_layered";
        const bool layered_after_flip_mode =
            construction == "conditioned_robust_tuned_syndrome_flip_layered";
        const bool full_hz_mode = construction == "conditioned_robust_fullhz";
        const bool full_hx_bp_mode = construction == "conditioned_fullhx_bp";
        const char* full_hx_iters_setting = std::getenv("GF256_FULLHX_BP_ITERS");
        const int full_hx_bp_iterations = full_hx_iters_setting
            ? std::stoi(full_hx_iters_setting) : 20;
        if (full_hx_bp_mode && (full_hx_bp_iterations < 1 || full_hx_bp_iterations > 200))
            throw std::invalid_argument("GF256_FULLHX_BP_ITERS must be in 1..200");
        const bool robust_mode = construction == "conditioned_robust" ||
                                 construction == "conditioned_robust_softx" ||
                                 construction == "conditioned_robust_structured_softx" ||
                                 construction == "conditioned_robust_augmented" ||
                                 construction == "conditioned_robust_layered" ||
                                 construction == "conditioned_robust_overcomplete" ||
                                 construction == "conditioned_robust_si" ||
                                 construction == "conditioned_robust_joint_projection_safe" ||
                                 construction == "conditioned_robust_refresh" || full_hz_mode ||
                                 construction == "conditioned_robust_joint_projection" ||
                                 construction == "conditioned_robust_joint_list" ||
                                 construction == "conditioned_robust_outer_osd" ||
                                 construction == "conditioned_robust_joint_osd" ||
                                 construction == "conditioned_robust_osd" ||
                                 construction == "conditioned_robust_osd2_weak" ||
                                 construction == "conditioned_robust_long" ||
                                      construction == "conditioned_robust_tuned" ||
                                      construction == "conditioned_robust_tuned_restart_belief_global" ||
                                      construction == "conditioned_robust_tuned_syndrome_flip" ||
                                      construction == "conditioned_robust_tuned_syndrome_flip_layered" ||
                                      construction == "conditioned_robust_adaptive" ||
                                      construction == "conditioned_robust_guess" ||
                                      construction == "conditioned_robust_bpgd" ||
                                 construction == "conditioned_robust_restart_belief" ||
                                      construction == "conditioned_robust_restart_belief_global" ||
                                      construction == "conditioned_robust_tuned_restart_belief_global";
        const bool long_turbo_mode = construction == "conditioned_robust_long";
        const bool tuned_turbo_mode = construction == "conditioned_robust_tuned" ||
                                      construction == "conditioned_robust_tuned_restart_belief_global" ||
                                      construction == "conditioned_robust_tuned_syndrome_flip" ||
                                      construction == "conditioned_robust_tuned_syndrome_flip_layered" ||
                                      construction == "conditioned_robust_softx" ||
                                      construction == "conditioned_robust_structured_softx" ||
                                      construction == "conditioned_robust_augmented" ||
                                      construction == "conditioned_robust_layered" ||
                                      construction == "conditioned_robust_overcomplete" ||
                                      construction == "conditioned_robust_si" ||
                                      construction == "conditioned_robust_joint_projection_safe";
        const bool stabilizer_inactivation_mode =
            construction == "conditioned_robust_si";
        const bool soft_x_mode = construction == "conditioned_robust_softx";
        const bool structured_soft_x_mode =
            construction == "conditioned_robust_structured_softx";
        const bool augmented_decoder_mode =
            construction == "conditioned_robust_augmented";
        const bool layered_outer_schedule_mode =
            construction == "conditioned_robust_layered";
        const bool overcomplete_decoder_mode =
            construction == "conditioned_robust_overcomplete";
        const bool adaptive_turbo_mode = construction == "conditioned_robust_adaptive" ||
                                      construction == "conditioned_robust_guess" ||
                                      construction == "conditioned_robust_osd2_weak" ||
                                      construction == "conditioned_robust_restart_belief" ||
                                      construction == "conditioned_robust_restart_belief_global";
        const bool guided_guess_mode = construction == "conditioned_robust_guess";
        const bool guided_bpgd_mode = construction == "conditioned_robust_bpgd";
        const bool restart_belief_mode =
            construction == "conditioned_robust_restart_belief" ||
            restart_belief_global_mode;
        const bool independent_projection_mode =
            construction == "conditioned_projection_independent";
        const char* independent_iters_setting =
            std::getenv("GF256_INDEPENDENT_OUTER_ITERS");
        const int independent_outer_iterations = independent_iters_setting
            ? std::stoi(independent_iters_setting) : 20;
        if (independent_projection_mode &&
            (independent_outer_iterations < 1 || independent_outer_iterations > 200))
            throw std::invalid_argument("GF256_INDEPENDENT_OUTER_ITERS must be in 1..200");
        const bool weak_osd2_fallback_mode =
            construction == "conditioned_robust_osd2_weak";
        const bool osd_fallback_mode = construction == "conditioned_robust_osd";
        const bool outer_osd_mode = construction == "conditioned_robust_outer_osd";
        const bool joint_osd_mode = construction == "conditioned_robust_joint_osd";
        const bool projection_mode = outer_osd_mode || joint_osd_mode;
        const bool refresh_outer_mode = construction == "conditioned_robust_refresh";
        const bool joint_list_mode = construction == "conditioned_robust_joint_list";
        const bool preserve_valid_joint_projection_mode =
            construction == "conditioned_robust_joint_projection_safe";
        const bool joint_projection_mode = construction == "conditioned_robust_joint_projection" ||
                                           preserve_valid_joint_projection_mode ||
                                           joint_list_mode;
        const bool projection_risk_audit =
            std::getenv("GF256_CASCADE_PROJECTION_RISK_AUDIT") != nullptr;
        const bool projection_diagnostic_mode = projection_mode ||
            stabilizer_inactivation_mode ||
            preserve_valid_joint_projection_mode ||
            projection_risk_audit ||
            (conditioned_mode &&
             std::getenv("GF256_CASCADE_PROJECTION_DIAG") != nullptr);
        const bool joint_factorized_mode = construction == "joint_factorized";
        const bool cascade_mode = construction == "cascade" || joint_factorized_mode ||
            construction == "structured" || conditioned_mode;
        if (z_bias != 1.0 && !adaptive_turbo_mode &&
            !independent_projection_mode)
            throw invalid_argument("GF256_Z_BIAS supports adaptive or independent conditioned decoding only");
        srand(seed);
        cerr << "seed=" << seed << " construction=" << construction << '\n';
        const vector<int> field_poly = field_degree == 3 ? vector<int>{1,0,1,1}
            : field_degree == 8 ? vector<int>{1,0,0,0,1,1,1,0,1}
                                : vector<int>{1,0,0,0,0,0,0,1,0,0,1};
        GF2p gf(field_degree, field_poly);
        NB_BP_Decoder decoder(gf, 50);
        FieldBasis basis(gf);

        const string field_tag = "GF" + to_string(gf.getSize());
        const string hx_prefix = matrix_dir + "/HX" +
            ((construction == "direct" || construction == "tit") ? "" :
             (cascade_mode ? "_Cascade" : "_Video"));
        const string hz_prefix = matrix_dir + "/HZ" +
            ((construction == "direct" || construction == "tit") ? "" :
             (cascade_mode ? "_Cascade" : "_Video"));
        const string suffix = (construction == "direct" || construction == "tit")
            ? "_Safe_" + field_tag + ".txt" : "_" + field_tag + ".txt";
        const auto hx_logical = loadMatrix(hx_prefix + suffix, gf);
        const auto hz_logical = loadMatrix(hz_prefix + suffix, gf);
        SparseMatrix<int, RowMajor> cascade_outer_hz_check;
        SparseMatrix<int, RowMajor> cascade_outer_hx_check;
        SparseMatrix<int, RowMajor> cascade_inner_h2;
        SparseMatrix<int, RowMajor> cascade_inner_g2;
        std::vector<RedundantOuterCheck> redundant_outer_checks;
        std::optional<FieldRowSpace> cascade_x_stabilizer_space;
        std::optional<FieldRowSpace> cascade_z_stabilizer_space;
        if (cascade_mode) {
            cascade_outer_hz_check = loadMatrix(matrix_dir + "/HX_TIT_Outer_For_HZ_" + field_tag + ".txt", gf);
            cascade_outer_hx_check = loadMatrix(matrix_dir + "/HZ_TIT_Outer_For_HX_" + field_tag + ".txt", gf);
            const string h2_file=matrix_dir + "/H2_Inner_" + field_tag + ".txt";
            const string g2_file=matrix_dir + "/G2_Inner_" + field_tag + ".txt";
            if(ifstream(h2_file)) cascade_inner_h2=loadMatrix(h2_file,gf);
            if(ifstream(g2_file)) cascade_inner_g2=loadMatrix(g2_file,gf);
            if(static_cast<bool>(cascade_inner_h2.rows()) !=
               static_cast<bool>(cascade_inner_g2.rows()))
                throw runtime_error("Cascade inner H2 and G2 must be supplied together");
            cascade_x_stabilizer_space.emplace(cascade_outer_hx_check,gf,
                nullptr,cascade_inner_g2.rows() ? &cascade_inner_g2 : nullptr);
            cascade_z_stabilizer_space.emplace(cascade_outer_hz_check,gf,
                cascade_inner_h2.rows() ? &cascade_inner_h2 : nullptr,nullptr);
            if(overcomplete_decoder_mode) {
                redundant_outer_checks=buildPairwiseRedundantOuterChecks(
                    cascade_outer_hx_check,gf,
                    static_cast<size_t>(cascade_outer_hx_check.rows()));
                cerr << "[overcomplete-outer] original_rows="
                     << cascade_outer_hx_check.rows()
                     << " added_rows=" << redundant_outer_checks.size() << '\n';
            }
        }
        vector<int> inner_scales(inner_length,1);
        if (cascade_mode) {
            if(cascade_inner_h2.rows()) {
                if(cascade_inner_h2.cols()!=inner_length)
                    throw runtime_error("Saved H2 has the wrong inner length");
                if (cascade_inner_h2.rows()==1)
                    for (int j=0;j<inner_length;++j) {
                        if (cascade_inner_h2.coeff(0,j)==0)
                            throw runtime_error("One-row H2 has a zero column");
                        inner_scales[j]=gf.inv(cascade_inner_h2.coeff(0,j));
                    }
            }
            if(cascade_inner_g2.rows() &&
               (cascade_inner_g2.cols()!=inner_length ||
                cascade_inner_g2.rows()!=inner_length-cascade_inner_h2.rows()))
                throw runtime_error("Saved G2 dimensions do not match H2");
        }
        const auto hx_physical = loadBinaryMatrix(construction == "tit" ? matrix_dir + "/HX_TIT_Physical_Binary.txt" : hx_prefix + "_Physical_Binary.txt");
        const auto hz_physical = loadBinaryMatrix(construction == "tit" ? matrix_dir + "/HZ_TIT_Physical_Binary.txt" : hz_prefix + "_Physical_Binary.txt");

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
        const int cascade_inner_dimension = cascade_mode
            ? hx_logical.rows() / cascade_outer_hx_check.rows() : 0;
        const int physical_bits = symbol_count * field_degree;
        const vector<double> error_rates = argc >= 4 ? vector<double>{stod(argv[3])}
                                                   : vector<double>{0.001, 0.005, 0.01, 0.02, 0.04};
        cout << (z_bias == 1.0 ? "p_D" : "p_total")
                << "\tLogicalFailureRate\tTotal\tFailures\tSyndromeFailures"
                "\tXSyndromeFailures\tZSyndromeFailures"
                "\tMeanPhysicalErrors\tMeanResidualSymbols\tUndetectedLogicalFailures"
                "\tExactRecoveryFrames\tExactXFrames\tCascadeTrailingRecoveries"
                "\tCascadeStagedAttempts\tCascadeStagedValid"
                "\tXFallbackAttempts\tXFallbackSuccesses"
                "\tZSoftFallbackAttempts\tZSoftFallbackSuccesses"
                "\tZOsdFallbackAttempts\tZOsdFallbackSuccesses"
                "\tZTunedFallbackAttempts\tZTunedFallbackSuccesses"
                "\tRobustAttempts\tRobustSelections"
                "\tXFullHzBpAttempts\tXFullHzBpSuccesses"
                "\tXSoftPosteriorAttempts\tXSoftPosteriorSuccesses"
                "\tZAugmentedAttempts\tZAugmentedSelections\n";

        for (double p : error_rates) {
            pauli_channel = PauliChannel::zBiased(p, z_bias);
            const int evaluated_frames = target_trial >= 0 ? 1 : trials;
            int failures = 0;
            int syndrome_failures = 0;
            int x_syndrome_failures = 0;
            int z_syndrome_failures = 0;
            int undetected_logical_failures = 0;
            int exact_recovery_frames = 0;
            int exact_x_frames = 0;
            int cascade_trailing_recoveries = 0;
            int cascade_staged_attempts = 0;
            int cascade_staged_successes = 0;
            int x_fallback_attempts = 0;
            int x_fallback_successes = 0;
            int x_full_hz_bp_attempts = 0;
            int x_full_hz_bp_successes = 0;
            int x_soft_posterior_attempts = 0;
            int x_soft_posterior_successes = 0;
            int z_augmented_attempts = 0;
            int z_augmented_selections = 0;
            int z_soft_fallback_attempts = 0;
            int z_soft_fallback_successes = 0;
            int z_osd_fallback_attempts = 0;
            int z_osd_fallback_successes = 0;
            int z_tuned_fallback_attempts = 0;
            int z_tuned_fallback_successes = 0;
            int robust_attempts = 0;
            int robust_selections = 0;
            int physical_errors = 0;
            int residual_symbols = 0;
            for (int trial = 0; trial < trials; ++trial) {
                vector<int> error_x(physical_bits, 0), error_z(physical_bits, 0);
                samplePauli(error_x, error_z);
                if (target_trial >= 0 && trial != target_trial) continue;
                int frame_physical_errors = 0;
                for (int bit = 0; bit < physical_bits; ++bit) {
                    frame_physical_errors += error_x[bit] || error_z[bit];
                }
                physical_errors += frame_physical_errors;

                const vector<int> syndrome_x = getBinarySyndrome(hz_physical, error_x);
                const vector<int> syndrome_z = getBinarySyndrome(hx_physical, error_z);
                const vector<int> logical_syndrome_x =
                    binaryToFieldSymbols(syndrome_x, field_degree, basis, true);
                const vector<int> logical_syndrome_z =
                    binaryToFieldSymbols(syndrome_z, field_degree, basis, false);
                vector<int> corrected_x;
                vector<int> corrected_z;
                vector<int> alternate_z;
                ProjectionDiagnostics initial_projection_diagnostics;
                ProjectionDiagnostics retry_projection_diagnostics;
                ProjectionDiagnostics contender_projection_diagnostics;
                std::vector<std::pair<int,ProjectionDiagnostics>> layered_audit_diagnostics;
                bool projection_retry_attempted=false;
                double current_joint_score=-std::numeric_limits<double>::infinity();
                double retry_joint_score=-std::numeric_limits<double>::infinity();
                double contender_joint_score=-std::numeric_limits<double>::infinity();
                bool x_staged_valid = false;
                bool x_prior_valid = false;
                std::vector<SoftVec> soft_z_channel_prior;
                if (cascade_mode) {
                    if (joint_factorized_mode) {
                        // Both CSS factor graphs update in each iteration and
                        // exchange q-valued messages through the depolarizing
                        // channel factor. Unlike the staged decoder, neither
                        // component is hard-decided before decoding the other.
                        const char* iteration_setting = std::getenv("GF256_JOINT_ITERS");
                        decoder.max_iters = iteration_setting
                            ? std::stoi(iteration_setting) : 12;
                        if (decoder.max_iters <= 0 || decoder.max_iters > 200)
                            throw std::invalid_argument("GF256_JOINT_ITERS must be in 1..200");
                        auto joint = decoder.joint_decode_factorized(
                            hx_logical, hz_logical, logical_syndrome_x,
                            logical_syndrome_z, p);
                        corrected_x = std::move(joint.first);
                        corrected_z = std::move(joint.second);
                    } else if (full_hz_mode) {
                        // Experimental comparison path: every frame uses BP
                        // on the complete H_Z factor graph, so upper and lower
                        // checks exchange messages in the same iteration.
                        ++x_full_hz_bp_attempts;
                        vector<int> candidate = decoder.decode_component(
                            hz_logical, logical_syndrome_x, pauli_channel.marginalX(), true);
                        if (decoder.check_x_syndrome(
                                hz_logical,candidate,logical_syndrome_x)) {
                            corrected_x=std::move(candidate);
                            x_prior_valid=true;
                            ++x_full_hz_bp_successes;
                        } else {
                            corrected_x=std::move(candidate);
                        }
                    } else {
                        // Decode the five upper-selector slices, then recover
                        // the sixth symbol from the lower H_2 check.
                        ++cascade_staged_attempts;
                        vector<int> staged_x;
                        x_staged_valid = decoder.decode_cascade_x(
                            cascade_outer_hz_check, hz_logical, inner_length,
                            logical_syndrome_x, pauli_channel.marginalX(), staged_x);
                        if (x_staged_valid) {
                            ++cascade_staged_successes;
                            ++cascade_trailing_recoveries;
                            corrected_x = std::move(staged_x);
                            x_prior_valid = true;
                        } else {
                            corrected_x.assign(symbol_count, 0);
                            if (conditioned_mode) {
                                ++x_fallback_attempts;
                                vector<int> candidate = decoder.decode_component(
                                    hz_logical, logical_syndrome_x, pauli_channel.marginalX(), true);
                                if (decoder.check_x_syndrome(
                                        hz_logical,candidate,logical_syndrome_x)) {
                                    corrected_x=std::move(candidate);
                                    x_prior_valid=true;
                                    ++x_fallback_successes;
                                }
                            }
                        }
                    }
                    if(soft_x_mode || structured_soft_x_mode) {
                        ++x_soft_posterior_attempts;
                        const int q=gf.getSize();
                        std::vector<SoftVec> x_channel_priors(
                            symbol_count,SoftVec(q,0.0));
                        const double px=pauli_channel.marginalX();
                        for(int v=0;v<symbol_count;++v)
                            for(int x=0;x<q;++x) {
                                const int bits=basis.toDual(x);
                                double probability=1.0;
                                for(int bit=0;bit<field_degree;++bit)
                                    probability*=((bits>>bit)&1)?px:1.0-px;
                                x_channel_priors[v][x]=probability;
                            }
                        const auto x_soft=structured_soft_x_mode
                            ? runCascadeXSoftBP(decoder,cascade_outer_hz_check,
                                hz_logical,logical_syndrome_x,x_channel_priors,
                                20,20,10,gf)
                            : runOuterSoftBP(decoder,hz_logical,
                                logical_syndrome_x,x_channel_priors,20,gf);
                        if(x_soft.valid) {
                            corrected_x=x_soft.hard;
                            x_prior_valid=true;
                            soft_z_channel_prior=softXPosteriorToZPrior(
                                x_soft.posterior,basis,gf,pauli_channel);
                            ++x_soft_posterior_successes;
                        }
                    }
                    // The conditioned mode transfers the decoded physical X
                    // bits through the depolarizing-channel prior for Z.
                    if (construction == "structured" || conditioned_mode) {
                        if (full_hx_bp_mode) {
                            // Controlled ablation: keep the same HZ/X stage and
                            // conditional channel prior, but decode Z directly
                            // on the complete HX Tanner graph, with no outer
                            // projection decoder or structured fallback.
                            std::vector<SoftVec> priors(
                                symbol_count,SoftVec(gf.getSize()));
                            for (int v=0;v<symbol_count;++v) {
                                const int xb=x_prior_valid
                                    ? basis.toDual(corrected_x[v]) : 0;
                                for (int z=0;z<gf.getSize();++z) {
                                    double probability=1.0;
                                    for (int bit=0;bit<field_degree;++bit) {
                                        const double pz=x_prior_valid
                                            ? pauli_channel.zGivenX((xb>>bit)&1)
                                            : pauli_channel.marginalZ();
                                        probability*=((z>>bit)&1)
                                            ? pz : 1.0-pz;
                                    }
                                    priors[v][z]=probability;
                                }
                                normalizeLocal(priors[v]);
                            }
                            corrected_z=runOuterSoftBP(
                                decoder,hx_logical,logical_syndrome_z,
                                priors,full_hx_bp_iterations,gf).hard;
                        } else if (independent_projection_mode) {
                            // Decode each outer projection once from its
                            // channel marginal, then reconstruct the shared
                            // physical symbols without Turbo message exchange.
                            corrected_z = decodeHxCascadeTurbo(decoder,
                                cascade_outer_hx_check,hx_logical,logical_syndrome_z,
                                p,x_prior_valid ? &corrected_x : nullptr,basis,gf,
                                1,independent_outer_iterations,0.7,true,
                                nullptr,false,nullptr,false,
                                nullptr,false,false,true);
                        } else {
                            const auto* z_soft_prior=soft_z_channel_prior.empty()
                                ? nullptr : &soft_z_channel_prior;
                            const auto* x_hard_prior=(z_soft_prior==nullptr &&
                                conditioned_mode && x_prior_valid)
                                ? &corrected_x : nullptr;
                            corrected_z = decodeHxCascadeTurbo(decoder,
                                cascade_outer_hx_check, hx_logical, logical_syndrome_z,
                                p, x_hard_prior, basis, gf,
                                10,
                                tuned_turbo_mode ? 20 : 10,
                                tuned_turbo_mode ? 0.2 : 0.7,
                                true,nullptr,projection_mode,
                                projection_diagnostic_mode ? &initial_projection_diagnostics : nullptr,
                                joint_osd_mode,joint_osd_mode ? &current_joint_score : nullptr,
                                refresh_outer_mode,false,false,nullptr,false,z_soft_prior,
                                nullptr,layered_outer_schedule_mode,
                                overcomplete_decoder_mode ? &redundant_outer_checks : nullptr,
                                stabilizer_inactivation_mode ? &cascade_outer_hz_check : nullptr,
                                stabilizer_inactivation_mode ? 12 : 0,
                                stabilizer_inactivation_mode ? &initial_projection_diagnostics : nullptr);
                            if (conditioned_mode && !decoder.check_z_syndrome(
                                    hx_logical,corrected_z,logical_syndrome_z)) {
                                ++z_soft_fallback_attempts;
                                projection_retry_attempted=true;
                                vector<int> retry = decodeHxCascadeTurbo(decoder,
                                    cascade_outer_hx_check,hx_logical,logical_syndrome_z,
                                    p,x_hard_prior,basis,gf,
                                    long_turbo_mode ? 60 : (tuned_turbo_mode ? 20 : 20),
                                    20,
                                    long_turbo_mode ? 0.4 : (tuned_turbo_mode ? 0.2 : 0.7),
                                    true,nullptr,projection_mode,
                                    projection_diagnostic_mode ? &retry_projection_diagnostics : nullptr,
                                    joint_osd_mode,joint_osd_mode ? &retry_joint_score : nullptr,
                                    refresh_outer_mode,false,false,nullptr,false,z_soft_prior,
                                    nullptr,layered_outer_schedule_mode,
                                    overcomplete_decoder_mode ? &redundant_outer_checks : nullptr,
                                    stabilizer_inactivation_mode ? &cascade_outer_hz_check : nullptr,
                                    stabilizer_inactivation_mode ? 12 : 0,
                                    stabilizer_inactivation_mode ? &retry_projection_diagnostics : nullptr);
                                if (decoder.check_z_syndrome(
                                        hx_logical,retry,logical_syndrome_z)) {
                                    corrected_z=std::move(retry);
                                    if(joint_osd_mode) current_joint_score=retry_joint_score;
                                    ++z_soft_fallback_successes;
                                }
                            }
                        }
                    } else {
                        corrected_z = decoder.decode_component(
                            hx_logical, logical_syndrome_z, 2.0 * p / 3.0, false);
                    }
                } else if (construction == "tit") {
                    corrected_x = decoder.decode_component(
                        hz_logical, logical_syndrome_x, 2.0 * p / 3.0, true);
                    corrected_z = decoder.decode_component(
                        hx_logical, logical_syndrome_z, 2.0 * p / 3.0, false);
                } else {
                    const auto [estimate_x, estimate_z] = decoder.joint_decode(
                        hx_logical, hz_logical, syndrome_x, syndrome_z, p);
                    corrected_x = estimate_x;
                    corrected_z = estimate_z;
                }
                if (construction != "tit" && construction != "structured" &&
                    !joint_factorized_mode && !conditioned_mode &&
                    !decoder.check_z_syndrome(hx_logical, corrected_z, logical_syndrome_z)) {
                    ++z_soft_fallback_attempts;
                    const vector<int> soft_z = decoder.decode_component(
                        hx_logical, logical_syndrome_z, 2.0 * p / 3.0, false);
                    if (decoder.check_z_syndrome(hx_logical, soft_z, logical_syndrome_z)) {
                        corrected_z = soft_z;
                        ++z_soft_fallback_successes;
                    }
                }
                if (robust_mode) {
                    ++robust_attempts;
                    const auto contender=decodeHxCascadeTurbo(decoder,
                        cascade_outer_hx_check,hx_logical,logical_syndrome_z,p,
                        x_prior_valid ? &corrected_x : nullptr,basis,gf,
                        10,10,0.3,false,nullptr,projection_mode,
                        projection_diagnostic_mode ? &contender_projection_diagnostics : nullptr,
                        joint_osd_mode,joint_osd_mode ? &contender_joint_score : nullptr,
                        refresh_outer_mode);
                    const bool contender_valid=decoder.check_z_syndrome(
                        hx_logical,contender,logical_syndrome_z);
                    const bool current_valid=decoder.check_z_syndrome(
                        hx_logical,corrected_z,logical_syndrome_z);
                    auto score=[&](const vector<int>& candidate) {
                        double value=0.0;
                        for (int j=0;j<symbol_count;++j) {
                            const int xb=x_prior_valid ? basis.toDual(corrected_x[j]) : 0;
                            for (int bit=0;bit<field_degree;++bit) {
                                const double prior=x_prior_valid
                                    ? pauli_channel.zGivenX((xb>>bit)&1)
                                    : pauli_channel.marginalZ();
                                value+=std::log(((candidate[j]>>bit)&1)
                                    ? prior : 1.0-prior);
                            }
                        }
                        return value;
                    };
                    const double score_gain=contender_valid && current_valid
                        ? (joint_osd_mode
                            ? contender_joint_score-current_joint_score
                            : score(contender)-score(corrected_z)) : 0.0;
                    const bool should_replace = preserve_valid_joint_projection_mode
                        ? !current_valid
                        : (!current_valid || score_gain>1e-9);
                    if (contender_valid && should_replace) {
                        cerr << "[robust-selected] trial=" << trial
                             << " prior_syndrome=" << current_valid
                             << " score_gain=" << score_gain << '\n';
                        corrected_z=contender;
                        if(joint_osd_mode) current_joint_score=contender_joint_score;
                        ++robust_selections;
                    }
                }
                if (joint_projection_mode && !decoder.check_z_syndrome(
                        hx_logical,corrected_z,logical_syndrome_z)) {
                    ProjectionDiagnostics joint_diagnostics;
                    double joint_score=0.0;
                    auto joint_candidate=decodeHxCascadeTurbo(decoder,
                        cascade_outer_hx_check,hx_logical,logical_syndrome_z,p,
                        x_prior_valid ? &corrected_x : nullptr,basis,gf,
                        10,10,0.3,false,nullptr,true,&joint_diagnostics,
                        true,&joint_score,false,joint_list_mode);
                    const bool joint_valid=decoder.check_z_syndrome(
                        hx_logical,joint_candidate,logical_syndrome_z);
                    cerr << "[joint-projection] trial=" << trial
                         << " syndrome=" << joint_valid
                         << " bp_valid=" << joint_diagnostics.bp_valid
                         << " osd_repairs=" << joint_diagnostics.repaired
                         << " candidates=" << joint_diagnostics.joint_candidates_tested
                         << " updates=" << joint_diagnostics.joint_updates
                         << " score_gain=" << joint_diagnostics.joint_log_gain
                         << " joint_score=" << joint_score << '\n';
                    if (joint_valid) corrected_z=std::move(joint_candidate);
                }
                if (joint_osd_mode) {
                    // A final coupled candidate over the complete tensor
                    // constraints. This is algebraic bounded OSD-1, not BP on the
                    // short-cycle graph; it uses only syndrome and channel
                    // priors, never the sampled error.
                    std::vector<SoftVec> physical_prior(
                        symbol_count,SoftVec(gf.getSize()));
                    for(int v=0;v<symbol_count;++v) {
                        const int xb=x_prior_valid ? basis.toDual(corrected_x[v]) : 0;
                        for(int z=0;z<gf.getSize();++z) {
                            double probability=1.0;
                            for(int bit=0;bit<field_degree;++bit) {
                                const double pz=x_prior_valid
                                    ? (((xb>>bit)&1) ? 0.5 :
                                       (p/3.0)/(1.0-2.0*p/3.0))
                                    : 2.0*p/3.0;
                                probability*=((z>>bit)&1) ? pz : 1.0-pz;
                            }
                            physical_prior[v][z]=probability;
                        }
                        normalizeLocal(physical_prior[v]);
                    }
                    vector<int> channel_hard(symbol_count);
                    for(int v=0;v<symbol_count;++v)
                        channel_hard[v]=std::distance(
                            physical_prior[v].begin(),
                            std::max_element(physical_prior[v].begin(),
                                             physical_prior[v].end()));
                    const auto full_candidate=fullMatrixOsd1(
                        hx_logical,logical_syndrome_z,channel_hard,
                        physical_prior,gf,16,1,6,true);
                    auto physical_score=[&](const vector<int>& word) {
                        double value=0.0;
                        for(int v=0;v<symbol_count;++v)
                            value+=std::log(std::max(1e-300,physical_prior[v][word[v]]));
                        return value;
                    };
                    const bool previous_valid=decoder.check_z_syndrome(
                        hx_logical,corrected_z,logical_syndrome_z);
                    const double gain=full_candidate.syndrome_ok
                        ? full_candidate.log_posterior-physical_score(corrected_z)
                        : -std::numeric_limits<double>::infinity();
                    const bool selected=full_candidate.syndrome_ok &&
                        (!previous_valid || gain>1e-9);
                    if(selected) corrected_z=full_candidate.estimate;
                    cerr << "[joint-full-osd] trial=" << trial
                         << " valid=" << full_candidate.syndrome_ok
                         << " selected=" << selected
                         << " log_gain=" << gain << '\n';
                }
                // Conservative adaptive fallback: retain the baseline result
                // whenever it already satisfies the syndrome. Only if all
                // baseline schedules fail, try the lower-damping Turbo
                // schedule that rescued the diagnosed hard frame. This can
                // turn a baseline failure into a success, but cannot replace
                // a baseline-valid result.
                if (adaptive_turbo_mode && !decoder.check_z_syndrome(
                        hx_logical,corrected_z,logical_syndrome_z)) {
                    ++z_tuned_fallback_attempts;
                    auto tuned=decodeHxCascadeTurbo(decoder,
                        cascade_outer_hx_check,hx_logical,logical_syndrome_z,p,
                        x_prior_valid ? &corrected_x : nullptr,basis,gf,
                        10,20,0.2,false);
                    bool tuned_valid=decoder.check_z_syndrome(
                        hx_logical,tuned,logical_syndrome_z);
                    if (!tuned_valid) {
                        ++z_tuned_fallback_attempts;
                        tuned=decodeHxCascadeTurbo(decoder,
                            cascade_outer_hx_check,hx_logical,logical_syndrome_z,p,
                            x_prior_valid ? &corrected_x : nullptr,basis,gf,
                            20,20,0.2,false);
                        tuned_valid=decoder.check_z_syndrome(
                            hx_logical,tuned,logical_syndrome_z);
                    }
                    if (tuned_valid) {
                        corrected_z=std::move(tuned);
                        ++z_tuned_fallback_successes;
                    }
                }
                if (weak_osd2_fallback_mode && !decoder.check_z_syndrome(
                        hx_logical,corrected_z,logical_syndrome_z)) {
                    // Preserve all baseline adaptive schedules first. Only
                    // then try weak-pivot OSD-2, so this experiment cannot
                    // replace a baseline syndrome-valid correction.
                    ++z_osd_fallback_attempts;
                    std::vector<SoftVec> turbo_posterior;
                    const auto turbo_hard=decodeHxCascadeTurbo(
                        decoder,cascade_outer_hx_check,hx_logical,
                        logical_syndrome_z,p,
                        x_prior_valid ? &corrected_x : nullptr,
                        basis,gf,20,20,0.7,false,&turbo_posterior);
                    const auto osd=fullMatrixOsd1(
                        hx_logical,logical_syndrome_z,turbo_hard,
                        turbo_posterior,gf,12,2,8,true);
                    const bool previous_valid=decoder.check_z_syndrome(
                        hx_logical,corrected_z,logical_syndrome_z);
                    if (osd.syndrome_ok && !previous_valid) {
                        corrected_z=osd.estimate;
                        ++z_osd_fallback_successes;
                    }
                    cerr << "[Z-OSD2-weak-fallback] trial=" << trial
                         << " syndrome=" << osd.syndrome_ok
                         << " selected=" << (osd.syndrome_ok && !previous_valid)
                         << " candidates=" << osd.candidates_tested
                         << " log_posterior=" << osd.log_posterior << '\n';
                }
                if (guided_guess_mode && !decoder.check_z_syndrome(
                        hx_logical,corrected_z,logical_syndrome_z)) {
                    // Bounded, syndrome-guided guessing. This is an adaptation
                    // inspired by BPGD, not a reproduction of its full algorithm.
                    std::vector<SoftVec> marginals;
                    auto seed_candidate=decodeHxCascadeTurbo(decoder,
                        cascade_outer_hx_check,hx_logical,logical_syndrome_z,p,
                        x_prior_valid ? &corrected_x : nullptr,basis,gf,
                        10,20,0.2,false,&marginals,false,nullptr,false,nullptr,
                        false,false,false,nullptr,true);
                    std::vector<int> touched(symbol_count,0);
                    for(int r=0;r<hx_logical.rows();++r) {
                        int sum=0;
                        for(SparseMatrix<int,RowMajor>::InnerIterator it(hx_logical,r);it;++it)
                            sum=gf.add(sum,gf.mul(it.value(),seed_candidate[it.col()]));
                        if(sum!=logical_syndrome_z[r])
                            for(SparseMatrix<int,RowMajor>::InnerIterator it(hx_logical,r);it;++it)
                                ++touched[it.col()];
                    }
                    std::vector<std::tuple<double,int,int,double>> choices;
                    for(int v=0;v<symbol_count;++v) if(touched[v])
                        for(int b=0;b<field_degree;++b) {
                            double one=0;
                            for(int z=0;z<gf.getSize();++z)
                                if((z>>b)&1) one+=marginals[v][z];
                            choices.emplace_back(std::abs(one-0.5),v,b,one);
                        }
                    std::sort(choices.begin(),choices.end());
                    bool accepted=false;
                    // Four uncertain physical bits, two branches each.
                    for(size_t c=0;c<std::min<size_t>(4,choices.size()) && !accepted;++c) {
                        const auto [uncertainty,v,b,one]=choices[c];
                        for(int branch=0;branch<2 && !accepted;++branch) {
                            const int bit_value=(one>=0.5 ? 1 : 0)^branch;
                            std::vector<std::pair<int,int>> constraints(symbol_count,{0,0});
                            constraints[v]={1<<b,bit_value<<b};
                            auto candidate=decodeHxCascadeTurbo(decoder,
                                cascade_outer_hx_check,hx_logical,logical_syndrome_z,p,
                                x_prior_valid ? &corrected_x : nullptr,basis,gf,
                                10,20,0.2,false,nullptr,false,nullptr,false,nullptr,
                                false,false,false,&constraints,true);
                            accepted=decoder.check_z_syndrome(
                                hx_logical,candidate,logical_syndrome_z);
                            cerr << "[guided-guess] trial=" << trial << " variable=" << v
                                 << " bit=" << b << " value=" << bit_value
                                 << " posterior_one=" << one << " valid=" << accepted << '\n';
                            if(accepted) corrected_z=std::move(candidate);
                        }
                    }
                }
                if (guided_bpgd_mode && !decoder.check_z_syndrome(
                        hx_logical,corrected_z,logical_syndrome_z)) {
                    // Sequential bit fixing, following the BPGD principle:
                    // fix the most reliable unfixed bit touching an
                    // unsatisfied check, rerun the coupled decoder, then
                    // recompute reliabilities from its new marginals.
                    std::vector<std::pair<int,int>> constraints(symbol_count,{0,0});
                    std::vector<SoftVec> marginals;
                    auto candidate=decodeHxCascadeTurbo(decoder,
                        cascade_outer_hx_check,hx_logical,logical_syndrome_z,p,
                        x_prior_valid ? &corrected_x : nullptr,basis,gf,
                        10,20,0.2,false,&marginals,false,nullptr,false,nullptr,
                        false,false,false,nullptr,true);
                    bool accepted=decoder.check_z_syndrome(
                        hx_logical,candidate,logical_syndrome_z);
                    for(int step=0;step<8 && !accepted;++step) {
                        std::vector<int> touched(symbol_count,0);
                        for(int r=0;r<hx_logical.rows();++r) {
                            int sum=0;
                            for(SparseMatrix<int,RowMajor>::InnerIterator it(hx_logical,r);it;++it)
                                sum=gf.add(sum,gf.mul(it.value(),candidate[it.col()]));
                            if(sum!=logical_syndrome_z[r])
                                for(SparseMatrix<int,RowMajor>::InnerIterator it(hx_logical,r);it;++it)
                                    ++touched[it.col()];
                        }
                        double best_reliability=-1;
                        int best_v=-1,best_b=-1,best_value=0;
                        for(int v=0;v<symbol_count;++v) if(touched[v])
                            for(int b=0;b<field_degree;++b) {
                                if(constraints[v].first & (1<<b)) continue;
                                double probability_one=0;
                                for(int z=0;z<gf.getSize();++z)
                                    if((z>>b)&1) probability_one+=marginals[v][z];
                                const double reliability=std::abs(probability_one-0.5);
                                if(reliability>best_reliability) {
                                    best_reliability=reliability; best_v=v; best_b=b;
                                    best_value=probability_one>=0.5 ? 1 : 0;
                                }
                            }
                        if(best_v<0) break;
                        constraints[best_v]={1<<best_b,best_value<<best_b};
                        candidate=decodeHxCascadeTurbo(decoder,
                            cascade_outer_hx_check,hx_logical,logical_syndrome_z,p,
                            x_prior_valid ? &corrected_x : nullptr,basis,gf,
                            10,20,0.2,false,&marginals,false,nullptr,false,nullptr,
                            false,false,false,&constraints,true);
                        accepted=decoder.check_z_syndrome(
                            hx_logical,candidate,logical_syndrome_z);
                        cerr << "[guided-decimation] trial=" << trial << " step=" << step
                             << " variable=" << best_v << " bit=" << best_b
                             << " value=" << best_value << " reliability="
                             << best_reliability << " valid=" << accepted << '\n';
                    }
                    if(accepted) corrected_z=std::move(candidate);
                }
                if (restart_belief_mode && !decoder.check_z_syndrome(
                        hx_logical,corrected_z,logical_syndrome_z)) {
                    // Restart-Belief-inspired outer fallback. Unlike the
                    // existing BPGD experiment above, this ranks candidate
                    // error bits by ascending channel-output LLR and forces
                    // them to one, restarting the coupled Turbo decoder after
                    // each accumulated set of guesses. It is experimental:
                    // the paper's distance-based stopping rule is not used
                    // because this cascade code's distance is not established.
                    const char* branch_env=std::getenv(
                        "GF256_CASCADE_RB_BRANCHES");
                    const char* depth_env=std::getenv(
                        "GF256_CASCADE_RB_MAX_STEPS");
                    const int branch_limit=branch_env ? std::stoi(branch_env) : 4;
                    const int max_steps=depth_env ? std::stoi(depth_env) : 4;
                    if(branch_limit<1 || branch_limit>64 ||
                       max_steps<1 || max_steps>32)
                        throw std::invalid_argument(
                            "Restart-Belief branches must be 1..64 and steps 1..32");
                    std::vector<SoftVec> root_marginals;
                    auto root_candidate=decodeHxCascadeTurbo(decoder,
                        cascade_outer_hx_check,hx_logical,logical_syndrome_z,p,
                        x_prior_valid ? &corrected_x : nullptr,basis,gf,
                        10,20,0.2,false,&root_marginals,false,nullptr,false,nullptr,
                        false,false,false,nullptr,true);
                    int branches_started=0, bp_calls=0, valid_candidates=0;
                    double best_score=-std::numeric_limits<double>::infinity();
                    std::vector<int> best_candidate;
                    auto candidate_score=[&](const std::vector<int>& word) {
                        double score=0.0;
                        for(int v=0;v<symbol_count;++v) {
                            const int xb=x_prior_valid
                                ? basis.toDual(corrected_x[v]) : 0;
                            const int physical=gf.mul(
                                word[v],gf.inv(inner_scales[v%inner_length]));
                            for(int bit=0;bit<field_degree;++bit) {
                                const double pz=x_prior_valid
                                    ? pauli_channel.zGivenX((xb>>bit)&1)
                                    : pauli_channel.marginalZ();
                                score+=std::log(((physical>>bit)&1)
                                    ? pz : 1.0-pz);
                            }
                        }
                        return score;
                    };
                    if(decoder.check_z_syndrome(
                            hx_logical,root_candidate,logical_syndrome_z)) {
                        best_candidate=root_candidate;
                        best_score=candidate_score(root_candidate);
                        ++valid_candidates;
                    }
                    auto unsatisfied_symbols=[&](const std::vector<int>& word) {
                        std::vector<int> touched(symbol_count,0);
                        for(int row=0;row<hx_logical.rows();++row) {
                            int value=0;
                            for(SparseMatrix<int,RowMajor>::InnerIterator it(
                                    hx_logical,row);it;++it)
                                value=gf.add(value,gf.mul(it.value(),word[it.col()]));
                            if(value!=logical_syndrome_z[row])
                                for(SparseMatrix<int,RowMajor>::InnerIterator it(
                                        hx_logical,row);it;++it)
                                    ++touched[it.col()];
                        }
                        return touched;
                    };
                    struct RbChoice { double llr; int variable; int bit; };
                    auto least_llr_choice=[&](const std::vector<int>& word,
                            const std::vector<SoftVec>& marginals,
                            const std::vector<int>* touched,
                            const std::vector<std::pair<int,int>>& fixed) {
                        std::vector<RbChoice> choices;
                        for(int v=0;v<symbol_count;++v)
                            if(!touched || (*touched)[v])
                            for(int bit=0;bit<field_degree;++bit) {
                                if(fixed[v].first&(1<<bit)) continue;
                                double p_one=0.0;
                                for(int z=0;z<gf.getSize();++z) {
                                    const int physical=gf.mul(
                                        z,gf.inv(inner_scales[v%inner_length]));
                                    if((physical>>bit)&1) p_one+=marginals[v][z];
                                }
                                p_one=std::clamp(p_one,1e-12,1.0-1e-12);
                                choices.push_back({std::log((1.0-p_one)/p_one),v,bit});
                            }
                        std::stable_sort(choices.begin(),choices.end(),
                            [](const RbChoice& a,const RbChoice& b) {
                                if(a.llr!=b.llr) return a.llr<b.llr;
                                if(a.variable!=b.variable) return a.variable<b.variable;
                                return a.bit<b.bit;
                            });
                        (void)word;
                        return choices;
                        };
                    std::vector<std::pair<int,int>> no_constraints(
                        symbol_count,{0,0});
                    const auto root_touched=unsatisfied_symbols(root_candidate);
                    const auto roots=least_llr_choice(root_candidate,root_marginals,
                        restart_belief_global_mode ? nullptr : &root_touched,
                        no_constraints);
                    for(size_t root=0;root<std::min<size_t>(
                            static_cast<size_t>(branch_limit),roots.size());++root) {
                        ++branches_started;
                        std::vector<std::pair<int,int>> constraints(
                            symbol_count,{0,0});
                        const auto first=roots[root];
                        constraints[first.variable]={1<<first.bit,1<<first.bit};
                        for(int step=0;step<max_steps;++step) {
                            std::vector<SoftVec> marginals;
                            auto candidate=decodeHxCascadeTurbo(decoder,
                                cascade_outer_hx_check,hx_logical,
                                logical_syndrome_z,p,
                                x_prior_valid ? &corrected_x : nullptr,basis,gf,
                                10,20,0.2,false,&marginals,false,nullptr,false,nullptr,
                                false,false,false,&constraints,true);
                            ++bp_calls;
                            bool respects=true;
                            for(int v=0;v<symbol_count && respects;++v) {
                                const auto [mask,value]=constraints[v];
                                const int physical=gf.mul(candidate[v],
                                    gf.inv(inner_scales[v%inner_length]));
                                respects=((physical^value)&mask)==0;
                            }
                            if(!respects) break;
                            if(decoder.check_z_syndrome(
                                    hx_logical,candidate,logical_syndrome_z)) {
                                ++valid_candidates;
                                const double score=candidate_score(candidate);
                                if(score>best_score) {
                                    best_score=score;
                                    best_candidate=std::move(candidate);
                                }
                                break;
                            }
                            const auto touched=unsatisfied_symbols(candidate);
                            const auto choices=least_llr_choice(candidate,marginals,
                                restart_belief_global_mode ? nullptr : &touched,
                                constraints);
                            if(choices.empty()) break;
                            const auto next=choices.front();
                            constraints[next.variable].first|=1<<next.bit;
                            constraints[next.variable].second|=1<<next.bit;
                        }
                    }
                    const bool rb_selected=!best_candidate.empty();
                    if(rb_selected) corrected_z=std::move(best_candidate);
                    cerr << "[restart-belief-fallback] trial=" << trial
                         << " selection="
                         << (restart_belief_global_mode ? "global" : "unsatisfied")
                         << " branches=" << branches_started
                         << " bp_calls=" << bp_calls
                         << " syndrome_valid_candidates=" << valid_candidates
                         << " selected=" << rb_selected << '\n';
                }
                if (syndrome_flip_mode && !decoder.check_z_syndrome(
                        hx_logical,corrected_z,logical_syndrome_z)) {
                    // Chase-style GF(256) syndrome-flip fallback. A candidate
                    // physical-bit pattern t is mapped to the field-symbol
                    // coordinates; the shifted syndrome is s + H_X t. If BP
                    // returns e' for that syndrome, e'+t satisfies the original
                    // syndrome because the field has characteristic two.
                    // This is reliability-guided, low-weight BP-SF-inspired
                    // decoding; it does not reproduce the paper's oscillation
                    // ranking or its binary decoder.
                    const char* count_env=std::getenv("GF256_CASCADE_SF_BITS");
                    const int max_candidates=count_env ? std::stoi(count_env) : 8;
                    if(max_candidates<1 || max_candidates>256)
                        throw std::invalid_argument(
                            "GF256_CASCADE_SF_BITS must be in 1..256");
                    const char* weight_env=std::getenv(
                        "GF256_CASCADE_SF_MAX_WEIGHT");
                    const int max_weight=weight_env ? std::stoi(weight_env) : 1;
                    if(max_weight<1 || max_weight>2)
                        throw std::invalid_argument(
                            "GF256_CASCADE_SF_MAX_WEIGHT must be 1 or 2");
                    std::vector<SoftVec> root_marginals;
                    const auto root=decodeHxCascadeTurbo(decoder,
                        cascade_outer_hx_check,hx_logical,logical_syndrome_z,p,
                        x_prior_valid ? &corrected_x : nullptr,basis,gf,
                        10,20,0.2,false,&root_marginals);
                    struct SfChoice { double uncertainty; int variable; int bit; };
                    std::vector<SfChoice> choices;
                    choices.reserve(static_cast<size_t>(symbol_count*field_degree));
                    for(int v=0;v<symbol_count;++v) {
                        for(int bit=0;bit<field_degree;++bit) {
                            double probability_one=0.0;
                            for(int symbol=0;symbol<gf.getSize();++symbol) {
                                const int physical=gf.mul(symbol,
                                    gf.inv(inner_scales[v%inner_length]));
                                if((basis.toDual(physical)>>bit)&1)
                                    probability_one+=root_marginals[v][symbol];
                            }
                            choices.push_back({std::abs(probability_one-0.5),v,bit});
                        }
                    }
                    std::stable_sort(choices.begin(),choices.end(),
                        [](const SfChoice& a,const SfChoice& b) {
                            if(a.uncertainty!=b.uncertainty)
                                return a.uncertainty<b.uncertainty;
                            if(a.variable!=b.variable) return a.variable<b.variable;
                            return a.bit<b.bit;
                        });
                    auto physical_log_score=[&](const std::vector<int>& word) {
                        double score=0.0;
                        for(int v=0;v<symbol_count;++v) {
                            const int physical=gf.mul(word[v],
                                gf.inv(inner_scales[v%inner_length]));
                            const int zbits=basis.toDual(physical);
                            const int xbits=x_prior_valid
                                ? basis.toDual(corrected_x[v]) : 0;
                            for(int bit=0;bit<field_degree;++bit) {
                                const double pz=x_prior_valid
                                    ? pauli_channel.zGivenX((xbits>>bit)&1)
                                    : pauli_channel.marginalZ();
                                score+=std::log(((zbits>>bit)&1)
                                    ? pz : 1.0-pz);
                            }
                        }
                        return score;
                    };
                    int attempts=0,valid_candidates=0,selected_variable=-1,
                        selected_bit=-1;
                    double best_score=-std::numeric_limits<double>::infinity();
                    std::vector<int> best_candidate;
                    const size_t limit=std::min<size_t>(
                        static_cast<size_t>(max_candidates),choices.size());
                    std::vector<std::vector<SfChoice>> patterns;
                    patterns.reserve(limit+limit*(limit-1)/2);
                    for(size_t i=0;i<limit;++i)
                        patterns.push_back({choices[i]});
                    if(max_weight>=2)
                        for(size_t i=0;i<limit;++i)
                            for(size_t j=i+1;j<limit;++j)
                                patterns.push_back({choices[i],choices[j]});
                    for(const auto& pattern:patterns) {
                        std::vector<int> test_error(symbol_count,0);
                        for(const auto& choice:pattern) {
                            const int v=choice.variable;
                            const int bit=choice.bit;
                            const int delta=gf.mul(basis.fromDual(1<<bit),
                                inner_scales[v%inner_length]);
                            test_error[v]=gf.add(test_error[v],delta);
                        }
                        std::vector<int> shifted_syndrome=logical_syndrome_z;
                        for(int row=0;row<hx_logical.rows();++row) {
                            int shift=0;
                            for(SparseMatrix<int,RowMajor>::InnerIterator it(
                                    hx_logical,row);it;++it)
                                shift=gf.add(shift,gf.mul(it.value(),
                                    test_error[it.col()]));
                            shifted_syndrome[row]=gf.add(
                                shifted_syndrome[row],shift);
                        }
                        std::vector<SoftVec> shifted_marginals;
                        auto shifted=decodeHxCascadeTurbo(decoder,
                            cascade_outer_hx_check,hx_logical,shifted_syndrome,p,
                            x_prior_valid ? &corrected_x : nullptr,basis,gf,
                            10,20,0.2,false,&shifted_marginals);
                        ++attempts;
                        for(int j=0;j<symbol_count;++j)
                            shifted[j]=gf.add(shifted[j],test_error[j]);
                        if(!decoder.check_z_syndrome(
                                hx_logical,shifted,logical_syndrome_z))
                            continue;
                        ++valid_candidates;
                        const double score=physical_log_score(shifted);
                        if(score>best_score) {
                            best_score=score;
                            best_candidate=std::move(shifted);
                            selected_variable=pattern.front().variable;
                            selected_bit=pattern.front().bit;
                        }
                    }
                    const bool selected=!best_candidate.empty();
                    if(selected) corrected_z=std::move(best_candidate);
                    cerr << "[syndrome-flip] trial=" << trial
                         << " attempts=" << attempts
                         << " max_weight=" << max_weight
                         << " valid_candidates=" << valid_candidates
                         << " selected=" << selected
                         << " variable=" << selected_variable
                         << " bit=" << selected_bit << '\n';
                }
                if (layered_after_flip_mode && !decoder.check_z_syndrome(
                        hx_logical,corrected_z,logical_syndrome_z)) {
                    // Fallback uses the same tuned channel prior and outer
                    // layered schedule as conditioned_robust_layered. It is
                    // invoked only after tuned and syndrome-flip fail the
                    // full syndrome; selection never uses the sampled error.
                    int calls=0;
                    for(int rounds : {10,20}) {
                        ++calls;
                        ProjectionDiagnostics layered_diagnostics;
                        auto candidate=decodeHxCascadeTurbo(decoder,
                            cascade_outer_hx_check,hx_logical,logical_syndrome_z,
                            p,x_prior_valid ? &corrected_x : nullptr,
                            basis,gf,rounds,20,0.2,true,
                            nullptr,false,projection_diagnostic_mode ? &layered_diagnostics : nullptr,false,nullptr,false,
                            false,false,nullptr,false,nullptr,nullptr,true);
                        if(projection_diagnostic_mode)
                            layered_audit_diagnostics.emplace_back(rounds,std::move(layered_diagnostics));
                        if(decoder.check_z_syndrome(
                                hx_logical,candidate,logical_syndrome_z)) {
                            corrected_z=std::move(candidate);
                            break;
                        }
                    }
                    cerr << "[layered-after-flip] trial=" << trial
                         << " calls=" << calls
                         << " selected=" << decoder.check_z_syndrome(
                                hx_logical,corrected_z,logical_syndrome_z)
                         << '\n';
                }
                if (osd_fallback_mode && !decoder.check_z_syndrome(
                        hx_logical, corrected_z, logical_syndrome_z)) {
                    ++z_osd_fallback_attempts;
                    std::vector<SoftVec> conditional_z_priors(
                        symbol_count, SoftVec(gf.getSize(), 0.0));
                    for (int j = 0; j < symbol_count; ++j) {
                        const int xbits = x_prior_valid
                            ? basis.toDual(corrected_x[j]) : 0;
                        for (int z = 0; z < gf.getSize(); ++z) {
                            double probability = 1.0;
                            for (int bit = 0; bit < field_degree; ++bit) {
                                const double pz = x_prior_valid
                                    ? (((xbits >> bit) & 1) ? 0.5
                                       : (p / 3.0) / (1.0 - 2.0 * p / 3.0))
                                    : 2.0 * p / 3.0;
                                probability *= ((z >> bit) & 1)
                                    ? pz : 1.0 - pz;
                            }
                            conditional_z_priors[j][z] = probability;
                        }
                        normalizeLocal(conditional_z_priors[j]);
                    }
                    std::vector<SoftVec> turbo_posterior;
                    const auto turbo_candidate = decodeHxCascadeTurbo(
                        decoder, cascade_outer_hx_check, hx_logical,
                        logical_syndrome_z, p,
                        x_prior_valid ? &corrected_x : nullptr,
                        basis, gf, 20, 20, 0.7, false, &turbo_posterior);
                    const auto osd = fullMatrixOsd1(
                        hx_logical, logical_syndrome_z, turbo_candidate,
                        turbo_posterior, gf,
                        weak_osd2_fallback_mode ? 12 : 256,
                        weak_osd2_fallback_mode ? 2 : 1, 8,
                        weak_osd2_fallback_mode);
                    if (osd.syndrome_ok) {
                        corrected_z = osd.estimate;
                        ++z_osd_fallback_successes;
                    }
                    cerr << "[Z-OSD-fallback] trial=" << trial
                         << " weak_pivot_osd2=" << weak_osd2_fallback_mode
                         << " syndrome=" << osd.syndrome_ok
                         << " candidates=" << osd.candidates_tested
                         << " log_posterior=" << osd.log_posterior << '\n';
                }
                if(augmented_decoder_mode && !decoder.check_z_syndrome(
                        hx_logical,corrected_z,logical_syndrome_z)) {
                    const int m=cascade_outer_hx_check.rows();
                    const int duplicate_count=std::max(1,(m+9)/10);
                    const int checks_per_outer = m > 0 && hx_logical.rows() % m == 0
                        ? hx_logical.rows() / m : 0;
                    std::vector<double> row_weights(static_cast<size_t>(m), 1.0);
                    if (checks_per_outer > 0) {
                        for (int outer_row = 0; outer_row < m; ++outer_row) {
                            int violated_projections = 0;
                            for (int projection = 0;
                                 projection < checks_per_outer; ++projection) {
                                const int row = outer_row * checks_per_outer + projection;
                                int value = 0;
                                for (SparseMatrix<int,RowMajor>::InnerIterator it(
                                         hx_logical,row); it; ++it)
                                    value = gf.add(value,gf.mul(
                                        it.value(),corrected_z[it.col()]));
                                violated_projections +=
                                    value != logical_syndrome_z[row];
                            }
                            // Prioritize outer checks associated with the current
                            // syndrome residual while retaining nonzero mass for
                            // exploratory restarts. This uses decoder-visible
                            // syndrome information only, never the true error.
                            row_weights[outer_row] += 4.0 * violated_projections;
                        }
                    }
                    double best_score=-std::numeric_limits<double>::infinity();
                    std::vector<int> best_candidate;
                    for(int attempt=0;attempt<3;++attempt) {
                        ++z_augmented_attempts;
                        std::mt19937 augmentation_rng(
                            seed ^ (static_cast<unsigned int>(trial+1)*0x9e3779b9u) ^
                            (static_cast<unsigned int>(attempt+1)*0x85ebca6bu));
                        auto weights = row_weights;
                        std::vector<int> rows;
                        rows.reserve(static_cast<size_t>(duplicate_count));
                        for (int selected=0; selected<std::min(duplicate_count,m);
                             ++selected) {
                            std::discrete_distribution<int> choose(
                                weights.begin(),weights.end());
                            const int row=choose(augmentation_rng);
                            rows.push_back(row);
                            weights[row]=0.0;
                        }
                        auto candidate=decodeHxCascadeTurbo(
                            decoder,cascade_outer_hx_check,hx_logical,
                            logical_syndrome_z,p,
                            x_prior_valid ? &corrected_x : nullptr,
                            basis,gf,10,20,0.2,false,nullptr,false,nullptr,
                            false,nullptr,false,false,false,nullptr,false,nullptr,
                            &rows);
                        if(!decoder.check_z_syndrome(
                                hx_logical,candidate,logical_syndrome_z)) continue;
                        double score=0.0;
                        for(int v=0;v<symbol_count;++v) {
                            const int xb=x_prior_valid ? basis.toDual(corrected_x[v]) : 0;
                            for(int bit=0;bit<field_degree;++bit) {
                                const double pz=x_prior_valid
                                    ? pauli_channel.zGivenX((xb>>bit)&1)
                                    : pauli_channel.marginalZ();
                                score+=std::log(((candidate[v]>>bit)&1)
                                    ? pz : 1.0-pz);
                            }
                        }
                        if(score>best_score) {
                            best_score=score;
                            best_candidate=std::move(candidate);
                        }
                    }
                    const bool selected_augmented=!best_candidate.empty();
                    if(selected_augmented) {
                        corrected_z=std::move(best_candidate);
                        ++z_augmented_selections;
                    }
                    cerr << "[Z-augmented-retry] trial=" << trial
                         << " attempts=" << z_augmented_attempts
                         << " selected=" << selected_augmented << '\n';
                }
                if (target_diagnostics && conditioned_mode) {
                    const bool usable_x = decoder.check_x_syndrome(
                        hz_logical,corrected_x,logical_syndrome_x);
                    alternate_z = decodeHxCascadeTurbo(decoder,
                        cascade_outer_hx_check,hx_logical,logical_syndrome_z,
                        p,usable_x ? &corrected_x : nullptr,basis,gf,
                        20,20,0.7);
                }
                const bool x_syndrome_ok = decoder.check_x_syndrome(
                    hz_logical, corrected_x, logical_syndrome_x);
                const bool z_syndrome_ok = decoder.check_z_syndrome(
                    hx_logical, corrected_z, logical_syndrome_z);
                if (conditioned_mode && !z_syndrome_ok) {
                    int unsatisfied = 0;
                    for (int row = 0; row < hx_logical.rows(); ++row) {
                        int value = 0;
                        for (SparseMatrix<int, RowMajor>::InnerIterator it(
                                 hx_logical, row); it; ++it) {
                            value = gf.add(value, gf.mul(
                                it.value(), corrected_z[it.col()]));
                        }
                        unsatisfied += value != logical_syndrome_z[row];
                    }
                    cerr << "[Z-syndrome-failure] trial=" << trial
                         << " unsatisfied_checks=" << unsatisfied
                         << " robust=" << robust_mode << '\n';
                }
                const bool syndrome_ok = x_syndrome_ok && z_syndrome_ok;
                if (!x_syndrome_ok) ++x_syndrome_failures;
                if (!z_syndrome_ok) ++z_syndrome_failures;
                if (!syndrome_ok) {
                    ++syndrome_failures;
                }

                const vector<int> actual_x = binaryToFieldSymbols(error_x, field_degree, basis, true);
                const vector<int> actual_z = binaryToFieldSymbols(error_z, field_degree, basis, false);
                if (joint_osd_mode) {
                    auto physical_log_likelihood=[&](const vector<int>& z) {
                        double value=0.0;
                        for(int v=0;v<symbol_count;++v) {
                            const int xb=x_prior_valid ? basis.toDual(corrected_x[v]) : 0;
                            for(int bit=0;bit<field_degree;++bit) {
                                const double pz=x_prior_valid
                                    ? (((xb>>bit)&1) ? 0.5 :
                                       (p/3.0)/(1.0-2.0*p/3.0))
                                    : 2.0*p/3.0;
                                value+=std::log(((z[v]>>bit)&1) ? pz : 1.0-pz);
                            }
                        }
                        return value;
                    };
                    // Truth is read only after the decoder has selected its
                    // candidate, for diagnosis; it never influences decoding.
                    cerr << "[joint-oracle-diagnostic] trial=" << trial
                         << " candidate_minus_truth_loglik="
                         << physical_log_likelihood(corrected_z)-
                            physical_log_likelihood(actual_z) << '\n';
                }
                if (projection_diagnostic_mode) {
                    const int outer_blocks=cascade_outer_hx_check.cols();
                    const int outer_checks=cascade_outer_hx_check.rows();
                    if(outer_checks<=0 || hx_logical.rows()%outer_checks!=0)
                        throw std::runtime_error(
                            "Cannot infer G2 dimensions for projection diagnostics");
                    const int projection_count=hx_logical.rows()/outer_checks;
                    int reference_row=-1,reference_block=-1,reference_value=0;
                    for(int r=0;r<outer_checks && reference_row<0;++r)
                        for(SparseMatrix<int,RowMajor>::InnerIterator it(
                                cascade_outer_hx_check,r);it;++it) {
                            reference_row=r;
                            reference_block=it.col();
                            reference_value=it.value();
                            break;
                        }
                    if(reference_row<0)
                        throw std::runtime_error(
                            "Cannot infer G2 from zero outer matrix in diagnostics");
                    const int inverse_reference=gf.inv(reference_value);
                    vector<int> diagnostic_g2(
                        static_cast<size_t>(projection_count*inner_length),0);
                    for(int a=0;a<projection_count;++a)
                        for(int c=0;c<inner_length;++c)
                            diagnostic_g2[a*inner_length+c]=gf.mul(
                                inverse_reference,hx_logical.coeff(
                                    reference_row*projection_count+a,
                                    reference_block*inner_length+c));
                    for(int r=0;r<outer_checks;++r)
                        for(int block=0;block<outer_blocks;++block)
                            for(int a=0;a<projection_count;++a)
                                for(int c=0;c<inner_length;++c)
                                    if(hx_logical.coeff(
                                            r*projection_count+a,
                                            block*inner_length+c)!=gf.mul(
                                                cascade_outer_hx_check.coeff(r,block),
                                                diagnostic_g2[a*inner_length+c]))
                                        throw std::runtime_error(
                                            "Projection diagnostics require HX=H_outer kron G2");
                    for(int a=0;a<projection_count;++a) {
                        vector<int> true_projection(outer_blocks);
                        vector<int> decoded_projection(outer_blocks);
                        vector<int> projected_syndrome(outer_checks);
                        int differences=0;
                        for(int j=0;j<outer_blocks;++j) {
                            for(int c=0;c<inner_length;++c) {
                                const int coefficient=
                                    diagnostic_g2[a*inner_length+c];
                                true_projection[j]=gf.add(true_projection[j],
                                    gf.mul(coefficient,
                                        actual_z[j*inner_length+c]));
                                decoded_projection[j]=gf.add(decoded_projection[j],
                                    gf.mul(coefficient,
                                        corrected_z[j*inner_length+c]));
                            }
                            differences+=true_projection[j]!=decoded_projection[j];
                        }
                        for(int r=0;r<outer_checks;++r)
                            projected_syndrome[r]=
                                logical_syndrome_z[r*projection_count+a];
                        cerr << "[outer-projection] trial=" << trial
                             << " basis=matrix_G2"
                             << " coordinate=" << a
                             << " g2_row_weight=" << std::count_if(
                                    diagnostic_g2.begin()+a*inner_length,
                                    diagnostic_g2.begin()+(a+1)*inner_length,
                                    [](int value){return value!=0;})
                             << " true_syndrome=" << decoder.check_z_syndrome(
                                    cascade_outer_hx_check,true_projection,projected_syndrome)
                             << " decoded_syndrome=" << decoder.check_z_syndrome(
                                    cascade_outer_hx_check,decoded_projection,projected_syndrome)
                             << " symbol_differences=" << differences << '\n';
                    }
                    auto print_projection_diagnostics=[&](const char* stage,
                                                          const ProjectionDiagnostics& d) {
                        cerr << "[outer-decoder] trial=" << trial
                             << " stage=" << stage
                             << " bp_valid=" << d.bp_valid
                             << " raw_valid=" << d.raw_valid
                             << " final_valid=" << d.final_valid
                             << " outer_refreshed=" << d.outer_refreshed
                             << " repaired=" << d.repaired
                             << " repair_failed=" << d.repair_failed
                             << " hard_disagreements=" << d.hard_disagreements
                             << " joint_candidates_tested=" << d.joint_candidates_tested
                             << " joint_updates=" << d.joint_updates
                             << " si_attempts=" << d.inactivation_attempts
                             << " si_successes=" << d.inactivation_successes
                             << " si_selected_checks=" << d.inactivated_checks
                             << " joint_log_gain=" << d.joint_log_gain << '\n';
                        if(d.pair_audit_left!=-1)
                            cerr << "[pair-candidate-audit] trial=" << trial
                                 << " stage=" << stage
                                 << " left=" << d.pair_audit_left
                                 << " physical_candidates=" << d.pair_audit_candidates << '\n';
                        for(int variant=0;variant<3;++variant) {
                          const auto& audit_list=variant==2 ? d.pair_candidate_list :
                              variant==1 ? d.neighbor_candidate_list : d.outer_candidate_list;
                          const std::string audit_stage=std::string(stage)+
                              (variant==2 ? "_pair" : variant==1 ? "_neighbor" : "");
                          if(audit_list.size()==static_cast<size_t>(projection_count)) {
                            for(int a=0;a<projection_count;++a) {
                                vector<int> truth(outer_blocks,0);
                                vector<int> projected_syndrome(outer_checks);
                                for(int j=0;j<outer_blocks;++j)
                                    for(int c=0;c<inner_length;++c)
                                        truth[j]=gf.add(truth[j],gf.mul(
                                            diagnostic_g2[a*inner_length+c],
                                            actual_z[j*inner_length+c]));
                                for(int r=0;r<outer_checks;++r)
                                    projected_syndrome[r]=
                                        logical_syndrome_z[r*projection_count+a];
                                int min_symbol_distance=outer_blocks+1;
                                bool true_projection_covered=false;
                                bool true_outer_coset_covered=false;
                                const auto truth_coset=cascade_z_stabilizer_space->reduce(truth);
                                for(const auto& candidate:audit_list[a]) {
                                    int distance=0;
                                    for(int j=0;j<outer_blocks;++j)
                                        distance+=candidate[j]!=truth[j];
                                    min_symbol_distance=std::min(
                                        min_symbol_distance,distance);
                                    true_projection_covered|=distance==0;
                                    if(!true_outer_coset_covered)
                                        true_outer_coset_covered=
                                            cascade_z_stabilizer_space->reduce(candidate)==truth_coset;
                                }
                                cerr << "[outer-candidate-audit] trial=" << trial
                                     << " stage=" << audit_stage
                                     << " coordinate=" << a
                                     << " candidates="
                                     << audit_list[a].size()
                                     << " true_projection_covered="
                                     << true_projection_covered
                                     << " min_symbol_distance="
                                     << (audit_list[a].empty()
                                            ? -1 : min_symbol_distance)
                                     << " true_outer_coset_covered="
                                     << true_outer_coset_covered << '\n';
                            }
                        }
                        }
                    };
                    print_projection_diagnostics("initial",initial_projection_diagnostics);
                    if (projection_retry_attempted)
                        print_projection_diagnostics("retry",retry_projection_diagnostics);
                    print_projection_diagnostics("contender",contender_projection_diagnostics);
                    for(const auto& item:layered_audit_diagnostics) {
                        const auto stage="layered"+std::to_string(item.first);
                        print_projection_diagnostics(stage.c_str(),item.second);
                    }
                }
                if (target_diagnostics && conditioned_mode) {
                    auto channel_score = [&](const vector<int>& candidate) {
                        double score=0.0;
                        for (int j=0; j<symbol_count; ++j) {
                            const int xb=basis.toDual(corrected_x[j]);
                            for (int bit=0; bit<field_degree; ++bit) {
                                const double prior=((xb>>bit)&1) ? 0.5 :
                                    (p/3.0)/(1.0-2.0*p/3.0);
                                score+=std::log(((candidate[j]>>bit)&1) ? prior : 1.0-prior);
                            }
                        }
                        return score;
                    };
                    auto report_z_candidate = [&](const string& label,
                                                  const vector<int>& candidate) {
                        vector<int> residual(symbol_count);
                        int differing_symbols=0;
                        for (int j=0;j<symbol_count;++j) {
                            residual[j]=gf.add(actual_z[j],candidate[j]);
                            differing_symbols+=candidate[j]!=actual_z[j];
                        }
                        int unsatisfied_checks=0;
                        for (int row=0;row<hx_logical.rows();++row) {
                            int value=0;
                            for (SparseMatrix<int,RowMajor>::InnerIterator it(
                                     hx_logical,row);it;++it)
                                value=gf.add(value,gf.mul(
                                    it.value(),candidate[it.col()]));
                            unsatisfied_checks+=value!=logical_syndrome_z[row];
                        }
                        cerr << "[Z-candidate] trial=" << trial
                             << " name=" << label
                             << " syndrome=" << decoder.check_z_syndrome(
                                    hx_logical,candidate,logical_syndrome_z)
                             << " unsatisfied_checks=" << unsatisfied_checks
                             << " logical=" << !isCascadeZStabilizer(
                                    residual,*cascade_z_stabilizer_space,gf,inner_length,
                                    cascade_inner_dimension,inner_scales)
                             << " exact=" << (candidate==actual_z)
                             << " differing_symbols=" << differing_symbols << '\n';
                    };
                    const vector<int> full_hx_bp = decoder.decode_component(
                        hx_logical,logical_syndrome_z,2.0*p/3.0,false);
                    const vector<int> oracle_x_turbo = decodeHxCascadeTurbo(
                        decoder,cascade_outer_hx_check,hx_logical,logical_syndrome_z,
                        p,&actual_x,basis,gf);
                    std::vector<SoftVec> conditional_z_priors(
                        symbol_count,SoftVec(gf.getSize(),0.0));
                    for (int j=0;j<symbol_count;++j) {
                        const int xbits=basis.toDual(corrected_x[j]);
                        for (int z=0;z<gf.getSize();++z) {
                            double probability=1.0;
                            for (int bit=0;bit<field_degree;++bit) {
                                const double pz=((xbits>>bit)&1) ? 0.5 :
                                    (p/3.0)/(1.0-2.0*p/3.0);
                                probability*=((z>>bit)&1) ? pz : 1.0-pz;
                            }
                            conditional_z_priors[j][z]=probability;
                        }
                        normalizeLocal(conditional_z_priors[j]);
                    }
                    const auto full_hx_conditioned = runOuterSoftBP(
                        decoder,hx_logical,logical_syndrome_z,
                        conditional_z_priors,200,gf);
                    report_z_candidate("current_conditioned_robust",corrected_z);
                    report_z_candidate("full_HX_BP_marginal_prior",full_hx_bp);
                    report_z_candidate("turbo_oracle_true_X_prior",oracle_x_turbo);
                    report_z_candidate("full_HX_BP_decoded_X_prior",
                                       full_hx_conditioned.hard);
                    const auto repaired_z = repairSyndromeBeam(
                        hx_logical,logical_syndrome_z,corrected_z,
                        conditional_z_priors,gf,128,12);
                    report_z_candidate("syndrome_beam_repair",repaired_z.estimate);
                    std::vector<SoftVec> full_bp_posterior(
                        symbol_count,SoftVec(gf.getSize(),0.0));
                    for(int j=0;j<symbol_count;++j) {
                        full_bp_posterior[j]=conditional_z_priors[j];
                        const SoftVec& ext=full_hx_conditioned.extrinsic[j][0];
                        for(int z=0;z<gf.getSize();++z)
                            full_bp_posterior[j][z]*=ext[z];
                        normalizeLocal(full_bp_posterior[j]);
                    }
                    if (std::getenv("GF256_CASCADE_ENABLE_GIBBS_AUDIT") != nullptr) {
                        int gibbs_sweeps=30;
                        if (const char* value=std::getenv(
                                "GF256_CASCADE_GIBBS_AUDIT_SWEEPS"))
                            gibbs_sweeps=std::stoi(value);
                        if (gibbs_sweeps<12 || gibbs_sweeps>2000)
                            throw std::invalid_argument(
                                "Gibbs audit sweeps must be in [12,2000]");
                        const int gibbs_burn=std::max(10,gibbs_sweeps/4);
                        int gibbs_thin=std::max(1,gibbs_sweeps/50);
                        if (const char* value=std::getenv(
                                "GF256_CASCADE_GIBBS_AUDIT_SAMPLE_EVERY"))
                            gibbs_thin=std::stoi(value);
                        if (gibbs_thin<=0)
                            throw std::invalid_argument(
                                "Gibbs sample interval must be positive");
                        uint32_t gibbs_seed=20261110;
                        if (const char* value=std::getenv(
                                "GF256_CASCADE_GIBBS_AUDIT_SEED"))
                            gibbs_seed=static_cast<uint32_t>(std::stoul(value));
                        int gibbs_block_updates=0;
                        if (const char* value=std::getenv(
                                "GF256_CASCADE_GIBBS_AUDIT_BLOCK_UPDATES"))
                            gibbs_block_updates=std::stoi(value);
                        if (gibbs_block_updates<0 || gibbs_block_updates>4096)
                            throw std::invalid_argument(
                                "Gibbs block updates must be in [0,4096]");
                        int gibbs_logical_updates=0;
                        if (const char* value=std::getenv(
                                "GF256_CASCADE_GIBBS_LOGICAL_UPDATES"))
                            gibbs_logical_updates=std::stoi(value);
                        if(gibbs_logical_updates<0 || gibbs_logical_updates>4096)
                            throw std::invalid_argument(
                                "Gibbs logical updates must be in [0,4096]");
                        int gibbs_stabilizer_updates=0;
                        if (const char* value=std::getenv(
                                "GF256_CASCADE_GIBBS_STABILIZER_UPDATES"))
                            gibbs_stabilizer_updates=std::stoi(value);
                        if(gibbs_stabilizer_updates<0 || gibbs_stabilizer_updates>4096)
                            throw std::invalid_argument(
                                "Gibbs stabilizer updates must be in [0,4096]");
                        std::vector<double> gibbs_betas={1.0,0.5,0.2};
                        if (const char* value=std::getenv(
                                "GF256_CASCADE_GIBBS_BETAS")) {
                            gibbs_betas.clear();
                            std::stringstream beta_stream(value);
                            std::string beta_token;
                            while(std::getline(beta_stream,beta_token,','))
                                gibbs_betas.push_back(std::stod(beta_token));
                        }
                        std::vector<int> gibbs_initial=full_hx_conditioned.hard;
                        bool logical_search_osd_seed=false;
                        int logical_search_osd_seed_candidates=0;
                        if (std::getenv(
                                "GF256_CASCADE_GIBBS_LOGICAL_LIST_TOPK") != nullptr) {
                            const auto seed_osd=fullMatrixOsd1(
                                hx_logical,logical_syndrome_z,
                                full_hx_conditioned.hard,conditional_z_priors,gf,256);
                            logical_search_osd_seed_candidates=
                                seed_osd.candidates_tested;
                            if(seed_osd.syndrome_ok) {
                                gibbs_initial=seed_osd.estimate;
                                logical_search_osd_seed=true;
                            }
                        }
                        const auto gibbs=sampleGfSyndromeCosets(
                            hx_logical,logical_syndrome_z,
                            gibbs_initial,conditional_z_priors,
                            *cascade_z_stabilizer_space,inner_length,
                            cascade_inner_dimension,inner_scales,gf,
                            gibbs_sweeps,gibbs_burn,gibbs_thin,gibbs_seed,
                            gibbs_block_updates,gibbs_logical_updates,gibbs_betas,
                            cascade_mode ? &hz_logical : nullptr,
                            gibbs_stabilizer_updates);
                        report_z_candidate("syndrome_space_gibbs_coset_mode",
                                           gibbs.estimate);
                        std::vector<double> logical_list_scores;
                        std::vector<std::string> logical_list_keys;
                        int logical_list_evaluated=0;
                        int logical_list_syndrome_valid=0;
                        int logical_list_direction_count=0;
                        int logical_descent_passes=0;
                        int logical_descent_moves=0;
                        long long logical_descent_evaluated=0;
                        std::vector<int> logical_descent_word;
                        double logical_descent_score=-
                            std::numeric_limits<double>::infinity();
                        std::map<std::string,double> logical_beam_class_log_mass;
                        int logical_beam_bp_calls=0;
                        int logical_beam_valid_candidates=0;
                        int logical_beam_constraint_consistent=0;
                        long long logical_beam_constraint_violations=0;
                        int logical_beam_retained_nodes=0;
                        int logical_beam_leaf_repairs=0;
                        int logical_beam_leaf_repair_successes=0;
                        int logical_beam_width_used=0;
                        int logical_beam_depth_used=0;
                        int beam_turbo_rounds=3, beam_outer_iterations=10;
                        double beam_damping=0.3;
                        double beam_constraint_penalty=25.0;
                        if(const char* value=std::getenv(
                                "GF256_CASCADE_CONSTRAINT_PENALTY_EXPONENT"))
                            beam_constraint_penalty=std::stod(value);
                        if (std::getenv(
                                "GF256_CASCADE_GIBBS_LOGICAL_LIST_TOPK") != nullptr) {
                            if (!gibbs.syndrome_ok || gibbs.logical_directions.empty())
                                throw std::runtime_error(
                                    "Logical candidate search needs a syndrome-valid Gibbs base and logical directions");
                            int direction_limit=128;
                            if (const char* value=std::getenv(
                                    "GF256_CASCADE_GIBBS_LOGICAL_LIST_DIRECTIONS"))
                                direction_limit=std::stoi(value);
                            int top_k=256;
                            if (const char* value=std::getenv(
                                    "GF256_CASCADE_GIBBS_LOGICAL_LIST_TOPK"))
                                top_k=std::stoi(value);
                            if (direction_limit<=0 || direction_limit>
                                    static_cast<int>(gibbs.logical_directions.size()) ||
                                top_k<=0 || top_k>100000)
                                throw std::invalid_argument(
                                    "Logical list directions/top-k are outside supported bounds");
                            logical_list_direction_count=direction_limit;
                            const int q=gf.getSize();
                            auto symbol_log_score=[&](int variable,int value) {
                                return std::log(std::max(1e-300,
                                    conditional_z_priors[variable][value]));
                            };
                            struct DirectionQuality { double gain; size_t index; };
                            std::vector<DirectionQuality> quality;
                            quality.reserve(gibbs.logical_directions.size());
                            for(size_t d=0;d<gibbs.logical_directions.size();++d) {
                                const auto& direction=gibbs.logical_directions[d];
                                double best_gain=-std::numeric_limits<double>::infinity();
                                for(int bit=0;bit<field_degree;++bit) {
                                    const int alpha=1<<bit;
                                    double gain=symbol_log_score(direction.variable,
                                        gf.add(gibbs.estimate[direction.variable],alpha))-
                                        symbol_log_score(direction.variable,
                                            gibbs.estimate[direction.variable]);
                                    for(const auto& [variable,coefficient]:
                                            direction.pivot_updates) {
                                        const int old_value=gibbs.estimate[variable];
                                        const int new_value=gf.add(old_value,
                                            gf.mul(coefficient,alpha));
                                        gain+=symbol_log_score(variable,new_value)-
                                              symbol_log_score(variable,old_value);
                                    }
                                    best_gain=std::max(best_gain,gain);
                                }
                                quality.push_back({best_gain,d});
                            }
                            std::partial_sort(quality.begin(),
                                quality.begin()+direction_limit,quality.end(),
                                [](const auto& a,const auto& b) {
                                    return a.gain>b.gain;
                                });
                            struct LogicalCandidate {
                                double score;
                                int direction;
                                int alpha;
                            };
                            struct MinScoreFirst {
                                bool operator()(const LogicalCandidate& a,
                                                const LogicalCandidate& b) const {
                                    return a.score>b.score;
                                }
                            };
                            std::priority_queue<LogicalCandidate,
                                std::vector<LogicalCandidate>,MinScoreFirst> best;
                            double base_score=0.0;
                            for(int v=0;v<symbol_count;++v)
                                base_score+=symbol_log_score(v,gibbs.estimate[v]);
                            best.push({base_score,-1,0});
                            for(int selected=0;selected<direction_limit;++selected) {
                                const auto& direction=gibbs.logical_directions[
                                    quality[selected].index];
                                for(int alpha=1;alpha<q;++alpha) {
                                    double score=base_score;
                                    const int old_free=gibbs.estimate[direction.variable];
                                    const int new_free=gf.add(old_free,alpha);
                                    score+=symbol_log_score(direction.variable,new_free)-
                                           symbol_log_score(direction.variable,old_free);
                                    for(const auto& [variable,coefficient]:
                                            direction.pivot_updates) {
                                        const int old_value=gibbs.estimate[variable];
                                        const int new_value=gf.add(old_value,
                                            gf.mul(coefficient,alpha));
                                        score+=symbol_log_score(variable,new_value)-
                                               symbol_log_score(variable,old_value);
                                    }
                                    ++logical_list_evaluated;
                                    if(static_cast<int>(best.size())<top_k)
                                        best.push({score,static_cast<int>(quality[selected].index),alpha});
                                    else if(score>best.top().score) {
                                        best.pop();
                                        best.push({score,
                                            static_cast<int>(quality[selected].index),alpha});
                                    }
                                }
                            }
                            std::vector<LogicalCandidate> ranked;
                            ranked.reserve(best.size());
                            while(!best.empty()) {
                                ranked.push_back(best.top());
                                best.pop();
                            }
                            std::sort(ranked.begin(),ranked.end(),
                                [](const auto& a,const auto& b) {
                                    return a.score>b.score;
                                });
                            logical_list_scores.reserve(ranked.size());
                            logical_list_keys.reserve(ranked.size());
                            for(const auto& entry:ranked) {
                                std::vector<int> candidate=gibbs.estimate;
                                if(entry.direction>=0) {
                                    const auto& direction=gibbs.logical_directions[
                                        entry.direction];
                                    candidate[direction.variable]=gf.add(
                                        candidate[direction.variable],entry.alpha);
                                    for(const auto& [variable,coefficient]:
                                            direction.pivot_updates)
                                        candidate[variable]=gf.add(candidate[variable],
                                            gf.mul(coefficient,entry.alpha));
                                }
                                if(!decoder.check_z_syndrome(hx_logical,candidate,
                                        logical_syndrome_z))
                                    continue;
                                ++logical_list_syndrome_valid;
                                logical_list_scores.push_back(entry.score);
                                logical_list_keys.push_back(cascadeZCosetKey(candidate,
                                    *cascade_z_stabilizer_space,gf,inner_length,
                                    cascade_inner_dimension,inner_scales));
                            }
                            // Greedy coordinate ascent composes multiple
                            // independent logical moves. Every move is chosen
                            // from the channel posterior only and preserves the
                            // full syndrome by construction.
                            std::vector<int> descent=gibbs.estimate;
                            double descent_score=base_score;
                            for(int pass=0;pass<3;++pass) {
                                int moves_this_pass=0;
                                for(int selected=0;selected<direction_limit;++selected) {
                                    const auto& direction=gibbs.logical_directions[
                                        quality[selected].index];
                                    const int old_free=descent[direction.variable];
                                    int best_value=old_free;
                                    double best_score=descent_score;
                                    for(int value=0;value<q;++value) {
                                        ++logical_descent_evaluated;
                                        const int delta=gf.add(old_free,value);
                                        double score=symbol_log_score(direction.variable,value);
                                        for(const auto& [variable,coefficient]:
                                                direction.pivot_updates) {
                                            const int old_value=descent[variable];
                                            const int new_value=gf.add(old_value,
                                                gf.mul(coefficient,delta));
                                            score+=symbol_log_score(variable,new_value);
                                            score-=symbol_log_score(variable,old_value);
                                        }
                                        score+=descent_score-
                                            symbol_log_score(direction.variable,old_free);
                                        if(score>best_score+1e-12) {
                                            best_score=score;
                                            best_value=value;
                                        }
                                    }
                                    if(best_value!=old_free) {
                                        const int delta=gf.add(old_free,best_value);
                                        descent[direction.variable]=best_value;
                                        for(const auto& [variable,coefficient]:
                                                direction.pivot_updates)
                                            descent[variable]=gf.add(descent[variable],
                                                gf.mul(coefficient,delta));
                                        descent_score=best_score;
                                        ++moves_this_pass;
                                        ++logical_descent_moves;
                                    }
                                }
                                logical_descent_passes=pass+1;
                                if(moves_this_pass==0) break;
                            }
                            if(!decoder.check_z_syndrome(hx_logical,descent,
                                    logical_syndrome_z))
                                throw std::runtime_error(
                                    "Logical coordinate ascent changed the syndrome");
                            logical_descent_word=descent;
                            logical_descent_score=descent_score;
                            ++logical_list_syndrome_valid;
                            logical_list_scores.push_back(descent_score);
                            logical_list_keys.push_back(cascadeZCosetKey(descent,
                                *cascade_z_stabilizer_space,gf,inner_length,
                                cascade_inner_dimension,inner_scales));
                        }
                        if(std::getenv(
                                "GF256_CASCADE_ENABLE_LOGICAL_BEAM_AUDIT") != nullptr) {
                            if(hx_logical.rows()!=cascade_outer_hx_check.rows()*
                                    (inner_length-1))
                                throw std::invalid_argument(
                                    "Logical bit beam currently requires the path inner-check layout");
                            int beam_width=4, max_depth=4;
                            if(const char* value=std::getenv(
                                    "GF256_CASCADE_LOGICAL_BEAM_WIDTH"))
                                beam_width=std::stoi(value);
                            if(const char* value=std::getenv(
                                    "GF256_CASCADE_LOGICAL_BEAM_DEPTH"))
                                max_depth=std::stoi(value);
                            if(const char* value=std::getenv(
                                    "GF256_CASCADE_LOGICAL_BEAM_TURBO_ROUNDS"))
                                beam_turbo_rounds=std::stoi(value);
                            if(const char* value=std::getenv(
                                    "GF256_CASCADE_LOGICAL_BEAM_OUTER_ITERATIONS"))
                                beam_outer_iterations=std::stoi(value);
                            if(const char* value=std::getenv(
                                    "GF256_CASCADE_LOGICAL_BEAM_DAMPING"))
                                beam_damping=std::stod(value);
                            if(beam_width<=0 || beam_width>64 || max_depth<=0 ||
                               max_depth>16 || beam_turbo_rounds<=0 ||
                               beam_outer_iterations<=0 || beam_damping<=0.0 ||
                               beam_damping>1.0)
                                throw std::invalid_argument(
                                    "Logical beam width/depth are outside supported bounds");
                            logical_beam_width_used=beam_width;
                            logical_beam_depth_used=max_depth;
                            struct BeamNode {
                                std::vector<std::pair<int,int>> constraints;
                                std::vector<int> word;
                                std::vector<SoftVec> marginals;
                                double path_log_probability=0.0;
                            };
                            auto run_masked_turbo=[&](
                                    const std::vector<std::pair<int,int>>& constraints,
                                    std::vector<SoftVec>& marginals) {
                                ++logical_beam_bp_calls;
                                return decodeHxCascadeTurbo(decoder,
                                    cascade_outer_hx_check,hx_logical,
                                    logical_syndrome_z,p,
                                    x_prior_valid ? &corrected_x : nullptr,
                                    basis,gf,beam_turbo_rounds,
                                    beam_outer_iterations,beam_damping,false,&marginals,
                                    false,nullptr,false,nullptr,false,false,false,
                                    &constraints,true);
                            };
                            auto constraint_violations=[&](
                                    const std::vector<int>& word,
                                    const std::vector<std::pair<int,int>>& constraints) {
                                int violations=0;
                                for(int v=0;v<symbol_count;++v) {
                                    const auto [mask,value]=constraints[v];
                                    if(!mask) continue;
                                    const int physical=gf.mul(word[v],
                                        gf.inv(inner_scales[v%inner_length]));
                                    violations+=__builtin_popcount(static_cast<unsigned>(
                                        (physical^value)&mask));
                                }
                                logical_beam_constraint_violations+=violations;
                                if(violations==0) ++logical_beam_constraint_consistent;
                                return violations;
                            };
                            std::vector<std::pair<int,int>> root_constraints(
                                symbol_count,{0,0});
                            std::vector<SoftVec> root_marginals;
                            auto root_word=run_masked_turbo(root_constraints,
                                                           root_marginals);
                            constraint_violations(root_word,root_constraints);
                            std::vector<BeamNode> beam;
                            beam.push_back({root_constraints,root_word,
                                            std::move(root_marginals),0.0});
                            std::unordered_set<std::string> unique_words;
                            std::vector<std::vector<int>> valid_words;
                            auto word_key=[](const std::vector<int>& word) {
                                std::string key(word.size()*2,'\0');
                                for(size_t i=0;i<word.size();++i) {
                                    key[2*i]=static_cast<char>(word[i]&0xff);
                                    key[2*i+1]=static_cast<char>((word[i]>>8)&0xff);
                                }
                                return key;
                            };
                            auto record_valid=[&](const std::vector<int>& word) {
                                if(decoder.check_z_syndrome(hx_logical,word,
                                        logical_syndrome_z) &&
                                   unique_words.insert(word_key(word)).second)
                                    valid_words.push_back(word);
                            };
                            record_valid(root_word);
                            for(int depth=0;depth<max_depth && !beam.empty();++depth) {
                                std::vector<BeamNode> next;
                                for(const auto& node:beam) {
                                    std::vector<char> touched(symbol_count,0);
                                    bool any_unsatisfied=false;
                                    for(int r=0;r<hx_logical.rows();++r) {
                                        int value=0;
                                        for(SparseMatrix<int,RowMajor>::InnerIterator it(
                                                hx_logical,r);it;++it)
                                            value=gf.add(value,gf.mul(it.value(),
                                                node.word[it.col()]));
                                        if(value!=logical_syndrome_z[r]) {
                                            any_unsatisfied=true;
                                            for(SparseMatrix<int,RowMajor>::InnerIterator it(
                                                    hx_logical,r);it;++it)
                                                touched[it.col()]=1;
                                        }
                                    }
                                    if(!any_unsatisfied)
                                        std::fill(touched.begin(),touched.end(),1);
                                    double least_reliability=std::numeric_limits<double>::infinity();
                                    int branch_variable=-1,branch_bit=-1;
                                    double branch_probability_one=0.5;
                                    for(int v=0;v<symbol_count;++v) if(touched[v]) {
                                        for(int bit=0;bit<field_degree;++bit) {
                                            const int mask=1<<bit;
                                            if(node.constraints[v].first&mask) continue;
                                            double probability_one=0.0;
                                            for(int z=0;z<gf.getSize();++z)
                                                if((z>>bit)&1)
                                                    probability_one+=node.marginals[v][z];
                                            const double reliability=std::abs(
                                                probability_one-0.5);
                                            if(reliability<least_reliability) {
                                                least_reliability=reliability;
                                                branch_variable=v;
                                                branch_bit=bit;
                                                branch_probability_one=probability_one;
                                            }
                                        }
                                    }
                                    if(branch_variable<0) continue;
                                    for(int bit_value=0;bit_value<=1;++bit_value) {
                                        BeamNode child;
                                        child.constraints=node.constraints;
                                        child.constraints[branch_variable].first|=
                                            1<<branch_bit;
                                        child.constraints[branch_variable].second=
                                            (child.constraints[branch_variable].second&
                                                ~(1<<branch_bit))|
                                            (bit_value<<branch_bit);
                                        const double probability=bit_value
                                            ? branch_probability_one
                                            : 1.0-branch_probability_one;
                                        child.path_log_probability=
                                            node.path_log_probability+
                                            std::log(std::max(1e-300,probability));
                                        child.word=run_masked_turbo(
                                            child.constraints,child.marginals);
                                        if(constraint_violations(child.word,
                                                child.constraints)!=0)
                                            continue;
                                        record_valid(child.word);
                                        next.push_back(std::move(child));
                                    }
                                }
                                std::sort(next.begin(),next.end(),
                                    [](const BeamNode& a,const BeamNode& b) {
                                        return a.path_log_probability>
                                               b.path_log_probability;
                                    });
                                std::unordered_set<std::string> unique_paths;
                                std::vector<BeamNode> pruned;
                                pruned.reserve(std::min<size_t>(beam_width,next.size()));
                                for(auto& node:next) {
                                    std::string key;
                                    key.reserve(node.constraints.size()*2);
                                    for(const auto& [mask,value]:node.constraints) {
                                        key.push_back(static_cast<char>(mask&0xff));
                                        key.push_back(static_cast<char>(value&0xff));
                                    }
                                    if(!unique_paths.insert(std::move(key)).second)
                                        continue;
                                    pruned.push_back(std::move(node));
                                    if(static_cast<int>(pruned.size())>=beam_width)
                                        break;
                                }
                                beam=std::move(pruned);
                                logical_beam_retained_nodes+=
                                    static_cast<int>(beam.size());
                            }
                            for(const auto& node:beam) {
                                ++logical_beam_leaf_repairs;
                                const auto repaired=repairSyndromeBeam(
                                    hx_logical,logical_syndrome_z,node.word,
                                    conditional_z_priors,gf,16,8,
                                    &node.constraints,&inner_scales,inner_length);
                                if(repaired.valid) {
                                    ++logical_beam_leaf_repair_successes;
                                    record_valid(repaired.estimate);
                                }
                            }
                            logical_beam_valid_candidates=
                                static_cast<int>(valid_words.size());
                            for(const auto& word:valid_words) {
                                double score=0.0;
                                for(int v=0;v<symbol_count;++v)
                                    score+=std::log(std::max(1e-300,
                                        conditional_z_priors[v][word[v]]));
                                const auto key=cascadeZCosetKey(word,
                                    *cascade_z_stabilizer_space,gf,inner_length,
                                    cascade_inner_dimension,inner_scales);
                                auto found=logical_beam_class_log_mass.find(key);
                                if(found==logical_beam_class_log_mass.end())
                                    logical_beam_class_log_mass.emplace(key,score);
                                else {
                                    const double high=std::max(found->second,score);
                                    found->second=high+std::log(
                                        std::exp(found->second-high)+
                                        std::exp(score-high));
                                }
                            }
                        }
                        const std::string true_key=cascadeZCosetKey(actual_z,
                            *cascade_z_stabilizer_space,gf,inner_length,
                            cascade_inner_dimension,inner_scales);
                        const int true_class_samples=std::count(
                            gibbs.sampled_class_keys.begin(),
                            gibbs.sampled_class_keys.end(),true_key);
                        if(std::getenv(
                                "GF256_CASCADE_ENABLE_LOGICAL_BEAM_AUDIT") != nullptr) {
                            std::vector<std::pair<double,std::string>> ranked;
                            for(const auto& [key,mass]:logical_beam_class_log_mass)
                                ranked.emplace_back(mass,key);
                            std::sort(ranked.begin(),ranked.end(),
                                [](const auto& a,const auto& b) {
                                    return a.first>b.first;
                                });
                            int true_rank=0;
                            for(size_t i=0;i<ranked.size();++i)
                                if(ranked[i].second==true_key) {
                                    true_rank=static_cast<int>(i)+1;
                                    break;
                                }
                            cerr << "[Gibbs-logical-beam-audit] trial=" << trial
                                 << " beam_width=" << logical_beam_width_used
                                 << " max_depth=" << logical_beam_depth_used
                                 << " turbo_rounds=" << beam_turbo_rounds
                                 << " outer_iterations=" << beam_outer_iterations
                                 << " damping=" << beam_damping
                                 << " constraint_penalty_exponent="
                                 << beam_constraint_penalty
                                 << " masked_turbo_calls=" << logical_beam_bp_calls
                                 << " constraint_consistent_calls="
                                 << logical_beam_constraint_consistent
                                 << " constraint_bit_violations="
                                 << logical_beam_constraint_violations
                                 << " retained_beam_nodes="
                                 << logical_beam_retained_nodes
                                 << " leaf_repair_calls="
                                 << logical_beam_leaf_repairs
                                 << " leaf_repair_successes="
                                 << logical_beam_leaf_repair_successes
                                 << " syndrome_valid_words="
                                 << logical_beam_valid_candidates
                                 << " logical_classes="
                                 << logical_beam_class_log_mass.size()
                                 << " true_class_covered=" << (true_rank>0)
                                 << " true_class_rank=" << true_rank
                                 << " selected_class_correct="
                                 << (!ranked.empty() && ranked.front().second==true_key)
                                 << " best_class_log_mass="
                                 << (ranked.empty() ?
                                     -std::numeric_limits<double>::infinity() :
                                     ranked.front().first) << '\n';
                        }
                        if(logical_list_evaluated>0) {
                            int true_rank=0;
                            for(size_t i=0;i<logical_list_keys.size();++i)
                                if(logical_list_keys[i]==true_key) {
                                    true_rank=static_cast<int>(i)+1;
                                    break;
                                }
                            auto logical_coordinate_weight=[&](
                                    const std::vector<int>& word) {
                                std::vector<int> free_delta(gibbs.free_variables.size());
                                for(size_t f=0;f<gibbs.free_variables.size();++f)
                                    free_delta[f]=gf.add(word[gibbs.free_variables[f]],
                                        gibbs.estimate[gibbs.free_variables[f]]);
                                for(size_t pivot=0;pivot<
                                        gibbs.stabilizer_echelon.size();++pivot) {
                                    if(gibbs.stabilizer_echelon[pivot].empty()) continue;
                                    const int factor=free_delta[pivot];
                                    if(!factor) continue;
                                    for(const auto& [column,value]:
                                            gibbs.stabilizer_echelon[pivot])
                                        free_delta[column]=gf.add(free_delta[column],
                                            gf.mul(factor,value));
                                }
                                int weight=0;
                                for(size_t f=0;f<free_delta.size();++f)
                                    if(gibbs.stabilizer_echelon[f].empty() &&
                                       free_delta[f]!=0) ++weight;
                                return weight;
                            };
                            const int true_logical_weight=
                                logical_coordinate_weight(actual_z);
                            const int descent_logical_weight=
                                logical_coordinate_weight(logical_descent_word);
                            double true_log_posterior=0.0;
                            for(int v=0;v<symbol_count;++v)
                                true_log_posterior+=std::log(std::max(1e-300,
                                    conditional_z_priors[v][actual_z[v]]));
                            double base_log_posterior=0.0;
                            for(int v=0;v<symbol_count;++v)
                                base_log_posterior+=std::log(std::max(1e-300,
                                    conditional_z_priors[v][gibbs.estimate[v]]));
                            cerr << "[Gibbs-logical-candidate-audit] trial=" << trial
                                 << " selected_directions="
                                 << logical_list_direction_count
                                 << " candidates_scored=" << logical_list_evaluated
                                 << " coordinate_ascent_candidates_scored="
                                 << logical_descent_evaluated
                                 << " osd_seed_used=" << logical_search_osd_seed
                                 << " osd_seed_candidates="
                                 << logical_search_osd_seed_candidates
                                 << " coordinate_ascent_passes="
                                 << logical_descent_passes
                                 << " coordinate_ascent_moves="
                                 << logical_descent_moves
                                 << " top_k=" << logical_list_keys.size()
                                 << " syndrome_valid=" << logical_list_syndrome_valid
                                 << " distinct_classes=" << std::set<std::string>(
                                     logical_list_keys.begin(),logical_list_keys.end()).size()
                                 << " true_class_covered=" << (true_rank>0)
                                 << " true_class_rank=" << true_rank
                                 << " true_class_logical_weight_from_base="
                                 << true_logical_weight
                                 << " coordinate_ascent_class_logical_weight="
                                 << descent_logical_weight
                                 << " coordinate_ascent_conditional_log_prior="
                                 << logical_descent_score
                                 << " gibbs_base_conditional_log_prior="
                                 << base_log_posterior
                                 << " truth_conditional_channel_log_prior_audit="
                                 << true_log_posterior
                                 << " best_list_conditional_log_prior="
                                 << (logical_list_scores.empty() ?
                                     -std::numeric_limits<double>::infinity() :
                                     logical_list_scores.front()) << '\n';
                        }
                        cerr << "[Gibbs-coset-audit] trial=" << trial
                             << " syndrome_valid=" << gibbs.syndrome_ok
                             << " samples=" << gibbs.samples
                             << " kernel_dimension=" << gibbs.kernel_dimension
                             << " stabilizer_rank=" << gibbs.stabilizer_rank
                             << " logical_dimension=" << gibbs.logical_dimension
                             << " distinct_classes=" << gibbs.distinct_classes
                             << " winning_class_samples=" << gibbs.winning_class_samples
                             << " updates=" << gibbs.updates
                             << " targeted_logical_updates=" << gibbs.logical_updates
                             << " changed_updates=" << gibbs.changed_updates
                             << " exchange_attempts=" << gibbs.exchange_attempts
                             << " exchange_accepts=" << gibbs.exchange_accepts
                             << " cold_class_persistence=" << gibbs.class_persistence
                             << " cold_score_lag1=" << gibbs.score_lag1_correlation
                             << " replica_change_rates=";
                        for(size_t replica=0;replica<gibbs.replica_change_rates.size();++replica)
                            cerr << (replica ? "," : "")
                                 << gibbs_betas[replica] << ":"
                                 << gibbs.replica_change_rates[replica];
                        cerr
                             << " true_class_sampled=" << (true_class_samples>0)
                             << " true_class_samples=" << true_class_samples
                             << " seed=" << gibbs_seed
                             << " sweeps=" << gibbs_sweeps
                             << " burn_in=" << gibbs_burn
                             << " thin=" << gibbs_thin
                             << " block_updates_per_sweep="
                             << gibbs_block_updates
                             << " logical_updates_per_sweep="
                             << gibbs_logical_updates
                             << " stabilizer_updates_per_sweep="
                             << gibbs_stabilizer_updates
                             << " inverse_temperatures=";
                        for(size_t replica=0;replica<gibbs_betas.size();++replica)
                            cerr << (replica ? "," : "") << gibbs_betas[replica];
                        cerr << '\n';
                        if (std::getenv("GF256_CASCADE_GIBBS_AUDIT_ONLY") != nullptr)
                            return 0;
                    }
                    auto osd_z=localizedOsd(
                        hx_logical,logical_syndrome_z,
                        full_hx_conditioned.hard,full_bp_posterior,gf,64);
                    int osd_passes=1;
                    report_z_candidate("localized_OSD_pass1",osd_z.estimate);
                    for(int pass=2;pass<=6 && !osd_z.valid;++pass) {
                        auto next_osd=localizedOsd(
                            hx_logical,logical_syndrome_z,osd_z.estimate,
                            full_bp_posterior,gf,64);
                        if(next_osd.unsatisfied>=osd_z.unsatisfied) break;
                        osd_z=std::move(next_osd); osd_passes=pass;
                        report_z_candidate("localized_OSD_pass"+std::to_string(pass),
                                           osd_z.estimate);
                    }
                    cerr << "[Z-full-HX-conditioned] trial=" << trial
                         << " max_iterations=200 syndrome=" << full_hx_conditioned.valid
                         << " x_exact=" << (corrected_x==actual_x) << '\n';
                    cerr << "[Z-syndrome-beam] trial=" << trial
                         << " valid=" << repaired_z.valid
                         << " depth=" << repaired_z.depth
                         << " unsatisfied=" << repaired_z.unsatisfied
                         << " prior_cost_delta=" << repaired_z.prior_cost
                         << " beam_width=128 max_depth=12\n";
                    cerr << "[Z-localized-OSD] trial=" << trial
                         << " valid=" << osd_z.valid
                         << " candidates_tested=" << osd_z.depth
                         << " unsatisfied=" << osd_z.unsatisfied
                         << " negative_log_posterior=" << osd_z.prior_cost
                         << " free_symbols=64 order=1 passes=" << osd_passes << '\n';
                    const auto full_osd_z=fullMatrixOsd1(
                        hx_logical,logical_syndrome_z,
                        full_hx_conditioned.hard,full_bp_posterior,gf,256);
                    report_z_candidate("full_matrix_OSD_order1",full_osd_z.estimate);
                    cerr << "[Z-full-matrix-OSD] trial=" << trial
                         << " valid=" << full_osd_z.syndrome_ok
                         << " rank=" << full_osd_z.rank
                         << " candidates_tested=" << full_osd_z.candidates_tested
                         << " log_posterior=" << full_osd_z.log_posterior
                         << " weakest_free_symbols=256 order=1\n";
                    std::vector<std::vector<int>> osd2_candidates;
                    const auto full_osd2_z=fullMatrixOsd1(
                        hx_logical,logical_syndrome_z,
                        full_hx_conditioned.hard,full_bp_posterior,gf,12,2,8,
                        false,&osd2_candidates);
                    report_z_candidate("full_matrix_OSD_order2_local",full_osd2_z.estimate);
                    cerr << "[Z-full-matrix-OSD2] trial=" << trial
                         << " valid=" << full_osd2_z.syndrome_ok
                         << " rank=" << full_osd2_z.rank
                         << " candidates_tested=" << full_osd2_z.candidates_tested
                         << " log_posterior=" << full_osd2_z.log_posterior
                         << " weakest_free_symbols=12 alternatives=8 order=2\n";
                    // Diagnostic only: compare the best individual OSD word
                    // with a posterior-mass ranking over logical cosets. The
                    // true error is read only after candidate generation and
                    // is used solely to measure list coverage, never to choose.
                    struct CosetMass {
                        double log_mass=-std::numeric_limits<double>::infinity();
                        int words=0;
                    };
                    std::map<std::string,CosetMass> coset_masses;
                    int signature_stabilizer_mismatches=0;
                    if (!osd2_candidates.empty()) {
                        const auto reference_key=cascadeZCosetKey(
                            osd2_candidates.front(),*cascade_z_stabilizer_space,
                            gf,inner_length,cascade_inner_dimension,inner_scales);
                        const int checks=std::min<int>(64,osd2_candidates.size());
                        for (int i=0; i<checks; ++i) {
                            std::vector<int> difference(symbol_count);
                            for (int v=0; v<symbol_count; ++v)
                                difference[v]=gf.add(osd2_candidates.front()[v],
                                                     osd2_candidates[i][v]);
                            const bool same_key=reference_key==cascadeZCosetKey(
                                osd2_candidates[i],*cascade_z_stabilizer_space,
                                gf,inner_length,cascade_inner_dimension,inner_scales);
                            const bool stabilizer=isCascadeZStabilizer(
                                difference,*cascade_z_stabilizer_space,gf,
                                inner_length,cascade_inner_dimension,inner_scales);
                            signature_stabilizer_mismatches+=same_key!=stabilizer;
                        }
                    }
                    for (const auto& candidate_word:osd2_candidates) {
                        double word_score=0.0;
                        for (int v=0; v<symbol_count; ++v)
                            word_score+=std::log(std::max(1e-300,
                                full_bp_posterior[v][candidate_word[v]]));
                        auto label=cascadeZCosetKey(candidate_word,
                            *cascade_z_stabilizer_space,gf,inner_length,
                            cascade_inner_dimension,inner_scales);
                        auto& bucket=coset_masses[label];
                        const double high=std::max(bucket.log_mass,word_score);
                        bucket.log_mass=high+std::log(
                            std::exp(bucket.log_mass-high)+std::exp(word_score-high));
                        ++bucket.words;
                    }
                    if (!coset_masses.empty()) {
                        const auto true_label=cascadeZCosetKey(actual_z,
                            *cascade_z_stabilizer_space,gf,inner_length,
                            cascade_inner_dimension,inner_scales);
                        auto selected=coset_masses.begin();
                        for (auto it=std::next(coset_masses.begin());
                             it!=coset_masses.end(); ++it)
                            if (it->second.log_mass>selected->second.log_mass)
                                selected=it;
                        std::vector<std::pair<double,std::string>> ranked_classes;
                        ranked_classes.reserve(coset_masses.size());
                        for (const auto& [label,mass]:coset_masses)
                            ranked_classes.emplace_back(mass.log_mass,label);
                        std::sort(ranked_classes.begin(),ranked_classes.end(),
                            [](const auto& a,const auto& b){return a.first>b.first;});
                        int true_class_rank=0;
                        for (size_t rank=0; rank<ranked_classes.size(); ++rank)
                            if (ranked_classes[rank].second==true_label) {
                                true_class_rank=static_cast<int>(rank)+1;
                                break;
                            }
                        const bool true_class_covered=true_class_rank!=0;
                        const bool class_selected_correctly=selected->first==true_label;
                        cerr << "[OSD2-logical-class-audit] trial=" << trial
                             << " valid_words=" << osd2_candidates.size()
                             << " classes=" << coset_masses.size()
                             << " true_class_covered=" << true_class_covered
                             << " true_class_rank=" << true_class_rank
                             << " mass_selected_class_correct=" << class_selected_correctly
                             << " selected_class_words=" << selected->second.words
                             << " selected_class_log_mass=" << selected->second.log_mass
                             << " best_single_word_log_score="
                             << full_osd2_z.log_posterior
                             << " signature_stabilizer_mismatches="
                             << signature_stabilizer_mismatches << "/"
                             << std::min<size_t>(64,osd2_candidates.size()) << '\n';
                    }
                    osd2_candidates.clear();
                    osd2_candidates.shrink_to_fit();
                    // Change only the deterministic pivot ordering among
                    // exactly tied reliabilities. This keeps the same BP
                    // posterior and syndrome while testing whether one OSD
                    // basis alone caused the logical-coset coverage miss.
                    std::vector<std::vector<int>> tied_basis_candidates;
                    const uint32_t tied_basis_seed=20261102;
                    const auto tied_basis_osd=fullMatrixOsd1(
                        hx_logical,logical_syndrome_z,
                        full_hx_conditioned.hard,full_bp_posterior,gf,12,2,8,
                        false,&tied_basis_candidates,tied_basis_seed);
                    std::set<std::string> tied_basis_classes;
                    for (const auto& candidate_word:tied_basis_candidates)
                        tied_basis_classes.insert(cascadeZCosetKey(
                            candidate_word,*cascade_z_stabilizer_space,gf,
                            inner_length,cascade_inner_dimension,inner_scales));
                    const auto tied_true_label=cascadeZCosetKey(actual_z,
                        *cascade_z_stabilizer_space,gf,inner_length,
                        cascade_inner_dimension,inner_scales);
                    cerr << "[OSD2-tied-basis-audit] trial=" << trial
                         << " seed=" << tied_basis_seed
                         << " syndrome_valid=" << tied_basis_osd.syndrome_ok
                         << " valid_words=" << tied_basis_candidates.size()
                         << " classes=" << tied_basis_classes.size()
                         << " true_class_covered="
                         << tied_basis_classes.count(tied_true_label)
                         << " candidates_tested=" << tied_basis_osd.candidates_tested
                         << " log_posterior=" << tied_basis_osd.log_posterior << '\n';
                    tied_basis_candidates.clear();
                    tied_basis_candidates.shrink_to_fit();
                    std::vector<std::vector<int>> weak_pivot_candidates;
                    const auto weak_pivot_osd=fullMatrixOsd1(
                        hx_logical,logical_syndrome_z,
                        full_hx_conditioned.hard,full_bp_posterior,gf,12,2,8,
                        true,&weak_pivot_candidates);
                    report_z_candidate("full_matrix_OSD2_weak_pivots",
                                       weak_pivot_osd.estimate);
                    std::set<std::string> weak_pivot_classes;
                    for (const auto& candidate_word:weak_pivot_candidates)
                        weak_pivot_classes.insert(cascadeZCosetKey(
                            candidate_word,*cascade_z_stabilizer_space,gf,
                            inner_length,cascade_inner_dimension,inner_scales));
                    cerr << "[OSD2-weak-pivot-audit] trial=" << trial
                         << " syndrome_valid=" << weak_pivot_osd.syndrome_ok
                         << " valid_words=" << weak_pivot_candidates.size()
                         << " classes=" << weak_pivot_classes.size()
                         << " true_class_covered="
                         << weak_pivot_classes.count(tied_true_label)
                         << " candidates_tested=" << weak_pivot_osd.candidates_tested
                         << " log_posterior=" << weak_pivot_osd.log_posterior << '\n';
                    weak_pivot_candidates.clear();
                    weak_pivot_candidates.shrink_to_fit();
                    std::vector<std::vector<int>> weak_pivot16_candidates;
                    const auto weak_pivot16_osd=fullMatrixOsd1(
                        hx_logical,logical_syndrome_z,
                        full_hx_conditioned.hard,full_bp_posterior,gf,16,2,8,
                        true,&weak_pivot16_candidates);
                    report_z_candidate("full_matrix_OSD2_weak_pivots_free16",
                                       weak_pivot16_osd.estimate);
                    std::set<std::string> weak_pivot16_classes;
                    for (const auto& candidate_word:weak_pivot16_candidates)
                        weak_pivot16_classes.insert(cascadeZCosetKey(
                            candidate_word,*cascade_z_stabilizer_space,gf,
                            inner_length,cascade_inner_dimension,inner_scales));
                    cerr << "[OSD2-weak-pivot-free16-audit] trial=" << trial
                         << " syndrome_valid=" << weak_pivot16_osd.syndrome_ok
                         << " valid_words=" << weak_pivot16_candidates.size()
                         << " classes=" << weak_pivot16_classes.size()
                         << " true_class_covered="
                         << weak_pivot16_classes.count(tied_true_label)
                         << " candidates_tested=" << weak_pivot16_osd.candidates_tested
                         << " log_posterior=" << weak_pivot16_osd.log_posterior << '\n';
                    weak_pivot16_candidates.clear();
                    weak_pivot16_candidates.shrink_to_fit();
                    // Expanded list sweeps are resource-intensive and are
                    // opt-in. Enable only for targeted coverage audits.
                    if (std::getenv("GF256_CASCADE_ENABLE_EXPANDED_OSD_AUDIT") != nullptr) {
                    std::vector<std::vector<int>> weak_pivot256_candidates;
                    int true_free_deviations=-1, true_selected_free_deviations=-1;
                    const auto weak_pivot256_osd=fullMatrixOsd1(
                        hx_logical,logical_syndrome_z,
                        full_hx_conditioned.hard,full_bp_posterior,gf,256,1,8,
                        true,&weak_pivot256_candidates,0,0,&actual_z,
                        &true_free_deviations,&true_selected_free_deviations);
                    report_z_candidate("full_matrix_OSD1_weak_pivots_free256",
                                       weak_pivot256_osd.estimate);
                    std::set<std::string> weak_pivot256_classes;
                    for (const auto& candidate_word:weak_pivot256_candidates)
                        weak_pivot256_classes.insert(cascadeZCosetKey(
                            candidate_word,*cascade_z_stabilizer_space,gf,
                            inner_length,cascade_inner_dimension,inner_scales));
                    cerr << "[OSD1-weak-pivot-free256-audit] trial=" << trial
                         << " syndrome_valid=" << weak_pivot256_osd.syndrome_ok
                         << " valid_words=" << weak_pivot256_candidates.size()
                         << " classes=" << weak_pivot256_classes.size()
                         << " true_class_covered="
                         << weak_pivot256_classes.count(tied_true_label)
                         << " candidates_tested=" << weak_pivot256_osd.candidates_tested
                         << " truth_free_deviations=" << true_free_deviations
                         << " truth_selected_free_deviations="
                         << true_selected_free_deviations
                         << " log_posterior=" << weak_pivot256_osd.log_posterior << '\n';
                    weak_pivot256_candidates.clear();
                    weak_pivot256_candidates.shrink_to_fit();
                    std::vector<std::vector<int>> weak_pivot512_candidates;
                    int true_free512_deviations=-1, true_selected_free512_deviations=-1;
                    const auto weak_pivot512_osd=fullMatrixOsd1(
                        hx_logical,logical_syndrome_z,
                        full_hx_conditioned.hard,full_bp_posterior,gf,512,1,8,
                        true,&weak_pivot512_candidates,0,0,&actual_z,
                        &true_free512_deviations,
                        &true_selected_free512_deviations);
                    report_z_candidate("full_matrix_OSD1_weak_pivots_free512",
                                       weak_pivot512_osd.estimate);
                    std::set<std::string> weak_pivot512_classes;
                    for (const auto& candidate_word:weak_pivot512_candidates)
                        weak_pivot512_classes.insert(cascadeZCosetKey(
                            candidate_word,*cascade_z_stabilizer_space,gf,
                            inner_length,cascade_inner_dimension,inner_scales));
                    cerr << "[OSD1-weak-pivot-free512-audit] trial=" << trial
                         << " syndrome_valid=" << weak_pivot512_osd.syndrome_ok
                         << " valid_words=" << weak_pivot512_candidates.size()
                         << " classes=" << weak_pivot512_classes.size()
                         << " true_class_covered="
                         << weak_pivot512_classes.count(tied_true_label)
                         << " candidates_tested=" << weak_pivot512_osd.candidates_tested
                         << " truth_free_deviations=" << true_free512_deviations
                         << " truth_selected_free_deviations="
                         << true_selected_free512_deviations
                         << " log_posterior=" << weak_pivot512_osd.log_posterior << '\n';
                    weak_pivot512_candidates.clear();
                    weak_pivot512_candidates.shrink_to_fit();
                    std::vector<std::vector<int>> weak_pivot32_order2_candidates;
                    const auto weak_pivot32_order2=fullMatrixOsd1(
                        hx_logical,logical_syndrome_z,
                        full_hx_conditioned.hard,full_bp_posterior,gf,32,2,2,
                        true,&weak_pivot32_order2_candidates);
                    report_z_candidate("full_matrix_OSD2_weak_pivots_free32_alt2",
                                       weak_pivot32_order2.estimate);
                    std::set<std::string> weak_pivot32_order2_classes;
                    for (const auto& candidate_word:weak_pivot32_order2_candidates)
                        weak_pivot32_order2_classes.insert(cascadeZCosetKey(
                            candidate_word,*cascade_z_stabilizer_space,gf,
                            inner_length,cascade_inner_dimension,inner_scales));
                    cerr << "[OSD2-weak-pivot-free32-alt2-audit] trial=" << trial
                         << " syndrome_valid=" << weak_pivot32_order2.syndrome_ok
                         << " valid_words=" << weak_pivot32_order2_candidates.size()
                         << " classes=" << weak_pivot32_order2_classes.size()
                         << " true_class_covered="
                         << weak_pivot32_order2_classes.count(tied_true_label)
                         << " candidates_tested=" << weak_pivot32_order2.candidates_tested
                         << " log_posterior=" << weak_pivot32_order2.log_posterior << '\n';
                    weak_pivot32_order2_candidates.clear();
                    weak_pivot32_order2_candidates.shrink_to_fit();
                    std::vector<std::vector<int>> weak_pivot64_order2_candidates;
                    const auto weak_pivot64_order2=fullMatrixOsd1(
                        hx_logical,logical_syndrome_z,
                        full_hx_conditioned.hard,full_bp_posterior,gf,64,2,2,
                        true,&weak_pivot64_order2_candidates);
                    report_z_candidate("full_matrix_OSD2_weak_pivots_free64_alt2",
                                       weak_pivot64_order2.estimate);
                    std::set<std::string> weak_pivot64_order2_classes;
                    for (const auto& candidate_word:weak_pivot64_order2_candidates)
                        weak_pivot64_order2_classes.insert(cascadeZCosetKey(
                            candidate_word,*cascade_z_stabilizer_space,gf,
                            inner_length,cascade_inner_dimension,inner_scales));
                    cerr << "[OSD2-weak-pivot-free64-alt2-audit] trial=" << trial
                         << " syndrome_valid=" << weak_pivot64_order2.syndrome_ok
                         << " valid_words=" << weak_pivot64_order2_candidates.size()
                         << " classes=" << weak_pivot64_order2_classes.size()
                         << " true_class_covered="
                         << weak_pivot64_order2_classes.count(tied_true_label)
                         << " candidates_tested=" << weak_pivot64_order2.candidates_tested
                         << " log_posterior=" << weak_pivot64_order2.log_posterior << '\n';
                    weak_pivot64_order2_candidates.clear();
                    weak_pivot64_order2_candidates.shrink_to_fit();
                    std::vector<std::vector<int>> weak_pivot128_order2_candidates;
                    const auto weak_pivot128_order2=fullMatrixOsd1(
                        hx_logical,logical_syndrome_z,
                        full_hx_conditioned.hard,full_bp_posterior,gf,128,2,1,
                        true,&weak_pivot128_order2_candidates);
                    report_z_candidate("full_matrix_OSD2_weak_pivots_free128_alt1",
                                       weak_pivot128_order2.estimate);
                    std::set<std::string> weak_pivot128_order2_classes;
                    for (const auto& candidate_word:weak_pivot128_order2_candidates)
                        weak_pivot128_order2_classes.insert(cascadeZCosetKey(
                            candidate_word,*cascade_z_stabilizer_space,gf,
                            inner_length,cascade_inner_dimension,inner_scales));
                    cerr << "[OSD2-weak-pivot-free128-alt1-audit] trial=" << trial
                         << " syndrome_valid=" << weak_pivot128_order2.syndrome_ok
                         << " valid_words=" << weak_pivot128_order2_candidates.size()
                         << " classes=" << weak_pivot128_order2_classes.size()
                         << " true_class_covered="
                         << weak_pivot128_order2_classes.count(tied_true_label)
                         << " candidates_tested=" << weak_pivot128_order2.candidates_tested
                         << " log_posterior=" << weak_pivot128_order2.log_posterior << '\n';
                    weak_pivot128_order2_candidates.clear();
                    weak_pivot128_order2_candidates.shrink_to_fit();
                    }
                    std::set<std::string> union_perturbed_classes;
                    for (const auto& [seed,window]:
                         std::vector<std::pair<uint32_t,int>>{
                             {20261103,8},{20261104,64},
                             {20261105,256},{20261106,1024}}) {
                        std::vector<std::vector<int>> perturbed_candidates;
                        const auto perturbed_osd=fullMatrixOsd1(
                            hx_logical,logical_syndrome_z,
                            full_hx_conditioned.hard,full_bp_posterior,gf,
                            12,2,8,false,&perturbed_candidates,seed,window);
                        std::set<std::string> local_classes;
                        for (const auto& candidate_word:perturbed_candidates) {
                            auto key=cascadeZCosetKey(candidate_word,
                                *cascade_z_stabilizer_space,gf,inner_length,
                                cascade_inner_dimension,inner_scales);
                            local_classes.insert(key);
                            union_perturbed_classes.insert(std::move(key));
                        }
                        cerr << "[OSD2-perturbed-basis-audit] trial=" << trial
                             << " seed=" << seed << " window=" << window
                             << " syndrome_valid=" << perturbed_osd.syndrome_ok
                             << " valid_words=" << perturbed_candidates.size()
                             << " classes=" << local_classes.size()
                             << " true_class_covered="
                             << local_classes.count(tied_true_label)
                             << " candidates_tested=" << perturbed_osd.candidates_tested
                             << " log_posterior=" << perturbed_osd.log_posterior
                             << '\n';
                        perturbed_candidates.clear();
                        perturbed_candidates.shrink_to_fit();
                    }
                    cerr << "[OSD2-perturbed-basis-union] trial=" << trial
                         << " windows=4 unique_classes="
                         << union_perturbed_classes.size()
                         << " true_class_covered="
                         << union_perturbed_classes.count(tied_true_label) << '\n';
                    int residual_bits=0;
                    vector<int> affected_blocks;
                    for (int j=0; j<symbol_count; ++j) {
                        const int delta=gf.add(actual_z[j],corrected_z[j]);
                        residual_bits+=__builtin_popcount(static_cast<unsigned>(delta));
                        if (delta && (affected_blocks.empty() || affected_blocks.back()!=j/inner_length))
                            affected_blocks.push_back(j/inner_length);
                    }
                    cerr << "[MAP-diagnostic] trial=" << trial
                         << " candidate_minus_truth_loglik="
                         << channel_score(corrected_z)-channel_score(actual_z)
                         << " residual_bits=" << residual_bits
                         << " affected_blocks=" << affected_blocks.size() << '\n';
                    if (std::getenv("GF256_CASCADE_SKIP_TARGET_SCHEDULES")==nullptr)
                    for (const auto& [rounds, iterations, blend] :
                         std::vector<std::tuple<int,int,double>>{{3,10,0.3},
                             {5,10,0.3},{10,10,0.3},{20,10,0.3},
                             {10,20,0.1},{10,20,0.2},{10,20,0.4},
                             {10,20,0.5},{10,20,0.6},{10,20,0.8},
                             {20,20,0.2},{20,20,0.5},{20,20,0.8}}) {
                        const auto candidate=decodeHxCascadeTurbo(decoder,
                            cascade_outer_hx_check,hx_logical,logical_syndrome_z,
                            p,&corrected_x,basis,gf,rounds,iterations,blend,false);
                        vector<int> residual(symbol_count);
                        for (int j=0;j<symbol_count;++j)
                            residual[j]=gf.add(actual_z[j],candidate[j]);
                        cerr << "[schedule] trial=" << trial
                             << " rounds=" << rounds << " outer=" << iterations
                             << " blend=" << blend
                             << " syndrome=" << decoder.check_z_syndrome(
                                 hx_logical,candidate,logical_syndrome_z)
                             << " logical=" << !isCascadeZStabilizer(
                                 residual,*cascade_z_stabilizer_space,gf,inner_length,
                                 cascade_inner_dimension,inner_scales)
                             << " delta_loglik="
                             << channel_score(candidate)-channel_score(corrected_z)
                             << '\n';
                    }
                }
                if (!alternate_z.empty()) {
                    vector<int> alternate_residual(symbol_count);
                    int differences=0;
                    for (int j=0; j<symbol_count; ++j) {
                        alternate_residual[j]=gf.add(actual_z[j],alternate_z[j]);
                        differences+=alternate_z[j]!=corrected_z[j];
                    }
                    cerr << "[alternate] trial=" << trial
                         << " syndrome=" << decoder.check_z_syndrome(
                                hx_logical,alternate_z,logical_syndrome_z)
                         << " logical=" << !isCascadeZStabilizer(
                                alternate_residual,*cascade_z_stabilizer_space,gf,inner_length,
                                cascade_inner_dimension,inner_scales)
                         << " candidate_differences=" << differences << '\n';
                }
                exact_x_frames += corrected_x == actual_x;
                vector<int> residual_x(symbol_count), residual_z(symbol_count);
                bool exact_recovery = true;
                for (int symbol = 0; symbol < symbol_count; ++symbol) {
                    residual_x[symbol] = gf.add(actual_x[symbol], corrected_x[symbol]);
                    residual_z[symbol] = gf.add(actual_z[symbol], corrected_z[symbol]);
                    residual_symbols += residual_x[symbol] != 0 || residual_z[symbol] != 0;
                    exact_recovery &= residual_x[symbol] == 0 && residual_z[symbol] == 0;
                }
                exact_recovery_frames += exact_recovery;

                bool logical_failure = !syndrome_ok;
                if (syndrome_ok && !exact_recovery) {
                    if (cascade_x_stabilizer_space) {
                        logical_failure = !isCascadeXStabilizer(
                            residual_x,*cascade_x_stabilizer_space,gf,inner_length,
                            cascade_inner_dimension,inner_scales)
                            || !isCascadeZStabilizer(
                                residual_z,*cascade_z_stabilizer_space,gf,inner_length,
                                cascade_inner_dimension,inner_scales);
                    } else {
                        logical_failure = !isInStabilizerRowSpace(hx_logical, residual_x, gf)
                            || !isInStabilizerRowSpace(hz_logical, residual_z, gf);
                    }
                }
                if (target_trial >= 0) {
                    const char* export_prefix = std::getenv("GF256_FRAME_EXPORT_PREFIX");
                    if (export_prefix != nullptr) {
                        const std::string prefix(export_prefix);
                        auto write_named_vector = [](std::ostream& out,
                                                     const std::string& name,
                                                     const std::vector<int>& values) {
                            out << name;
                            for (int value : values) out << '\t' << value;
                            out << '\n';
                        };
                        std::vector<int> x_hat_bits(physical_bits, 0);
                        std::vector<int> z_hat_bits(physical_bits, 0);
                        std::vector<int> x_truth_bits(physical_bits, 0);
                        std::vector<int> z_truth_bits(physical_bits, 0);
                        for (int symbol = 0; symbol < symbol_count; ++symbol) {
                            const int x_hat = basis.toDual(corrected_x[symbol]);
                            const int x_truth = basis.toDual(actual_x[symbol]);
                            for (int bit = 0; bit < field_degree; ++bit) {
                                const int index = symbol * field_degree + bit;
                                x_hat_bits[index] = (x_hat >> bit) & 1;
                                z_hat_bits[index] = (corrected_z[symbol] >> bit) & 1;
                                x_truth_bits[index] = (x_truth >> bit) & 1;
                                z_truth_bits[index] = (actual_z[symbol] >> bit) & 1;
                            }
                        }
                        std::ofstream decoder_input(prefix + ".decoder_input.tsv");
                        if (!decoder_input)
                            throw std::runtime_error("Could not write decoder input export");
                        decoder_input << "p_D\t" << p << "\nz_bias\t" << z_bias
                                      << "\nfield_degree\t" << field_degree
                                      << "\ninner_length\t" << inner_length << '\n';
                        write_named_vector(decoder_input, "syndrome_z", syndrome_z);
                        write_named_vector(decoder_input, "x_hat_bits", x_hat_bits);
                        write_named_vector(decoder_input, "current_z_bits", z_hat_bits);

                        std::ofstream truth(prefix + ".truth.tsv");
                        if (!truth)
                            throw std::runtime_error("Could not write frame truth export");
                        write_named_vector(truth, "error_x_bits", error_x);
                        write_named_vector(truth, "error_z_bits", error_z);

                        std::ofstream baseline(prefix + ".baseline.tsv");
                        if (!baseline)
                            throw std::runtime_error("Could not write baseline export");
                        baseline << "logical_failure\t" << logical_failure
                                 << "\nx_syndrome_ok\t" << x_syndrome_ok
                                 << "\nz_syndrome_ok\t" << z_syndrome_ok
                                 << "\nexact_x\t" << (corrected_x == actual_x)
                                 << "\nexact_z\t" << (corrected_z == actual_z) << '\n';
                    }
                }
                if (std::getenv("GF256_FRAME_RESULTS") != nullptr) {
                    cerr << "[frame-result] trial=" << trial
                         << " logical_failure=" << logical_failure
                         << " x_syndrome_ok=" << x_syndrome_ok
                         << " z_syndrome_ok=" << z_syndrome_ok
                         << " exact_x=" << (corrected_x == actual_x)
                         << " physical_errors=" << frame_physical_errors << '\n';
                }
                if (projection_risk_audit) {
                    cerr << "[projection-risk-audit] trial=" << trial
                         << " risk_hard_disagreements="
                         << initial_projection_diagnostics.hard_disagreements
                         << " bp_valid_projections="
                         << initial_projection_diagnostics.bp_valid
                         << " raw_valid_projections="
                         << initial_projection_diagnostics.raw_valid
                         << " final_valid_projections="
                         << initial_projection_diagnostics.final_valid
                         << " syndrome_ok=" << syndrome_ok
                         << " logical_failure=" << logical_failure << '\n';
                }
                if (logical_failure) {
                    ++failures;
                    if (syndrome_ok) ++undetected_logical_failures;
                    if (conditioned_mode)
                        cerr << "[failure] trial=" << trial
                             << " Xstage=" << x_staged_valid
                             << " Xexact=" << (corrected_x == actual_x)
                             << " Xsyn=" << x_syndrome_ok
                             << " Zsyn=" << z_syndrome_ok
                             << " residual_symbols=" << std::count_if(
                                 residual_z.begin(), residual_z.end(),
                                 [](int value){return value!=0;}) << '\n';
                }
                if (target_trial < 0 && (trial + 1) % 50 == 0)
                    cerr << "[progress] frames=" << trial + 1
                         << " failures=" << failures
                         << " syndrome_failures=" << syndrome_failures << '\n';
            }
            cout << p << '\t' << static_cast<double>(failures) / evaluated_frames << '\t'
                 << evaluated_frames << '\t' << failures << '\t' << syndrome_failures << '\t'
                 << x_syndrome_failures << '\t' << z_syndrome_failures << '\t'
                 << static_cast<double>(physical_errors) / evaluated_frames << '\t'
                 << static_cast<double>(residual_symbols) / evaluated_frames << '\t'
                 << undetected_logical_failures << '\t' << exact_recovery_frames << '\t'
                 << exact_x_frames << '\t'
                 << cascade_trailing_recoveries << '\t' << cascade_staged_attempts << '\t'
                 << cascade_staged_successes << '\t'
                 << x_fallback_attempts << '\t' << x_fallback_successes << '\t'
                 << z_soft_fallback_attempts << '\t'
                 << z_soft_fallback_successes << '\t'
                 << z_osd_fallback_attempts << '\t'
                 << z_osd_fallback_successes << '\t'
                 << z_tuned_fallback_attempts << '\t'
                 << z_tuned_fallback_successes << '\t'
                 << robust_attempts << '\t' << robust_selections << '\t'
                 << x_full_hz_bp_attempts << '\t'
                 << x_full_hz_bp_successes << '\t'
                 << x_soft_posterior_attempts << '\t'
                 << x_soft_posterior_successes << '\t'
                 << z_augmented_attempts << '\t'
                 << z_augmented_selections << '\n';
        }
    } catch (const exception& error) {
        cerr << "Simulation error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}

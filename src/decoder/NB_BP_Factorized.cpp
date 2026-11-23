#include "NB_BP_Decoder.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {
using Matrix = Eigen::SparseMatrix<int, Eigen::RowMajor>;

struct Edge {
    int variable;
    int coefficient;
    std::vector<double> v2c, c2v;
};

struct Graph {
    const Matrix& check;
    std::vector<std::vector<Edge>> rows;
    std::vector<std::vector<std::pair<int, int>>> variable_edges;

    Graph(const Matrix& matrix, const std::vector<double>& initial)
        : check(matrix), rows(matrix.rows()), variable_edges(matrix.cols()) {
        for (int row = 0; row < matrix.rows(); ++row) {
            for (Matrix::InnerIterator entry(matrix, row); entry; ++entry) {
                const int index = static_cast<int>(rows[row].size());
                rows[row].push_back({static_cast<int>(entry.col()), entry.value(), initial,
                                     std::vector<double>(initial.size(), 1.0 / initial.size())});
                variable_edges[entry.col()].emplace_back(row, index);
            }
        }
    }
};

void normalize(std::vector<double>& values) {
    constexpr double floor = 1e-300;
    double total = 0;
    for (double& value : values) {
        value = std::max(floor, value);
        total += value;
    }
    if (!std::isfinite(total) || total <= 0)
        throw std::runtime_error("Invalid factorized BP message");
    for (double& value : values) value /= total;
}

void fwht(std::vector<double>& values, bool inverse) {
    const int q = static_cast<int>(values.size());
    for (int length = 1; length < q; length <<= 1) {
        for (int i = 0; i < q; i += 2 * length) {
            for (int offset = 0; offset < length; ++offset) {
                const double a = values[i + offset], b = values[i + length + offset];
                values[i + offset] = a + b;
                values[i + length + offset] = a - b;
            }
        }
    }
    if (inverse)
        for (double& value : values) value /= q;
}

void check_step(Graph& graph, const std::vector<int>& syndromes,
                const GF2p& gf, int q) {
    for (int row = 0; row < graph.check.rows(); ++row) {
        auto& edges = graph.rows[row];
        const size_t degree = edges.size();
        if (!degree) continue;
        std::vector<std::vector<double>> transformed(degree, std::vector<double>(q));
        for (size_t index = 0; index < degree; ++index) {
            for (int value = 0; value < q; ++value)
                transformed[index][gf.mul(edges[index].coefficient, value)] +=
                    edges[index].v2c[value];
            fwht(transformed[index], false);
        }
        std::vector<std::vector<double>> prefix(degree + 1, std::vector<double>(q, 1.0));
        std::vector<std::vector<double>> suffix(degree + 1, std::vector<double>(q, 1.0));
        for (size_t index = 0; index < degree; ++index)
            for (int value = 0; value < q; ++value)
                prefix[index + 1][value] = prefix[index][value] * transformed[index][value];
        for (size_t index = degree; index > 0; --index)
            for (int value = 0; value < q; ++value)
                suffix[index - 1][value] = suffix[index][value] * transformed[index - 1][value];
        for (size_t index = 0; index < degree; ++index) {
            std::vector<double> other(q);
            for (int value = 0; value < q; ++value)
                other[value] = prefix[index][value] * suffix[index + 1][value];
            fwht(other, true);
            auto& edge = edges[index];
            for (int value = 0; value < q; ++value) {
                const int target = gf.sub(syndromes[row], gf.mul(edge.coefficient, value));
                edge.c2v[value] = std::max(0.0, other[target]);
            }
            normalize(edge.c2v);
        }
    }
}

// The depolarizing channel prior is a tensor product of this 2-by-2 kernel,
// with rows indexed by physical X bits and columns by physical Z bits.
std::vector<double> channel_transform(std::vector<double> values,
                                      int degree, double p_depolar) {
    const double a = 1.0 - p_depolar, b = p_depolar / 3.0;
    const int q = static_cast<int>(values.size());
    for (int bit = 0; bit < degree; ++bit) {
        const int stride = 1 << bit;
        for (int start = 0; start < q; start += 2 * stride) {
            for (int offset = 0; offset < stride; ++offset) {
                const int zero = start + offset, one = zero + stride;
                const double v0 = values[zero], v1 = values[one];
                values[zero] = a * v0 + b * v1;
                values[one] = b * v0 + b * v1;
            }
        }
    }
    return values;
}

bool valid_syndrome(const Graph& graph, const std::vector<int>& estimate,
                    const std::vector<int>& expected, const GF2p& gf) {
    for (int row = 0; row < graph.check.rows(); ++row) {
        int total = 0;
        for (const Edge& edge : graph.rows[row])
            total = gf.add(total, gf.mul(edge.coefficient, estimate[edge.variable]));
        if (total != expected[row]) return false;
    }
    return true;
}
} // namespace

std::pair<std::vector<int>, std::vector<int>> NB_BP_Decoder::joint_decode_factorized(
    const Matrix& HX_nb, const Matrix& HZ_nb,
    const std::vector<int>& logical_syn_x, const std::vector<int>& logical_syn_z,
    double p_depolar) {
    if (HX_nb.cols() != HZ_nb.cols() ||
        logical_syn_x.size() != static_cast<size_t>(HZ_nb.rows()) ||
        logical_syn_z.size() != static_cast<size_t>(HX_nb.rows()) ||
        p_depolar < 0 || p_depolar > 1)
        throw std::invalid_argument("Invalid factorized joint BP dimensions or channel rate");
    const double fm = 2.0 * p_depolar / 3.0;
    std::vector<double> initial_x(q), initial_z(q);
    for (int symbol = 0; symbol < q; ++symbol) {
        const int x_weight = __builtin_popcount(static_cast<unsigned>(basis.toDual(symbol)));
        const int z_weight = __builtin_popcount(static_cast<unsigned>(symbol));
        initial_x[symbol] = std::pow(fm, x_weight) * std::pow(1.0 - fm, e - x_weight);
        initial_z[symbol] = std::pow(fm, z_weight) * std::pow(1.0 - fm, e - z_weight);
    }
    Graph x_graph(HZ_nb, initial_x), z_graph(HX_nb, initial_z);
    std::vector<int> x_dual(q);
    for (int symbol = 0; symbol < q; ++symbol) x_dual[symbol] = basis.toDual(symbol);
    std::vector<int> estimate_x(HX_nb.cols(), 0), estimate_z(HX_nb.cols(), 0);
    for (int iteration = 0; iteration < max_iters; ++iteration) {
        check_step(x_graph, logical_syn_x, gf, q);
        check_step(z_graph, logical_syn_z, gf, q);
        for (int variable = 0; variable < HX_nb.cols(); ++variable) {
            std::vector<double> x_likelihood(q, 1.0), z_likelihood(q, 1.0);
            const auto& x_edges = x_graph.variable_edges[variable];
            const auto& z_edges = z_graph.variable_edges[variable];
            for (const auto& [row, index] : x_edges) {
                const auto& message = x_graph.rows[row][index].c2v;
                for (int value = 0; value < q; ++value)
                    x_likelihood[value] *= message[value];
                normalize(x_likelihood);
            }
            for (const auto& [row, index] : z_edges) {
                const auto& message = z_graph.rows[row][index].c2v;
                for (int value = 0; value < q; ++value)
                    z_likelihood[value] *= message[value];
                normalize(z_likelihood);
            }
            std::vector<double> z_to_x = channel_transform(z_likelihood, e, p_depolar);
            std::vector<double> x_dual_likelihood(q);
            for (int value = 0; value < q; ++value)
                x_dual_likelihood[x_dual[value]] = x_likelihood[value];
            std::vector<double> x_to_z = channel_transform(x_dual_likelihood, e, p_depolar);

            // The channel transforms are the exact marginal sums over the
            // opposite component. Their products with check likelihoods give
            // the two BP beliefs in O(q), without a q-by-q joint MAP scan.
            double best_x = -1.0, best_z = -1.0;
            for (int x = 0; x < q; ++x) {
                const double belief = x_likelihood[x] * z_to_x[x_dual[x]];
                if (belief > best_x) {
                    best_x = belief;
                    estimate_x[variable] = x;
                }
            }
            for (int z = 0; z < q; ++z) {
                const double belief = z_likelihood[z] * x_to_z[z];
                if (belief > best_z) {
                    best_z = belief;
                    estimate_z[variable] = z;
                }
            }
            for (const auto& [target_row, target_index] : x_edges) {
                auto& outgoing = x_graph.rows[target_row][target_index].v2c;
                std::fill(outgoing.begin(), outgoing.end(), 1.0);
                for (const auto& [row, index] : x_edges) {
                    if (row == target_row && index == target_index) continue;
                    const auto& message = x_graph.rows[row][index].c2v;
                    for (int value = 0; value < q; ++value)
                        outgoing[value] *= message[value];
                    normalize(outgoing);
                }
                for (int value = 0; value < q; ++value)
                    outgoing[value] *= z_to_x[x_dual[value]];
                normalize(outgoing);
            }
            for (const auto& [target_row, target_index] : z_edges) {
                auto& outgoing = z_graph.rows[target_row][target_index].v2c;
                std::fill(outgoing.begin(), outgoing.end(), 1.0);
                for (const auto& [row, index] : z_edges) {
                    if (row == target_row && index == target_index) continue;
                    const auto& message = z_graph.rows[row][index].c2v;
                    for (int value = 0; value < q; ++value)
                        outgoing[value] *= message[value];
                    normalize(outgoing);
                }
                for (int value = 0; value < q; ++value)
                    outgoing[value] *= x_to_z[value];
                normalize(outgoing);
            }
        }
        if (valid_syndrome(x_graph, estimate_x, logical_syn_x, gf) &&
            valid_syndrome(z_graph, estimate_z, logical_syn_z, gf)) break;
    }
    return {estimate_x, estimate_z};
}

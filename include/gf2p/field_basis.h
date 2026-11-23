#pragma once

#include <stdexcept>
#include <vector>

#include "gf2p/gf2p.h"

// The binary X check blocks are transposed multiplication matrices. Physical
// X bits therefore use the trace-dual coordinates of a GF(q) element, while
// physical Z bits use its polynomial-basis coordinates.
class FieldBasis {
public:
    explicit FieldBasis(const GF2p& gf) : to_dual_(gf.getSize()), from_dual_(gf.getSize(), -1) {
        const int degree = gf.getPower();
        if (gf.getSize() != (1 << degree)) {
            throw std::invalid_argument("Field basis requires a binary extension field");
        }
        for (int value = 0; value < gf.getSize(); ++value) {
            int dual = 0;
            for (int bit = 0; bit < degree; ++bit) {
                const int basis = 1 << bit;
                int power = gf.mul(value, basis);
                int trace = 0;
                for (int k = 0; k < degree; ++k) {
                    trace = gf.add(trace, power);
                    power = gf.mul(power, power);
                }
                if (trace != 0 && trace != 1) {
                    throw std::runtime_error("Invalid field trace");
                }
                dual |= trace << bit;
            }
            if (from_dual_[dual] != -1) {
                throw std::runtime_error("Trace-dual coordinate map is not invertible");
            }
            to_dual_[value] = dual;
            from_dual_[dual] = value;
        }
    }

    int toDual(int polynomial_value) const { return to_dual_.at(polynomial_value); }
    int fromDual(int dual_bits) const { return from_dual_.at(dual_bits); }

private:
    std::vector<int> to_dual_;
    std::vector<int> from_dual_;
};

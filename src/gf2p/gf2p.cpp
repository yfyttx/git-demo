#include "gf2p/gf2p.h"

GF2p::GF2p(int p_val, const std::vector<int>& irreducible_poly_vec) : p(p_val) {
    if (p <= 0 || p > 16) throw std::invalid_argument("p must be between 1 and 16");
    size = 1 << p;
    mod = size - 1;
    irreducible_poly_val = 0;
    for (int bit : irreducible_poly_vec) irreducible_poly_val = (irreducible_poly_val << 1) | bit;
    alpha = findPrimitiveRoot();
    alpha_to_poly.resize(mod);
    poly_to_log.resize(size, -1);
    int current_poly = 1;
    for (int i = 0; i < mod; ++i) {
        alpha_to_poly[i] = current_poly;
        poly_to_log[current_poly] = i;
        current_poly = polyMul(current_poly, alpha);
    }
}
int GF2p::polyAdd(int a, int b) const { return a ^ b; }
int GF2p::polyMul(int a, int b) const {
    int res = 0;
    while (b > 0) {
        if (b & 1) res ^= a;
        a <<= 1;
        if (a & (1 << p)) a ^= irreducible_poly_val;
        b >>= 1;
    }
    return res;
}
bool GF2p::isPrimitiveRoot(int a) const {
    if (a == 0 || a == 1) return false;
    int current = 1;
    for (int i = 0; i < mod - 1; ++i) {
        current = polyMul(current, a);
        if (current == 1) return false;
    }
    return polyMul(current, a) == 1;
}
int GF2p::findPrimitiveRoot() {
    for (int a = 2; a < size; ++a) if (isPrimitiveRoot(a)) return a;
    throw std::runtime_error("No primitive root found");
}
int GF2p::add(int a, int b) const { return polyAdd(a, b); }
int GF2p::mul(int a, int b) const {
    if (a == 0 || b == 0) return 0;
    return alpha_to_poly.at((poly_to_log.at(a) + poly_to_log.at(b)) % mod);
}
int GF2p::inv(int a) const {
    if (a == 0) throw std::runtime_error("No inverse for 0");
    return alpha_to_poly.at((mod - poly_to_log.at(a)) % mod);
}
int GF2p::log(int a) const {
    if (a == 0) throw std::runtime_error("Log of 0 is undefined");
    return poly_to_log.at(a);
}
int GF2p::exp(int x) const { return alpha_to_poly.at((x % mod + mod) % mod); }
int GF2p::getSize() const { return size; }
int GF2p::getMod() const { return mod; }
int GF2p::getAlpha() const { return alpha; }
std::string GF2p::toString(int a) const {
    std::string s;
    for (int i = p - 1; i >= 0; --i) s += ((a >> i) & 1) ? '1' : '0';
    return s;
}
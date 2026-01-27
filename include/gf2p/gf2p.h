#ifndef GF2P_H
#define GF2P_H

#include <stdexcept>
#include <string>
#include <vector>

class GF2p {
  public:
    // Constructor
    GF2p(int p, const std::vector<int> &irreducible_poly_vec);

    // Public API
    int add(int a, int b) const;
    int mul(int a, int b) const;
    int inv(int a) const;
    int log(int a) const;
    int exp(int x) const;
    int getSize() const;
    int getMod() const;
    int getAlpha() const;
    std::string toString(int a) const;
    int getPower() const {
        return p;
    }

  private:
    // Member variables
    int p;
    int size;
    int mod;
    int irreducible_poly_val;       // For polynomial multiplication
    std::vector<int> alpha_to_poly; // Anti-log table: index -> polynomial form
    std::vector<int> poly_to_log;   // Log table: polynomial form -> index
    int alpha;                      // The primitive root (generator)

    // Private helper functions
    int polyAdd(int a, int b) const;
    int polyMul(int a, int b) const;
    bool isPrimitiveRoot(int a) const;
    int findPrimitiveRoot(); // 确认这里没有 const
};

#endif // GF2P_H

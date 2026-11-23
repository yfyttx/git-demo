"""Compare a two-block cascade with exact physical-error MAP enumeration."""

from itertools import product


def probability(error, flip_probabilities):
    value = 1.0
    for bit, p in zip(error, flip_probabilities):
        value *= p if bit else (1.0 - p)
    return value


def syndrome(error):
    # H1=[1,1]; G2 rows=[1 at col 0, 1 at col a+1].
    return tuple(error[a + 1] ^ error[0] ^ error[a + 7] ^ error[6]
                 for a in range(5))


flip_probabilities = (.031, .17, .043, .11, .07, .23,
                      .15, .035, .29, .12, .051, .19)
direct_mass = {}
direct_max = {}
for error in product(range(2), repeat=12):
    syn = syndrome(error)
    likelihood = probability(error, flip_probabilities)
    direct_mass[syn] = direct_mass.get(syn, 0.0) + likelihood
    direct_max[syn] = max(direct_max.get(syn, 0.0), likelihood)

local_mass = []
local_max = []
for block in range(2):
    masses = {}
    maxima = {}
    p = flip_probabilities[6*block:6*block+6]
    for projected in product(range(2), repeat=5):
        values = []
        for tail in range(2):
            error = (tail,) + tuple(tail ^ u for u in projected)
            values.append(probability(error, p))
        masses[projected] = sum(values)
        maxima[projected] = max(values)
    local_mass.append(masses)
    local_max.append(maxima)

largest_mass_difference = 0.0
largest_map_difference = 0.0
for syn in product(range(2), repeat=5):
    mass = 0.0
    best = 0.0
    for left in product(range(2), repeat=5):
        right = tuple(a ^ b for a, b in zip(left, syn))
        mass += local_mass[0][left] * local_mass[1][right]
        best = max(best, local_max[0][left] * local_max[1][right])
    largest_mass_difference = max(largest_mass_difference,
                                  abs(mass - direct_mass[syn]))
    largest_map_difference = max(largest_map_difference,
                                 abs(best - direct_max[syn]))

assert largest_mass_difference < 1e-14, largest_mass_difference
assert largest_map_difference < 1e-14, largest_map_difference
print("exact MAP audit: 32/32 syndromes agree")
print(f"maximum syndrome probability difference: {largest_mass_difference:.3g}")
print(f"maximum MAP likelihood difference: {largest_map_difference:.3g}")

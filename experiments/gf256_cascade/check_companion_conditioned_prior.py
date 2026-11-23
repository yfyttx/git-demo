#!/usr/bin/env python3
"""Exhaustive local verification of the nine-bit conditional G2 prior."""
import json
import numpy as np
from run_companion_staged_xz import conditioned_z_prior
from run_companion_z_bp import OuterBP

pd = .07125
patterns = (np.arange(512)[:,None] >> np.arange(9)) & 1
actual = conditioned_z_prior(patterns, pd)
u = (patterns[:,:8] ^ patterns[:,8,None]) @ (1 << np.arange(8))
p = np.where(patterns == 1, .5, pd/(3-2*pd))
expected = np.zeros_like(actual)
for j in range(512):
    weights = np.prod(np.where(patterns, p[j], 1-p[j]), axis=1)
    np.add.at(expected[j], u, weights)
error = float(np.max(np.abs(actual-expected)))
assert error < 1e-14
bp = OuterBP([[(0,1),(1,2)]], 2)
prior = np.arange(1,257,dtype=float)
prior /= prior.sum()
a = bp.decode(np.array([3]), prior, 5)
b = bp.decode(np.array([3]), np.tile(prior,(2,1)), 5)
assert np.array_equal(a[0],b[0]) and a[1:] == b[1:]
print(json.dumps(dict(p_D=pd, x_patterns=512, z_patterns_per_x=512,
                     pushforward_max_absolute_error=error, shared_vs_tiled_prior_bp_equal=True),indent=2))

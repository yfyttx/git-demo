"""Layered GF(256) sum-product BP on row groups with disjoint variables."""

import numpy as np

from run_companion_z_bp import fwht


def row_layers(columns):
    """Greedily place checks into parallel groups with no shared variable."""
    groups = []
    occupied = []
    for row, variables in enumerate(columns):
        used = set(map(int, variables))
        if len(used) != len(variables):
            raise ValueError('Duplicate variable in one check')
        for group, seen in zip(groups, occupied):
            if not used & seen:
                group.append(row)
                seen.update(used)
                break
        else:
            groups.append([row])
            occupied.append(used)
    return [np.asarray(group) for group in groups]


def normalize_log(values):
    weights = np.exp(values - values.max(axis=-1, keepdims=True))
    return weights / weights.sum(axis=-1, keepdims=True)


class LayeredBP:
    def __init__(self, bp):
        self.bp = bp
        self.layers = [(rows, bp.cols[rows], bp.invperm[rows], bp.products[rows])
                       for rows in row_layers(bp.cols)]

    def decode(self, syndrome, prior, max_iters, *, damping=0.15, reverse=False):
        """Return a syndrome-valid estimate or the final marginal decision."""
        if prior.shape != (self.bp.n, 256) or max_iters < 1 or not 0 <= damping < 1:
            raise ValueError('Invalid layered BP parameters')
        post = np.log(np.maximum(prior, 1e-300))
        log_messages = np.full((*self.bp.cols.shape, 256), -np.log(256.0))
        layers = self.layers[::-1] if reverse else self.layers
        lookup = [products ^ syndrome[rows, None, None] for rows, _, _, products in layers]
        for iteration in range(1, max_iters+1):
            for (rows, cols, invperm, _), indices in zip(layers, lookup):
                old = log_messages[rows]
                extrinsic = normalize_log(post[cols] - old)
                weighted = np.take_along_axis(extrinsic, invperm, axis=-1)
                transformed = fwht(weighted)
                unit = np.ones_like(transformed[:, :1])
                left = np.concatenate((unit, np.cumprod(transformed[:, :-1], axis=1)), axis=1)
                right = np.concatenate((np.cumprod(transformed[:, :0:-1], axis=1)[:, ::-1], unit), axis=1)
                conv = fwht(left*right) / 256
                messages = np.maximum(np.take_along_axis(conv, indices, axis=-1), 1e-300)
                messages /= messages.sum(axis=-1, keepdims=True)
                if damping:
                    messages = (1-damping)*messages + damping*np.exp(old)
                new = np.log(np.maximum(messages, 1e-300))
                post[cols] += new - old
                log_messages[rows] = new
            estimate = np.argmax(post, axis=1)
            if np.array_equal(self.bp.syndrome(estimate), syndrome):
                return estimate, iteration, True, post
        return estimate, max_iters, False, post

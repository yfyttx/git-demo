"""Bounded reliability-ordered syndrome decoding over GF(256)."""

import numpy as np

from construct_companion_cascade import mul


class SyndromeOSD:
    """Order 0 or 1 reprocessing with weak independent columns as pivots."""

    def __init__(self, bp):
        self.bp = bp
        nrows = len(bp.cols)
        self.matrix = np.zeros((nrows, bp.n), dtype=np.uint8)
        self.matrix[np.arange(nrows)[:, None], bp.cols] = bp.h
        self.product = np.array([[mul(a, b) for b in range(256)] for a in range(256)],
                                dtype=np.uint8)
        self.inverse = np.argmax(self.product == 1, axis=1)

    def reduce(self, syndrome, column_order):
        """Reduce the augmented system with deterministic column tie breaking."""
        nrows, ncols = self.matrix.shape
        augmented = np.empty((nrows, ncols+1), dtype=np.uint8)
        augmented[:, :-1] = self.matrix[:, column_order]
        augmented[:, -1] = syndrome
        pivots = []
        for col in range(ncols):
            rank = len(pivots)
            rows = np.flatnonzero(augmented[rank:, col])
            if not rows.size:
                continue
            pivot = rank + rows[0]
            augmented[[rank, pivot]] = augmented[[pivot, rank]]
            augmented[rank, col:] = self.product[self.inverse[augmented[rank, col]],
                                                  augmented[rank, col:]]
            targets = np.flatnonzero(augmented[:, col])
            targets = targets[targets != rank]
            augmented[targets, col:] ^= self.product[augmented[targets, col, None],
                                                     augmented[rank, col:]]
            pivots.append(col)
            if len(pivots) == nrows:
                break
        if np.any(augmented[len(pivots):, -1]):
            raise ValueError('Inconsistent outer syndrome')
        return augmented[:len(pivots)], np.asarray(pivots, dtype=int)

    def decode(self, syndrome, hard, log_posterior, prior, *, order):
        """Use BP posterior for search and channel likelihood for ranking."""
        ranked_values = np.argsort(-log_posterior, axis=1, kind='stable')[:, :5]
        confidence = (log_posterior[np.arange(self.bp.n), ranked_values[:, 0]] -
                      log_posterior[np.arange(self.bp.n), ranked_values[:, 1]])
        column_order = np.argsort(confidence, kind='stable')
        reduced, pivot_positions = self.reduce(syndrome, column_order)
        free_positions = np.setdiff1d(np.arange(self.bp.n), pivot_positions)
        pivots = column_order[pivot_positions]
        free = column_order[free_positions]
        base = hard.copy()
        rhs = np.bitwise_xor.reduce(self.product[reduced[:, free_positions], hard[free]], axis=1)
        base[pivots] = reduced[:, -1] ^ rhs

        candidates = [base]
        if order == 1:
            for position, symbol in zip(free_positions[:32], free[:32]):
                for alternative in ranked_values[symbol, 1:]:
                    candidate = base.copy()
                    candidate[symbol] = alternative
                    candidate[pivots] ^= self.product[reduced[:, position], hard[symbol] ^ alternative]
                    candidates.append(candidate)

        candidates = np.asarray(candidates)
        # For each block this is the likelihood of the reconstructed head and
        # tail up to the known parity normalization constant.
        scores = np.log(np.maximum(prior, 1e-300))[np.arange(self.bp.n), candidates].sum(axis=1)
        estimate = candidates[int(np.argmax(scores))]
        if not np.array_equal(self.bp.syndrome(estimate), syndrome):
            raise ValueError('OSD estimate did not satisfy the outer syndrome')
        return estimate, len(candidates)

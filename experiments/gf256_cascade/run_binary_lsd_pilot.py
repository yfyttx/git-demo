#!/usr/bin/env python3
"""Run a paired Z-side binary BP+LSD decode on one exported simulator frame.

The C++ simulator writes decoder inputs and simulation truth to separate files.
Only syndrome and the already-decoded X estimate are passed to BP+LSD; truth is
read afterward solely to score syndrome and cascade stabilizer equivalence.

Requires ldpc==2.1.0, numpy, and scipy.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.sparse import csr_matrix
from ldpc.bplsd_decoder import BpLsdDecoder


FIELD_POLY = 0x11D


def read_tsv(path: Path) -> tuple[dict[str, str], dict[str, np.ndarray]]:
    scalars: dict[str, str] = {}
    vectors: dict[str, np.ndarray] = {}
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) == 2:
            scalars[fields[0]] = fields[1]
        elif len(fields) > 2:
            vectors[fields[0]] = np.asarray([int(x) for x in fields[1:]], dtype=np.uint8)
        else:
            raise ValueError(f"Malformed export line in {path}: {line!r}")
    return scalars, vectors


def read_alist(path: Path, binary: bool) -> csr_matrix:
    with path.open() as stream:
        n, m, q = map(int, stream.readline().split())
        if (binary and q != 2) or (not binary and q != 256):
            raise ValueError(f"Unexpected field size in {path}: {q}")
        stream.readline()  # variable degrees
        stream.readline()  # check degrees
        rows, cols, values = [], [], []
        for row in range(m):
            entries = list(map(int, stream.readline().split()))
            if len(entries) % 2:
                raise ValueError(f"Malformed sparse row {row} in {path}")
            for i in range(0, len(entries), 2):
                col, value = entries[i] - 1, entries[i + 1]
                if not 0 <= col < n:
                    raise ValueError(f"Column index outside matrix in {path}")
                actual_value = value if binary else gf_exp(value)
                if actual_value:
                    rows.append(row)
                    cols.append(col)
                    values.append(actual_value)
        return csr_matrix((np.asarray(values, dtype=np.int32), (rows, cols)), shape=(m, n))


def gf_mul(a: int, b: int) -> int:
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 0x100:
            a ^= FIELD_POLY
    return result


def gf_exp(log_value: int) -> int:
    value = 1
    for _ in range(log_value % 255):
        value = gf_mul(value, 2)
    return value


def rowspace_basis(matrix: csr_matrix) -> tuple[list[np.ndarray], list[int]]:
    dense = matrix.toarray().astype(np.int32)
    pivots: list[int] = []
    rank = 0
    for col in range(dense.shape[1]):
        pivot = next((r for r in range(rank, dense.shape[0]) if dense[r, col]), None)
        if pivot is None:
            continue
        dense[[rank, pivot]] = dense[[pivot, rank]]
        inv = next(x for x in range(1, 256) if gf_mul(int(dense[rank, col]), x) == 1)
        dense[rank, col:] = [gf_mul(int(x), inv) for x in dense[rank, col:]]
        for r in range(rank + 1, dense.shape[0]):
            factor = int(dense[r, col])
            if factor:
                dense[r, col:] ^= np.asarray(
                    [gf_mul(factor, int(x)) for x in dense[rank, col:]], dtype=np.int32
                )
        pivots.append(col)
        rank += 1
        if rank == dense.shape[0]:
            break
    return [dense[r].copy() for r in range(rank)], pivots


def in_rowspace(vector: np.ndarray, basis: list[np.ndarray], pivots: list[int]) -> bool:
    value = vector.astype(np.int32, copy=True)
    for row, pivot in zip(basis, pivots):
        factor = int(value[pivot])
        if factor:
            value[pivot:] ^= np.asarray(
                [gf_mul(factor, int(x)) for x in row[pivot:]], dtype=np.int32
            )
    return not np.any(value)


def bits_to_symbols(bits: np.ndarray, degree: int) -> np.ndarray:
    if bits.size % degree:
        raise ValueError("Bit vector length is not divisible by field degree")
    return bits.reshape(-1, degree).astype(np.int32) @ (1 << np.arange(degree))


def cascade_z_stabilizer(residual: np.ndarray, outer_basis: list[np.ndarray],
                         outer_pivots: list[int]) -> bool:
    n2 = 12
    if residual.size % n2:
        raise ValueError("Unexpected symbol count for n2=12 cascade")
    blocks = residual.size // n2
    for group in range(2):
        for offset in range(5):
            coordinate = np.asarray(
                [residual[j * n2 + 5 * group + offset] ^ residual[j * n2 + 10 + group]
                 for j in range(blocks)], dtype=np.int32
            )
            if not in_rowspace(coordinate, outer_basis, outer_pivots):
                return False
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", type=Path, required=True,
                        help="C++ export prefix (without .decoder_input.tsv)")
    parser.add_argument("--matrix-dir", type=Path, required=True)
    parser.add_argument("--max-iter", type=int, default=30)
    parser.add_argument("--lsd-method", choices=["LSD_0", "LSD_E", "LSD_CS"], default="LSD_0")
    parser.add_argument("--lsd-order", type=int, default=0)
    args = parser.parse_args()

    scalars, inputs = read_tsv(Path(f"{args.prefix}.decoder_input.tsv"))
    _, truth = read_tsv(Path(f"{args.prefix}.truth.tsv"))
    _, baseline = read_tsv(Path(f"{args.prefix}.baseline.tsv"))
    degree = int(scalars["field_degree"])
    n2 = int(scalars["inner_length"])
    if degree != 8 or n2 != 12:
        raise ValueError("This pilot is currently defined for GF(256), inner_length=12")

    load_started = perf_counter()
    hx_binary = read_alist(args.matrix_dir / "HX_Cascade_Physical_Binary.txt", binary=True).astype(np.uint8)
    outer_hx = read_alist(args.matrix_dir / "HX_TIT_Outer_For_HZ_GF256.txt", binary=False)
    matrix_load_seconds = perf_counter() - load_started
    syndrome = inputs["syndrome_z"]
    x_hat = inputs["x_hat_bits"]
    actual_z = truth["error_z_bits"]
    baseline_z = inputs["current_z_bits"]
    if hx_binary.shape != (syndrome.size, x_hat.size):
        raise ValueError(f"Matrix/input shape mismatch: H={hx_binary.shape}, "
                         f"syndrome={syndrome.size}, bits={x_hat.size}")

    p_d = float(scalars["p_D"])
    z_bias = float(scalars["z_bias"])
    p_low = p_d / (z_bias + 2.0)
    p_z = np.where(x_hat != 0, 0.5, (z_bias * p_low) / (1.0 - p_d + z_bias * p_low))
    decode_started = perf_counter()
    decoder = BpLsdDecoder(
        hx_binary,
        error_channel=p_z.astype(np.float64),
        max_iter=args.max_iter,
        bp_method="minimum_sum",
        ms_scaling_factor=0.625,
        schedule="parallel",
        lsd_method=args.lsd_method,
        lsd_order=args.lsd_order,
        bits_per_step=1,
    )
    estimate = np.asarray(decoder.decode(syndrome), dtype=np.uint8).reshape(-1)
    decode_seconds = perf_counter() - decode_started
    syndrome_ok = np.array_equal((hx_binary @ estimate) % 2, syndrome)
    outer_basis, outer_pivots = rowspace_basis(outer_hx)

    def logical_success(candidate: np.ndarray) -> bool:
        if not np.array_equal((hx_binary @ candidate) % 2, syndrome):
            return False
        residual = bits_to_symbols(actual_z ^ candidate, degree)
        return cascade_z_stabilizer(residual, outer_basis, outer_pivots)

    baseline_success = logical_success(baseline_z)
    lsd_success = logical_success(estimate)
    baseline_syndrome = np.array_equal((hx_binary @ baseline_z) % 2, syndrome)
    result = {
        "matrix_shape": list(hx_binary.shape),
        "p_D": p_d,
        "conditional_z_error_rate_min": float(p_z.min()),
        "conditional_z_error_rate_max": float(p_z.max()),
        "bp_lsd_syndrome_ok": bool(syndrome_ok),
        "bp_lsd_logical_success_z": bool(lsd_success),
        "baseline_syndrome_ok_z": bool(baseline_syndrome),
        "baseline_logical_success_z": bool(baseline_success),
        "bp_lsd_converged": bool(decoder.converge),
        "bp_iterations": int(decoder.iter),
        "matrix_load_seconds": matrix_load_seconds,
        "decode_seconds": decode_seconds,
        "lsd_method": args.lsd_method,
        "lsd_order": args.lsd_order,
    }
    output_path = Path(f"{args.prefix}.binary_lsd.json")
    import json
    output_path.write_text(json.dumps(result, indent=2) + "\n")
    estimate_path = Path(f"{args.prefix}.binary_lsd_bits.txt")
    estimate_path.write_text(" ".join(map(str, estimate.tolist())) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

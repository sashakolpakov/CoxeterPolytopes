"""Temporary exact-modular audit of the D0 high-rank involutions."""

from __future__ import annotations

import ast
from collections import Counter
from pathlib import Path

import numpy as np
from numba import njit, prange


SOURCE = Path("/private/tmp/compdataChams.txt")
N = 22
P = 1009


def extract_assignment(text: str, name: str) -> str:
    marker = f"{name}:="
    start = text.index(marker) + len(marker)
    while text[start].isspace():
        start += 1
    depth = 0
    for end in range(start, len(text)):
        if text[end] == "[":
            depth += 1
        elif text[end] == "]":
            depth -= 1
            if depth == 0:
                return text[start : end + 1]
    raise ValueError(name)


@njit(cache=True)
def inverse_mod(value: int, prime: int) -> int:
    result = 1
    base = value % prime
    exponent = prime - 2
    while exponent:
        if exponent & 1:
            result = (result * base) % prime
        base = (base * base) % prime
        exponent >>= 1
    return result


@njit(cache=True)
def rank_mod(matrix: np.ndarray, prime: int) -> int:
    a = matrix.copy()
    rows, cols = a.shape
    rank = 0
    for col in range(cols):
        pivot = -1
        for row in range(rank, rows):
            if a[row, col] % prime:
                pivot = row
                break
        if pivot < 0:
            continue
        if pivot != rank:
            for j in range(col, cols):
                temp = a[rank, j]
                a[rank, j] = a[pivot, j]
                a[pivot, j] = temp
        inv = inverse_mod(a[rank, col], prime)
        for j in range(col, cols):
            a[rank, j] = (a[rank, j] * inv) % prime
        for row in range(rows):
            if row == rank:
                continue
            factor = a[row, col] % prime
            if factor:
                for j in range(col, cols):
                    a[row, j] = (a[row, j] - factor * a[rank, j]) % prime
        rank += 1
        if rank == cols:
            break
    return rank


@njit(cache=True)
def pair_fixed_ranks(involutions: np.ndarray, prime: int) -> np.ndarray:
    count = len(involutions)
    out = np.empty(count * (count - 1) // 2, dtype=np.int16)
    identity = np.eye(N, dtype=np.int64)
    stacked = np.empty((2 * N, N), dtype=np.int64)
    cursor = 0
    for left in range(count):
        for right in range(left + 1, count):
            stacked[:N] = involutions[left] - identity
            stacked[N:] = involutions[right] - identity
            out[cursor] = N - rank_mod(stacked, prime)
            cursor += 1
    return out


@njit(cache=True, parallel=True)
def fixed_rank_upper_bounds(products: np.ndarray, prime: int) -> np.ndarray:
    count = len(products)
    out = np.empty(count, dtype=np.int16)
    for index in prange(count):
        matrix = products[index].copy()
        for diagonal in range(N):
            matrix[diagonal, diagonal] -= 1
        out[index] = N - rank_mod(matrix, prime)
    return out


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    ss = np.asarray(ast.literal_eval(extract_assignment(text, "SAutXD[0]")), dtype=np.int64)
    ts = np.asarray(ast.literal_eval(extract_assignment(text, "TAutXD[0]")), dtype=np.int64)
    identity = np.eye(N, dtype=np.int64)
    selectors: list[tuple[int, int]] = []
    matrices: list[np.ndarray] = []
    for s_index, s in enumerate(ss):
        products = s @ ts
        for t_index in np.flatnonzero(np.trace(products, axis1=1, axis2=2) == 6):
            matrix = products[t_index]
            if np.array_equal(matrix @ matrix, identity):
                selectors.append((s_index, int(t_index)))
                matrices.append(matrix)
    involutions = np.asarray(matrices, dtype=np.int64)
    if len({matrix.tobytes() for matrix in involutions}) != len(involutions):
        raise ArithmeticError("duplicate involutions")
    ranks = pair_fixed_ranks(involutions, P)
    print("factors", len(ss), len(ts), "group", len(ss) * len(ts))
    print("rank14 involutions", len(involutions), "first selectors", selectors[:10])
    print("pair fixed ranks mod", P, dict(sorted(Counter(map(int, ranks)).items())))
    print("pair max", int(ranks.max()), "pair min", int(ranks.min()))

    modular_histogram: Counter[int] = Counter()
    unresolved: list[tuple[int, int, int]] = []
    for start in range(0, len(ss), 50):
        products = (ss[start : start + 50, None] @ ts[None]).reshape(-1, N, N)
        fixed_bounds = fixed_rank_upper_bounds(products, P)
        modular_histogram.update(map(int, fixed_bounds))
        for offset in np.flatnonzero(fixed_bounds > 10):
            s_index = start + int(offset) // len(ts)
            t_index = int(offset) % len(ts)
            unresolved.append((s_index, t_index, int(fixed_bounds[offset])))
    print("all group modular fixed upper-bound histogram", dict(sorted(modular_histogram.items())))
    print("unresolved >10", len(unresolved), Counter(bound for _, _, bound in unresolved))


if __name__ == "__main__":
    main()

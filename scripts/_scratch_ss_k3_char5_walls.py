"""Temporary wall-orbit and strict-facet exploration for the char-5 stock."""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp

from _scratch_ss_k3_char5_d0_pair_rank import extract_assignment


ROOT = Path(__file__).resolve().parents[1]
BASIC = ROOT / "data/external/katsura_kondo_shimada_char5/compdataB.txt"
CHAMS = ROOT / "data/external/katsura_kondo_shimada_char5/compdataChams.txt"
N = 22
REPRESENTATIVES = {
    0: ((0, 206),),
    1: ((0, 42), (1, 9)),
    2: ((0, 0), (0, 1), (0, 26), (12, 38)),  # last is flat index 614 for |T|=48
}


def cycles_and_normals(walls: np.ndarray, matrix: np.ndarray):
    lookup = {tuple(map(int, wall)): index for index, wall in enumerate(walls)}
    permutation = []
    # The archived wall rows are dual covectors, whereas SAutXD/TAutXD act
    # on primal row vectors.  Thus the dual right action is inverse-transpose.
    dual_action = np.asarray(sp.Matrix(matrix.tolist()).inv().T.tolist(), dtype=np.int64)
    for wall in walls:
        image = tuple(map(int, wall @ dual_action))
        if image not in lookup:
            raise ArithmeticError((image, tuple(map(int, wall))))
        permutation.append(lookup[image])
    unseen = set(range(len(walls)))
    cycles = []
    normals = []
    while unseen:
        left = min(unseen)
        right = permutation[left]
        cycle = tuple(sorted({left, right}))
        if permutation[right] != left:
            raise ArithmeticError("not an involution")
        unseen.difference_update(cycle)
        cycles.append(cycle)
        normals.append(walls[list(cycle)].sum(axis=0))
    return cycles, np.asarray(normals, dtype=np.int64)


def find_bad_strict_pair(normals: np.ndarray, base: np.ndarray, dual: np.ndarray):
    q = np.einsum("ij,jk,ik->i", normals, dual, normals)
    gram_rows = normals @ dual
    strict = []
    certificates = {}
    for index, normal in enumerate(normals):
        square = int(q[index])
        if square >= 0:
            continue
        p = int(base @ dual @ normal)
        witness = (-square) * base + p * normal
        witness_square = int(witness @ dual @ witness)
        pairings = gram_rows @ witness
        if pairings[index] != 0 or witness_square <= 0:
            continue
        if np.all(np.delete(pairings, index) > 0):
            certificates[index] = (witness_square, int(pairings.min()), p)
            for other in strict:
                product = int(normal @ dual @ normals[other])
                ratio = Fraction(product * product, int(q[index] * q[other]))
                if 0 < ratio < 1 and ratio not in {
                    Fraction(1, 4), Fraction(1, 2), Fraction(3, 4)
                }:
                    return strict + [index], (other, index, product, int(q[other]), int(q[index]), ratio), certificates
            strict.append(index)
    return strict, None, certificates


def main() -> None:
    basic_text = BASIC.read_text(encoding="utf-8")
    chams_text = CHAMS.read_text(encoding="utf-8")
    gram = sp.Matrix(ast.literal_eval(extract_assignment(basic_text, "GramSX")))
    dual_sp = 5 * gram.inv()
    assert all(value.q == 1 for value in dual_sp)
    dual = np.asarray(dual_sp.tolist(), dtype=np.int64)
    for chamber, selectors in REPRESENTATIVES.items():
        nested = ast.literal_eval(extract_assignment(chams_text, f"walls[{chamber}]"))
        sizes = [len(block) for block in nested]
        walls = np.asarray([wall for block in nested for wall in block], dtype=np.int64)
        ss = np.asarray(ast.literal_eval(extract_assignment(chams_text, f"SAutXD[{chamber}]")), dtype=np.int64)
        ts = np.asarray(ast.literal_eval(extract_assignment(chams_text, f"TAutXD[{chamber}]")), dtype=np.int64)
        base = np.asarray(ast.literal_eval(extract_assignment(chams_text, f"w[{chamber}]"))[:N], dtype=np.int64)
        print("chamber", chamber, "blocks", sizes, "group", len(ss) * len(ts))
        for selector in selectors:
            matrix = ss[selector[0]] @ ts[selector[1]]
            assert np.array_equal(matrix @ matrix, np.eye(N, dtype=np.int64))
            cycles, normals = cycles_and_normals(walls, matrix)
            fixed_by_block = []
            offset = 0
            for size in sizes:
                fixed_by_block.append(sum(len(cycle) == 1 and offset <= cycle[0] < offset + size for cycle in cycles))
                offset += size
            q = np.einsum("ij,jk,ik->i", normals, dual, normals)
            strict, bad, certificates = find_bad_strict_pair(normals, base, dual)
            print(" selector", selector, "cycles", len(cycles), "singletons", sum(len(c)==1 for c in cycles), "fixed blocks", fixed_by_block)
            print(" q", dict(sorted(Counter(map(int, q)).items())))
            print(" strict scanned", len(strict), "bad", bad, "certs", None if bad is None else {bad[0]:certificates[bad[0]],bad[1]:certificates[bad[1]]})


if __name__ == "__main__":
    main()

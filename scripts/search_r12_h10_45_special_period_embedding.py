#!/usr/bin/env python3
"""Search primitive E7(-2) embeddings for the R12 H10/45 special-period lift.

The complement is T=U(2)^2+D6(-1).  A candidate is seven norm-minus-four
vectors with the E7(-2) Gram matrix.  Orthogonal sign reversal is tested
integrally, and its action on A_T is filtered by the fixed-q2 histogram
forced by the archived Picard transposition.  The search is deterministic
for a supplied seed; it is a discovery tool, not the eventual certificate.
"""

from __future__ import annotations

import argparse
import itertools
import json
import time

import numpy as np
import sympy as sp


def gf2_rank(matrix: np.ndarray) -> int:
    work = (matrix.copy() & 1).astype(np.uint8)
    rank = 0
    for column in range(work.shape[1]):
        pivots = np.flatnonzero(work[rank:, column])
        if not len(pivots):
            continue
        pivot = rank + int(pivots[0])
        work[[rank, pivot]] = work[[pivot, rank]]
        for row in range(work.shape[0]):
            if row != rank and work[row, column]:
                work[row] ^= work[rank]
        rank += 1
    return rank


def stock() -> dict[str, np.ndarray]:
    simple = np.zeros((6, 6), dtype=int)
    for index in range(5):
        simple[index, index] = 1
        simple[index, index + 1] = -1
    simple[5, 4] = simple[5, 5] = 1

    gram = np.zeros((10, 10), dtype=int)
    gram[0, 1] = gram[1, 0] = 2
    gram[2, 3] = gram[3, 2] = 2
    gram[4:, 4:] = -simple @ simple.T

    e7_cartan = np.array(
        [
            [2, -1, 0, 0, 0, 0, 0],
            [-1, 2, -1, 0, 0, 0, 0],
            [0, -1, 2, -1, 0, 0, -1],
            [0, 0, -1, 2, -1, 0, 0],
            [0, 0, 0, -1, 2, -1, 0],
            [0, 0, 0, 0, -1, 2, 0],
            [0, 0, -1, 0, 0, 0, 2],
        ],
        dtype=int,
    )
    twice_inverse = np.array(
        (2 * sp.Matrix(e7_cartan).inv()).tolist(), dtype=int
    )

    discriminant_basis = np.array(
        [
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 1, 0, 1, 0],
            [0, 0, 0, 0, 1, 0, 1, 0, 0, 1],
        ],
        dtype=int,
    )
    bits = np.array(list(itertools.product((0, 1), repeat=6)), dtype=int)
    bit_coordinates = {
        tuple((row @ discriminant_basis) % 2): row for row in bits
    }
    q2 = np.array(
        [
            int((row @ discriminant_basis) @ gram @ (row @ discriminant_basis) // 2)
            % 4
            for row in bits
        ],
        dtype=int,
    )
    return {
        "simple": simple,
        "gram": gram,
        "twice_inverse": twice_inverse,
        "discriminant_basis": discriminant_basis,
        "bits": bits,
        "q2": q2,
        "bit_coordinates": bit_coordinates,
    }


def norm_minus_four_candidates(simple: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    rows = []
    for row in itertools.product(range(-1, 2), repeat=10):
        if sum(row[4:]) % 2:
            continue
        square = (
            4 * row[0] * row[1]
            + 4 * row[2] * row[3]
            - sum(value * value for value in row[4:])
        )
        if square == -4:
            rows.append(row)
    standard = np.array(rows, dtype=np.int16)
    inverse = np.array(sp.Matrix(simple).inv().tolist(), dtype=float)
    lattice = np.hstack(
        [
            standard[:, :4],
            np.rint(standard[:, 4:] @ inverse).astype(np.int16),
        ]
    )
    if not np.array_equal(lattice[:, 4:] @ simple, standard[:, 4:]):
        raise ArithmeticError("D6 coordinate conversion failed")
    return standard, lattice


def search(
    trials: int, seed: int, random_central: bool, central_batch: int
) -> dict[str, object]:
    data = stock()
    gram = data["gram"]
    standard, lattice = norm_minus_four_candidates(data["simple"])
    pairing = np.zeros((10, 10), dtype=np.int16)
    pairing[0, 1] = pairing[1, 0] = 2
    pairing[2, 3] = pairing[3, 2] = 2
    for index in range(4, 10):
        pairing[index, index] = -1

    edges = {(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (2, 6)}

    def required(left: int, right: int) -> int:
        return 2 if tuple(sorted((left, right))) in edges else 0

    fixed_standard = np.array(
        [0, 0, 0, 0, 1, 1, 1, 1, 0, 0], dtype=np.int16
    )
    fixed_index = int(np.flatnonzero(np.all(standard == fixed_standard, axis=1))[0])
    order = [1, 3, 6, 0, 4, 5]
    rng = np.random.default_rng(seed)
    counts = {"complete": 0, "primitive": 0, "integral": 0}
    histograms: dict[str, int] = {}
    started = time.time()
    relation_cache: dict[int, tuple[np.ndarray, np.ndarray]] = {}

    def relation_masks(candidate_index: int) -> tuple[np.ndarray, np.ndarray]:
        cached = relation_cache.get(candidate_index)
        if cached is None:
            products = standard @ pairing @ standard[candidate_index]
            cached = (products == 0, products == 2)
            relation_cache[candidate_index] = cached
        return cached

    fixed_masks = relation_masks(fixed_index)
    fixed_domains = {
        node: np.flatnonzero(fixed_masks[required(node, 2) // 2])
        for node in order
    }

    central_index = fixed_index
    random_domains: dict[int, np.ndarray] | None = None
    for trial in range(trials):
        if random_central and trial % central_batch == 0:
            central_index = int(rng.integers(len(standard)))
            central_masks = relation_masks(central_index)
            random_domains = {
                node: np.flatnonzero(central_masks[required(node, 2) // 2])
                for node in order
            }
        central = standard[central_index]
        domains = (
            random_domains.copy()  # type: ignore[union-attr]
            if random_central
            else fixed_domains.copy()
        )
        selected: dict[int, int] = {}
        good = True
        for position, node in enumerate(order):
            choices = domains[node]
            if not len(choices):
                good = False
                break
            choice = int(rng.choice(choices))
            selected[node] = choice
            choice_masks = relation_masks(choice)
            for later in order[position + 1 :]:
                previous = domains[later]
                mask = choice_masks[required(later, node) // 2]
                domains[later] = previous[mask[previous]]
                if not len(domains[later]):
                    good = False
                    break
            if not good:
                break
        if not good:
            continue

        counts["complete"] += 1
        embedding = np.array(
            [
                lattice[central_index] if node == 2 else lattice[selected[node]]
                for node in range(7)
            ],
            dtype=int,
        )
        if gf2_rank(embedding) != 7:
            continue
        counts["primitive"] += 1

        numerator = gram @ embedding.T @ data["twice_inverse"] @ embedding
        if np.any(numerator % 2):
            continue
        counts["integral"] += 1
        involution = np.eye(10, dtype=int) + numerator // 2
        discr_action = np.array(
            [
                data["bit_coordinates"][tuple((row @ involution) % 2)]
                for row in data["discriminant_basis"]
            ],
            dtype=int,
        )
        fixed = np.all((data["bits"] @ discr_action) % 2 == data["bits"], axis=1)
        histogram = tuple(
            int(value)
            for value in np.bincount(data["q2"][fixed], minlength=4)
        )
        key = ",".join(map(str, histogram))
        histograms[key] = histograms.get(key, 0) + 1
        if histogram == (8, 8, 8, 8):
            return {
                "found": True,
                "trial": trial,
                "seed": seed,
                "random_central": random_central,
                "central_batch": central_batch,
                "counts": counts,
                "histograms": histograms,
                "elapsed_seconds": time.time() - started,
                "embedding": embedding.tolist(),
                "involution": involution.tolist(),
                "discriminant_action": discr_action.tolist(),
            }

    return {
        "found": False,
        "trials": trials,
        "seed": seed,
        "random_central": random_central,
        "central_batch": central_batch,
        "counts": counts,
        "histograms": histograms,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=20260807)
    parser.add_argument("--random-central", action="store_true")
    parser.add_argument("--central-batch", type=int, default=1_000)
    args = parser.parse_args()
    print(
        json.dumps(
            search(
                args.trials,
                args.seed,
                args.random_central,
                args.central_batch,
            ),
            indent=2,
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exactly reject one H9 cap construction, independently of its numeric search.

All base roots have square two.  In their first-nine-root basis, a cap is
(v, sqrt(T)), where v has its prescribed nine side pairings and
T = 2 - (v,v).  Exact Q(phi) square comparisons eliminate every radical
from the Coxeter-angle and strict-witness checks.  No NumPy/SciPy or search
transcript is imported by this verifier.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.exact.matrix import inertia, rank


ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / "data/enumerations/bugaenko_h9_caps_20261013.json"
ZERO, ONE, TWO, PHI = Golden(0), Golden(1), Golden(2), Golden(0, 1)
LABELS = {2: ZERO, 3: -ONE, 5: -PHI}


def dot(left, right):
    return sum((a * b for a, b in zip(left, right)), ZERO)


def row_times_matrix(vector, matrix):
    return [dot(vector, column) for column in zip(*matrix)]


def inverse(matrix):
    size = len(matrix)
    work = [
        list(row) + [ONE if i == j else ZERO for j in range(size)]
        for i, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(i for i in range(column, size) if work[i][column] != ZERO)
        work[pivot], work[column] = work[column], work[pivot]
        divisor = work[column][column]
        work[column] = [entry / divisor for entry in work[column]]
        for i in range(size):
            if i != column:
                multiplier = work[i][column]
                work[i] = [
                    x - multiplier * y for x, y in zip(work[i], work[column])
                ]
    return [row[size:] for row in work]


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def verify(path: Path) -> dict:
    certificate = json.loads(path.read_text())
    source = ROOT / certificate["source"]
    require(
        sha256(source.read_bytes()).hexdigest() == certificate["source_sha256"],
        "the calibrated source changed",
    )
    data = json.loads(source.read_text())
    gram = [[Golden.from_json(entry) for entry in row] for row in data["wall_gram"]]
    base = [row[:9] for row in gram[:9]]
    signature = inertia(base)
    require((signature.positive, signature.negative, signature.zero) == (8, 1, 0),
            "the base is not Lorentzian of rank nine")
    require(rank(gram) == 9, "the eleven-root Gram has changed rank")
    base_inverse = inverse(base)
    inverse_cache = {}
    caps = []
    for recipe in certificate["cap_recipes"]:
        omitted = tuple(recipe["omitted"])
        chosen = [i for i in range(11) if i not in omitted]
        require(len(chosen) == len(recipe["labels"]) == 9, "invalid cap recipe")
        if omitted not in inverse_cache:
            inverse_cache[omitted] = inverse([[gram[i][j] for j in chosen] for i in chosen])
        prescribed = [LABELS[label] for label in recipe["labels"]]
        coefficients = row_times_matrix(prescribed, inverse_cache[omitted])
        pairings = row_times_matrix(coefficients, [gram[i] for i in chosen])
        vector = row_times_matrix(pairings[:9], base_inverse)
        square = dot(vector, pairings[:9])
        require(square.sign() < 0, "a cap projection is not timelike")
        require(all(c in LABELS.values() or (c + TWO).sign() <= 0 for c in pairings),
                "a cap has a forbidden angle with a side wall")
        caps.append((vector, pairings, TWO - square))

    angles = Counter()
    for i, j in combinations(range(len(caps)), 2):
        vector, _, left_t = caps[i]
        _, right_pairings, right_t = caps[j]
        rational_part = dot(vector, right_pairings[:9])
        radical_square = left_t * right_t
        matched = False
        for label, target in LABELS.items():
            difference = target - rational_part
            if difference.sign() > 0 and difference * difference == radical_square:
                angles[str(label)] += 1
                matched = True
                break
        if not matched:
            shifted = rational_part + TWO
            require(shifted.sign() < 0 and (shifted * shifted - radical_square).sign() > 0,
                    f"caps {i}, {j} have a forbidden or unchecked angle")
            angles["ultraparallel"] += 1

    witnesses = {}
    for key, expected_sign in (("strict_hyperideal_witness", 1), ("strict_timelike_witness", -1)):
        encoded = certificate[key]
        denominator = encoded["denominator"]
        vector = [Golden(Fraction(x, denominator)) for x in encoded["basis_coefficients_numerators"]]
        z = Golden(Fraction(encoded["extra_coordinate_numerator"], denominator))
        require(z.sign() > 0, "a witness violates the base wall")
        side_pairings = row_times_matrix(vector, gram[:9])
        require(all(value.sign() < 0 for value in side_pairings), "a witness is not strict at a side")
        cap_margins = []
        for _, pairings, t in caps:
            value = dot(vector, pairings[:9])
            margin = value * value - z * z * t
            require(value.sign() < 0 and margin.sign() > 0, "a witness violates a cap")
            cap_margins.append(margin)
        norm = dot(vector, row_times_matrix(vector, base)) + z * z
        require(norm.sign() == expected_sign, "a witness has the wrong Lorentz sign")
        witnesses[key] = {
            "square": norm.to_json(),
            "minimum_side_margin": min(-value for value in side_pairings).to_json(),
            "minimum_squared_cap_margin": min(cap_margins).to_json(),
            "all_34_inequalities_strict": True,
        }
    return {
        "status": "PASS_EXACT_INFINITE_VOLUME",
        "dimension": 9,
        "wall_count": 12 + len(caps),
        "cap_pair_angles": dict(sorted(angles.items())),
        "witnesses": witnesses,
        "proof": "The strict timelike-to-spacelike segment meets the light cone strictly inside all halfspaces, giving an open ideal-boundary patch and infinite hyperbolic volume.",
        "scope": "This one candidate only; no claim that H9 or the cap-construction method is impossible.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", nargs="?", type=Path, default=DEFAULT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify(args.certificate)
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()

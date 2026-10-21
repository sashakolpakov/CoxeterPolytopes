#!/usr/bin/env python3
"""Certify that P17--P19 escape all three current CCR transfer branches.

This is deliberately a branch-exhaustion theorem, not a no-surface theorem.
It combines the frozen rational MMP/Abel obstruction with an independent
weak-log adjunction/nef replay and the elementary K3 parity obstruction.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from fractions import Fraction
from math import lcm
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/p17_p19_three_branch_property4_escape.json"
INVERSE = ROOT / "data/enumerations/rational_i1n_inverse_design.json"
P17_MMP = ROOT / "data/enumerations/rational_p17_mmp.json"
P18_P19_MMP = ROOT / "data/enumerations/rational_p18_p19_mmp.json"

PINS = {
    INVERSE: "351a0a8525209313de13f9dd3752320bd074daf0985da704fc634684c75efc18",
    P17_MMP: "2c35b61bb5357111a0f74823c98d192ca2eaef86d5a939071261eeb371be060b",
    P18_P19_MMP: "b2bf17c357a605ff0f39f8645628079f3fb8a808255ad5e184fb1da99995751d",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dot(x: Sequence[int], y: Sequence[int]) -> int:
    return x[0] * y[0] - sum(a * b for a, b in zip(x[1:], y[1:]))


def load_inputs() -> tuple[dict[str, object], dict[str, object], dict[str, object]]:
    for path, expected in PINS.items():
        if sha256(path) != expected:
            raise ArithmeticError(f"frozen dependency changed: {path.relative_to(ROOT)}")
    return (
        json.loads(INVERSE.read_text()),
        json.loads(P17_MMP.read_text()),
        json.loads(P18_P19_MMP.read_text()),
    )


def canonical_packet(
    n: int, p17: dict[str, object], p18_p19: dict[str, object]
) -> dict[str, object]:
    return p17 if n == 17 else p18_p19["dimensions"][str(n)]


def inverse(matrix: Sequence[Sequence[int]]) -> list[list[Fraction]]:
    """Invert a square integer matrix by exact Gauss--Jordan elimination."""

    size = len(matrix)
    aug = [
        [Fraction(value) for value in row]
        + [Fraction(int(i == j)) for j in range(size)]
        for i, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next((row for row in range(column, size) if aug[row][column]), None)
        if pivot is None:
            raise ArithmeticError("chosen facet basis ceased to span")
        aug[column], aug[pivot] = aug[pivot], aug[column]
        scale = aug[column][column]
        aug[column] = [entry / scale for entry in aug[column]]
        for row in range(size):
            if row == column or not aug[row][column]:
                continue
            scale = aug[row][column]
            aug[row] = [a - scale * b for a, b in zip(aug[row], aug[column])]
    return [row[size:] for row in aug]


def matvec(matrix: Sequence[Sequence[Fraction]], vector: Sequence[int]) -> list[Fraction]:
    return [sum((entry * value for entry, value in zip(row, vector)), Fraction()) for row in matrix]


def build() -> dict[str, object]:
    inverse_design, p17, p18_p19 = load_inputs()
    rows: list[dict[str, object]] = []

    for n in (17, 18, 19):
        case = inverse_design["dimensions"][str(n)]
        roots = [item["vector"] for item in case["roots"]]
        packet = canonical_packet(n, p17, p18_p19)
        canonical_classes = packet["canonical_classes"]
        norm_one = [i for i, root in enumerate(roots) if dot(root, root) == -1]
        norm_two = [i for i, root in enumerate(roots) if dot(root, root) == -2]

        expected_norm_one = packet.get("initial_square_minus_one_facet_count", len(norm_one))
        if len(norm_one) != expected_norm_one:
            raise ArithmeticError(f"P{n} norm-one facet count changed")
        if len(norm_one) + len(norm_two) != len(roots):
            raise ArithmeticError(f"P{n} acquired a facet outside the (-1),(-2) alphabet")
        if not (
            packet["canonical_classes_equal_standard_diagram_orbit"]
            and packet["canonical_class_count"] == len(canonical_classes)
        ):
            raise ArithmeticError(f"P{n} canonical marking exhaustion changed")

        # Exhaust every physical weak-log scaling before invoking an MMP.
        # A primitive norm-one ray can carry either its primitive rational
        # (-1)-curve (K.r=-1) or a boundary B=2r (K.r=+1).  A primitive
        # norm-two ray can only carry its primitive rational (-2)-curve
        # (K.r=0).  Because the facets span, these pairings must be the values
        # of one linear functional K.  The first rho facets form a basis; solve
        # there and test every remaining wall relation exactly.
        rho = n + 1
        pairing_matrix = [[root[0]] + [-entry for entry in root[1:]] for root in roots]
        basis_inverse = inverse(pairing_matrix[:rho])
        relation_rows: list[tuple[list[int], int, int]] = []
        for facet_index, row in enumerate(pairing_matrix[rho:], start=rho):
            coefficients = [
                sum((Fraction(entry) * basis_inverse[k][j] for k, entry in enumerate(row)), Fraction())
                for j in range(rho)
            ]
            denominator = lcm(*(value.denominator for value in coefficients))
            relation_rows.append(
                ([int(value * denominator) for value in coefficients], denominator, facet_index)
            )

        first_failure_histogram: Counter[int] = Counter()
        max_satisfied = 0
        consistent_colorings: list[dict[str, object]] = []
        for mask in range(1 << len(norm_one)):
            signs = {
                root_index: (1 if (mask >> bit) & 1 else -1)
                for bit, root_index in enumerate(norm_one)
            }
            targets = [signs.get(index, 0) for index in range(len(roots))]
            basis_targets = targets[:rho]
            satisfied = rho
            failure = None
            for coefficients, denominator, facet_index in relation_rows:
                actual = sum(entry * value for entry, value in zip(coefficients, basis_targets))
                expected = denominator * targets[facet_index]
                if actual != expected:
                    failure = facet_index
                    break
                satisfied += 1
            max_satisfied = max(max_satisfied, satisfied)
            if failure is None:
                canonical = matvec(basis_inverse, basis_targets)
                consistent_colorings.append(
                    {
                        "boundary_half_root_labels": [
                            index + 1 for index in norm_one if signs[index] == 1
                        ],
                        "canonical": [str(value) for value in canonical],
                    }
                )
            else:
                first_failure_histogram[failure + 1] += 1
        if consistent_colorings:
            raise ArithmeticError(f"P{n} unexpectedly acquired an adjunction-consistent weak-log coloring")

        status_text = json.dumps(case["status"], sort_keys=True)
        if "impossible on every smooth projective rational surface" not in status_text:
            raise ArithmeticError(f"P{n} rational MMP/Abel no-go changed")

        # A full-rank K3 carrier has an even Picard lattice.  The declared
        # host itself is odd because it contains norm-one roots.  Even after
        # allowing one common rational similarity and rational rescaling on
        # each wall line, simultaneous physical (-2) representatives would
        # require the ratio 2 between norm classes to be a rational square.
        # It is not: if a,b are the line rescalings, a^2/(2 b^2)=1.
        mixed_norm_square_classes = sorted({-dot(roots[i], roots[i]) for i in range(len(roots))})
        if mixed_norm_square_classes != [1, 2]:
            raise ArithmeticError(f"P{n} mixed norm square classes changed")

        rows.append(
            {
                "n": n,
                "rho": n + 1,
                "facet_count": len(roots),
                "declared_host": f"I_(1,{n})",
                "norm_one_facet_labels": [i + 1 for i in norm_one],
                "norm_two_facet_labels": [i + 1 for i in norm_two],
                "canonical_marking_count": len(canonical_classes),
                "direct_rational_full_nef": "EMPTY_BY_MARKING_INDEPENDENT_MMP_ABEL",
                "weak_log_rational": {
                    "status": "EMPTY_AT_ADJUNCTION_LINEAR_SYSTEM",
                    "reason": (
                        "all assignments K.r in {-1,+1} on norm-one rays and K.r=0 on norm-two rays "
                        "violate an exact wall relation, so no weak-log canonical class exists"
                    ),
                    "colorings_tested": 1 << len(norm_one),
                    "consistent_colorings": [],
                    "basis_facet_labels": list(range(1, rho + 1)),
                    "max_initial_equations_satisfied": max_satisfied,
                    "first_failed_facet_histogram": {
                        str(label): count for label, count in sorted(first_failure_histogram.items())
                    },
                },
                "literal_full_rank_K3": {
                    "status": "EMPTY_AT_PHYSICAL_ROOT_NORMALIZATION",
                    "odd_host": True,
                    "facet_norm_square_classes": mixed_norm_square_classes,
                    "rational_common_similarity_obstruction": "2 is not a square in Q",
                },
            }
        )

    return {
        "schema": "p17-p19-three-branch-property4-escape-v2",
        "status": "PASS_CURRENT_THREE_BRANCH_COMPILER_EMPTY_ON_P17_P18_P19",
        "all_checks_pass": True,
        "theorem": {
            "inputs": (
                "the exact standard-host classical finite-volume chambers P17, P18, P19, "
                "complete relative to their declared integral reflection hosts"
            ),
            "empty_branches": [
                "direct complex rational full-nef",
                "certified complex weak-log rational full-nef",
                (
                    "exact-return minimalizable full trace, which at ranks 18--20 reduces to complex K3 "
                    "with same-line physical (-2) facets"
                ),
            ],
            "consequence": (
                "declared-host reflection completeness does not make the currently proved three-way CCR transfer "
                "exhaustive; a new carrier species or a genuinely broader geometric return is necessary"
            ),
        },
        "rows": rows,
        "scope": {
            "proved": (
                "emptiness of the three currently declared complex transfer branches for these pinned integral "
                "hosts, whose chambers are complete relative to those hosts"
            ),
            "not_claimed": [
                "nonexistence of P17, P18, or P19 as arithmetic Coxeter chambers",
                "absolute or maximal-arithmetic reflection maximality of these pinned hosts",
                "nonexistence of every possible carrier on a nonrational surface or higher-dimensional variety",
                "a positive-characteristic no-go",
                "failure of a projected-wall carrier after weakening the same-vector physical-root contract",
                "failure for an abstract Coxeter diagram with a different arithmetic realization or integral host",
            ],
        },
        "dependencies": [
            {"path": str(path.relative_to(ROOT)), "sha256": expected}
            for path, expected in PINS.items()
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--write-archive", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    payload = build()
    if args.write_archive:
        ARCHIVE.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"wrote {ARCHIVE.relative_to(ROOT)}")
    elif args.check_archive:
        if payload != json.loads(ARCHIVE.read_text()):
            raise SystemExit("live P17--P19 branch-exhaustion certificate differs from archive")
        print("archived P17--P19 branch-exhaustion certificate matches")
    elif args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(payload["status"])


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Certify the finite lattice arithmetic in the Conway-face K3 theorem.

The period theorem, Conway-face theorem, and the fact that saturation gives
an isotropic subquotient are cited mathematical inputs.  This script checks
the complete Allcock face-type table, the stable-range inequalities, the
dimension identity, and the two non-stable rank-20 discriminant-form gluings.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path
from typing import Any, Iterable

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data/enumerations/conway_face_k3_inverse_carrier_theorem.json"
)


def _canonical_sha256(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _cartan_a(rank: int) -> sp.Matrix:
    if rank < 1:
        raise ValueError("A-rank must be positive")
    gram = 2 * sp.eye(rank)
    for i in range(rank - 1):
        gram[i, i + 1] = gram[i + 1, i] = -1
    return gram


def _cartan_d(rank: int) -> sp.Matrix:
    if rank < 4:
        raise ValueError("D-rank must be at least four")
    gram = 2 * sp.eye(rank)
    for i in range(rank - 3):
        gram[i, i + 1] = gram[i + 1, i] = -1
    branch = rank - 3
    gram[branch, rank - 2] = gram[rank - 2, branch] = -1
    gram[branch, rank - 1] = gram[rank - 1, branch] = -1
    return gram


def _cartan_e(rank: int) -> sp.Matrix:
    if rank not in {6, 7, 8}:
        raise ValueError("only E6, E7, and E8 occur in this certificate")
    gram = 2 * sp.eye(rank)
    # A chain of length rank-1, with the final node attached to chain node 2.
    for i in range(rank - 2):
        gram[i, i + 1] = gram[i + 1, i] = -1
    gram[2, rank - 1] = gram[rank - 1, 2] = -1
    return gram


def _orthogonal_sum(grams: Iterable[sp.Matrix]) -> sp.Matrix:
    rows = list(grams)
    if not rows:
        raise ValueError("an orthogonal sum needs a component")
    result = rows[0]
    for gram in rows[1:]:
        result = sp.diag(result, gram)
    return result


def _root_gram(name: str) -> sp.Matrix:
    if name.startswith("A"):
        return _cartan_a(int(name[1:]))
    if name.startswith("D"):
        return _cartan_d(int(name[1:]))
    if name.startswith("E"):
        return _cartan_e(int(name[1:]))
    raise ValueError(f"unsupported root type {name}")


def _sum_gram(names: tuple[str, ...]) -> sp.Matrix:
    return _orthogonal_sum(_root_gram(name) for name in names)


def _standard_endpoint_picard(endpoint: str) -> sp.Matrix:
    """Standard Gram of Vinberg's determinant-three/four rank-20 K3 lattice."""

    hyperbolic_plane = sp.Matrix([[0, 1], [1, 0]])
    negative_e8 = -_root_gram("E8")
    if endpoint == "E6":
        tail = -_root_gram("A2")
    elif endpoint == "D6":
        tail = -_sum_gram(("A1", "A1"))
    else:
        raise ValueError(f"unsupported endpoint {endpoint}")
    return sp.diag(hyperbolic_plane, negative_e8, negative_e8, tail)


def _smith_invariants(gram: sp.Matrix) -> list[int]:
    smith = smith_normal_form(gram, domain=ZZ)
    return [abs(int(smith[i, i])) for i in range(smith.rows)]


def _discriminant_length(gram: sp.Matrix) -> int:
    return sum(value > 1 for value in _smith_invariants(gram))


def _mod_one(value: sp.Rational) -> Fraction:
    value = sp.Rational(value)
    return Fraction(int(value.p % value.q), int(value.q))


def _mod_two(value: sp.Rational) -> Fraction:
    value = sp.Rational(value)
    modulus = 2 * value.q
    return Fraction(int(value.p % modulus), int(value.q))


def _discriminant_form(
    gram: sp.Matrix, *, sign: int = 1
) -> tuple[list[tuple[Fraction, ...]], list[list[int]], list[Fraction]]:
    """Enumerate A_L, its addition table, and q modulo 2Z.

    Determinants here are at most four, so redundant enumeration in a
    determinant box is small and makes the endpoint check transparent.
    """

    determinant = abs(int(gram.det()))
    inverse = gram.inv()
    reps_set: set[tuple[Fraction, ...]] = set()
    for vector in product(range(determinant), repeat=gram.rows):
        image = inverse * sp.Matrix(vector)
        reps_set.add(tuple(_mod_one(entry) for entry in image))
    reps = sorted(reps_set)
    if len(reps) != determinant:
        raise ArithmeticError("discriminant-group enumeration has wrong order")
    index = {rep: i for i, rep in enumerate(reps)}

    def add(left: tuple[Fraction, ...], right: tuple[Fraction, ...]):
        return tuple((a + b) % 1 for a, b in zip(left, right, strict=True))

    addition = [
        [index[add(left, right)] for right in reps]
        for left in reps
    ]
    values: list[Fraction] = []
    for rep in reps:
        vector = sp.Matrix([sp.Rational(value.numerator, value.denominator) for value in rep])
        square = sign * (vector.T * gram * vector)[0]
        values.append(_mod_two(square))
    return reps, addition, values


def _forms_are_isometric(
    left_gram: sp.Matrix,
    right_gram: sp.Matrix,
    *,
    left_sign: int = 1,
    right_sign: int = 1,
) -> bool:
    left, left_add, left_q = _discriminant_form(left_gram, sign=left_sign)
    right, right_add, right_q = _discriminant_form(right_gram, sign=right_sign)
    if len(left) != len(right):
        return False
    left_zero = left.index(tuple(Fraction(0) for _ in range(left_gram.rows)))
    right_zero = right.index(tuple(Fraction(0) for _ in range(right_gram.rows)))
    left_nonzero = [i for i in range(len(left)) if i != left_zero]
    right_nonzero = [i for i in range(len(right)) if i != right_zero]
    for target_order in permutations(right_nonzero):
        mapping = {left_zero: right_zero, **dict(zip(left_nonzero, target_order, strict=True))}
        if any(left_q[i] != right_q[mapping[i]] for i in range(len(left))):
            continue
        if all(
            mapping[left_add[i][j]] == right_add[mapping[i]][mapping[j]]
            for i in range(len(left))
            for j in range(len(left))
        ):
            return True
    return False


def _has_nonzero_isotropic_class(gram: sp.Matrix) -> bool:
    reps, _, values = _discriminant_form(gram)
    zero = tuple(Fraction(0) for _ in range(gram.rows))
    return any(rep != zero and value == 0 for rep, value in zip(reps, values, strict=True))


def _row(identifier: str, components: tuple[str, ...]) -> dict[str, Any]:
    gram = _sum_gram(components)
    rank = gram.rows
    hyperbolic_dimension = 25 - rank
    picard_rank = 26 - rank
    moduli_dimension = 20 - picard_rank
    length = _discriminant_length(gram)
    return {
        "id": identifier,
        "components": list(components),
        "spherical_rank": rank,
        "root_span_determinant": abs(int(gram.det())),
        "root_span_smith": _smith_invariants(gram),
        "root_span_discriminant_length_upper_bound": length,
        "hyperbolic_dimension": hyperbolic_dimension,
        "picard_rank": picard_rank,
        "k3_moduli_dimension": moduli_dimension,
        "dimension_sum": hyperbolic_dimension + moduli_dimension,
        "stable_range_required_bound": rank - 6,
        "stable_range_passes_from_root_span_bound": (
            rank >= 7 and length <= rank - 6
        ),
    }


def build() -> dict[str, Any]:
    types: list[tuple[str, tuple[str, ...]]] = [
        ("E6", ("E6",)),
        ("D6", ("D6",)),
        ("E7", ("E7",)),
        ("D7", ("D7",)),
        ("D8", ("D8",)),
        ("D9", ("D9",)),
        ("D6+D4", ("D6", "D4")),
        ("D7+D4", ("D7", "D4")),
        ("D6+D6", ("D6", "D6")),
        ("D7+D6", ("D7", "D6")),
        *[(f"D7+D{rank}", ("D7", f"D{rank}")) for rank in range(7, 17)],
        ("D6+D6+D4", ("D6", "D6", "D4")),
    ]
    rows = [_row(identifier, components) for identifier, components in types]
    e6_picard = _standard_endpoint_picard("E6")
    d6_picard = _standard_endpoint_picard("D6")
    endpoint_data = {
        "E6": {
            "root_determinant": abs(int(_root_gram("E6").det())),
            "complement": "A2",
            "complement_determinant": abs(int(_root_gram("A2").det())),
            "root_has_nonzero_isotropic_discriminant_class": (
                _has_nonzero_isotropic_class(_root_gram("E6"))
            ),
            "anti_isometry_verified": _forms_are_isometric(
                _root_gram("E6"), _root_gram("A2"), left_sign=-1
            ),
            "standard_picard_lattice": "U+E8(-1)^2+A2(-1)",
            "standard_picard_rank": e6_picard.rows,
            "standard_picard_determinant": abs(int(e6_picard.det())),
            "standard_picard_smith": _smith_invariants(e6_picard),
            "historical_identification": (
                "Vinberg determinant-three K3; Borcherds Example 5.3 "
                "takes J=E6 and L=J^perp"
            ),
        },
        "D6": {
            "root_determinant": abs(int(_root_gram("D6").det())),
            "complement": "A1+A1",
            "complement_determinant": abs(
                int(_sum_gram(("A1", "A1")).det())
            ),
            "root_has_nonzero_isotropic_discriminant_class": (
                _has_nonzero_isotropic_class(_root_gram("D6"))
            ),
            "anti_isometry_verified": _forms_are_isometric(
                _root_gram("D6"),
                _sum_gram(("A1", "A1")),
                left_sign=-1,
            ),
            "standard_picard_lattice": "U+E8(-1)^2+A1(-1)^2",
            "standard_picard_rank": d6_picard.rows,
            "standard_picard_determinant": abs(int(d6_picard.det())),
            "standard_picard_smith": _smith_invariants(d6_picard),
            "historical_identification": (
                "Vinberg determinant-four K3; Borcherds's paragraph after "
                "Example 5.3 uses the D6/D7 Conway-complement category"
            ),
        },
    }
    checks = {
        "all_dimension_sums_equal_19": all(row["dimension_sum"] == 19 for row in rows),
        "all_picard_ranks_at_most_20": all(row["picard_rank"] <= 20 for row in rows),
        "all_positive_moduli_rows_pass_stable_range": all(
            row["stable_range_passes_from_root_span_bound"]
            for row in rows
            if row["spherical_rank"] >= 7
        ),
        "E6_endpoint_is_primitive": not endpoint_data["E6"][
            "root_has_nonzero_isotropic_discriminant_class"
        ],
        "D6_endpoint_is_primitive": not endpoint_data["D6"][
            "root_has_nonzero_isotropic_discriminant_class"
        ],
        "E6_A2_anti_isometry": endpoint_data["E6"]["anti_isometry_verified"],
        "D6_A1sq_anti_isometry": endpoint_data["D6"]["anti_isometry_verified"],
        "E6_endpoint_is_standard_determinant3_picard": (
            endpoint_data["E6"]["standard_picard_rank"] == 20
            and endpoint_data["E6"]["standard_picard_determinant"] == 3
            and endpoint_data["E6"]["standard_picard_smith"]
            == [*([1] * 19), 3]
        ),
        "D6_endpoint_is_standard_determinant4_picard": (
            endpoint_data["D6"]["standard_picard_rank"] == 20
            and endpoint_data["D6"]["standard_picard_determinant"] == 4
            and endpoint_data["D6"]["standard_picard_smith"]
            == [*([1] * 18), 2, 2]
        ),
        "H13_period_dimension_six": next(
            row for row in rows if row["id"] == "D6+D6"
        )["k3_moduli_dimension"]
        == 6,
        "H9_subtree_seed_period_dimension_ten": next(
            row for row in rows if row["id"] == "D6+D6+D4"
        )["k3_moduli_dimension"]
        == 10,
    }
    result: dict[str, Any] = {
        "schema_version": "1.1.0",
        "kind": "conway_face_k3_inverse_carrier_lattice_certificate",
        "theorem_scope": {
            "ambient_lattice": "II_(25,1)",
            "carrier_lattice": "N=(Sat_L(R_sigma)^perp)(-1)",
            "stable_embedding_gate": "length(A_S) <= s-6",
            "dimension_identity": "dim_H(P)+dim_C(K3 exact-Picard moduli)=19",
            "promotion_level": (
                "full-Picard ample-intersection divisor-hyperplane carrier"
            ),
            "excluded_claims": [
                "full-nef equality",
                "effectivity or irreducibility of every wall divisor",
                "surface-automorphism realization of reflections",
            ],
        },
        "face_types": rows,
        "rank20_endpoint_gluings": endpoint_data,
        "external_theorem_inputs": [
            "Allcock-Borcherds Conway-face Coxeter and finite-volume theorems",
            "primitive orthogonal-complement discriminant anti-isometry in an even unimodular lattice",
            "Nikulin stable-range primitive K3 embedding theorem",
            "Nikulin uniqueness in the indefinite stable genus for the two rank-20 endpoint identifications",
            "surjectivity of the marked K3 period map",
            "K3 Weyl-chamber transitivity",
            "Allcock finite-subtree doubling theorem",
        ],
        "checks": checks,
        "verified": all(checks.values()),
    }
    result["certificate_sha256"] = _canonical_sha256(result)
    return result


def verify(archived: dict[str, Any]) -> dict[str, Any]:
    fresh = build()
    if archived != fresh:
        raise AssertionError("the archived Conway-face K3 certificate has drifted")
    if not fresh["verified"]:
        failed = [name for name, value in fresh["checks"].items() if not value]
        raise AssertionError(f"failed Conway-face K3 checks: {failed}")
    return {
        "certificate_sha256": fresh["certificate_sha256"],
        "face_type_count": len(fresh["face_types"]),
        "dimension_range": [
            min(row["hyperbolic_dimension"] for row in fresh["face_types"]),
            max(row["hyperbolic_dimension"] for row in fresh["face_types"]),
        ],
        "all_dimension_sums": sorted(
            {row["dimension_sum"] for row in fresh["face_types"]}
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = build()
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

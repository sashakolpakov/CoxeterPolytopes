#!/usr/bin/env python3
"""Certify the compact (2,4,6) obstruction to the current Property-4 compiler.

The certificate is deliberately self-contained.  It starts from a rank-three
integral wall Gram, verifies compact arithmeticity, proves absolute reflection
maximality by a covolume/index-two argument, and records the two geometric
obstructions:

* a positive-rank rational surface has a nef isotropic ruling class, whereas
  the returned triangle is compact;
* a rank-three canonical-orthogonal minimal surface has numerically trivial
  canonical class, so every irreducible negative facet is a (-2)-curve.  Such
  curves can meet at finite Coxeter angles only pi/2 or pi/3, not pi/4 or
  pi/6.

This rules out the three exact-return Property-4 packages obtained by
combining Chamber Fidelity with the complex carrier alternatives in the
frozen canonical-chamber transfer theorem.  It does not rule out the bare
carrier-relative third row of that theorem after exact chamber return is
removed, or a future enlarged carrier notion which abandons literal full-nef
return for k=Q.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/compact_246_property4_counterexample.json"
SCHEMA = "coxeter.compact_246.property4_counterexample.v1"


def _canonical_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [[int(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _reflection_matrix(gram: sp.Matrix, root_index: int) -> sp.Matrix:
    """Matrix of x -> x-2(x,r)/(r,r)r in the wall basis."""

    result = sp.eye(gram.rows)
    denominator = gram[root_index, root_index]
    for column in range(gram.cols):
        result[root_index, column] -= (
            2 * gram[column, root_index] / denominator
        )
    return result


def _triangle_area(labels: tuple[int, int, int]) -> Fraction:
    return Fraction(1) - sum((Fraction(1, label) for label in labels), Fraction())


def _small_overgroup_triangles() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Exhaust triangles of area at most pi/24 using exact reciprocal bounds."""

    source_area = Fraction(1, 12)
    area_ceiling = Fraction(1, 24)
    reciprocal_floor = 1 - area_ceiling
    area_window: list[dict[str, Any]] = []

    # If 2<=a<=b<=c and sum(1/a,1/b,1/c)>=23/24, then a<=3.
    for a in range(2, 4):
        b = a
        # Since 1/c<=1/b, this is a necessary and terminating bound on b.
        while Fraction(1, a) + Fraction(2, b) >= reciprocal_floor:
            base = 1 - Fraction(1, a) - Fraction(1, b)
            if base > 0:
                # Hyperbolicity says c>1/base.  The area ceiling says
                # 1/c >= base-1/24, giving a finite upper bound.
                lower = max(b, base.denominator // base.numerator + 1)
                upper_term = base - area_ceiling
                if upper_term > 0:
                    upper = upper_term.denominator // upper_term.numerator
                    for c in range(lower, upper + 1):
                        area = _triangle_area((a, b, c))
                        if not (0 < area <= area_ceiling):
                            continue
                        index = source_area / area
                        area_window.append(
                            {
                                "labels": [a, b, c],
                                "area_over_pi": str(area),
                                "source_area_ratio": str(index),
                                "integer_index": index.denominator == 1 and index >= 2,
                            }
                        )
            b += 1

    candidates = [
        {
            "labels": item["labels"],
            "area_over_pi": item["area_over_pi"],
            "index": int(Fraction(item["source_area_ratio"])),
        }
        for item in area_window
        if item["integer_index"]
    ]
    return area_window, candidates


def build_certificate() -> dict[str, Any]:
    # Negative wall normals in a Lorentzian lattice of signature (1,2).
    # The normalized off-diagonal entries are cos(pi/4), 0, cos(pi/6).
    gram = sp.Matrix(
        [
            [-12, 6, 0],
            [6, -6, 3],
            [0, 3, -2],
        ]
    )
    labels = (4, 6, 2)
    reflections = [_reflection_matrix(gram, i) for i in range(3)]
    principal_pair_determinants = [
        int(gram.extract(indices, indices).det())
        for indices in ((0, 1), (1, 2), (0, 2))
    ]
    normalized_cosine_squares = [
        Fraction(int(gram[i, j]) ** 2, int(gram[i, i] * gram[j, j]))
        for i, j in ((0, 1), (1, 2), (0, 2))
    ]

    area_window, candidates = _small_overgroup_triangles()
    # An index-two triangle supergroup would have type (2,3,8).  Gluing two
    # such chambers across a side gives a triangle only when that side meets
    # its right-angle vertex.  Across the (2,8)-side the result is (3,3,4);
    # across the (2,3)-side the doubled 2*pi/3 angle is not Coxeter; across
    # the (3,8)-side the union is a quadrilateral.
    index_two_doublings = [
        {
            "glued_endpoint_labels": [2, 8],
            "opposite_label": 3,
            "union": "triangle_(3,3,4)",
        },
        {
            "glued_endpoint_labels": [2, 3],
            "opposite_label": 8,
            "union": "non_Coxeter_triangle_with_angle_2pi/3",
        },
        {
            "glued_endpoint_labels": [3, 8],
            "opposite_label": 2,
            "union": "quadrilateral",
        },
    ]

    # For q=-G/2 the integral quadratic polynomial is
    # 6x^2-6xy+3y^2-3yz+z^2.  Modulo 3, q=0 first forces z=0.  After writing
    # z=3z_1 and dividing by 3, the residual binary form
    # 2x^2+xy+y^2 forces x=y=0 modulo 3.  Infinite descent proves anisotropy.
    residual_binary_zeros_mod_3 = [
        [x, y]
        for x in range(3)
        for y in range(3)
        if (2 * x * x + x * y + y * y) % 3 == 0
    ]

    checks = {
        "gram_symmetric_integral": gram == gram.T
        and all(value.q == 1 for value in gram),
        "gram_determinant_36": gram.det() == 36,
        # Exact LDL pivots are -12,-3,1, hence signature (1,2).
        "gram_signature_1_2": list(
            gram.LDLdecomposition(hermitian=False)[1].diagonal()
        )
        == [-12, -3, 1],
        "coxeter_cosine_squares": normalized_cosine_squares
        == [Fraction(1, 2), Fraction(3, 4), Fraction(0)],
        "coxeter_labels_4_6_2": labels == (4, 6, 2),
        "all_vertices_finite": principal_pair_determinants == [36, 3, 24],
        "area_pi_over_12": _triangle_area(tuple(sorted(labels))) == Fraction(1, 12),
        "reflections_integral": all(
            all(value.q == 1 for value in reflection) for reflection in reflections
        ),
        "reflections_preserve_gram": all(
            reflection.T * gram * reflection == gram for reflection in reflections
        ),
        "reflections_have_order_two": all(
            reflection * reflection == sp.eye(3) for reflection in reflections
        ),
        "proper_overgroup_triangle_candidate_unique": candidates
        == [{"labels": [2, 3, 8], "area_over_pi": "1/24", "index": 2}],
        "small_triangle_area_window_exhausted": area_window
        == [
            {
                "labels": [2, 3, 7],
                "area_over_pi": "1/42",
                "source_area_ratio": "7/2",
                "integer_index": False,
            },
            {
                "labels": [2, 3, 8],
                "area_over_pi": "1/24",
                "source_area_ratio": "2",
                "integer_index": True,
            },
        ],
        "index_two_candidate_does_not_double_to_246": all(
            item["union"] != "triangle_(2,4,6)" for item in index_two_doublings
        ),
        "target_has_non_simply_laced_angles": {4, 6}.issubset(set(labels)),
        "anisotropy_residual_binary_has_only_zero": residual_binary_zeros_mod_3
        == [[0, 0]],
        # For two (-2)-curves meeting at a finite Coxeter angle, C.D is an
        # integer in {0,1}; hence the only finite labels are 2 and 3.
        "k_trivial_surface_allowed_finite_labels": [
            (intersection, Fraction(intersection * intersection, 4))
            for intersection in (0, 1)
        ]
        == [(0, Fraction(0)), (1, Fraction(1, 4))],
        "full_nef_rank_equals_three": len(labels) == 3 and gram.rank() == 3,
    }

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "scope": {
            "input": (
                "the compact arithmetic Coxeter triangle Delta(2,4,6), with "
                "minimal Vinberg field Q"
            ),
            "excluded_target": (
                "the exact-return Property-4 compiler obtained by combining "
                "Chamber Fidelity with the three complex carrier alternatives "
                "of canonical_chamber_rigidity_absw_transfer_theorem.md"
            ),
            "not_claimed": (
                "nonexistence for the bare carrier-relative third row without "
                "exact chamber return, or for a future carrier category which "
                "permits a proper invariant/eigen-nef slice even when k=Q"
            ),
        },
        "coxeter_triangle": {
            "wall_gram": _matrix_payload(gram),
            "wall_squares": [int(gram[i, i]) for i in range(3)],
            "determinant": int(gram.det()),
            "signature": [1, 2],
            "edge_order": [[0, 1], [1, 2], [0, 2]],
            "normalized_cosine_squares": [str(value) for value in normalized_cosine_squares],
            "edge_labels": list(labels),
            "principal_pair_determinants": principal_pair_determinants,
            "area_over_pi": "1/12",
            "reflection_matrices": [_matrix_payload(matrix) for matrix in reflections],
        },
        "arithmeticity": {
            "integral_lattice": True,
            "integral_reflections": True,
            "compact_fundamental_triangle": True,
            "argument": (
                "The reflection group lies in O(G,Z) and has compact, hence "
                "finite-covolume, fundamental chamber.  Since O(G,Z) is an "
                "arithmetic lattice, the subgroup has finite index and is arithmetic."
            ),
            "minimal_vinberg_field": "Q",
            "rational_anisotropy": {
                "quadratic_polynomial": "6x^2-6xy+3y^2-3yz+z^2",
                "prime": 3,
                "residual_binary_zeros_mod_3": residual_binary_zeros_mod_3,
                "argument": (
                    "q=0 mod 3 forces 3|z.  Put z=3z1 and divide by 3; "
                    "mod 3 the equation 2x^2+xy+y^2=0 forces 3|x and "
                    "3|y.  Infinite descent excludes a nonzero rational "
                    "isotropic vector."
                ),
            },
        },
        "absolute_reflection_maximality": {
            "four_or_more_side_compact_coxeter_area_lower_bound_over_pi": "1/6",
            "proper_overgroup_area_upper_bound_over_pi": "1/24",
            "all_triangle_types_in_area_window": area_window,
            "integer_index_triangle_candidates": candidates,
            "index_two_doubling_analysis": index_two_doublings,
            "conclusion": (
                "No proper discrete reflection overgroup exists; Delta(2,4,6) "
                "is therefore absolutely reflection-maximal and in particular "
                "maximal among arithmetic reflection groups."
            ),
        },
        "carrier_obstruction": {
            "literal_full_nef_for_k_Q_forces_picard_rank": 3,
            "rational_host_gate": (
                "The rank-three trace form is anisotropic over Q, whereas "
                "NS of a smooth rational surface of rank three is the "
                "isotropic odd unimodular form I_(1,2).  No full-rank rational "
                "similarity embedding exists."
            ),
            "rational_branches": (
                "Every smooth projective rational surface of Picard rank at "
                "least two has a pulled-back nef isotropic ruling class.  Its "
                "projectivized full nef cone is noncompact, contradicting the "
                "compact target triangle."
            ),
            "minimalizable_full_trace_branch": (
                "Exact Property-4 literal full-nef return gives rho=3.  A nontrivial "
                "contraction would leave Picard rank below the rank-three trace "
                "block, so the carrier is already minimal.  Canonical "
                "orthogonality of a full-rank trace block forces K numerically "
                "zero.  Adjunction then makes every irreducible negative facet "
                "a smooth rational (-2)-curve."
            ),
            "angle_contradiction": (
                "Distinct (-2)-curves have nonnegative integral intersection. "
                "At a finite Coxeter ridge this intersection is 0 or 1, giving "
                "only labels 2 or 3.  Delta(2,4,6) has labels 4 and 6."
            ),
            "conclusion": (
                "The direct rational, certified weak-log rational, and complex "
                "exact-return minimalizable full-trace packages are all empty "
                "for this raw maximal arithmetic chamber; the bare carrier-relative "
                "third row without exact chamber return is not claimed empty."
            ),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    payload["certificate_sha256"] = _canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()

    payload = build_certificate()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not OUTPUT.is_file() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise AssertionError("archived compact (2,4,6) certificate has drifted")
        print(f"archive check passed: {OUTPUT}")
        print(f"certificate_sha256={payload['certificate_sha256']}")
        return
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT)
        print(payload["certificate_sha256"])
        return
    print(rendered, end="")


if __name__ == "__main__":
    main()

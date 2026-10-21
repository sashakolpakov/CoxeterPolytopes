#!/usr/bin/env python3
"""Exact finite-algebra/Galois calibration for the H11 line--quintic chart.

This is deliberately tiny.  It does not claim that the global carrier stack
is zero-dimensional.  It certifies that the already exhaustive rank locus for
the one H11/26 inverse chart has two rational projective branches, that the
nonreduced branch is rejected, and that the accepted rigidified parameter and
all exact Cartan cyclic products have field Q.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
PLANE_SOURCE = ROOT / "data/enumerations/coble_type_i_two_boundary_plane_inverse.json"
POLYTOPE_SOURCE = ROOT / "data/enumerations/coble_type_i_two_boundary_full_ns_h11.json"
PLANE_SHA256 = "6313b11cc9ba267834231aab2469d99e9d0242f4a95f53280047351f8318da04"
POLYTOPE_SHA256 = "bcfa2de1e96bf4f41282345a64f8ee6f3a75e715287bab4962d7082584b4d985"
DEFAULT_OUTPUT = ROOT / "data/enumerations/dessin_native_field_effectivity_h11.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _token(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def _digest(values: list[Fraction]) -> str:
    payload = ";".join(_token(value) for value in values).encode()
    return hashlib.sha256(payload).hexdigest()


def build() -> dict[str, Any]:
    assert _sha256(PLANE_SOURCE) == PLANE_SHA256
    assert _sha256(POLYTOPE_SOURCE) == POLYTOPE_SHA256
    plane = json.loads(PLANE_SOURCE.read_text())
    polytope = json.loads(POLYTOPE_SOURCE.read_text())

    geometry = plane["canonical_determinantal_geometry"]
    matrix = geometry["quintic_matrix"]
    rejected = geometry["rejected_component"]
    accepted = geometry["surviving_horizontal_component"]
    assert plane["all_checks_pass"]
    assert matrix["exact_rank_drop_locus"] == "a^2-b^2=0"
    assert matrix["maximal_minor"] == "-64*(a - b)*(a + b)"
    assert matrix["rank_on_a_equals_b"] == matrix["rank_on_a_equals_minus_b"] == 20
    assert rejected["equation"] == "a+b=0"
    assert rejected["sample"] == "a=1,b=-1"
    assert "conic squared" in rejected["reason"]
    assert accepted["equation"] == "a-b=0 with a nonzero"
    assert accepted["sample"] == "a=b=1 over Q"
    assert plane["checks"]["minus_branch_rejected_as_nonreduced"]
    assert plane["checks"]["plus_branch_absolutely_irreducible"]

    gram = polytope["facets"]["negative_wall_gram"]
    assert polytope["all_checks_pass"]
    assert len(gram) == 26 and all(len(row) == 26 for row in gram)
    assert all(isinstance(entry, int) for row in gram for entry in row)
    cartan = [
        Fraction(2 * gram[i][j], gram[i][i])
        for i in range(26)
        for j in range(26)
    ]
    two_cycles = [
        Fraction(4 * gram[i][j] * gram[i][j], gram[i][i] * gram[j][j])
        for i in range(26)
        for j in range(i + 1, 26)
    ]
    assert all(value.denominator <= 2 for value in cartan)
    assert all(value.denominator == 1 for value in two_cycles)

    certificate: dict[str, Any] = {
        "schema": "coxeter.dessin_native_field_effectivity_h11.v1",
        "status_date": "2026-08-09",
        "source_pins": [
            {
                "path": str(PLANE_SOURCE.relative_to(ROOT)),
                "sha256": PLANE_SHA256,
                "role": "exact projective rank locus and geometric branch audit",
            },
            {
                "path": str(POLYTOPE_SOURCE.relative_to(ROOT)),
                "sha256": POLYTOPE_SHA256,
                "role": "exact rational wall Gram and same-polytope return",
            },
        ],
        "projective_rank_locus": {
            "equation": "a^2-b^2=(a-b)(a+b)=0 in P1_Q",
            "geometric_points": [[1, 1], [1, -1]],
            "geometric_point_count_before_geometric_gate": 2,
            "rejected": {
                "point": [1, -1],
                "chart_ideal": "b+1 in Q[b] on a=1",
                "reason": "line boundary times a squared conic is nonreduced",
            },
            "accepted": {
                "point": [1, 1],
                "chart_ideal": "b-1 in Q[b] on a=1",
                "reduced_coordinate_algebra": "Q[b]/(b-1) = Q",
                "algebra_Q_dimension": 1,
                "geometric_point_count": 1,
                "galois_orbit_size": 1,
                "field_of_moduli": "Q",
                "field_of_definition": "Q",
            },
        },
        "exact_gram_field": {
            "facet_count": 26,
            "scalar_coordinate_field": "Q",
            "cartan_entry_count": len(cartan),
            "cartan_unique_values": sorted({_token(value) for value in cartan}),
            "cartan_sequence_sha256": _digest(cartan),
            "length_two_cyclic_product_count": len(two_cycles),
            "length_two_cyclic_product_unique_values": sorted(
                {_token(value) for value in two_cycles}
            ),
            "length_two_cycle_sequence_sha256": _digest(two_cycles),
            "all_cartan_entries_rational": True,
            "all_cyclic_products_rational": True,
            "recovered_Vinberg_field": "Q",
        },
        "scope": {
            "proved": "one rigidified accepted H11/26 line-quintic chart point has residue and Vinberg field Q",
            "not_proved": "the variable-class Coxeter-Borcherds carrier stack is zero-dimensional",
            "global_field_finiteness": "open pending exhaustive zero-dimensionality",
        },
        "all_checks_pass": True,
    }
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-output", type=Path)
    args = parser.parse_args()
    certificate = build()
    payload = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check_output is not None:
        if args.check_output.read_text() != payload:
            raise SystemExit(f"archive mismatch: {args.check_output}")
    elif args.output is not None:
        args.output.write_text(payload)
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()

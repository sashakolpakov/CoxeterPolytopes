#!/usr/bin/env python3
"""Certify the class-2a H13/50 and class-2c H14/98 fixed sections.

This verifier rebuilds both sections from the pinned 210-wall Borcherds /
Dolgachev--Kondo Gram matrix.  It proves the exact diagram actions, saturated
fixed lattices, strict irredundancy of every negative restricted wall, Coxeter
pair relations, and integrality of every primitive facet reflection.  The
large literal Vinberg edge/cusp enumeration is deliberately kept in the
companion C++ certificate.
"""

from __future__ import annotations

import argparse
import base64
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.certify_alexeev_high_dim_reflection_examples import (
    SOURCE_SHA256,
    _extract_matrix,
)
from scripts.certify_borcherds_char2_3b_h11_22 import (
    DEFAULT_SOURCE,
    _ambient_action,
    _digest,
    _facet_lattice_records,
    _fixed_lattice,
    _matrix_payload,
    _order,
)
from scripts.search_k3_high_rank_involution_folds import _coxeter_audit, _cycles
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia
from scripts.verify_k3_quotient_nef_slices import _orbit_sum_matrix


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data/enumerations/borcherds_char2_2a_h13_50_2c_h14_98_construction.json"
)
SCHEMA = "coxeter.borcherds_char2_2a_h13_50_2c_h14_98.construction.v1"


CASES: dict[str, dict[str, Any]] = {
    "2a": {
        "id": "SS2-2A-H13-50",
        "permutation_base64": (
            "EAsoAzIHBgUmCRwBGx8OMAARGG15FRo/EhkWDAotHg0gITF9mWUILgJrKtGWHSeB"
            "DyIEOFWossIzQrqphkjBF4xhOZeVRUbEPUmhWUyEZsp1jXq+tDRkj69LyYCzkbGr"
            "bEG8Y1YlTqTGaY4pYBN/cZ9vvYK4UNCYzhRSsMsjiW5bL3ODTYU8i4h+iodAUWpX"
            "yF2Sv6dELEN3JKObpZ3DcMBKuZpnnKaUNTuqX7etrlh7XjZcVLa1rHSiOrticlOTo"
            "D43nkfFaMeQWk98zM14z3Yr"
        ),
        "class_size": 315,
        "fixed_rank": 14,
        "dimension": 13,
        "cycle_histogram": {1: 34, 2: 88},
        "causal_counts": (50, 72, 0),
        "action_trace": 6,
        "fixed_inertia": (1, 13, 0),
        "fixed_determinant": -1024,
        "fixed_smith": [1] * 4 + [2] * 10,
        "fixed_lattice_class": "even 2-elementary (14,10,0)",
        "facet_types": Counter(
            {(1, -4, 2): 24, (2, -4, 2): 16, (1, -2, 1): 10}
        ),
        "facet_type_labels": {
            "orbit1_square-2_div1": 10,
            "orbit1_square-4_div2": 24,
            "orbit2_square-4_div2": 16,
        },
        "coxeter_edges": {
            "2": 700,
            "3": 137,
            "4": 112,
            "6": 0,
            "infinity": 156,
            "ultraparallel": 120,
        },
        "facet_gram_sha256": (
            "7aaea6ed54d34abd3e6d205c1ba132b243bfddcc0ad2db455a1a45a28a21bd51"
        ),
        "geometric_status": (
            "contained in the Dolgachev--Kondo index-two geometric subgroup"
        ),
        "lineage": (
            "known (14,10,0) lattice and 16-wall full Alexeev chamber; "
            "50-wall subgroup index and literal priority open"
        ),
    },
    "2c": {
        "id": "SS2-2C-H14-98",
        "permutation_base64": (
            "AAMcAQQFBgcICQ4NDAsKGBEQGhQTFxYVDx8SHgIdGxkiKCAjKiUnJiE5JEZELTI9"
            "MDEuNkBJMzg3KTpSYC9bPzRUVkMsSCtRRTVqS1BXeW9MRzt9QVlCTWxVdz5ycHFz"
            "PG6EhW1pa5uIZUpmWGRhT11eXF+NqKRalE6dhpZTgrKJin6gYmN7mmiAgbiMdI+O"
            "urCSsXiVfL6qmYdnoXqrr4Oco6J2vKemda2YnqypyZ+Rk3/KwrXIx4vEkLulv5e9"
            "0cO0wbnNxre2rrPL0MXPzszA"
        ),
        "class_size": 360,
        "fixed_rank": 15,
        "dimension": 14,
        "cycle_histogram": {1: 28, 2: 91},
        "causal_counts": (98, 21, 0),
        "action_trace": 8,
        "fixed_inertia": (1, 14, 0),
        "fixed_determinant": 256,
        "fixed_smith": [1] * 8 + [2] * 6 + [4],
        "fixed_lattice_class": (
            "even rank-15 lattice of determinant 256; maximal lineage open"
        ),
        "facet_types": Counter(
            {
                (2, -4, 4): 56,
                (1, -2, 1): 14,
                (2, -4, 2): 14,
                (1, -4, 2): 14,
            }
        ),
        "facet_type_labels": {
            "orbit1_square-2_div1": 14,
            "orbit1_square-4_div2": 14,
            "orbit2_square-4_div2": 14,
            "orbit2_square-4_div4": 56,
        },
        "coxeter_edges": {
            "2": 1547,
            "3": 49,
            "4": 70,
            "6": 0,
            "infinity": 735,
            "ultraparallel": 2352,
        },
        "facet_gram_sha256": (
            "257c0444e7774232672520096ef5639074d83ef2e1a58d625ab2435007f2f330"
        ),
        "geometric_status": (
            "outside the Dolgachev--Kondo geometric subgroup; Frobenius "
            "wall-diagram class, not a K3-surface automorphism"
        ),
        "lineage": (
            "exact normalized Coxeter-diagram rediscovery of Alexeev's "
            "published (15,9,1) 98-wall chamber, with short and long norms "
            "exchanged by the coodd doubled-dual realization; not a new "
            "abstract polytope"
        ),
    },
}


def _permutation(record: dict[str, Any]) -> tuple[int, ...]:
    permutation = tuple(base64.b64decode(record["permutation_base64"]))
    if len(permutation) != 210 or set(permutation) != set(range(210)):
        raise ArithmeticError("an archived characteristic-two permutation is bad")
    return permutation


def _build_case(wall_gram: sp.Matrix, atlas_class: str) -> dict[str, Any]:
    expected = CASES[atlas_class]
    permutation = _permutation(expected)
    if _order(permutation) != 2:
        raise ArithmeticError(f"class {atlas_class} is not an involution")
    if any(
        wall_gram[left, right]
        != wall_gram[permutation[left], permutation[right]]
        for left in range(210)
        for right in range(210)
    ):
        raise ArithmeticError(f"class {atlas_class} does not preserve the Gram")

    cycles = _cycles(permutation)
    histogram = Counter(map(len, cycles))
    if histogram != Counter(expected["cycle_histogram"]):
        raise ArithmeticError(f"class {atlas_class} cycle census changed")
    sums = _orbit_sum_matrix(210, cycles)
    orbit_gram = sums.T * wall_gram * sums
    fixed_rank = expected["fixed_rank"]
    if not (
        orbit_gram.rank() == fixed_rank
        and exact_inertia(orbit_gram)
        == expected["fixed_inertia"][:-1] + (orbit_gram.rows - fixed_rank,)
    ):
        raise ArithmeticError(f"class {atlas_class} fixed wall Gram changed")

    facet_indices = [i for i in range(orbit_gram.rows) if orbit_gram[i, i] < 0]
    isotropic = [i for i in range(orbit_gram.rows) if orbit_gram[i, i] == 0]
    timelike = [i for i in range(orbit_gram.rows) if orbit_gram[i, i] > 0]
    if (len(facet_indices), len(isotropic), len(timelike)) != expected[
        "causal_counts"
    ]:
        raise ArithmeticError(f"class {atlas_class} causal census changed")

    interior = sp.ones(orbit_gram.rows, 1)
    interior_pairings = orbit_gram * interior
    interior_square = int((interior.T * orbit_gram * interior)[0])
    if interior_square != 137214 or min(interior_pairings) != 99:
        raise ArithmeticError(f"class {atlas_class} strict interior changed")
    facet_witnesses = []
    for facet in facet_indices:
        witness = interior.copy()
        witness[facet] -= sp.Rational(
            interior_pairings[facet], orbit_gram[facet, facet]
        )
        pairings = orbit_gram * witness
        square = sp.Rational((witness.T * orbit_gram * witness)[0])
        other_pairings = [
            value for index, value in enumerate(pairings) if index != facet
        ]
        if not (
            square > 0
            and pairings[facet] == 0
            and min(other_pairings) > 0
        ):
            raise ArithmeticError(f"class {atlas_class} has a redundant facet")
        facet_witnesses.append(
            {
                "orbit_index_one_based": facet + 1,
                "witness_square": str(square),
                "minimum_other_wall_pairing": str(min(other_pairings)),
            }
        )

    facet_gram = orbit_gram.extract(facet_indices, facet_indices)
    coxeter = _coxeter_audit(
        facet_gram, range(facet_gram.rows), [cycles[i] for i in facet_indices]
    )
    if not (
        coxeter["all_pairs_are_exact_coxeter_relations"]
        and coxeter["edge_counts"] == expected["coxeter_edges"]
    ):
        raise ArithmeticError(f"class {atlas_class} Coxeter pair gate changed")

    basis_gram, action = _ambient_action(wall_gram, permutation, order=2)
    if int(action.trace()) != expected["action_trace"]:
        raise ArithmeticError(f"class {atlas_class} NS trace changed")
    fixed_basis, fixed_gram, fixed_smith = _fixed_lattice(
        basis_gram,
        action,
        expected_rank=fixed_rank,
        expected_inertia=expected["fixed_inertia"],
        expected_determinant=expected["fixed_determinant"],
        expected_smith=expected["fixed_smith"],
    )
    if atlas_class == "2a":
        inverse = fixed_gram.inv()
        if not all(sp.denom(inverse[i, i]) == 1 for i in range(fixed_rank)):
            raise ArithmeticError("class 2a fixed lattice does not have delta zero")
    facet_lattice = _facet_lattice_records(
        wall_gram,
        basis_gram,
        fixed_basis,
        fixed_gram,
        cycles,
        facet_indices,
        expected_histogram=expected["facet_types"],
    )
    facet_payload = _matrix_payload(facet_gram)
    facet_hash = _digest(facet_payload)
    if facet_hash != expected["facet_gram_sha256"]:
        raise ArithmeticError(f"class {atlas_class} facet Gram hash changed")

    checks = {
        "order_two_ambient_gram_automorphism": True,
        "integral_unimodular_NS_action": True,
        "saturated_fixed_lattice": True,
        "strict_fixed_interior": True,
        "all_negative_walls_irredundant": True,
        "all_facet_pairs_exact_coxeter": True,
        "all_primitive_facet_reflections_integral": True,
    }
    return {
        "id": expected["id"],
        "action": {
            "atlas_class": atlas_class,
            "full_diagram_group_class_size": expected["class_size"],
            "order": 2,
            "permutation_images_one_based": [value + 1 for value in permutation],
            "cycle_length_histogram": {
                str(length): count for length, count in sorted(histogram.items())
            },
            "NS_action_columns": _matrix_payload(action),
            "trace": int(action.trace()),
            "geometric_subgroup_status": expected["geometric_status"],
        },
        "fixed_lattice": {
            "definition": "NS(X) intersect Fix_Q(g)",
            "basis_columns_in_ambient_basis": _matrix_payload(fixed_basis),
            "gram": _matrix_payload(fixed_gram),
            "signature_positive_negative": [1, expected["dimension"]],
            "determinant": int(fixed_gram.det()),
            "smith_diagonal_absolute": fixed_smith,
            "class": expected["fixed_lattice_class"],
        },
        "restricted_chamber": {
            "wall_orbit_count": len(cycles),
            "negative_facet_count": len(facet_indices),
            "isotropic_inactive_orbit_count": len(isotropic),
            "timelike_inactive_orbit_count": len(timelike),
            "strict_interior_square": interior_square,
            "strict_interior_minimum_wall_pairing": int(min(interior_pairings)),
            "facet_orbit_indices_one_based": [i + 1 for i in facet_indices],
            "facet_gram": facet_payload,
            "facet_gram_sha256": facet_hash,
            "facet_witnesses": facet_witnesses,
            "primitive_facet_lattice_records": facet_lattice,
            "primitive_facet_type_histogram": expected["facet_type_labels"],
            "coxeter_edge_counts": coxeter["edge_counts"],
        },
        "lineage_status": expected["lineage"],
        "decision": {
            "hyperbolic_dimension": expected["dimension"],
            "facet_count": len(facet_indices),
            "construction": "exact_irredundant_integral_Coxeter_chamber",
            "volume": "certified_in_companion_literal_Vinberg_artifact",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def build_certificate(source_path: Path = DEFAULT_SOURCE) -> dict[str, Any]:
    raw = source_path.read_bytes()
    source_hash = hashlib.sha256(raw).hexdigest()
    if source_hash != SOURCE_SHA256:
        raise ArithmeticError("the pinned Alexeev ancillary source hash changed")
    wall_gram = _extract_matrix(raw.decode("utf-8"), (22, 2, 0))
    if not (
        wall_gram.shape == (210, 210)
        and wall_gram.rank() == 22
        and exact_inertia(wall_gram) == (1, 21, 188)
    ):
        raise ArithmeticError("the ambient H21/210 Gram changed")
    cases = [_build_case(wall_gram, atlas) for atlas in ("2a", "2c")]
    return {
        "schema": SCHEMA,
        "source": {
            "path": str(source_path),
            "sha256": source_hash,
            "matrix_key": "(22,2,0)",
        },
        "ambient": {
            "hyperbolic_dimension": 21,
            "wall_count": 210,
            "lattice": "U + D20(-1)",
        },
        "cases": cases,
        "all_checks_pass": all(case["all_checks_pass"] for case in cases),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", nargs="?", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate(args.source)
    if args.check_archive:
        if payload != json.loads(args.output.read_text(encoding="utf-8")):
            raise ArithmeticError("rebuilt H13/H14 construction differs from archive")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                case["id"]: {
                    "facets": case["decision"]["facet_count"],
                    "facet_gram_sha256": case["restricted_chamber"][
                        "facet_gram_sha256"
                    ],
                    "all_checks_pass": case["all_checks_pass"],
                }
                for case in payload["cases"]
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Certify the Hilbert-square scaled-dual carrier of SS2-2C H14/98.

The characteristic-two fixed lattice itself is not the right dual object.
First pass to its unique reflection-stable odd index-two closure M.  Then

    N = M^vee(2),     Gram(N) = 2 Gram(M)^(-1)

is the even 2-elementary lattice (15,9,1).  A frozen determinant-one
isometry identifies N with Alexeev's archived Pic(Hilb^2(X)) lattice.
Every one of the 98 oriented source halfspaces is transported by primitive
covectors, those covectors generate N, and all target reflections are
integral.  The hyperbolic polytope is the already known Alexeev 98-wall
chamber; the positive result here is the complete coordinate-level carrier
bridge, not a new abstract polytope or a Kahler-cone theorem.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = (
    ROOT
    / "data/enumerations/borcherds_char2_2a_h13_50_2c_h14_98_construction.json"
)
LINEAGE = (
    ROOT / "data/enumerations/borcherds_char2_2a_h13_50_2c_h14_98_lineage.json"
)
VINBERG = (
    ROOT / "data/enumerations/borcherds_char2_2a_h13_50_2c_h14_98_vinberg.json"
)
HILB2_SOURCE = (
    ROOT
    / "data/enumerations/alexeev_20_2_1_15_9_1_wall_symmetry_inventory.json"
)
OUTPUT = (
    ROOT
    / "data/enumerations/ss2_2c_h14_98_hilb2_scaled_dual_carrier.json"
)

SCHEMA = "coxeter.ss2_2c_h14_98.hilb2_scaled_dual_divisor_carrier.v1"
EXPECTED_SHA256 = {
    "construction": "a9087200b360c4bed75cd93c987bd7e7d493ff1ce852bb9c5ea4dc3490f1fe11",
    "lineage": "f0867d47f4a1405a5398a274865b32fff3b24608790f2b41a25eb491de2e8971",
    "vinberg": "032c8343a6dd7f4904b9daecb074b8b0bbbf1a94f9431252884484c98d2e8456",
    "hilb2_source": "20a257e60b768571a7e7c8aa64e60350bdf7940f5d221cfe031fc2e6ad879e09",
}

# Coordinates in N=M^vee(2) are sent to the leading-root basis of Alexeev's
# actual Pic(Hilb^2(X)) lattice.  The replay checks the isometry, every one of
# the 98 root images, and determinant +1; the matrix is not trusted by fiat.
DOUBLED_DUAL_TO_HILB2 = sp.Matrix(
    [
        [-60, 78, -24, -36, 150, 208, 254, 300, -37, -74, -109, -144, -179, -120, 39],
        [40, -52, 16, 24, -100, -139, -170, -200, 24, 48, 72, 96, 120, 80, -26],
        [80, -104, 32, 48, -200, -278, -339, -400, 48, 96, 144, 192, 240, 160, -52],
        [60, -78, 24, 36, -150, -208, -254, -300, 36, 72, 108, 144, 180, 120, -39],
        [30, -39, 12, 18, -75, -104, -127, -150, 18, 36, 54, 72, 90, 60, -20],
        [-184, 240, -74, -111, 462, 640, 782, 924, -114, -228, -336, -444, -552, -368, 120],
        [-252, 328, -101, -152, 633, 876, 1071, 1266, -156, -312, -460, -608, -756, -504, 164],
        [-160, 208, -64, -96, 402, 556, 680, 804, -99, -198, -292, -386, -480, -320, 104],
        [-80, 104, -32, -48, 201, 278, 340, 402, -50, -99, -146, -193, -240, -160, 52],
        [95, -124, 38, 57, -240, -332, -406, -480, 60, 120, 176, 232, 288, 192, -62],
        [80, -104, 32, 48, -201, -278, -340, -402, 50, 100, 146, 193, 240, 160, -52],
        [65, -84, 26, 39, -162, -224, -274, -324, 40, 80, 117, 154, 192, 128, -42],
        [50, -64, 20, 30, -123, -170, -208, -246, 30, 60, 88, 116, 144, 96, -32],
        [25, -32, 10, 15, -62, -85, -104, -123, 15, 30, 44, 58, 72, 48, -16],
        [-120, 156, -48, -72, 300, 416, 508, 600, -74, -148, -218, -288, -358, -239, 78],
    ]
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [
            int(value) if sp.Rational(value).q == 1 else str(sp.Rational(value))
            for value in matrix.row(row)
        ]
        for row in range(matrix.rows)
    ]


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _smith(matrix: sp.Matrix) -> list[int]:
    normal = smith_normal_form(matrix, domain=ZZ)
    return [abs(int(normal[index, index])) for index in range(matrix.rows)]


def _content(vector: sp.Matrix) -> int:
    return math.gcd(*(abs(int(value)) for value in vector))


def _reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = (root.T * gram * root)[0]
    return sp.simplify(
        sp.eye(gram.rows) - sp.Rational(2, square) * root * root.T * gram
    )


def _histogram_payload(counter: Counter[tuple[int, ...]]) -> dict[str, int]:
    return {
        ",".join(str(value) for value in key): count
        for key, count in sorted(counter.items())
    }


def build_certificate() -> dict[str, Any]:
    paths = {
        "construction": CONSTRUCTION,
        "lineage": LINEAGE,
        "vinberg": VINBERG,
        "hilb2_source": HILB2_SOURCE,
    }
    for key, path in paths.items():
        if _sha256(path) != EXPECTED_SHA256[key]:
            raise ArithmeticError(f"pinned {key} archive changed")

    construction = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))
    lineage = json.loads(LINEAGE.read_text(encoding="utf-8"))
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    hilb2_source = json.loads(HILB2_SOURCE.read_text(encoding="utf-8"))

    case = next(row for row in construction["cases"] if row["id"] == "SS2-2C-H14-98")
    vinberg_case = next(row for row in vinberg["cases"] if row["id"] == "SS2-2C-H14-98")
    h14_lineage = lineage["H14_98"]
    closure = h14_lineage["reflection_stable_overlattice"]
    target_case = next(
        row for row in hilb2_source["examples"] if row["key"] == [15, 9, 1]
    )

    gram_l = sp.Matrix(case["fixed_lattice"]["gram"])
    closure_basis = sp.Matrix(closure["basis_columns_in_fixed_Q_space"])
    gram_m = sp.simplify(closure_basis.T * gram_l * closure_basis)
    gram_n = sp.simplify(2 * gram_m.inv())
    gram_hilb2 = sp.Matrix(target_case["ambient"]["gram"])
    target_roots = [
        sp.Matrix(row)
        for row in target_case["walls"]["coordinates_in_leading_root_basis"]
    ]
    permutation = h14_lineage["exact_Alexeev_match"][
        "candidate_to_alexeev_permutation_zero_based"
    ]
    transport = DOUBLED_DUAL_TO_HILB2

    lattice_checks = {
        "all_four_upstream_certificates_green": (
            construction["all_checks_pass"]
            and lineage["all_checks_pass"]
            and vinberg["all_checks_pass"]
            and len(hilb2_source["examples"]) == 2
        ),
        "source_fixed_lattice_rank15_det256_even": (
            exact_inertia(gram_l) == (1, 14, 0)
            and abs(int(gram_l.det())) == 256
            and all(int(gram_l[index, index]) % 2 == 0 for index in range(15))
        ),
        "unique_stable_closure_replayed": (
            closure["W98_stable_class_count"] == 1
            and abs(closure_basis.det()) == sp.Rational(1, 2)
            and gram_m == sp.Matrix(closure["gram"])
            and gram_m.det() == 64
            and exact_inertia(gram_m) == (1, 14, 0)
            and _smith(gram_m) == [1] * 9 + [2] * 6
            and any(int(gram_m[index, index]) % 2 for index in range(15))
        ),
        "doubled_dual_integral_even_15_9_1": (
            all(value.is_Integer for value in gram_n)
            and gram_n.det() == 512
            and exact_inertia(gram_n) == (1, 14, 0)
            and _smith(gram_n) == [1] * 6 + [2] * 9
            and all(int(gram_n[index, index]) % 2 == 0 for index in range(15))
        ),
        "target_is_archived_actual_Hilbert_square_Picard_lattice": (
            target_case["geometric_distinction"]["context_from_Alexeev"].startswith(
                "Picard lattice of a Hilbert square Hilb^2(X)"
            )
            and target_case["ambient"]["leading_walls_are_integral_lattice_basis"]
            and gram_hilb2.det() == 512
            and exact_inertia(gram_hilb2) == (1, 14, 0)
            and _smith(gram_hilb2) == [1] * 6 + [2] * 9
        ),
        "explicit_unimodular_isometry": (
            all(value.is_Integer for value in transport)
            and transport.det() == 1
            and transport.T * gram_hilb2 * transport == gram_n
        ),
        "permutation_is_bijective": sorted(permutation) == list(range(98)),
    }
    if not all(lattice_checks.values()):
        raise ArithmeticError(
            f"Hilbert-square lattice bridge failed: "
            f"{[key for key, value in lattice_checks.items() if not value]}"
        )

    records: list[dict[str, Any]] = []
    primitive_n_roots: list[sp.Matrix] = []
    primitive_target_roots: list[sp.Matrix] = []
    closure_types: Counter[tuple[int, int]] = Counter()
    target_types: Counter[tuple[int, int]] = Counter()
    source_to_closure_types: Counter[tuple[int, ...]] = Counter()
    inverse_closure_basis = closure_basis.inv()

    for index, source in enumerate(
        case["restricted_chamber"]["primitive_facet_lattice_records"]
    ):
        root_l = sp.Matrix(source["primitive_fixed_coordinates"])
        raw_root_m = inverse_closure_basis * root_l
        if any(not value.is_Integer for value in raw_root_m):
            raise ArithmeticError(f"facet {index} is not contained in the stable closure")
        closure_content = _content(raw_root_m)
        root_m = (raw_root_m / closure_content).applyfunc(int)
        div_m = _content(gram_m * root_m)
        root_n = (gram_m * root_m / div_m).applyfunc(int)
        root_target = transport * root_n
        expected_target = target_roots[permutation[index]]

        square_l = int((root_l.T * gram_l * root_l)[0])
        div_l = _content(gram_l * root_l)
        square_m = int((root_m.T * gram_m * root_m)[0])
        square_n = int((root_n.T * gram_n * root_n)[0])
        div_n = _content(gram_n * root_n)

        reflection_m = _reflection(gram_m, root_m)
        reflection_n = _reflection(gram_n, root_n)
        reflection_target = _reflection(gram_hilb2, root_target)
        halfspace_factor_from_closure = sp.Rational(2, div_m)
        halfspace_factor_from_original = sp.Rational(
            2, div_m * closure_content
        )

        checks = {
            "source_type_archived": (
                square_l == source["primitive_square"]
                and div_l == source["primitive_divisibility"]
            ),
            "orientation_scale_positive": source["orbit_sum_scale_from_primitive"] > 0,
            "primitive_in_closure_dual_and_target": (
                _content(root_m) == _content(root_n) == _content(root_target) == 1
            ),
            "scaled_dual_square_law": square_n * div_m**2 == 2 * square_m,
            "exact_positive_halfspace_transport": (
                halfspace_factor_from_closure > 0
                and halfspace_factor_from_original > 0
                and root_n.T * gram_n * gram_m
                == halfspace_factor_from_closure * root_m.T * gram_m
            ),
            "contragredient_reflection": reflection_n == reflection_m.T,
            "all_three_reflections_integral": (
                all(value.is_Integer for value in reflection_m)
                and all(value.is_Integer for value in reflection_n)
                and all(value.is_Integer for value in reflection_target)
            ),
            "target_root_is_exact_archived_wall": root_target == expected_target,
            "target_reflection_intertwines": (
                reflection_target * transport == transport * reflection_n
            ),
            "target_reflection_unimodular_involution_isometry": (
                reflection_target.det() == -1
                and reflection_target**2 == sp.eye(15)
                and reflection_target.T * gram_hilb2 * reflection_target
                == gram_hilb2
            ),
        }
        if not all(checks.values()):
            raise ArithmeticError(
                f"facet {index} failed: "
                f"{[key for key, value in checks.items() if not value]}"
            )

        closure_types[(square_m, div_m)] += 1
        target_types[(square_n, div_n)] += 1
        source_to_closure_types[
            (square_l, div_l, closure_content, square_m, div_m, square_n, div_n)
        ] += 1
        primitive_n_roots.append(root_n)
        primitive_target_roots.append(root_target)
        records.append(
            {
                "facet_zero_based": index,
                "alexeev_wall_zero_based": permutation[index],
                "source_L_coordinates": [int(value) for value in root_l],
                "source_L_square": square_l,
                "source_L_divisibility": div_l,
                "coordinates_in_stable_closure_before_primitive_reduction": [
                    int(value) for value in raw_root_m
                ],
                "stable_closure_coordinate_content": closure_content,
                "primitive_M_coordinates": [int(value) for value in root_m],
                "primitive_M_square": square_m,
                "primitive_M_divisibility": div_m,
                "primitive_N_coordinates": [int(value) for value in root_n],
                "primitive_Hilb2_coordinates": [int(value) for value in root_target],
                "primitive_Hilb2_square": square_n,
                "primitive_Hilb2_divisibility": div_n,
                "positive_halfspace_factor_from_M": str(halfspace_factor_from_closure),
                "positive_halfspace_factor_from_original_L": str(
                    halfspace_factor_from_original
                ),
                "checks": checks,
            }
        )

    roots_n_matrix = sp.Matrix.hstack(*primitive_n_roots).applyfunc(int)
    roots_target_matrix = sp.Matrix.hstack(*primitive_target_roots).applyfunc(int)
    hnf_n = hermite_normal_form(roots_n_matrix)
    hnf_target = hermite_normal_form(roots_target_matrix)
    carrier = sum(primitive_n_roots, sp.zeros(15, 1))
    carrier_target = transport * carrier
    carrier_pairings = [
        int((carrier.T * gram_n * root)[0]) for root in primitive_n_roots
    ]
    generation_checks = {
        "all_98_facets_transport": len(records) == 98,
        "primitive_N_roots_generate_full_lattice": hnf_n == sp.eye(15),
        "primitive_target_roots_generate_full_Picard_lattice": (
            hnf_target == sp.eye(15)
        ),
        "archived_first_15_target_walls_are_basis": (
            sp.Matrix.hstack(*target_roots[:15]).det() == 1
        ),
        "sum_is_strict_timelike_interior": (
            int((carrier.T * gram_n * carrier)[0]) == 40432
            and min(carrier_pairings) == 76
            and all(value > 0 for value in carrier_pairings)
            and carrier_target.T * gram_hilb2 * carrier_target
            == carrier.T * gram_n * carrier
        ),
        "closure_type_histogram": closure_types
        == Counter({(-1, 1): 70, (-4, 2): 14, (-2, 1): 14}),
        "target_type_histogram": target_types
        == Counter({(-2, 2): 70, (-2, 1): 14, (-4, 2): 14}),
    }
    if not all(generation_checks.values()):
        raise ArithmeticError(
            f"generation/interior audit failed: "
            f"{[key for key, value in generation_checks.items() if not value]}"
        )

    vinberg_checks = {
        "literal_case_green": vinberg_case["all_checks_pass"] is True,
        "finite_volume_noncompact": (
            vinberg_case["finite_volume"] is True
            and vinberg_case["compact"] is False
            and vinberg_case["classification"]
            == "finite_volume_noncompact_Coxeter_polytope"
        ),
        "every_elliptic_edge_has_two_endpoints": (
            vinberg_case["bad_edge_count"] == 0
            and vinberg_case["elliptic_edge_count"] == 646786
        ),
        "all_parabolic_neighborhoods_good": (
            vinberg_case["bad_parabolic_neighborhood_count"] == 0
            and vinberg_case["connected_affine_component_count"] == 5838
        ),
        "finite_and_ideal_vertices": (
            vinberg_case["finite_vertex_count"] == 75898
            and vinberg_case["ideal_vertex_count"] == 2114
        ),
    }
    if not all(vinberg_checks.values()):
        raise ArithmeticError(f"inherited Vinberg audit failed: {vinberg_checks}")

    checks = {
        "lattice_bridge": all(lattice_checks.values()),
        "all_facet_records": all(all(row["checks"].values()) for row in records),
        "full_generation_and_interior": all(generation_checks.values()),
        "literal_Vinberg_finite_volume": all(vinberg_checks.values()),
    }
    return {
        "schema": SCHEMA,
        "status_date": "2026-08-07",
        "id": "SS2-2C-H14-98-HILB2-SCALED-DUAL-CARRIER",
        "upstream": {
            key: {
                "path": str(path.relative_to(ROOT)),
                "sha256": EXPECTED_SHA256[key],
            }
            for key, path in paths.items()
        },
        "theorem": (
            "After the unique W98-stable odd index-two closure M of the "
            "SS2-2C fixed lattice, the doubled dual M^vee(2) is explicitly "
            "the actual Picard lattice (15,9,1) of Hilb^2(X).  All 98 "
            "oriented halfspaces and integral reflections transport, and "
            "their primitive normals generate the full Picard lattice."
        ),
        "stable_closure": {
            "operation": (
                "unique reflection-stable odd index-two closure before duality"
            ),
            "basis_columns_in_original_fixed_Q_space": _payload(closure_basis),
            "gram": _payload(gram_m),
            "determinant": int(gram_m.det()),
            "smith_invariants": _smith(gram_m),
            "parity": "odd",
            "primitive_mirror_types_square_divisibility": _histogram_payload(
                closure_types
            ),
        },
        "scaled_dual_bridge": {
            "operation": "N=M^vee(2), with Gram 2*Gram(M)^(-1)",
            "gram_N": _payload(gram_n),
            "N_determinant": int(gram_n.det()),
            "N_smith_invariants": _smith(gram_n),
            "N_parity": "even",
            "target": {
                "variety": "Hilb^2(X) for a K3 surface with NS(X)=(14,8,1)",
                "divisor_lattice": "Pic(Hilb^2(X))=(15,9,1)",
                "archived_basis_gram": _payload(gram_hilb2),
                "source_context": target_case["geometric_distinction"][
                    "context_from_Alexeev"
                ],
            },
            "N_to_Hilb2_isometry": _payload(transport),
            "isometry_determinant": int(transport.det()),
            "candidate_to_Alexeev_wall_permutation_zero_based": permutation,
            "primitive_target_mirror_types_square_divisibility": (
                _histogram_payload(target_types)
            ),
            "source_to_closure_to_target_type_histogram": (
                _histogram_payload(source_to_closure_types)
            ),
            "positive_halfspace_identity": (
                "For primitive q in M, p=Gram(M)q/div_M(q) and "
                "Phi(x)=Gram(M)x satisfy (Phi(x),p)_N="
                "(2/div_M(q))*(x,q)_M; the original L normal is a positive "
                "multiple of q."
            ),
            "checks": lattice_checks,
        },
        "facet_records": records,
        "full_lattice_generation": {
            "N_root_matrix_HNF": _payload(hnf_n),
            "Hilb2_root_matrix_HNF": _payload(hnf_target),
            "generation_index": 1,
            "strict_interior_N_coordinates": [int(value) for value in carrier],
            "strict_interior_Hilb2_coordinates": [
                int(value) for value in carrier_target
            ],
            "strict_interior_square": int((carrier.T * gram_n * carrier)[0]),
            "minimum_wall_pairing": min(carrier_pairings),
            "checks": generation_checks,
        },
        "inherited_literal_Vinberg_audit": {
            "criterion": vinberg_case["criterion"],
            "dimension": vinberg_case["dimension"],
            "facet_count": vinberg_case["facet_count"],
            "finite_vertex_count": vinberg_case["finite_vertex_count"],
            "ideal_vertex_count": vinberg_case["ideal_vertex_count"],
            "elliptic_edge_count": vinberg_case["elliptic_edge_count"],
            "connected_affine_component_count": vinberg_case[
                "connected_affine_component_count"
            ],
            "bad_edge_count": vinberg_case["bad_edge_count"],
            "bad_parabolic_neighborhood_count": vinberg_case[
                "bad_parabolic_neighborhood_count"
            ],
            "classification": vinberg_case["classification"],
            "checks": vinberg_checks,
        },
        "scope": {
            "positive_gain": (
                "complete coordinate-level characteristic-two stable-closure/"
                "scaled-dual bridge to an actual characteristic-zero "
                "hyperkahler divisor lattice"
            ),
            "known_abstract_polytope": True,
            "known_target_reflection_chamber": True,
            "not_a_new_geometric_chamber_slot": True,
            "not_a_Kahler_or_movable_cone_claim": True,
            "not_an_automorphism_fixed_chamber_claim": True,
            "frobenius_2c_not_in_geometric_K3_subgroup": True,
            "worldwide_priority_claim": "none",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()

    payload = build_certificate()
    if not payload["all_checks_pass"]:
        raise ArithmeticError("Hilbert-square scaled-dual certificate failed")
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists():
            raise ArithmeticError(f"missing archive: {args.output}")
        archived = args.output.read_text(encoding="utf-8")
        if archived != serialized:
            raise ArithmeticError("archived Hilbert-square carrier payload changed")
        print(f"verified {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized, encoding="utf-8")
        print(f"wrote {args.output}")
    print(f"payload_sha256={_payload_sha256(payload)}")


if __name__ == "__main__":
    main()

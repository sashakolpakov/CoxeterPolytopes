#!/usr/bin/env python3
"""Attach the A18-C04 H10/28 chamber to an explicit elliptic K3 family.

This is a composition certificate.  It pins the independent nine-parameter
U+D8(-1)+A1(-1) Weierstrass-family theorem and the independent scaled-dual
A18-C04 carrier theorem.  A deterministic sequence of reflections in
Roulleau's 14 published curve walls reduces the carrier marking.  In the
geometric fiber basis all 28 primitive mirrors then have coefficients of
absolute value at most six.

The transported strict chamber interior is big and nef.  Adding Roulleau's
published ample class gives a square-1408 ample divisor lying strictly inside
all 28 halfspaces.  Thus the known finite-volume noncompact H10/28 chamber is
realized literally in the divisor lattice of every very-general marked member
of the equation family.  It is not asserted to be the surface's nef chamber or
to arise from surface automorphisms.
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
from sympy.matrices.normalforms import hermite_normal_form


ROOT = Path(__file__).resolve().parents[1]
FAMILY = ROOT / "data/enumerations/ud8a1_weierstrass_family.json"
CARRIER = ROOT / "data/enumerations/a18_c04_ud8a1_scaled_dual_carrier.json"
ROULLEAU = ROOT / "data/enumerations/roulleau_rank11_direct_graph_denominator.json"
OUTPUT = ROOT / "data/enumerations/a18_c04_ud8a1_weierstrass_carrier.json"

SCHEMA = "coxeter.a18_c04.ud8a1_weierstrass_divisor_carrier.v1"
EXPECTED_SHA256 = {
    "family": "6669585141296dcff27267566d7f9a7f21e205961ea4c2999626a3dc96b72041",
    "carrier": "52974f5036bd20c512d2e9218c6300079af266c1c924c7732c695689f1f13932",
    "roulleau": "d82e31de7dd9d1f567cba9cc9a98588168628714521842c87943a9c76ec20d43",
}

# Roulleau's published ample class in the archived curve basis
# [A1,A2,A3,A4,A5,A6,A7,A8,A9,A11,A12].
ROULLEAU_AMPLE = sp.Matrix([-23, -45, -35, -24, -12, 1, 15, 15, 9, -16, -31])
EXPECTED_AMPLE_WALL_PAIRINGS = [1, 1, 1, 1, 1, 1, 1, 18, 12, 2, 1, 1, 1, 18]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [
            int(value) if sp.Rational(value).q == 1 else str(sp.Rational(value))
            for value in matrix.row(row)
        ]
        for row in range(matrix.rows)
    ]


def _content(vector: sp.Matrix) -> int:
    return math.gcd(*(abs(int(value)) for value in vector))


def _reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = (root.T * gram * root)[0]
    return sp.simplify(
        sp.eye(gram.rows) - sp.Rational(2, square) * root * root.T * gram
    )


def _expression(coordinates: sp.Matrix, labels: list[str]) -> str:
    terms: list[str] = []
    for coefficient, label in zip(coordinates, labels):
        coefficient = int(coefficient)
        if coefficient == 0:
            continue
        if coefficient == 1:
            term = label
        elif coefficient == -1:
            term = f"-{label}"
        else:
            term = f"{coefficient}*{label}"
        terms.append(term)
    if not terms:
        return "0"
    return " + ".join(terms).replace("+ -", "- ")


def _reduce_to_nef(
    vector: sp.Matrix, gram: sp.Matrix, walls: list[sp.Matrix]
) -> tuple[sp.Matrix, sp.Matrix, list[int]]:
    """Reflect a positive vector into the archived nef chamber exactly.

    At every step choose the most negative wall pairing, breaking ties by
    the one-based wall number.  This makes the 226-step replay deterministic.
    """

    reduced = vector
    action = sp.eye(gram.rows)
    sequence: list[int] = []
    for _ in range(10_000):
        pairings = [int((reduced.T * gram * wall)[0]) for wall in walls]
        negative = [(pairing, index) for index, pairing in enumerate(pairings) if pairing < 0]
        if not negative:
            return action, reduced, sequence
        _, index = min(negative)
        reflection = _reflection(gram, walls[index])
        reduced = reflection * reduced
        action = reflection * action
        sequence.append(index + 1)
    raise ArithmeticError("Roulleau nef reduction did not terminate")


def build_certificate() -> dict[str, Any]:
    paths = {"family": FAMILY, "carrier": CARRIER, "roulleau": ROULLEAU}
    for key, path in paths.items():
        if _sha256(path) != EXPECTED_SHA256[key]:
            raise ArithmeticError(f"pinned {key} archive changed")

    family = json.loads(FAMILY.read_text(encoding="utf-8"))
    carrier = json.loads(CARRIER.read_text(encoding="utf-8"))
    roulleau = json.loads(ROULLEAU.read_text(encoding="utf-8"))
    source = next(row for row in roulleau["cases"] if row["id"] == "R11_UD8A1")
    if not (
        family["all_checks_pass"]
        and carrier["all_checks_pass"]
        and roulleau["all_checks_pass"]
        and source["all_checks_pass"]
    ):
        raise ArithmeticError("an upstream positive theorem is not green")

    gram_standard = sp.Matrix(family["very_general_picard_lattice"]["gram"])
    gram_roulleau = sp.Matrix(source["source"]["basis_gram"])
    gram_n = sp.Matrix(carrier["scaled_dual_bridge"]["gram_N"])
    # Rows give standard fiber-basis vectors in the Roulleau curve basis.
    standard_rows_in_roulleau = sp.Matrix(
        family["explicit_fiber_basis_bridge"]["standard_to_roulleau_rows_T"]
    )
    standard_columns_to_roulleau = standard_rows_in_roulleau.T
    n_to_roulleau = sp.Matrix(
        carrier["scaled_dual_bridge"]["N_to_K3_isometry"]
    )

    walls_roulleau = [
        sp.Matrix(row) for row in source["source"]["wall_coordinate_rows"]
    ]
    interior_n = sp.Matrix(
        carrier["full_lattice_generation"]["strict_interior_N_coordinates"]
    )
    initial_interior_roulleau = n_to_roulleau * interior_n
    weyl_action, nef_interior_roulleau, reflection_sequence = _reduce_to_nef(
        initial_interior_roulleau, gram_roulleau, walls_roulleau
    )
    reduced_n_to_roulleau = weyl_action * n_to_roulleau
    n_to_standard = standard_columns_to_roulleau.inv() * reduced_n_to_roulleau

    wall_squares = [int((wall.T * gram_roulleau * wall)[0]) for wall in walls_roulleau]
    nef_pairings = [
        int((nef_interior_roulleau.T * gram_roulleau * wall)[0])
        for wall in walls_roulleau
    ]
    reduction_checks = {
        "four_pinned_archives_green": (
            family["all_checks_pass"]
            and carrier["all_checks_pass"]
            and roulleau["all_checks_pass"]
            and source["all_checks_pass"]
        ),
        "all_14_published_walls_are_minus_two": wall_squares == [-2] * 14,
        "deterministic_reflection_count_226": len(reflection_sequence) == 226,
        "weyl_action_integral_unimodular": (
            all(value.is_Integer for value in weyl_action)
            and weyl_action.det() == 1
        ),
        "weyl_action_isometry": (
            weyl_action.T * gram_roulleau * weyl_action == gram_roulleau
        ),
        "reduced_interior_exact": list(map(int, nef_interior_roulleau))
        == [-41, -82, -63, -44, -22, 0, 22, 27, 15, -30, -60],
        "reduced_interior_square_508": (
            nef_interior_roulleau.T * gram_roulleau * nef_interior_roulleau
        )[0]
        == 508,
        "reduced_interior_is_nef": nef_pairings
        == [0, 0, 0, 3, 0, 0, 10, 20, 24, 0, 0, 8, 0, 30],
        "reduced_transport_integral_unimodular": (
            all(value.is_Integer for value in reduced_n_to_roulleau)
            and reduced_n_to_roulleau.det() == 1
        ),
        "fiber_transport_integral_unimodular": (
            all(value.is_Integer for value in n_to_standard)
            and n_to_standard.det() == -1
        ),
        "fiber_transport_gram_identity": (
            n_to_standard.T * gram_standard * n_to_standard == gram_n
        ),
        "two_step_transport_identity": (
            standard_columns_to_roulleau * n_to_standard
            == reduced_n_to_roulleau
        ),
    }
    if not all(reduction_checks.values()):
        raise ArithmeticError(
            "equation-basis reduction failed: "
            f"{[key for key, value in reduction_checks.items() if not value]}"
        )

    labels = [
        "F",
        "(F+O)",
        "D1",
        "D2",
        "D3",
        "D4",
        "D5",
        "D6",
        "D7",
        "D8",
        "C",
    ]
    roots_n: list[sp.Matrix] = []
    roots_standard: list[sp.Matrix] = []
    facet_records: list[dict[str, Any]] = []
    type_histogram: Counter[tuple[int, int]] = Counter()
    maximum_root_coordinate = 0

    for archived in carrier["facet_records"]:
        root_n = sp.Matrix(archived["primitive_oriented_N_coordinates"])
        archived_roulleau_root = sp.Matrix(
            archived["primitive_oriented_K3_coordinates"]
        )
        reduced_roulleau_root = reduced_n_to_roulleau * root_n
        root_standard = n_to_standard * root_n
        square = int((root_standard.T * gram_standard * root_standard)[0])
        divisibility = _content(gram_standard * root_standard)
        reflection_n = _reflection(gram_n, root_n)
        reflection_standard = _reflection(gram_standard, root_standard)
        maximum_root_coordinate = max(
            maximum_root_coordinate, *(abs(int(value)) for value in root_standard)
        )
        type_histogram[(square, divisibility)] += 1

        checks = {
            "archived_Roulleau_root_replayed": n_to_roulleau * root_n
            == archived_roulleau_root,
            "Weyl_reduction_replayed": weyl_action * archived_roulleau_root
            == reduced_roulleau_root,
            "fiber_basis_transport_replayed": (
                standard_columns_to_roulleau * root_standard
                == reduced_roulleau_root
            ),
            "primitive": _content(root_standard) == 1,
            "square_and_divisibility_preserved": (
                square == archived["primitive_K3_square"]
                and divisibility == archived["primitive_K3_divisibility"]
            ),
            "reflection_integral_unimodular_involution": (
                all(value.is_Integer for value in reflection_standard)
                and reflection_standard.det() == -1
                and reflection_standard**2 == sp.eye(11)
            ),
            "reflection_isometry": (
                reflection_standard.T
                * gram_standard
                * reflection_standard
                == gram_standard
            ),
            "reflection_intertwines": (
                reflection_standard * n_to_standard
                == n_to_standard * reflection_n
            ),
        }
        if not all(checks.values()):
            raise ArithmeticError(
                f"facet {archived['facet_zero_based']} failed equation attachment: "
                f"{[key for key, value in checks.items() if not value]}"
            )

        roots_n.append(root_n)
        roots_standard.append(root_standard)
        facet_records.append(
            {
                "facet_zero_based": archived["facet_zero_based"],
                "source_orbit_index_one_based": archived[
                    "source_orbit_index_one_based"
                ],
                "source_presentation_orientation": archived[
                    "source_presentation_orientation"
                ],
                "positive_halfspace_factor_from_original_L": archived[
                    "positive_halfspace_factor_from_original_L"
                ],
                "primitive_N_coordinates": [int(value) for value in root_n],
                "primitive_equation_fiber_basis_coordinates": [
                    int(value) for value in root_standard
                ],
                "divisor_class_expression": _expression(root_standard, labels),
                "square": square,
                "divisibility": divisibility,
                "reflection_matrix_in_equation_fiber_basis": _payload(
                    reflection_standard
                ),
                "checks": checks,
            }
        )

    roots_standard_matrix = sp.Matrix.hstack(*roots_standard).applyfunc(int)
    roots_n_matrix = sp.Matrix.hstack(*roots_n).applyfunc(int)
    gram_28_standard = roots_standard_matrix.T * gram_standard * roots_standard_matrix
    gram_28_n = roots_n_matrix.T * gram_n * roots_n_matrix
    interior_standard = n_to_standard * interior_n
    ample_standard = standard_columns_to_roulleau.inv() * ROULLEAU_AMPLE
    ample_inside_standard = interior_standard + ample_standard
    ample_inside_roulleau = nef_interior_roulleau + ROULLEAU_AMPLE
    ample_wall_pairings = [
        int((ROULLEAU_AMPLE.T * gram_roulleau * wall)[0])
        for wall in walls_roulleau
    ]
    imported_interior_pairings = [
        int((interior_standard.T * gram_standard * root)[0])
        for root in roots_standard
    ]
    ample_inside_imported_pairings = [
        int((ample_inside_standard.T * gram_standard * root)[0])
        for root in roots_standard
    ]
    ample_inside_nef_pairings = [
        int((ample_inside_roulleau.T * gram_roulleau * wall)[0])
        for wall in walls_roulleau
    ]
    generation_checks = {
        "all_28_roots_transported": len(roots_standard) == 28,
        "all_28_oriented_Gram_entries_preserved": gram_28_standard == gram_28_n,
        "all_roots_generate_full_equation_NS": (
            hermite_normal_form(roots_standard_matrix) == sp.eye(11)
        ),
        "root_type_histogram": type_histogram
        == Counter({(-2, 1): 16, (-2, 2): 8, (-4, 2): 4}),
        "maximum_absolute_root_coordinate_six": maximum_root_coordinate == 6,
        "interior_is_sum_of_28_oriented_roots": (
            sum(roots_standard, sp.zeros(11, 1)) == interior_standard
        ),
        "interior_square_508": (
            interior_standard.T * gram_standard * interior_standard
        )[0]
        == 508,
        "interior_strict_for_imported_chamber": (
            min(imported_interior_pairings) == 3
            and all(value > 0 for value in imported_interior_pairings)
        ),
        "published_ample_replayed": (
            (ROULLEAU_AMPLE.T * gram_roulleau * ROULLEAU_AMPLE)[0] == 208
            and ample_wall_pairings == EXPECTED_AMPLE_WALL_PAIRINGS
        ),
        "ample_inside_exact_coordinate_transport": (
            standard_columns_to_roulleau * ample_inside_standard
            == ample_inside_roulleau
        ),
        "ample_inside_square_1408": (
            ample_inside_standard.T * gram_standard * ample_inside_standard
        )[0]
        == 1408,
        "ample_inside_complete_nef_chamber": (
            min(ample_inside_nef_pairings) == 1
            and all(value > 0 for value in ample_inside_nef_pairings)
        ),
        "ample_inside_imported_chamber": (
            min(ample_inside_imported_pairings) == 7
            and all(value > 0 for value in ample_inside_imported_pairings)
        ),
    }
    if not all(generation_checks.values()):
        raise ArithmeticError(
            "equation-level generation/interior gate failed: "
            f"{[key for key, value in generation_checks.items() if not value]}"
        )

    vinberg = carrier["inherited_literal_Vinberg_audit"]
    vinberg_checks = {
        "classification": vinberg["classification"]
        == "finite_volume_noncompact_Coxeter_polytope",
        "dimension_facets": (vinberg["dimension"], vinberg["facet_count"])
        == (10, 28),
        "literal_edge_gate": (
            vinberg["elliptic_edge_count"] == 6956
            and vinberg["bad_edge_count"] == 0
        ),
        "literal_parabolic_gate": (
            vinberg["connected_affine_component_count"] == 204
            and vinberg["bad_parabolic_neighborhood_count"] == 0
        ),
        "all_archived_Vinberg_checks": all(vinberg["checks"].values()),
    }
    if not all(vinberg_checks.values()):
        raise ArithmeticError("inherited literal Vinberg gate changed")

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "id": "A18_4_0_C04_equation_level_UD8A1_H10_28",
        "theorem": (
            "Every very-general marked member of the displayed nine-parameter "
            "elliptic K3 family carries the known A18-C04 finite-volume "
            "noncompact H10/28 Coxeter chamber as 28 explicit primitive divisor "
            "mirrors.  They generate NS, their coefficients in the geometric "
            "fiber basis have absolute value at most six, and an explicit "
            "square-1408 ample divisor lies strictly inside all 28 halfspaces."
        ),
        "upstream": {
            key: {
                "path": str(path.relative_to(ROOT)),
                "sha256": EXPECTED_SHA256[key],
            }
            for key, path in paths.items()
        },
        "equation_family": {
            "parameters": family["weierstrass_family"]["parameters"],
            "P": family["weierstrass_family"]["P"],
            "Q": family["weierstrass_family"]["Q"],
            "equation": family["weierstrass_family"]["affine_equation"],
            "fiber_configuration": family["weierstrass_family"][
                "fiber_configuration"
            ],
            "nonempty_witness": family["weierstrass_family"][
                "nonempty_witness"
            ],
            "very_general_NS": family["very_general_picard_lattice"][
                "isometry_class"
            ],
            "very_general_MW": "trivial",
            "fiber_basis_labels": labels,
            "fiber_basis_geometry": family["explicit_fiber_basis_bridge"][
                "standard_basis"
            ],
        },
        "exact_composition": {
            "standard_rows_in_Roulleau_basis": _payload(
                standard_rows_in_roulleau
            ),
            "archived_N_to_Roulleau_basis": _payload(n_to_roulleau),
            "nef_reduction_wall_sequence_one_based": reflection_sequence,
            "nef_reduction_action_in_Roulleau_basis": _payload(weyl_action),
            "reduced_N_to_Roulleau_basis": _payload(reduced_n_to_roulleau),
            "reduced_N_to_equation_fiber_basis": _payload(n_to_standard),
            "determinant_N_to_equation_fiber_basis": int(n_to_standard.det()),
            "maximum_absolute_transport_entry": max(
                abs(int(value)) for value in n_to_standard
            ),
            "checks": reduction_checks,
        },
        "divisor_mirrors": {
            "count": 28,
            "type_histogram_square_divisibility": {
                f"{square},{divisibility}": count
                for (square, divisibility), count in sorted(type_histogram.items())
            },
            "maximum_absolute_fiber_basis_coordinate": maximum_root_coordinate,
            "oriented_28_by_28_Gram": _payload(gram_28_standard),
            "records": facet_records,
            "generation_index": 1,
        },
        "algebro_geometric_interior": {
            "transported_chamber_interior_fiber_coordinates": [
                int(value) for value in interior_standard
            ],
            "transported_chamber_interior_expression": _expression(
                interior_standard, labels
            ),
            "transported_chamber_interior_square": 508,
            "transported_chamber_minimum_pairing": min(
                imported_interior_pairings
            ),
            "nef_wall_pairings": nef_pairings,
            "status": "big_and_nef; strict_for_the_imported_28_wall_chamber",
            "published_ample_fiber_coordinates": [
                int(value) for value in ample_standard
            ],
            "published_ample_square": 208,
            "ample_inside_both_chambers_fiber_coordinates": [
                int(value) for value in ample_inside_standard
            ],
            "ample_inside_both_chambers_expression": _expression(
                ample_inside_standard, labels
            ),
            "ample_inside_both_chambers_square": 1408,
            "ample_inside_complete_nef_wall_pairings": ample_inside_nef_pairings,
            "ample_inside_imported_wall_pairings": ample_inside_imported_pairings,
            "minimum_complete_nef_wall_pairing": min(
                ample_inside_nef_pairings
            ),
            "minimum_imported_wall_pairing": min(
                ample_inside_imported_pairings
            ),
            "conclusion": (
                "the imported finite-volume Coxeter chamber has nonempty "
                "intersection with the K3 ample cone"
            ),
            "checks": generation_checks,
        },
        "literal_Vinberg_gate": {
            "classification": vinberg["classification"],
            "dimension": vinberg["dimension"],
            "facet_count": vinberg["facet_count"],
            "finite_vertex_count": vinberg["finite_vertex_count"],
            "ideal_vertex_count": vinberg["ideal_vertex_count"],
            "elliptic_edge_count": vinberg["elliptic_edge_count"],
            "connected_affine_component_count": vinberg[
                "connected_affine_component_count"
            ],
            "bad_edge_count": vinberg["bad_edge_count"],
            "bad_parabolic_neighborhood_count": vinberg[
                "bad_parabolic_neighborhood_count"
            ],
            "checks": vinberg_checks,
        },
        "lineage_and_scope": {
            "classical_surface_input": (
                "Roulleau published the U+D8(-1)+A1(-1) K3 lattice, fiber "
                "configuration, curve basis, ample class, and complete nef chamber"
            ),
            "known_abstract_chamber": (
                "A18-C04 H10/28 and its literal finite-volume Vinberg audit "
                "predate this composition"
            ),
            "new_positive_gain": (
                "a separately replayable equation-to-28-wall attachment with "
                "small divisor coordinates and an ample polarization inside the "
                "imported chamber"
            ),
            "not_the_complete_nef_chamber": True,
            "minus_four_mirrors_are_not_curve_reflections": True,
            "minus_two_classes_not_asserted_irreducible": True,
            "no_surface_automorphism_claim": True,
            "not_a_new_abstract_polytope": True,
            "worldwide_priority_claim": "none",
        },
        "checks": {
            "equation_family": family["all_checks_pass"],
            "scaled_dual_carrier": carrier["all_checks_pass"],
            "composition": all(reduction_checks.values()),
            "mirrors_and_ample_interior": all(generation_checks.values()),
            "literal_Vinberg": all(vinberg_checks.values()),
        },
        "all_checks_pass": True,
    }
    payload["payload_sha256"] = _payload_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    if args.check_archive:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if archived != payload:
            raise ArithmeticError("rebuilt equation-level H10/28 carrier changed")
        print(f"archive check passed: {OUTPUT}")
        return
    OUTPUT.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(OUTPUT)
    print(payload["payload_sha256"])


if __name__ == "__main__":
    main()

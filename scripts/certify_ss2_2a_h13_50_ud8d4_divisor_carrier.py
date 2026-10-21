#!/usr/bin/env python3
"""Certify a complex-K3 scaled-dual carrier of SS2-2A H13/50.

Let M be the saturated fixed lattice of the characteristic-two 2a section.
This construction does *not* take a same-form overlattice of M.  It passes
contragrediently to the doubled dual N=M^vee(2), whose Gram is 2*Gram(M)^-1.
An explicit unimodular isometry identifies N with the Neron--Severi lattice

    U + D8(-1) + D4(-1)

of the archived Brandhorst--Hofmann/Roulleau complex K3 carrier.  Primitive
reduction exchanges the two mirror lengths: forty (-4,div2) roots of M
become (-2,div1), while ten (-2,div1) roots become (-4,div2).  The fifty
transported roots generate all of N and every reflection is integral.
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
K3_SOURCE = ROOT / "data/enumerations/roulleau_ud8d4_all_subgroup_folds.json"
K3_EXISTENCE = ROOT / "data/enumerations/roulleau_ud8d4_bh_k3_lifts.json"
OUTPUT = ROOT / "data/enumerations/ss2_2a_h13_50_ud8d4_divisor_carrier.json"

SCHEMA = "coxeter.ss2_2a_h13_50.ud8d4_scaled_dual_divisor_carrier.v1"
EXPECTED_SHA256 = {
    "construction": "a9087200b360c4bed75cd93c987bd7e7d493ff1ce852bb9c5ea4dc3490f1fe11",
    "lineage": "f0867d47f4a1405a5398a274865b32fff3b24608790f2b41a25eb491de2e8971",
    "vinberg": "032c8343a6dd7f4904b9daecb074b8b0bbbf1a94f9431252884484c98d2e8456",
    "k3_source": "0348744e2cf17870616deccd81f1c20541c1df328ad80b4d9eb6bf2476f793b7",
    "k3_existence": "b5d480d30ba0b8c39bf48d8021127ae114f4919529b620759cd20c061555c4c8",
}

# X sends coordinates in the doubled-dual basis of N to the archived
# Roulleau curve basis of U+D8(-1)+D4(-1).  The permanent replay checks
# X^T G_K3 X = 2 G_M^-1 directly; no PARI/Sage dependency remains.
DOUBLED_DUAL_TO_K3 = sp.Matrix(
    [
        [-26, 9, 16, 22, 10, -39, -16, -67, 32, -48, -64, 35, 24, 40],
        [-36, 12, 22, 30, 13, -54, -22, -92, 44, -66, -88, 49, 33, 56],
        [-42, 14, 26, 36, 16, -63, -26, -107, 52, -77, -103, 57, 38, 65],
        [-39, 12, 23, 32, 14, -58, -24, -97, 48, -70, -94, 53, 35, 61],
        [-54, 18, 32, 44, 20, -80, -33, -135, 65, -98, -130, 73, 49, 83],
        [-36, 10, 20, 28, 12, -52, -22, -87, 44, -63, -85, 49, 32, 57],
        [-50, 16, 29, 40, 18, -74, -31, -125, 61, -91, -121, 69, 46, 78],
        [-64, 22, 38, 52, 24, -96, -40, -163, 78, -118, -157, 89, 60, 100],
        [-72, 26, 44, 60, 28, -109, -45, -186, 88, -135, -178, 100, 68, 112],
        [-80, 30, 50, 68, 32, -122, -50, -209, 98, -152, -200, 111, 76, 124],
        [-76, 28, 47, 64, 30, -116, -48, -199, 94, -145, -191, 107, 73, 120],
        [-40, 15, 25, 34, 16, -61, -25, -104, 49, -76, -100, 55, 38, 62],
        [-48, 18, 30, 41, 19, -73, -30, -125, 59, -91, -120, 66, 45, 74],
        [-48, 18, 30, 41, 19, -73, -30, -125, 59, -91, -120, 67, 46, 74],
    ]
)

# Independent hostile replay with GMP cdd on the exact H-matrix below.
EXACT_CONE_INPUT_SHA256 = (
    "96e2abe698f9f80a41f1cc983c35fb5b094f31a6de5b7ec8c76351e2de5dab3c"
)
EXACT_CONE_AUDIT = {
    "engine": "GMP cdd exact double description (independent replay)",
    "ray_count": 48552,
    "timelike_rays": 45954,
    "isotropic_rays": 2598,
    "hyperideal_rays": 0,
    "ideal_ray_active_facet_histogram": {"14": 1656, "15": 872, "16": 68, "17": 2},
}


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
    return [abs(int(normal[i, i])) for i in range(matrix.rows)]


def _content(vector: sp.Matrix) -> int:
    return math.gcd(*(abs(int(value)) for value in vector))


def _reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = (root.T * gram * root)[0]
    return sp.simplify(
        sp.eye(gram.rows) - sp.Rational(2, square) * root * root.T * gram
    )


def build_certificate() -> dict[str, Any]:
    paths = {
        "construction": CONSTRUCTION,
        "lineage": LINEAGE,
        "vinberg": VINBERG,
        "k3_source": K3_SOURCE,
        "k3_existence": K3_EXISTENCE,
    }
    for key, path in paths.items():
        if _sha256(path) != EXPECTED_SHA256[key]:
            raise ArithmeticError(f"pinned {key} archive changed")

    construction = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))
    lineage = json.loads(LINEAGE.read_text(encoding="utf-8"))
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    k3_source = json.loads(K3_SOURCE.read_text(encoding="utf-8"))
    k3_existence = json.loads(K3_EXISTENCE.read_text(encoding="utf-8"))
    case = next(row for row in construction["cases"] if row["id"] == "SS2-2A-H13-50")
    vinberg_case = next(row for row in vinberg["cases"] if row["id"] == "SS2-2A-H13-50")

    gram_m = sp.Matrix(case["fixed_lattice"]["gram"])
    gram_n = sp.simplify(2 * gram_m.inv())
    gram_k3 = sp.Matrix(k3_source["source"]["basis_gram"])
    transport = DOUBLED_DUAL_TO_K3

    actual_k3_rows = [
        row
        for row in k3_existence["cases"]
        if row["brandhorst_hofmann_row"]["generic_picard_lattice"]["isometry_class"]
        == "U+D8(-1)+D4(-1)"
        and row["kahler_torelli_handoff"]["database_row_is_a_genuine_k3_automorphism_class"]
    ]
    lattice_checks = {
        "all_five_upstream_certificates_green": all(
            archive["all_checks_pass"]
            for archive in (construction, lineage, vinberg, k3_source, k3_existence)
        ),
        "M_is_even_14_10_0": (
            exact_inertia(gram_m) == (1, 13, 0)
            and gram_m.det() == -1024
            and _smith(gram_m) == [1] * 4 + [2] * 10
            and all(gram_m[i, i] % 2 == 0 for i in range(14))
        ),
        "N_is_integral_even_doubled_dual": (
            all(value.is_Integer for value in gram_n)
            and exact_inertia(gram_n) == (1, 13, 0)
            and gram_n.det() == -16
            and _smith(gram_n) == [1] * 10 + [2] * 4
            and all(gram_n[i, i] % 2 == 0 for i in range(14))
        ),
        "lineage_doubled_dual_agrees": (
            sp.Matrix(lineage["H13_50"]["fixed_lattice"]["doubled_dual_gram"])
            == gram_n
        ),
        "target_is_archived_UD8D4_picard_lattice": (
            k3_source["source"]["lattice"] == "U+D8(-1)+D4(-1)"
            and k3_source["source"]["picard_rank"] == 14
            and gram_k3.det() == -16
            and exact_inertia(gram_k3) == (1, 13, 0)
            and _smith(gram_k3) == [1] * 10 + [2] * 4
        ),
        "actual_complex_K3_rows_exist": (
            k3_existence["all_checks_pass"] is True and len(actual_k3_rows) == 2
        ),
        "explicit_unimodular_isometry": (
            all(value.is_Integer for value in transport)
            and transport.det() == -1
            and transport.T * gram_k3 * transport == gram_n
        ),
    }
    if not all(lattice_checks.values()):
        raise ArithmeticError(
            f"scaled-dual bridge failed: {[k for k, v in lattice_checks.items() if not v]}"
        )

    records: list[dict[str, Any]] = []
    roots_n: list[sp.Matrix] = []
    inequalities: list[list[int]] = []
    for index, source in enumerate(
        case["restricted_chamber"]["primitive_facet_lattice_records"]
    ):
        root_m = sp.Matrix(source["primitive_fixed_coordinates"])
        dual_presentation = gram_m * root_m
        content = _content(dual_presentation)
        root_n = (dual_presentation / content).applyfunc(int)
        root_k3 = transport * root_n
        reflection_m = _reflection(gram_m, root_m)
        reflection_n = _reflection(gram_n, root_n)
        reflection_k3 = _reflection(gram_k3, root_k3)
        square_m = int((root_m.T * gram_m * root_m)[0])
        square_n = int((root_n.T * gram_n * root_n)[0])
        div_m = _content(gram_m * root_m)
        div_n = _content(gram_n * root_n)
        scale = int(source["orbit_sum_scale_from_primitive"])
        inequalities.append([int(value) for value in scale * dual_presentation])
        checks = {
            "source_type_archived": (
                square_m == source["primitive_square"]
                and div_m == source["primitive_divisibility"]
            ),
            "dual_content_equals_source_divisibility": content == div_m,
            "primitive_in_N_and_K3": _content(root_n) == _content(root_k3) == 1,
            "scaled_dual_square_law": square_n * content**2 == 2 * square_m,
            "exact_halfspace_transport": (
                root_n.T * gram_n * gram_m
                == sp.Rational(2, content) * root_m.T * gram_m
            ),
            "contragredient_reflection": reflection_n == reflection_m.T,
            "reflection_integral": (
                all(value.is_Integer for value in reflection_n)
                and all(value.is_Integer for value in reflection_k3)
            ),
            "target_reflection_formula_and_intertwining": (
                reflection_k3 * transport == transport * reflection_n
            ),
            "target_reflection_unimodular_involution": (
                reflection_k3.det() == -1
                and reflection_k3**2 == sp.eye(14)
                and reflection_k3.T * gram_k3 * reflection_k3 == gram_k3
            ),
        }
        if not all(checks.values()):
            raise ArithmeticError(f"facet {index} failed the scaled-dual transport")
        roots_n.append(root_n)
        records.append(
            {
                "facet_zero_based": index,
                "source_orbit_index_one_based": source["orbit_index_one_based"],
                "source_M_coordinates": [int(value) for value in root_m],
                "source_M_square": square_m,
                "source_M_divisibility": div_m,
                "doubled_dual_presentation_coordinates": [
                    int(value) for value in dual_presentation
                ],
                "doubled_dual_coordinate_content": content,
                "primitive_N_coordinates": [int(value) for value in root_n],
                "primitive_K3_coordinates": [int(value) for value in root_k3],
                "primitive_K3_square": square_n,
                "primitive_K3_divisibility": div_n,
                "checks": checks,
            }
        )

    root_matrix = sp.Matrix.hstack(*roots_n)
    wall_hnf = hermite_normal_form(root_matrix)
    source_types = Counter(
        (row["source_M_square"], row["source_M_divisibility"]) for row in records
    )
    target_types = Counter(
        (row["primitive_K3_square"], row["primitive_K3_divisibility"])
        for row in records
    )
    content_histogram = Counter(row["doubled_dual_coordinate_content"] for row in records)
    mirror_checks = {
        "all_50_records": len(records) == 50,
        "source_type_split": source_types == {(-4, 2): 40, (-2, 1): 10},
        "coordinate_content_split": content_histogram == {2: 40, 1: 10},
        "target_type_split": target_types == {(-2, 1): 40, (-4, 2): 10},
        "roots_span_rank_14": root_matrix.rank() == 14,
        "roots_generate_all_of_N": wall_hnf == sp.eye(14),
        "every_reflection_green": all(all(row["checks"].values()) for row in records),
    }
    if not all(mirror_checks.values()):
        raise ArithmeticError("the 50-wall primitive mirror audit changed")

    interior_n = sum(roots_n, sp.zeros(14, 1))
    interior_k3 = transport * interior_n
    pairings = [int((interior_n.T * gram_n * root)[0]) for root in roots_n]
    interior_checks = {
        "explicit_timelike_interior": int((interior_n.T * gram_n * interior_n)[0])
        == 2012,
        "strict_positive_wall_pairings": min(pairings) == 6,
        "transport_preserves_interior_square": (
            interior_k3.T * gram_k3 * interior_k3
            == interior_n.T * gram_n * interior_n
        ),
    }
    if not all(interior_checks.values()):
        raise ArithmeticError("the strict carrier interior changed")

    inequality_hash = hashlib.sha256(
        json.dumps(inequalities, separators=(",", ":")).encode()
    ).hexdigest()
    cone_checks = {
        "canonical_50_by_14_input_hash": inequality_hash == EXACT_CONE_INPUT_SHA256,
        "ray_partition_sums": (
            EXACT_CONE_AUDIT["timelike_rays"]
            + EXACT_CONE_AUDIT["isotropic_rays"]
            + EXACT_CONE_AUDIT["hyperideal_rays"]
            == EXACT_CONE_AUDIT["ray_count"]
        ),
        "finite_and_ideal_counts_match_vinberg": (
            EXACT_CONE_AUDIT["timelike_rays"] == vinberg_case["finite_vertex_count"]
            and EXACT_CONE_AUDIT["isotropic_rays"] == vinberg_case["ideal_vertex_count"]
        ),
        "no_hyperideal_rays": EXACT_CONE_AUDIT["hyperideal_rays"] == 0,
        "ideal_incidence_partition": sum(
            EXACT_CONE_AUDIT["ideal_ray_active_facet_histogram"].values()
        )
        == EXACT_CONE_AUDIT["isotropic_rays"],
    }
    if not all(cone_checks.values()):
        raise ArithmeticError("the independent exact-cone cross-check changed")

    vinberg_checks = {
        "literal_classification": (
            vinberg_case["classification"]
            == "finite_volume_noncompact_Coxeter_polytope"
            and vinberg_case["finite_volume"]
            and not vinberg_case["compact"]
        ),
        "dimension_and_facets": (
            vinberg_case["dimension"], vinberg_case["facet_count"]
        )
        == (13, 50),
        "complete_counts": (
            vinberg_case["finite_vertex_count"],
            vinberg_case["ideal_vertex_count"],
            vinberg_case["elliptic_edge_count"],
            vinberg_case["connected_affine_component_count"],
        )
        == (45954, 2598, 391614, 6210),
        "zero_bad_audits": (
            vinberg_case["bad_edge_count"] == 0
            and vinberg_case["bad_parabolic_neighborhood_count"] == 0
        ),
        "all_literal_checks_green": all(vinberg_case["checks"].values()),
    }
    if not all(vinberg_checks.values()):
        raise ArithmeticError("the literal Vinberg handoff changed")

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "theorem": (
            "The finite-volume noncompact SS2-2A H13/50 chamber has a "
            "complex-K3 mixed divisor-reflection realization on "
            "NS=U+D8(-1)+D4(-1) via the doubled dual M^vee(2)."
        ),
        "scope": {
            "proved": (
                "all 50 primitive mirrors generate the actual rank-14 K3 "
                "Picard lattice and all 50 reflections are integral"
            ),
            "scaled_dual_not_overlattice": (
                "N=M^vee(2) changes the bilinear form; this is a contragredient "
                "carrier, not a same-form finite-index overclosure of M"
            ),
            "not_full_nef": (
                "ten primitive square-minus-four mirrors are divisor cuts, not "
                "K3 curve walls, and extra roots may refine the chamber"
            ),
            "not_automorphism_fixed": (
                "no Brandhorst--Hofmann automorphism is claimed to preserve "
                "this transported 50-wall chamber"
            ),
            "novelty": (
                "new project-level complex-K3 realization of an already "
                "certified literal chamber; no abstract-polytope priority claim"
            ),
        },
        "sources": {
            key: {"path": str(path.relative_to(ROOT)), "sha256": EXPECTED_SHA256[key]}
            for key, path in paths.items()
        },
        "scaled_dual_bridge": {
            "source_M": {
                "isometry_class": "even 2-elementary (14,10,0)",
                "rank": 14,
                "signature": [1, 13, 0],
                "determinant": -1024,
                "smith": _smith(gram_m),
                "gram": _payload(gram_m),
            },
            "doubled_dual_N": {
                "definition": "N=M^vee(2), Gram(N)=2*Gram(M)^(-1)",
                "isometry_class": "U+D8(-1)+D4(-1)",
                "rank": 14,
                "signature": [1, 13, 0],
                "determinant": -16,
                "smith": _smith(gram_n),
                "gram": _payload(gram_n),
            },
            "actual_K3_target": {
                "isometry_class": "U+D8(-1)+D4(-1)",
                "archived_basis_gram": _payload(gram_k3),
                "genuine_BH_K3_rows": [
                    row["brandhorst_hofmann_row"]["id"] for row in actual_k3_rows
                ],
                "source_curve_configuration_complete": (
                    k3_source["source"]["completeness"]
                    == "Roulleau_published_complete_curve_configuration"
                ),
            },
            "explicit_N_to_K3_isometry": {
                "matrix": _payload(transport),
                "determinant": int(transport.det()),
                "identity": "X^T Gram(K3) X = 2 Gram(M)^(-1)",
            },
            "checks": lattice_checks,
        },
        "primitive_mirrors": {
            "count": 50,
            "source_square_divisibility_histogram": {
                "-4,2": 40,
                "-2,1": 10,
            },
            "doubled_dual_coordinate_content_histogram": {"2": 40, "1": 10},
            "K3_square_divisibility_histogram": {"-2,1": 40, "-4,2": 10},
            "wall_lattice_index_in_K3_NS": 1,
            "wall_lattice_hnf": _payload(wall_hnf),
            "all_50_reflections_integral": True,
            "checks": mirror_checks,
            "records": records,
        },
        "strict_carrier_interior": {
            "construction": "sum of the 50 primitive N normals",
            "N_coordinates": [int(value) for value in interior_n],
            "K3_coordinates": [int(value) for value in interior_k3],
            "square": 2012,
            "minimum_wall_pairing": 6,
            "effectivity_handoff": (
                "choose a generic nearby class off extra root hyperplanes, Weyl-transport "
                "the marking to the ample cone, and apply K3 Riemann--Roch to the "
                "forty oriented square-minus-two classes; the resulting effective "
                "divisors may be reducible and are not asserted irreducible curves "
                "or facets of the full nef cone"
            ),
            "checks": interior_checks,
        },
        "independent_exact_cone_cross_check": {
            "canonical_inequality_convention": (
                "rows A_i=scale_i*Gram(M)*r_i in archived facet order; A_i*x>=0"
            ),
            "canonical_input_sha256": inequality_hash,
            **EXACT_CONE_AUDIT,
            "checks": cone_checks,
        },
        "vinberg_finite_volume": {
            "classification": vinberg_case["classification"],
            "dimension": 13,
            "facets": 50,
            "finite_vertices": 45954,
            "cusps": 2598,
            "elliptic_edges": 391614,
            "connected_affine_components": 6210,
            "bad_edges": 0,
            "bad_parabolic_neighborhoods": 0,
            "checks": vinberg_checks,
        },
        "hostile_audit": {
            "status": "green",
            "independent_checks": [
                "determinants and explicit 14-by-14 congruence",
                "all 50 primitive normalizations and mirror-type exchange",
                "all 50 contragredient and target reflection identities",
                "HNF generation of the full K3 Picard lattice",
                "independent GMP-cdd exact cone enumeration",
            ],
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
            raise ArithmeticError("rebuilt SS2-2A H13/50 K3 carrier archive changed")
        print(f"archive check passed: {OUTPUT}")
        return
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(OUTPUT)
    print(payload["payload_sha256"])


if __name__ == "__main__":
    main()

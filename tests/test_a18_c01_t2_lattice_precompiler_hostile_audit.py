from __future__ import annotations

import ast
from collections import Counter
from functools import lru_cache
import json

from scripts.audit_a18_c01_t2_lattice_precompiler_hostile import (
    OUTPUT,
    ROOT,
    audit,
)


@lru_cache(maxsize=1)
def _payload() -> dict:
    return audit()


def test_hostile_archive_replays_from_the_single_pinned_gram() -> None:
    payload = _payload()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    assert payload["sole_input"]["construction_id"] == "A18_4_0_C01"
    assert payload["sole_input"]["field"] == (
        "section.restricted_mirrors.facet_gram"
    )
    assert payload["sole_input"]["facet_gram_sha256"] == (
        "6979b16de7583c20277077a142df863e9f9a7209c443e9252505234b15d66c56"
    )


def test_complete_scalings_and_two_sheet_compatible_lanes() -> None:
    payload = _payload()
    scaling = payload["scaling_enumeration"]
    assert scaling["complete_ray_count"] == 18
    assert scaling["one_two_sheet_compatible_lane_count"] == 2
    assert sum(row["one_two_sheet_compatible"] for row in scaling["rays"]) == 2

    lanes = {row["label"]: row for row in payload["lanes"]}
    assert lanes["displayed"]["normal_norms_absolute"] == [
        2, 4, 4, 4, 2, 4, 4, 4, 2, 4, 4, 4, 4, 2
    ]
    assert lanes["fully_reversed"]["normal_norms_absolute"] == [
        4, 2, 2, 2, 4, 2, 2, 2, 4, 2, 2, 2, 2, 4
    ]


def test_all_stable_even_overlattices_and_obstructions() -> None:
    payload = _payload()
    lanes = {row["label"]: row for row in payload["lanes"]}
    displayed = lanes["displayed"]["overlattices"]
    reversed_lane = lanes["fully_reversed"]["overlattices"]
    assert len(displayed) == 16
    assert len(reversed_lane) == 5
    assert Counter(row["glue_dimension"] for row in displayed) == {
        0: 1,
        1: 7,
        2: 7,
        3: 1,
    }
    assert Counter(row["glue_dimension"] for row in reversed_lane) == {
        0: 1,
        1: 3,
        2: 1,
    }
    assert displayed[0]["discriminant_length"] == 10
    assert reversed_lane[0]["radical_long_wall_classes_one_based"] == [
        1,
        5,
        9,
        14,
    ]
    assert all(
        pair["folded_pairing"] == 0
        for row in displayed + reversed_lane
        for pair in row["repeated_long_class_pairs"]
    )
    survivors = [
        row
        for row in displayed + reversed_lane
        if all(row["gates"].values())
    ]
    assert len(survivors) == 1
    assert survivors[0]["determinant_absolute"] == 256
    assert survivors[0]["smith_diagonal_absolute"] == [1, 1] + [2] * 8


def test_explicit_u_e8_minus2_and_rank18_fold_back() -> None:
    candidate = _payload()["unique_candidate"]
    assert candidate["discriminant_invariants"] == [10, 8, 0]
    assert candidate["long_class_span_rank"] == 6
    fixed = candidate["fixed_lattice"]
    assert fixed["isometry_class"] == "U+E8(-2)"
    assert fixed["splitting_determinant_absolute"] == 1
    assert fixed["E8_minus_one_even_unimodular_negative_definite"] is True

    glue = candidate["long_class_glue"]
    assert glue["class_span_rank"] == 6
    assert glue["fixed_long_classes_as_bits"] == glue["anti_long_classes_as_bits"]
    assert glue["glue_index"] == 64

    rank18 = candidate["rank18"]
    assert rank18["abstract_2_elementary_invariants"] == [18, 4, 0]
    assert rank18["determinant"] == -16
    assert rank18["smith_diagonal_absolute"] == [1] * 14 + [2] * 4
    assert rank18["signature"] == [1, 17]
    assert rank18["involution_order"] == 2
    assert rank18["fixed_rank"] == 10 and rank18["anti_rank"] == 8
    assert rank18["anti_lattice"] == "E8(-2)"
    assert rank18["anti_lattice_rootless"] is True

    fold = candidate["fold_back"]
    assert fold["source_root_count"] == 24
    assert fold["source_roots_span_rank18_with_index_one"] is True
    assert fold["one_sheet_orbit_count"] == 4
    assert fold["two_sheet_orbit_count"] == 10
    assert fold["source_root_gram_pairing_histogram"] == {"0": 244, "1": 32}
    assert fold["folded_gram_equals_pinned_facet_gram"] is True


def test_scope_stays_ic0_and_no_production_or_promotion_import() -> None:
    payload = _payload()
    assert payload["scope"] == {
        "grade": "IC0 finite lattice precompiler and exact marked fold-back",
        "class_wide_IC2_or_IC4_claim": False,
        "new_polytope_claim": False,
        "full_nef_or_surface_automorphism_claim": False,
    }
    bridge = payload["marked_bridge_gate"]
    assert bridge["existing_C01_strong_promotion_identified"] is False
    assert "explicit integral marked isometry" in bridge[
        "required_for_datum_specific_IC2"
    ]

    script = ROOT / "scripts/audit_a18_c01_t2_lattice_precompiler_hostile.py"
    source = script.read_text(encoding="utf-8")
    assert "a18_c01_c07_k3_nef_promotions.json" not in source
    tree = ast.parse(source)
    imported: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.append(node.module)
    assert not any(name == "scripts" or name.startswith("scripts.") for name in imported)
    assert not any(
        name == "coxeter_polytopes" or name.startswith("coxeter_polytopes.")
        for name in imported
    )

from __future__ import annotations

from collections import Counter
import json

from scripts.inventory_shimada_oriented_weighted_wall_symmetries import (
    CONFIGS,
    _default_output,
)
from scripts.screen_shimada_oriented_weighted_high_sections import (
    _default_case_directory,
    _default_manifest,
    _inventory_cases,
    _sha256,
    audit_archive,
)
from scripts.shimada_disc11_atlas import ShimadaDisc11Atlas


DISCRIMINANT = 11
INVENTORY = _default_output(DISCRIMINANT)
CASE_DIRECTORY = _default_case_directory(DISCRIMINANT)
MANIFEST = _default_manifest(DISCRIMINANT)


def _load(path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_complete_oriented_weighted_wall_denominator() -> None:
    inventory = _load(INVENTORY)
    assert inventory["complete"] is True
    assert inventory["source"] == {
        "ambient_rank": 20,
        "chamber_count": 1098,
        "discriminant": 11,
        "gram_determinant": -11,
        "path": "/private/tmp/threediscs/disc11.txt",
        "sha256": "fa39e5b59f4a4b2b1b7fe7632e5737d1734ee3eb7e7023ae26bf2b83387afb68",
    }
    summary = inventory["summary"]
    assert summary["completed_chamber_count"] == 1098
    assert summary["pending_chamber_count"] == 0
    assert summary["archived_AutD_order_histogram"] == {
        "1": 1027,
        "2": 60,
        "3": 2,
        "4": 3,
        "6": 5,
        "20": 1,
    }
    assert summary["oriented_weighted_group_order_histogram"] == {
        "1": 814,
        "2": 270,
        "3": 2,
        "4": 5,
        "6": 5,
        "8": 1,
        "20": 1,
    }
    assert summary["weighted_group_strict_extension_count"] == 216
    assert summary["zero_published_Weyl_wall_orbit_count"] == 0
    assert summary["weighted_primal_action_denominator_histogram"] == {
        "1": 1215,
        "2": 223,
    }


def test_294_high_fixed_sections_have_exact_group_and_denominator_types() -> None:
    inventory = _load(INVENTORY)
    cases = _inventory_cases(inventory)
    assert len(cases) == 294
    assert Counter(record["hyperbolic_dimension"] for _, record in cases) == {
        10: 3,
        11: 71,
        12: 2,
        17: 2,
        18: 216,
    }
    detail = Counter()
    for chamber, record in cases:
        actions = {
            tuple(action["permutation"]): action
            for action in chamber["oriented_weighted_wall_group"]["all_actions"]
        }
        denominators = tuple(
            sorted(
                {
                    actions[tuple(element)]["primal_action_denominator_lcm"]
                    for element in record["elements"]
                }
            )
        )
        detail[
            (
                record["hyperbolic_dimension"],
                record["order"],
                record["representative_contained_in_archived_AutD"],
                denominators,
            )
        ] += 1
    assert detail == {
        (18, 2, False, (1, 2)): 216,
        (17, 2, False, (1, 2)): 1,
        (17, 4, False, (1, 2)): 1,
        (12, 2, False, (1, 2)): 2,
        (11, 2, True, (1,)): 69,
        (11, 4, False, (1, 2)): 2,
        (10, 4, False, (1, 2)): 2,
        (10, 8, False, (1, 2)): 1,
    }


def test_all_294_classes_have_a_certified_forbidden_facet_pair() -> None:
    manifest = _load(MANIFEST)
    summary = manifest["summary"]
    assert summary["inventory_case_count"] == 294
    assert summary["completed_case_count"] == 294
    assert summary["pending_case_count"] == 0
    assert summary["classification_histogram"] == {"not_coxeter": 294}
    assert summary["coxeter_candidate_count_requiring_volume_audit"] == 0
    assert summary["total_certified_non_coxeter_pair_count"] == 3557
    assert "0" not in summary["bad_pair_count_histogram"]
    assert summary["bad_pair_count_histogram"]["1"] == 3
    assert sum(summary["bad_pair_count_histogram"].values()) == 294
    assert summary["fully_resolved_case_count"] == 136
    assert summary["case_count_with_unresolved_additional_facets"] == 158
    assert manifest["survivors_requiring_volume_audit"] == []


def test_fixed_divisor_lattice_determinants_and_reflection_separation() -> None:
    manifest = _load(MANIFEST)
    determinant_by_dimension = {
        10: 5632,
        11: -2816,
        12: 2816,
        17: -176,
        18: 44,
    }
    counts = Counter()
    for row in manifest["cases"]:
        case = _load(CASE_DIRECTORY / f"{row['id']}.json")
        dimension = case["case"]["hyperbolic_dimension"]
        counts[dimension] += 1
        assert (
            case["saturated_fixed_divisor_lattice"]["determinant"]
            == determinant_by_dimension[dimension]
        )
        assert (
            case["restricted_mirror_geometry"]
            ["saturated_fixed_divisor_lattice_is_stable_under_all_active_reflections"]
            is False
        )
    assert counts == {10: 3, 11: 71, 12: 2, 17: 2, 18: 216}


def test_all_no_go_grams_and_one_full_rebuild_per_dimension() -> None:
    config = CONFIGS[DISCRIMINANT]
    inventory = _load(INVENTORY)
    manifest = audit_archive(
        ShimadaDisc11Atlas(config["source"]),
        DISCRIMINANT,
        config["source_sha256"],
        INVENTORY,
        _sha256(INVENTORY),
        _inventory_cases(inventory),
        CASE_DIRECTORY,
        MANIFEST,
    )
    assert manifest["summary"]["classification_histogram"] == {
        "not_coxeter": 294
    }

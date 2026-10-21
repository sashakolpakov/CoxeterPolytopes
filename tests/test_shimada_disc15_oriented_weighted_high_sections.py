from __future__ import annotations

from collections import Counter
import json

from scripts.inventory_shimada_oriented_weighted_wall_symmetries import (
    _default_output,
)
from scripts.screen_shimada_oriented_weighted_high_sections import (
    _default_case_directory,
    _default_manifest,
    _inventory_cases,
    _sha256,
    audit_archive,
)
from scripts.inventory_shimada_oriented_weighted_wall_symmetries import CONFIGS
from scripts.shimada_disc11_atlas import ShimadaDisc11Atlas


DISCRIMINANT = 15
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
        "chamber_count": 2051,
        "discriminant": 15,
        "gram_determinant": -15,
        "path": "/private/tmp/threediscs/disc15.txt",
        "sha256": "f2fdcf80bd980c363fbdb68957a85878284bdb5e7d75c4dac92604aefefc4081",
    }
    summary = inventory["summary"]
    assert summary["completed_chamber_count"] == 2051
    assert summary["pending_chamber_count"] == 0
    assert summary["archived_AutD_order_histogram"] == {
        "1": 1998,
        "2": 44,
        "3": 2,
        "4": 2,
        "6": 1,
        "8": 1,
        "12": 2,
        "24": 1,
    }
    assert summary["oriented_weighted_group_order_histogram"] == {
        "1": 1840,
        "2": 196,
        "3": 2,
        "4": 6,
        "8": 2,
        "12": 1,
        "16": 1,
        "24": 2,
        "48": 1,
    }
    assert summary["weighted_group_strict_extension_count"] == 171
    assert summary["zero_published_Weyl_wall_orbit_count"] == 114
    assert summary["weighted_primal_action_denominator_histogram"] == {
        "1": 2258,
        "2": 144,
    }


def test_high_fixed_sections_split_into_three_exact_arithmetic_families() -> None:
    inventory = _load(INVENTORY)
    cases = _inventory_cases(inventory)
    assert len(cases) == 236
    assert Counter(record["hyperbolic_dimension"] for _, record in cases) == {
        11: 59,
        13: 33,
        18: 144,
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
        (11, 2, True, (1,)): 59,
        (13, 2, False, (1,)): 33,
        (18, 2, False, (1, 2)): 144,
    }


def test_all_236_classes_have_a_certified_forbidden_facet_pair() -> None:
    manifest = _load(MANIFEST)
    summary = manifest["summary"]
    assert summary["inventory_case_count"] == 236
    assert summary["completed_case_count"] == 236
    assert summary["pending_case_count"] == 0
    assert summary["classification_histogram"] == {"not_coxeter": 236}
    assert summary["coxeter_candidate_count_requiring_volume_audit"] == 0
    assert "0" not in summary["bad_pair_count_histogram"]
    assert summary["bad_pair_count_histogram"]["1"] == 12
    assert sum(summary["bad_pair_count_histogram"].values()) == 236
    assert manifest["survivors_requiring_volume_audit"] == []
    assert all(
        row["status"] == "complete"
        and row["classification"] == "not_coxeter"
        and row["non_coxeter_pair_count"] >= 1
        for row in manifest["cases"]
    )


def test_fixed_divisor_lattices_and_mirror_lattices_remain_distinct() -> None:
    manifest = _load(MANIFEST)
    determinant_by_dimension = {11: -3840, 13: -320, 18: 60}
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
    assert counts == {11: 59, 13: 33, 18: 144}


def test_all_no_go_grams_and_selected_full_arithmetic_rebuilds() -> None:
    config = CONFIGS[DISCRIMINANT]
    inventory = _load(INVENTORY)
    atlas = ShimadaDisc11Atlas(config["source"])
    manifest = audit_archive(
        atlas,
        DISCRIMINANT,
        config["source_sha256"],
        INVENTORY,
        _sha256(INVENTORY),
        _inventory_cases(inventory),
        CASE_DIRECTORY,
        MANIFEST,
    )
    assert manifest["summary"]["classification_histogram"] == {
        "not_coxeter": 236
    }

from __future__ import annotations

import ast
from collections import Counter
from functools import lru_cache
import json
from pathlib import Path

from scripts.audit_jacobian_kummer_h12_doubling_cube_hostile_audit import (
    EXPECTED_BY_MASK,
    EXPECTED_FACET_HISTOGRAM,
    OUTPUT,
    ROOT,
    build_certificate,
)


@lru_cache(maxsize=1)
def _rebuilt() -> dict:
    return build_certificate()


def test_independent_rebuild_matches_frozen_hostile_archive() -> None:
    rebuilt = _rebuilt()
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert rebuilt == archived
    assert rebuilt["all_checks_pass"]
    assert all(rebuilt["checks"].values())
    assert not rebuilt["independence"]["production_cube_code_imported_or_invoked"]
    assert not rebuilt["independence"]["production_cube_archive_used_as_input"]


def test_c2_four_action_and_all_subset_boundaries() -> None:
    payload = _rebuilt()
    seed = payload["seed"]
    action = payload["C2_4_doubling_action"]
    assert seed["selected_doubling_walls_one_based"] == [31, 58, 59, 68]
    assert seed["selected_root_gram"] == [
        [-4, 0, 0, 0],
        [0, -4, 0, 0],
        [0, 0, -4, 0],
        [0, 0, 0, -4],
    ]
    assert action["abstract_group"] == "C2^4"
    assert action["order"] == 16
    assert action["checks"]["C2_4_action_is_faithful"]
    assert action["checks"]["all_seed_Coxeter_relations_respect_C2_4_quotient"]

    descendants = payload["descendants"]
    assert len(descendants) == 15
    recovered = {
        int(row["subset_mask_binary"], 2): row["facet_count"]
        for row in descendants
    }
    assert recovered == EXPECTED_BY_MASK
    assert Counter(row["facet_count"] for row in descendants) == Counter(
        EXPECTED_FACET_HISTOGRAM
    )
    for row in descendants:
        facets = row["facet_count"]
        assert row["cube_chamber_count"] == 2 ** row["selected_wall_count"]
        assert row["integral_reflection_count"] == facets
        assert row["strict_facet_witnesses"]["count"] == facets
        assert row["strict_facet_witnesses"]["minimum_off_target_slack"] == 1
        assert row["exact_unordered_facet_pair_count"] == facets * (facets - 1) // 2
        assert sum(row["relation_histogram"].values()) == facets * (facets - 1) // 2
        assert row["finite_volume"] and not row["compact"]
        assert row["all_exact_checks_pass"]
        assert all(row["checks"].values())


def test_lower_bound_and_bounded_lineage_scope() -> None:
    payload = _rebuilt()
    census = payload["census"]
    lineage = payload["bounded_lineage"]
    assert census["certified_pairwise_nonisomorphic_lower_bound"] == 7
    assert census["additional_facet_count_types"] == [523, 559, 941, 1004, 1827]
    assert census["additional_types_beyond_singles_lower_bound"] == 5
    assert lineage["ordinary_record_count"] == 1214
    assert lineage["doubling_record_count"] == 2959
    assert all(row["absent_from_both_pinned_lanes"] for row in lineage["candidate_headers"])
    assert lineage["worldwide_priority"] == "open"
    assert "not a publication-complete denominator" in lineage["scope_warning"]


def test_hostile_script_has_no_production_module_dependency() -> None:
    script = (
        ROOT
        / "scripts/audit_jacobian_kummer_h12_doubling_cube_hostile_audit.py"
    )
    tree = ast.parse(script.read_text(encoding="utf-8"))
    imported = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.append(node.module)
    forbidden = "scripts.certify_jacobian_kummer_h12_doubling_cube"
    assert forbidden not in imported


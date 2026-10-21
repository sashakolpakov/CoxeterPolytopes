from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import pytest

from scripts.independent_h9_14_record30_a0_target_labels import (
    _finite_code,
    _lanner_codes,
    replay,
)


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "data/enumerations/h9_14_record30_a0_candidate.json"


def _induced_code(
    code: tuple[int, ...], vertices: tuple[int, ...]
) -> tuple[int, ...]:
    labels = dict(zip(combinations(range(5), 2), code, strict=True))
    return tuple(labels[pair] for pair in combinations(vertices, 2))


def test_complete_lanner_table_is_minimally_nonfinite() -> None:
    placements = _lanner_codes()
    assert len(placements) == 420
    for code in placements:
        assert not _finite_code(5, code)
        for omitted in range(5):
            kept = tuple(vertex for vertex in range(5) if vertex != omitted)
            assert _finite_code(4, _induced_code(code, kept))


def test_candidate_fixture_replays_exact_literal_and_label_counts() -> None:
    document = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    result = replay(document)
    assert result["compact_only"] is True
    assert result["imports_coxeter_polytopes"] is False
    assert result["complete_order_five_lanner_placement_count"] == 420
    assert result["target_orbit_count"] == 1
    assert result["terminal_sparse_state_count_histogram"] == {"0": 1}
    assert result["unrestricted_compact_label_no_go_target_count"] == 1
    assert result["sparse_label_survivor_count"] == 0
    assert result["finite_volume_noncompact_chamber_excluded"] is False
    audit = result["target_audits"][0]
    assert audit["minimal_transversal_histogram"] == {"5": 132}
    assert audit["actual_rank_nine_vertex_count"] == 132
    assert audit["facet_complement_sha256"] == (
        "878d579b2ef3045c8948828bac4912450598fa44f1635b2743d8a27e58e0d8e6"
    )
    assert [stage["passing_state_count"] for stage in audit["stages"]] == [
        420,
        96,
        6,
        0,
    ]
    assert audit["processed_support_indices"] == [0, 1, 2, 7]
    assert audit["unrestricted_compact_label_no_go"] is True


def test_connected_ridge_census_input_shape_is_consumed() -> None:
    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    census_shape = {
        "kind": "h9_14_record30_a0_connected_ridge_census",
        "target_orbit_count": 1,
        "target_orbit_representatives": [candidate["missing_A5s"]],
    }
    assert replay(census_shape)["target_audits"] == replay(candidate)[
        "target_audits"
    ]


def test_nontrivial_fixed_matching_relabeling_preserves_no_go() -> None:
    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    permutation = list(range(14))
    permutation[0], permutation[1] = 1, 0
    permutation[11], permutation[12] = 12, 11
    for source, target in zip(range(3, 11), (*range(4, 11), 3), strict=True):
        permutation[source] = target
    moved = {
        "target_orbit_representatives": [
            [
                sorted(permutation[vertex] for vertex in support)
                for support in candidate["missing_A5s"]
            ]
        ]
    }
    audit = replay(moved)["target_audits"][0]
    assert audit["actual_rank_nine_vertex_count"] == 132
    assert [stage["passing_state_count"] for stage in audit["stages"]] == [
        420,
        96,
        6,
        0,
    ]
    assert audit["unrestricted_compact_label_no_go"] is True


def test_two_distinct_132_vertex_support_fixtures_reach_zero() -> None:
    # These are deliberately not relabellings of the frozen multiwedge row.
    # They exercise two different dense common-triple overlap patterns.
    families = [
        [59, 61, 122, 124, 248, 1984, 3968, 14080],
        [59, 121, 122, 248, 1984, 3968, 6016, 7936],
    ]
    result = replay(families)
    assert result["target_orbit_count"] == 2
    assert result["terminal_sparse_state_count_histogram"] == {"0": 2}
    assert [
        [stage["passing_state_count"] for stage in audit["stages"]]
        for audit in result["target_audits"]
    ] == [[420, 96, 6, 0], [420, 126, 18, 0]]
    assert [
        audit["facet_complement_sha256"] for audit in result["target_audits"]
    ] == [
        "260f449de5d478463b9a0dc2d251fc8d10e319a27f716b01da988d1872fa49b7",
        "140ff6ec52d4338628c5e8557d25ba1e965a9c52a036b9adc9cda7fef9bfcbe8",
    ]


def test_malformed_target_is_rejected_before_label_search() -> None:
    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    malformed = {"target_orbit_representatives": [candidate["missing_A5s"][:7]]}
    with pytest.raises(ValueError, match="exactly eight"):
        replay(malformed)

    bad_count = {
        "target_orbit_count": 2,
        "target_orbit_representatives": [candidate["missing_A5s"]],
    }
    with pytest.raises(ValueError, match="declared target_orbit_count"):
        replay(bad_count)


def test_empty_census_does_not_emit_a_vacuous_label_theorem() -> None:
    result = replay(
        {"target_orbit_count": 0, "target_orbit_representatives": []}
    )
    assert result["target_orbit_count"] == 0
    assert result["unrestricted_compact_label_no_go_target_count"] == 0
    assert result[
        "all_nonempty_input_targets_have_unrestricted_compact_label_no_go"
    ] is False


def test_source_has_no_project_import() -> None:
    source = (
        ROOT / "scripts/independent_h9_14_record30_a0_target_labels.py"
    ).read_text(encoding="utf-8")
    assert "from coxeter_polytopes" not in source
    assert "import coxeter_polytopes" not in source

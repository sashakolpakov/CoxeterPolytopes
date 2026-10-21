"""Frozen checks for the independent hostile golden four-node replay."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "data/enumerations/golden_agcat_four_node_at_most_one_dotted_hostile_audit.json"
)
SOURCE = (
    ROOT
    / "scripts/independent_hostile_audit_golden_agcat_four_node_at_most_one_dotted.py"
)


def payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def test_independent_orbits_counts_and_digests() -> None:
    frozen = payload()
    replay = frozen["replay"]
    assert frozen["status"] == "PASS"
    assert frozen["verdict"]["four_node_local_tranche"] == "INDEPENDENT HOSTILE PASS"
    assert replay["alphabet"] == {
        "full_count": 6_875,
        "finite_count": 6,
        "dotted_count": 6_869,
        "dotted_b_major_sha256": (
            "75e964fae1b1f53981ea19a5e494c7924ebe34681d3e7d89cb74719c3ace542f"
        ),
    }
    assert replay["canonical_orbits"] == {
        "finite_only_connected": 2_380,
        "one_dotted_complements": 2_080,
        "one_dotted_expanded": 14_287_520,
    }
    counts = replay["counts"]
    assert (counts["finite_compatible"], counts["dotted_compatible"]) == (327, 532_734)
    assert (counts["finite_retained"], counts["dotted_retained"]) == (119, 180_402)
    assert (counts["finite_positive"], counts["finite_lorentzian"]) == (5, 114)
    assert counts["dotted_lorentzian"] == 180_402
    assert counts["dotted_nonlorentzian"] == 0
    assert replay["stream_digests"] == {
        "finite_only_cycle_compatible_sha256": (
            "687762835339ac1f56fed84bcc41975e51da471abba4d97a7ab4cafc64590e7a"
        ),
        "finite_only_two_place_retained_sha256": (
            "a8a04a852ce054ebfa1668156200998263939d7417d7c437662130b2320c50e8"
        ),
        "one_dotted_cycle_compatible_sha256": (
            "0c988e10b93d804457dcfb9f238a832c546fd749b0a66751e67c875462153cc5"
        ),
        "one_dotted_two_place_retained_sha256": (
            "7de03d51631e744ba34c248a23a933da855e841256f7f27a33c1cf09b8a0b0a1"
        ),
    }


def test_zero_split_strata_and_scope_are_explicit() -> None:
    frozen = payload()
    replay = frozen["replay"]
    counts = replay["counts"]
    assert (
        counts["finite_selected_leading_zero"],
        counts["finite_selected_full_zero"],
        counts["finite_conjugate_leading_zero"],
        counts["finite_conjugate_full_zero"],
    ) == (21, 3, 21, 3)
    assert (
        counts["dotted_selected_leading_zero"],
        counts["dotted_selected_full_zero"],
        counts["dotted_conjugate_leading_zero"],
        counts["dotted_conjugate_full_zero"],
    ) == (0, 0, 0, 0)
    assert len(replay["strata"]) == 48
    assert sum(row["count"] for row in replay["strata"]) == 327 + 532_734
    assert replay["mathematical_audit"]["local_scope_only"] is True
    assert frozen["verdict"]["new_compact_chamber"] is False
    assert frozen["verdict"]["golden_H22_H29_nonexistence"] == "OPEN"
    assert frozen["defect_history"]["repair_verified"].endswith(
        "the one-dotted split is 0/0"
    )


def test_global_path_family_and_source_independence() -> None:
    frozen = payload()
    rows = frozen["replay"]["global_path_rows"]
    assert [row["native_rank"] for row in rows] == list(range(23, 31))
    assert [row["admissible_dotted_value_count"] for row in rows] == [
        3_590, 3_584, 3_577, 3_571, 3_566, 3_561, 3_556, 3_552
    ]
    assert all(row["sample_dotted_product"] == [2, 2] for row in rows)
    assert all(
        row["sample_full_determinant"] == [2, -2 * (row["native_rank"] - 1)]
        for row in rows
    )
    source = SOURCE.read_text()
    assert "from coxeter_polytopes.classification.golden_four_node_packets" not in source
    assert "certify_golden_agcat_four_node_at_most_one_dotted" not in source
    assert "golden_triangle_packets" not in source
    assert "def element_squareclass" in source
    assert "def gram_determinant" in source
    assert "for local_permutation in permutations" in source

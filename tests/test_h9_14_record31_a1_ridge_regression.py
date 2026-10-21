from __future__ import annotations

from collections import Counter
import json
from math import comb
from pathlib import Path

from scripts.explore_h9_14_record31_a1_downstream import audit
from scripts.generate_h9_14_record31_a1_addition_engine import generate


ROOT = Path(__file__).resolve().parents[1]
INVALID = (
    ROOT
    / "data/enumerations/h9_14_record31_a1_invalid_ridge531_transcript.json"
)
CORRECTED = (
    ROOT / "data/enumerations/h9_14_record31_a1_orbit0_corrected.json"
)


def _independent(face: int) -> bool:
    edges = ((0, 11), (1, 12), (0, 13), (1, 13))
    return all(face & ((1 << left) | (1 << right)) != ((1 << left) | (1 << right))
               for left, right in edges)


def _a1_dehn_sommerville_vector() -> tuple[int, ...]:
    dimension = 9
    independent_counts = Counter(
        face.bit_count()
        for face in range(1, 1 << 14)
        if face.bit_count() <= 4 and _independent(face)
    )
    # f_{-1}, f_0, f_1, f_2, f_3.  The a=1 schedule removes one
    # independent four-set as its unique minimal missing A4.
    prefix = (
        1,
        independent_counts[1],
        independent_counts[2],
        independent_counts[3],
        independent_counts[4] - 1,
    )
    first_h = tuple(
        sum(
            (-1) ** (index - order)
            * comb(dimension - order, dimension - index)
            * prefix[order]
            for order in range(index + 1)
        )
        for index in range(5)
    )
    h = first_h + first_h[::-1]
    return tuple(
        sum(
            comb(dimension - index, dimension - order) * h[index]
            for index in range(order + 1)
        )
        for order in range(1, dimension + 1)
    )


def test_record31_a1_ridge_count_is_derived_as_540() -> None:
    f_vector = _a1_dehn_sommerville_vector()
    assert f_vector == (14, 87, 319, 770, 1281, 1477, 1142, 540, 120)
    assert f_vector[-2] == 540


def test_invalid_531_gate_transcript_preserves_the_caught_survivor() -> None:
    document = json.loads(INVALID.read_text(encoding="utf-8"))
    row = document["orbits"][0]
    assert row["pure_target_ridge_histogram"]["2:540"] == 1
    assert row["ridge_valid_target_orbit_count"] == 0
    assert row["ridge_valid_target_minimal_nonfaces"] == []


def test_generated_engine_uses_dehn_sommerville_ridge_derivation() -> None:
    source = generate()
    assert "DerivedScheduledRidgeCount" in source
    assert "face_prefix[4] -= kScheduledFours" in source
    assert "{{2, 531}}" not in source


def test_corrected_orbit0_has_one_target_and_zero_compact_label_states() -> None:
    document = json.loads(CORRECTED.read_text(encoding="utf-8"))
    row = document["orbits"][0]
    assert row["ridge_valid_pure_target_count"] == 1
    assert row["ridge_valid_target_orbit_count"] == 1
    result = audit(document)
    assert result["input_target_orbit_count"] == 1
    assert result["connected_mod2_homology_sphere_orbit_count"] == 1
    assert result["compact_label_survivor_orbit_count"] == 0
    target = result["target_orbits"][0]
    assert target["global_topology_audit"]["face_vector"][-2:] == [540, 120]
    assert target["sparse_compact_label_gate"]["terminal_sparse_state_count"] == 0

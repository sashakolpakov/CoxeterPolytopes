from __future__ import annotations

import json

from scripts.generate_h9_14_record31_a1_addition_engine import generate
from scripts.generate_h9_14_record31_a1_compact_key_engine import generate_compact
from scripts.verify_h9_14_record31_a1_compact_key_equivalence import (
    OUTPUT,
    verify_equivalence,
)


def test_compact_generator_changes_only_the_canonical_key_serialization_site() -> None:
    dense = generate()
    compact = generate_compact()
    assert dense != compact
    assert "compact canonical key failed dense-graph round trip" in compact
    assert "incidence_bits = row_count * kN" in compact
    assert "canonical.size() * sizeof(graph)" in dense
    assert "canonical.size() * sizeof(graph)" not in compact


def test_bounded_dense_and_compact_quotient_replays_are_exactly_equal() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    rebuilt = verify_equivalence()
    assert rebuilt == archived
    assert archived["dense_and_compact_transcripts_identical"] is True
    assert archived["compact_per_state_dense_round_trip_enabled"] is True
    assert archived["bounded_run"] == {
        "partition_depth_node_counts": [1, 56, 420, 8495, 195250, 27551],
        "positive_A4_mask": 60,
        "quotient_seen_counts": [1, 4, 47, 1393, 51035, 8396],
        "search_node_count": 231773,
        "source_orbit_index": 0,
        "target_orbit_count": 1,
    }

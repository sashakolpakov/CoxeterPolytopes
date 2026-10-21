from __future__ import annotations

from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path

import scripts.certify_jacobian_kummer_h12_five_vertex_fork_b_le_2_complete_census as combined
import scripts.certify_jacobian_kummer_h12_five_vertex_fork_streaming_census as fork
import scripts.certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census as p5


ROOT = Path(__file__).resolve().parents[1]
DENOMINATOR = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_pattern_denominator.json"
LOW_SUMMARY = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_b0_to_b1_streaming_census.json"
LOW_ROWS = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_b0_to_b1_streaming_census_rows.jsonl.gz"
B2_SUMMARY = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_b2_to_b2_streaming_census.json"
B2_ROWS = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_b2_to_b2_streaming_census_rows.jsonl.gz"
COMBINED = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_b_le_2_complete_census.json"


def _rows(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        return [json.loads(line) for line in stream]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_exact_twelve_pattern_denominator_partition() -> None:
    payload = json.loads(DENOMINATOR.read_text())
    assert payload["all_checks_pass"]
    assert len(payload["patterns"]) == 12
    assert payload["total_raw"] == 75_897_360
    assert payload["total_orbits"] == 199_362
    assert payload["orbit_count_by_B"] == {
        "0": 116,
        "1": 2_042,
        "2": 18_061,
        "3": 72_758,
        "4": 106_385,
    }
    for row in payload["patterns"].values():
        assert row["Burnside_fixed_configuration_sum"] // row["Burnside_denominator"] == row["orbit_count"]


def test_compact_transcripts_cover_and_classify_all_20219_rows() -> None:
    low_summary = json.loads(LOW_SUMMARY.read_text())
    b2_summary = json.loads(B2_SUMMARY.read_text())
    low = _rows(LOW_ROWS)
    b2 = _rows(B2_ROWS)
    assert len(low) == low_summary["compact_exact_transcript"]["row_count"] == 2_158
    assert len(b2) == b2_summary["compact_exact_transcript"]["row_count"] == 18_061
    assert _sha256(LOW_ROWS) == low_summary["compact_exact_transcript"]["deterministic_gzip_sha256"]
    assert _sha256(B2_ROWS) == b2_summary["compact_exact_transcript"]["deterministic_gzip_sha256"]
    rows = low + b2
    assert Counter(row["facet_count"] for row in rows) == Counter(
        {752: 116, 761: 2_042, 770: 18_061}
    )
    exact_keys = {
        (
            row["facet_count"],
            row["ordinary_degree_profile_multiset_sha256"],
            row.get("color_refinement_invariant_sha256", ""),
        )
        for row in rows
    }
    assert len(exact_keys) == 20_219
    assert sum(row["color_refinement_applied"] for row in low) == 6
    assert sum(row["color_refinement_applied"] for row in b2) == 914


def test_one_exact_boundary_per_completed_fork_pattern() -> None:
    rows = _rows(LOW_ROWS) + _rows(B2_ROWS)
    representatives = {}
    for row in rows:
        representatives.setdefault(row["fork_class_pattern"], row)
    assert set(representatives) == {
        "AA|AA",
        "BA|AA",
        "AB|AA",
        "AA|AB",
        "BB|AA",
        "BA|AB",
        "AB|AB",
        "AA|BB",
    }
    _carrier, gram, roots, interior, _bounded = p5.base._load_inputs()
    engine = fork.ExactForkBoundaryEngine(gram, roots, interior)
    for pattern, row in representatives.items():
        colors = tuple(
            value - p5.STOCK_START - 1 for value in row["edge_colors_one_based"]
        )
        boundary, data = engine.boundary(colors)
        assert len(boundary) == data["facet_count"] == 752 + 9 * fork._pattern_b_count(pattern)
        assert data["internal_mirror_count"] == 4
        assert data["minimum_strict_interior_facet_pairing"] > 0
        assert data["fixed_boundary_core_Smith_diagonal"] == [1] * p5.RANK


def test_combined_theorem_rebuild_and_claim_boundary() -> None:
    archived = json.loads(COMBINED.read_text())
    assert archived == combined.build()
    theorem = archived["theorem"]
    assert theorem["ordinary_Coxeter_type_count"] == 20_219
    assert theorem["type_facet_histogram"] == {"752": 116, "761": 2_042, "770": 18_061}
    assert theorem["finite_volume"] and not theorem["compact"]
    assert archived["strong_exact_five_vertex_union"]["pairwise_nonisomorphic_type_count"] == 237_068
    assert archived["full_fork_denominator"]["remaining_unclassified_orbit_count"] == 179_143
    assert archived["claim_boundary"]["star_vs_P5_B_le_2_separation_open"]
    assert not archived["claim_boundary"]["complete_fork_census_claimed"]

from __future__ import annotations

from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path

import scripts.certify_jacobian_kummer_h12_five_vertex_fork_b3_terminal_acceptance as accepted
import scripts.certify_jacobian_kummer_h12_five_vertex_fork_b_le_3_complete_census as combined
import scripts.certify_jacobian_kummer_h12_five_vertex_fork_streaming_census as fork
import scripts.certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census as p5


ROOT = Path(__file__).resolve().parents[1]
B3_SUMMARY = accepted.B3_SUMMARY
B3_ROWS = accepted.B3_ROWS
ACCEPTED = accepted.OUTPUT
COMBINED = combined.OUTPUT


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_terminal_acceptance_provenance_and_hostile_repairs() -> None:
    payload = json.loads(ACCEPTED.read_text())
    assert payload["all_checks_pass"]
    assert payload["production_provenance"]["exact_launch_time_code_sha256"] == (
        accepted.PINNED_CODE_SHA256
    )
    assert payload["production_provenance"][
        "exact_input_and_output_archive_sha256"
    ] == accepted.PINNED_ARCHIVE_SHA256
    alignment = payload["iterator_alignment_replay"]
    assert alignment["production_iterator_equals_archive_authenticated_iterator"]
    assert alignment["all_archived_colours_equal_regenerated_colours"]
    assert alignment["source_row_counts"] == {
        "K1,4_star": 6_117,
        "P5_ABBB": 36_247,
        "P5_BABB": 36_648,
    }
    assert alignment["total_same_facet_row_count"] == 79_012
    audit = payload["terminal_transcript_audit"]
    assert audit["color_refined_current_member_count"] == 6_569
    assert audit["corrected_color_refined_discrete_current_count"] == 6_446
    assert audit["production_summary_underreported_discrete_count_as"] == 3_975
    assert audit["diagnostic_only_correction_no_theorem_change"]


def test_compact_transcript_exact_counts_hashes_and_keys() -> None:
    summary = json.loads(B3_SUMMARY.read_text())
    count = refined = discrete = 0
    pattern_histogram: Counter[str] = Counter()
    facet_histogram: Counter[int] = Counter()
    exact_keys: set[tuple[int, str, str]] = set()
    transcript = hashlib.sha256()
    with gzip.open(B3_ROWS, "rb") as stream:
        for line in stream:
            transcript.update(line)
            row = json.loads(line)
            count += 1
            assert row["orbit_number"] == count
            colors = tuple(value - p5.STOCK_START - 1 for value in row["edge_colors_one_based"])
            assert sum(color >= p5.A_SIZE for color in colors) == 3
            pattern_histogram[row["fork_class_pattern"]] += 1
            facet_histogram[row["facet_count"]] += 1
            refined += int(row["color_refinement_applied"])
            discrete += int(row.get("color_refinement_discrete_canonical_form", False))
            exact_keys.add(
                (
                    row["facet_count"],
                    row["ordinary_degree_profile_multiset_sha256"],
                    row.get("color_refinement_invariant_sha256", ""),
                )
            )
    assert count == len(exact_keys) == 72_758
    assert pattern_histogram == Counter(
        {"BB|AB": 36_247, "BA|BB": 18_055, "AB|BB": 18_456}
    )
    assert facet_histogram == Counter({779: 72_758})
    assert refined == 6_569 and discrete == 6_446
    assert transcript.hexdigest() == summary["compact_exact_transcript"][
        "uncompressed_canonical_JSONL_sha256"
    ]
    assert _sha256(B3_ROWS) == summary["compact_exact_transcript"][
        "deterministic_gzip_sha256"
    ]


def test_one_exact_boundary_per_B3_fork_pattern() -> None:
    representatives: dict[str, dict] = {}
    with gzip.open(B3_ROWS, "rt", encoding="utf-8") as stream:
        for line in stream:
            row = json.loads(line)
            representatives.setdefault(row["fork_class_pattern"], row)
            if len(representatives) == 3:
                break
    assert set(representatives) == {"BB|AB", "BA|BB", "AB|BB"}
    _carrier, gram, roots, interior, _bounded = p5.base._load_inputs()
    engine = fork.ExactForkBoundaryEngine(gram, roots, interior)
    for row in representatives.values():
        colors = tuple(
            value - p5.STOCK_START - 1 for value in row["edge_colors_one_based"]
        )
        boundary, data = engine.boundary(colors)
        assert len(boundary) == data["facet_count"] == 779
        assert data["internal_mirror_count"] == 4
        assert data["minimum_strict_interior_facet_pairing"] > 0
        assert data["fixed_boundary_core_Smith_diagonal"] == [1] * p5.RANK
        assert data["marked_boundary_root_sha256"] == row["marked_boundary_root_sha256"]


def test_combined_B_le_3_theorem_rebuild_and_claim_boundary() -> None:
    archived = json.loads(COMBINED.read_text())
    assert archived == combined.build()
    theorem = archived["theorem"]
    assert theorem["ordinary_Coxeter_type_count"] == 92_977
    assert theorem["type_facet_histogram"] == {
        "752": 116,
        "761": 2_042,
        "770": 18_061,
        "779": 72_758,
    }
    assert theorem["finite_volume"] and not theorem["compact"]
    assert archived["strong_exact_five_vertex_union"][
        "pairwise_nonisomorphic_type_count"
    ] == 309_826
    assert archived["full_fork_denominator"][
        "remaining_unclassified_orbit_count"
    ] == 106_385
    claim = archived["claim_boundary"]
    assert claim["fork_B4_matrix_census_open_type_denominator"] == 106_385
    assert claim["omitted_star_B_le_2_row_count"] == 1_770
    assert claim["star_vs_P5_B_le_2_separation_open"]
    assert not claim["complete_fork_census_claimed"]
    assert claim["worldwide_priority"] == "open"
    assert not claim["worldwide_novelty_claim"]
    assert claim["weak_full_Picard_K3_carrier_not_full_nef_or_surface_action"]

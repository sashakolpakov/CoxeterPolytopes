from __future__ import annotations

from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path

import scripts.certify_jacobian_kummer_h12_five_vertex_fork_b4_terminal_acceptance as accepted
import scripts.certify_jacobian_kummer_h12_five_vertex_three_shape_complete_census as combined
import scripts.certify_jacobian_kummer_h12_five_vertex_fork_streaming_census as fork
import scripts.certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census as p5


ROOT = Path(__file__).resolve().parents[1]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def test_B4_terminal_acceptance_provenance_and_hostile_repairs() -> None:
    payload = json.loads(accepted.OUTPUT.read_text())
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
        "K1,4_star": 8_985,
        "P5_BBBB": 108_159,
    }
    assert alignment["total_same_facet_row_count"] == 117_144
    audit = payload["terminal_transcript_audit"]
    assert audit["color_refined_current_member_count"] == 14_096
    assert audit["corrected_color_refined_discrete_current_count"] == 13_945
    assert audit["production_summary_underreported_discrete_count_as"] == 9_300
    assert audit["diagnostic_only_correction_no_theorem_change"]


def test_B4_compact_transcript_exact_counts_hashes_and_keys() -> None:
    summary = json.loads(accepted.B4_SUMMARY.read_text())
    count = refined = discrete = 0
    patterns: Counter[str] = Counter()
    facets: Counter[int] = Counter()
    exact_keys: set[tuple[int, str, str]] = set()
    transcript = hashlib.sha256()
    with gzip.open(accepted.B4_ROWS, "rb") as stream:
        for line in stream:
            transcript.update(line)
            row = json.loads(line)
            count += 1
            assert row["orbit_number"] == count
            colors = tuple(
                value - p5.STOCK_START - 1 for value in row["edge_colors_one_based"]
            )
            assert sum(color >= p5.A_SIZE for color in colors) == 4
            patterns[row["fork_class_pattern"]] += 1
            facets[row["facet_count"]] += 1
            refined += int(row["color_refinement_applied"])
            discrete += int(
                row.get("color_refinement_discrete_canonical_form", False)
            )
            exact_keys.add(
                (
                    row["facet_count"],
                    row["ordinary_degree_profile_multiset_sha256"],
                    row.get("color_refinement_invariant_sha256", ""),
                )
            )
    assert count == len(exact_keys) == 106_385
    assert patterns == Counter({"BB|BB": 106_385})
    assert facets == Counter({788: 106_385})
    assert refined == 14_096 and discrete == 13_945
    assert transcript.hexdigest() == summary["compact_exact_transcript"][
        "uncompressed_canonical_JSONL_sha256"
    ]
    assert _sha256(accepted.B4_ROWS) == summary["compact_exact_transcript"][
        "deterministic_gzip_sha256"
    ]


def test_one_exact_B4_fork_boundary() -> None:
    with gzip.open(accepted.B4_ROWS, "rt", encoding="utf-8") as stream:
        row = json.loads(next(stream))
    _carrier, gram, roots, interior, _bounded = p5.base._load_inputs()
    engine = fork.ExactForkBoundaryEngine(gram, roots, interior)
    colors = tuple(
        value - p5.STOCK_START - 1 for value in row["edge_colors_one_based"]
    )
    boundary, data = engine.boundary(colors)
    assert len(boundary) == data["facet_count"] == 788
    assert data["internal_mirror_count"] == 4
    assert data["minimum_strict_interior_facet_pairing"] > 0
    assert data["fixed_boundary_core_Smith_diagonal"] == [1] * p5.RANK
    assert data["marked_boundary_root_sha256"] == row["marked_boundary_root_sha256"]


def test_complete_fork_and_three_shape_theorem() -> None:
    archived = json.loads(combined.OUTPUT.read_text())
    assert archived == combined.build()
    theorem = archived["theorem"]
    assert theorem["ordinary_Coxeter_type_count"] == 417_981
    assert theorem["pairwise_nonisomorphic"]
    assert theorem["shape_type_counts"] == {
        "path_P5": 201_747,
        "star_K1,4": 16_872,
        "degree_3_fork": 199_362,
    }
    assert theorem["type_facet_histogram"] == {
        "752": 277,
        "761": 4_302,
        "770": 38_103,
        "779": 151_770,
        "788": 223_529,
    }
    assert theorem["finite_volume"] and not theorem["compact"]
    complete_fork = archived["complete_fork_theorem"]
    assert complete_fork["ordinary_Coxeter_type_count"] == 199_362
    assert complete_fork["remaining_unclassified_orbit_count"] == 0
    assert complete_fork["type_count_by_B"] == {
        "0": 116,
        "1": 2_042,
        "2": 18_061,
        "3": 72_758,
        "4": 106_385,
    }
    claim = archived["claim_boundary"]
    assert claim["complete_fork_census_claimed"]
    assert claim["all_three_five_vertex_tree_shapes_complete"]
    assert claim["all_417981_five_vertex_denominator_rows_classified"]
    assert claim["worldwide_priority"] == "open"
    assert not claim["worldwide_novelty_claim"]


def test_complete_collision_accounting() -> None:
    archived = json.loads(combined.OUTPUT.read_text())
    collision = archived["exact_cross_shape_accounting"]
    assert collision["complete_P5_plus_star_type_count"] == 218_619
    assert collision["complete_P5_plus_star_plus_fork_type_count"] == 417_981
    assert collision["fork_internal_literal_degree_collision_group_count"] == 6_844
    assert collision["fork_internal_literal_degree_collision_member_count"] == 14_090
    assert collision["fork_cross_P5_star_literal_degree_match_count"] == 8_525
    assert collision["unresolved_pair_count"] == 0
    assert collision["all_equality_candidates_replayed_literally"]
    assert collision["all_literal_collisions_separated_by_exact_color_refinement"]

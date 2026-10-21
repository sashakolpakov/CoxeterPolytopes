from __future__ import annotations

import json

from scripts.certify_golden_lanner_exhaustive_census import OUTPUT, analyze


def _archive() -> dict[str, object]:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_exact_replay_matches_frozen_census() -> None:
    payload = analyze()
    assert payload == _archive()
    assert payload["all_checks_pass"] is True
    assert payload["summary"] == {
        "orbit_counts_by_dimension": {"2": 4, "3": 5, "4": 3},
        "total_orbits": 12,
        "previously_handpicked": 8,
        "newly_added": 4,
        "all_complete_norm_two_chambers": True,
        "all_compact": True,
        "all_stably_embed_in_K3": True,
        "all_return_the_same_polytope_as_eigen_nef_section": True,
        "subdivision_count": 0,
        "worldwide_new_abstract_polytope_claim": False,
    }


def test_unlabelled_diagram_enumeration_is_exhaustive_in_ranks_three_to_five() -> None:
    ranks = _archive()["finite_enumeration"]["ranks"]
    assert [
        (
            row["rank"],
            row["all_support_graphs"],
            row["prefilter_support_graphs"],
            row["decorated_labelled_grams_tested"],
            row["labelled_survivors"],
            row["unlabelled_orbits"],
        )
        for row in ranks
    ] == [
        (3, 8, 4, 20, 10, 4),
        (4, 64, 19, 176, 54, 5),
        (5, 1024, 137, 2384, 240, 3),
    ]
    assert _archive()["finite_enumeration"]["total_unlabelled_orbits"] == 12
    assert len({tuple(key) for row in ranks for key in row["canonical_keys"]}) == 12


def test_four_previously_missing_orbits_are_complete_compact_chambers() -> None:
    assert _archive()["newly_added_examples"] == [
        "star_5_3_3",
        "cycle4_one5",
        "cycle4_opposite5",
        "fork5_terminal5",
    ]
    assert _archive()["new_rows_exact_summary"] == {
        "star_5_3_3": {
            "dimension": 3,
            "facets": 4,
            "boundary_contacts": 26,
            "interior_cuts": 0,
            "trace_determinant": -16,
            "trace_discriminant_invariants": [2, 2, 2, 2],
        },
        "cycle4_one5": {
            "dimension": 3,
            "facets": 4,
            "boundary_contacts": 28,
            "interior_cuts": 0,
            "trace_determinant": -31,
            "trace_discriminant_invariants": [31],
        },
        "cycle4_opposite5": {
            "dimension": 3,
            "facets": 4,
            "boundary_contacts": 44,
            "interior_cuts": 0,
            "trace_determinant": -55,
            "trace_discriminant_invariants": [55],
        },
        "fork5_terminal5": {
            "dimension": 4,
            "facets": 5,
            "boundary_contacts": 108,
            "interior_cuts": 0,
            "trace_determinant": -16,
            "trace_discriminant_invariants": [4, 4],
        },
    }


def test_every_orbit_has_literal_vinberg_and_k3_return_certificates() -> None:
    examples = _archive()["examples"]
    assert len(examples) == 12
    for example in examples:
        assert all(example["checks"].values())
        assert example["vinberg_compactness"]["finite_vertices"] == (
            example["dimension"] + 1
        )
        assert example["vinberg_compactness"]["ideal_vertices"] == 0
        assert example["full_norm_two_arrangement"]["mixed_sign_roots"] == []
        assert example["full_norm_two_arrangement"]["interior_cut_count"] == 0
        assert example["full_norm_two_arrangement"]["seed_is_complete_chamber"]
        assert example["trace_k3_package"]["primitive_k3_embedding_automatic"]
        assert example["k3_eigen_nef_return"][
            "projectivized_eigen_nef_section_equals_closed_seed"
        ]
        assert example["lineage"]["subdivided"] is False
        assert example["lineage"]["worldwide_novelty_claim"] is False


def test_lanner_theorem_closes_all_higher_ranks() -> None:
    closure = _archive()["rank_ge_6_nonexistence"]
    assert closure["maximum_rank"] == 5
    assert closure["maximum_dimension"] == 4
    assert "no rank k>=6" in closure["application"]
    assert "Lanner" in closure["theorem"]

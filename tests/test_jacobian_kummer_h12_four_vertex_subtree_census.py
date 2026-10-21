from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from scripts.certify_jacobian_kummer_h12_four_vertex_subtree_census import (
    OUTPUT,
    _canonical_color_refinement,
    _canonical_refinement_colors,
    build,
)


@pytest.fixture(scope="module")
def payload() -> dict:
    rebuilt = build()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_complete_tree_coloring_denominator(payload: dict) -> None:
    denominator = payload["coloring_denominator"]
    assert denominator["P4"]["admissible_ordered_words"] == 1_379_952
    assert denominator["P4"]["unrooted_colored_paths_before_seed_action"] == 696_192
    assert denominator["P4"]["seed_action_orbits_after_reversal"] == 1_964
    assert denominator["P4"]["palindromic_outer_color_orbits"] == 47
    assert denominator["K1,3"]["unrooted_colored_stars_before_seed_action"] == 227_920
    assert denominator["K1,3"]["seed_action_orbits"] == 680
    assert denominator["total_unrooted_colorings_before_seed_action"] == 924_112
    assert denominator["total_seed_and_tree_automorphism_orbits"] == 2_644


def test_exact_facet_formula_and_census(payload: dict) -> None:
    theorem = payload["theorem"]
    assert theorem["total_pairwise_nonisomorphic_ordinary_Coxeter_polytopes"] == 2_644
    assert theorem["P4_ordinary_Coxeter_type_count"] == 1_964
    assert theorem["K13_ordinary_Coxeter_type_count"] == 680
    assert theorem["P4_and_K13_types_are_mutually_disjoint"]
    assert theorem["facet_formula"] == (
        "f(T)=607+9*B(T), where B counts edge-color occurrences in W77--W172"
    )
    assert theorem["type_facet_histogram"] == {
        "607": 23,
        "616": 170,
        "625": 844,
        "634": 1_607,
    }
    assert theorem["raw_unrooted_coloring_facet_histogram"] == {
        "607": 2_480,
        "616": 47_616,
        "625": 293_376,
        "634": 580_640,
    }
    assert theorem["reflection_index_and_volume_multiple"] == 4
    assert theorem["finite_volume"] and not theorem["compact"]


def test_ordinary_coxeter_nonisomorphism_and_carrier(payload: dict) -> None:
    certificate = payload["ordinary_Coxeter_nonisomorphism_certificate"]
    assert certificate["parallel_and_ultraparallel_both_labelled_infinity"]
    assert certificate["initial_exact_degree_profile_class_count"] == 2_535
    assert certificate["degree_profile_collision_group_count"] == 102
    assert certificate["degree_profile_collision_member_count"] == 211
    assert certificate["color_refined_member_count"] == 211
    assert certificate["color_refined_discrete_canonical_form_count"] == 195
    assert certificate["all_2644_combined_invariants_distinct"]
    audits = certificate["permutation_invariance_audit"]
    assert set(audits) == {"discrete", "nondiscrete"}
    for expected_discrete, name in ((True, "discrete"), (False, "nondiscrete")):
        audit = audits[name]
        assert audit["discrete"] is expected_discrete
        assert audit["tested_vertex_permutations"] == [
            "reverse",
            "even_vertices_then_odd_vertices",
        ]
        assert audit["canonical_numeric_color_IDs_equivariant"]
        assert audit["exact_final_invariant_word_unchanged"]

    representatives = payload["orbit_representatives"]
    assert len(representatives) == 2_644
    type_hashes = {
        row["ordinary_Coxeter_type_invariant_sha256"] for row in representatives
    }
    assert len(type_hashes) == 2_644
    for row in representatives:
        assert row["facet_count"] == 607 + 9 * row["B_edge_occurrence_count"]
        assert row["primitive_integral_reflection_root_count"] == row["facet_count"]
        assert row["fixed_boundary_core_Smith_diagonal"] == [1] * 13
        assert row["full_linear_rank"] == 13


@pytest.mark.parametrize(
    ("labels", "expected_discrete"),
    [
        (
            [
                [-1, 0, 1, 0, 1, 1],
                [0, -1, 0, 1, 0, 0],
                [1, 0, -1, 1, 0, 1],
                [0, 1, 1, -1, 1, 0],
                [1, 0, 0, 1, -1, 0],
                [1, 0, 1, 0, 0, -1],
            ],
            True,
        ),
        (
            [
                [-1, 1, 0, 0, 0],
                [1, -1, 1, 0, 0],
                [0, 1, -1, 1, 0],
                [0, 0, 1, -1, 1],
                [0, 0, 0, 1, -1],
            ],
            False,
        ),
    ],
)
def test_refinement_color_IDs_and_words_are_permutation_invariant(
    labels: list[list[int]], expected_discrete: bool
) -> None:
    matrix = np.array(labels, dtype=np.int8)
    colors, rounds = _canonical_refinement_colors(matrix)
    invariant, data = _canonical_color_refinement(matrix)
    assert data["color_refinement_discrete_canonical_form"] is expected_discrete
    permutations = (
        np.arange(len(matrix) - 1, -1, -1),
        np.concatenate((np.arange(0, len(matrix), 2), np.arange(1, len(matrix), 2))),
    )
    for permutation in permutations:
        permuted = matrix[np.ix_(permutation, permutation)]
        permuted_colors, permuted_rounds = _canonical_refinement_colors(permuted)
        permuted_invariant, _ = _canonical_color_refinement(permuted)
        assert np.array_equal(permuted_colors, colors[permutation])
        assert permuted_rounds == rounds
        assert permuted_invariant == invariant


def test_scope_and_lineage_are_bounded(payload: dict) -> None:
    assert payload["all_checks_pass"]
    assert all(payload["checks"].values())
    assert payload["allcock_and_carrier"]["same_full_Picard_K3_family"]
    assert payload["allcock_and_carrier"]["weak_divisor_reflection_carrier_only"]
    assert not payload["allcock_and_carrier"][
        "full_nef_effectivity_or_surface_action_claim"
    ]
    assert payload["bounded_lineage"]["headers_607_616_625_634_absent"]
    assert payload["bounded_lineage"]["evidence_grade"] == (
        "B1: bounded same-header absence only"
    )
    assert payload["bounded_lineage"]["worldwide_priority"] == "open"
    bounds = payload["arithmetic_runtime_bounds"]
    preflight = bounds["preflight_absolute_entry_bounds"]
    assert preflight["four_reflection_product_entry_absolute"] < bounds[
        "int64_maximum"
    ]
    assert preflight["transformed_seed_interior_entry_absolute"] < bounds[
        "int64_maximum"
    ]
    assert preflight["sum_of_four_transformed_interiors_entry_absolute"] < bounds[
        "int64_maximum"
    ]

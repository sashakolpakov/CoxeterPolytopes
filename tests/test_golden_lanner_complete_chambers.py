from __future__ import annotations

import json

from scripts.certify_golden_lanner_complete_chambers import OUTPUT, analyze


def _archive() -> dict[str, object]:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_exact_replay_matches_frozen_archive() -> None:
    payload = analyze()
    assert payload == _archive()
    assert payload["all_checks_pass"] is True
    assert payload["global_conclusion"] == {
        "all_eight_original_seeds_are_full_norm_two_chambers": True,
        "all_eight_are_compact": True,
        "all_eight_have_automatic_K3_eigen_nef_lifts": True,
        "period_dimensions": [14, 12, 10],
        "subdivision_count": 0,
        "worldwide_novelty_claim": False,
    }


def test_all_eight_full_root_chambers_are_literal_original_simplices() -> None:
    examples = {example["seed_id"]: example for example in _archive()["examples"]}
    expected = {
        "triangle_255": (2, 3, 9, 6, 12, 12),
        "triangle_335": (2, 3, 8, 6, 9, 9),
        "triangle_355": (2, 3, 10, 6, 12, 12),
        "triangle_555": (2, 3, 12, 6, 15, 15),
        "chain_353": (3, 4, 25, 25, 80, 63),
        "chain_535": (3, 4, 27, 25, 104, 72),
        "chain_5333": (4, 5, 64, 90, 2463, 2327),
        "chain_5335": (4, 5, 114, 90, 4992, 4633),
    }
    assert set(examples) == set(expected)
    for seed_id, values in expected.items():
        example = examples[seed_id]
        arrangement = example["full_norm_two_arrangement"]
        actual = (
            example["dimension"],
            example["simple_walls"]["count"],
            arrangement["boundary_contact_mirror_count"],
            arrangement["mixed_sign_patterns_up_to_global_sign"],
            arrangement["mixed_sign_assignments_tested"],
            arrangement["boundary_assignments_tested"],
        )
        assert actual == values
        assert arrangement["mixed_sign_norm_two_roots"] == []
        assert arrangement["interior_cut_count"] == 0
        assert arrangement["seed_is_complete_chamber"] is True
        assert example["lineage"]["same_as_input_seed"] is True
        assert example["lineage"]["subdivided_by_full_root_arrangement"] is False


def test_literal_vertex_and_vinberg_compactness_certificates() -> None:
    for example in _archive()["examples"]:
        dimension = example["dimension"]
        vertices = example["literal_vertices"]
        assert len(vertices) == dimension + 1
        assert example["field_identity_inertia"] == [dimension, 1, 0]
        assert example["field_conjugate_inertia"] == [dimension + 1, 0, 0]
        assert example["specified_interior_point"]["facet_pairings"] == [
            {"a": "-1", "b": "0"}
        ] * (dimension + 1)
        assert all(vertex["finite"] for vertex in vertices)
        assert all(
            vertex["vertex_gram_identity_inertia"] == [dimension, 0, 0]
            and vertex["vertex_gram_conjugate_inertia"] == [dimension, 0, 0]
            for vertex in vertices
        )
        assert example["vinberg_compactness"] == {
            "criterion": (
                "the full Gram has signature (n,1), and every deletion of "
                "one of the n+1 simplex walls is elliptic; hence all n+1 "
                "vertices are finite and the Coxeter simplex is compact"
            ),
            "finite_vertex_count": dimension + 1,
            "ideal_vertex_count": 0,
            "hyperideal_vertex_count": 0,
            "compact": True,
            "finite_volume": True,
        }


def test_trace_lattices_and_automatic_k3_promotions() -> None:
    expected = {
        "triangle_255": (-16, [4, 4], 14),
        "triangle_335": (-20, [2, 10], 14),
        "triangle_355": (-36, [6, 6], 14),
        "triangle_555": (-100, [10, 10], 14),
        "chain_353": (-11, [11], 12),
        "chain_535": (-19, [19], 12),
        "chain_5333": (-4, [2, 2], 10),
        "chain_5335": (-20, [2, 10], 10),
    }
    for example in _archive()["examples"]:
        trace = example["trace_k3_package"]
        assert (
            trace["determinant"],
            trace["discriminant_group_invariants"],
            trace["generic_k3_period_dimension"],
        ) == expected[example["seed_id"]]
        assert trace["even"] is True
        assert trace["primitive_k3_embedding_automatic"] is True
        assert trace["golden_action_self_adjoint"] is True
        assert trace["golden_action_minimal_polynomial"] is True
        assert trace["paired_roots_square_minus_two"] is True
        assert trace["same_wall_root_pairs_orthogonal"] is True
        assert trace["paired_off_diagonal_intersections_nonnegative"] is True
        assert trace["visible_k3_wall_shell_collapse"][
            "second_shell_maps_to_norm_two"
        ] is True
        promotion = example["k3_eigen_nef_promotion"]
        assert promotion["full_picard_projective_k3_exists"] is True
        assert promotion["projectivized_eigen_nef_section_equals_closed_seed"] is True
        assert promotion["facet_curve_classes"] == 2 * (example["dimension"] + 1)
        assert promotion["surface_automorphism_claim"] is False

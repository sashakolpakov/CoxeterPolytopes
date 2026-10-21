from __future__ import annotations

from functools import lru_cache
import json

import sympy as sp

from scripts.certify_coble_one_boundary_raw_k3_promotion import (
    OUTPUT,
    _all_discriminant_anti_isometries,
    build_certificate,
)


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


@lru_cache(maxsize=1)
def _live() -> dict:
    return build_certificate()


def test_live_exhaustive_replay_matches_archive() -> None:
    assert _live() == _archive()


def test_surface_free_denominator_and_non_tautological_iff() -> None:
    payload = _archive()
    raw = payload["surface_free_denominator"]
    assert raw["name"] == "D_1_split"
    assert payload["accepted_locus"]["name"] == "R_1_split"
    assert payload["finite_compiler"]["name"] == "C_1_split"
    assert "a pre-existing algebraic surface" in raw["not_a_hypothesis"]
    assert "equality between the input chamber and a nef cone" in raw["not_a_hypothesis"]
    assert "all marked wall systems" in raw["a_priori_scope"]
    assert "two-input exact calibration" in raw["a_priori_scope"]
    rendered_raw_definition = json.dumps(raw, sort_keys=True)
    assert "Kondo" not in rendered_raw_definition
    assert "Type I" not in rendered_raw_definition
    assert "passing branch" not in raw["definition"]
    assert "specialized arithmetic split censuses" in payload["finite_compiler"][
        "implementation_scope"
    ]
    assert "if and only if" in payload["promotion_theorem"]["intrinsic_iff"]
    assert "no compact P" in payload["scope"]["compactness"]
    for case in payload["cases"]:
        strong_coble = case["geometric_promotion"]["strong_coble_linear_systems"]
        assert "|-K_S| is empty" in strong_coble
        assert "|-2K_S|={B}" in strong_coble


def test_all_split_survivors_and_rank20_k3_embeddings_are_certified() -> None:
    cases = {case["raw_id"]: case for case in _archive()["cases"]}
    expected = {
        "raw_22": (27_648, 18, 1_024, -8, [1] * 18 + [2, 4], [[2, 0], [0, 4]], 8),
        "raw_30": (32_768, 16, 2_048, -16, [1] * 18 + [2, 8], [[2, 0], [0, 8]], 16),
    }
    for raw_id, (
        assignments,
        orbits,
        survivor_assignments,
        determinant,
        smith,
        complement,
        glue_index,
    ) in expected.items():
        case = cases[raw_id]
        split = case["exhaustive_split_compiler"]
        assert split["raw_assignment_count"] == assignments
        assert split["sheet_switch_orbit_count"] == orbits
        assert split["admissible_orbit_count"] == 1
        assert split["admissible_raw_assignment_count"] == survivor_assignments
        assert split["all_admissible_orbits_retained"] is True
        lattice = case["rank20_equivariant_lattice"]
        assert lattice["determinant"] == determinant
        assert lattice["smith_diagonal"] == smith
        assert all(lattice["checks"].values())
        embedding = case["primitive_equivariant_K3_embedding"]
        assert embedding["unique_complement_gram"] == complement
        assert embedding["graph_glue_index"] == glue_index
        assert embedding["fixed_lattice_nikulin_invariants"] == [11, 11, 1]
        assert all(embedding["checks"].values())

    assert [
        row["index_over_generated_L"]
        for row in cases["raw_22"]["primitive_equivariant_K3_embedding"][
            "even_overclosure_census"
        ]
    ] == [1]
    raw_30_closures = cases["raw_30"]["primitive_equivariant_K3_embedding"][
        "even_overclosure_census"
    ]
    assert [row["index_over_generated_L"] for row in raw_30_closures] == [1, 2]
    assert raw_30_closures[1]["determinant"] == -4
    rejected = raw_30_closures[1]
    complement = rejected["compatible_primitive_K3_complements"]
    assert len(complement) == 1
    assert complement[0]["gram"] == [[2, 0], [0, 2]]
    assert complement[0]["discriminant_anti_isometry_count"] == 2
    assert complement[0]["primitive_equivariant_K3_glue_count"] == 2
    assert complement[0]["marked_fixed_lattice_passing_glue_count"] == 0
    for glue in complement[0]["equivariant_glue_audits"]:
        assert glue["primitive_equivariant_K3_embedding_pass"] is True
        assert glue["marked_fixed_lattice_gate_pass"] is False
        assert glue["marked_M_primitive_smith"] == [1] * 10 + [2]
        assert glue["actual_fixed_lattice_nikulin_invariants"] == [11, 9, 1]
    assert rejected["primitive_equivariant_K3_embedding_gate_pass"] is True
    assert rejected["fixed_lattice_gate_pass"] is False
    assert rejected["degree_zero_root_count"] == 18
    assert rejected["anti_invariant_minus_two_root_count"] == 18
    root_audit = rejected["full_bounded_root_audit"]
    assert root_audit["enumeration"]["bound"] == 118
    assert root_audit["enumeration"]["signed_vector_count"] == 96
    assert root_audit["positive_smaller_root_degree_histogram"] == {"1": 39}
    assert root_audit["irreducibility_obstruction_count"] == 6
    assert sum(
        len(record["witnesses"])
        for record in root_audit["irreducibility_obstructions"]
    ) == 12
    assert [
        record["proposed_root_zero_based"]
        for record in root_audit["irreducibility_obstructions"]
    ] == [6, 7, 14, 15, 22, 23]
    assert rejected["period_ample_gate_pass"] is False


def test_order_four_discriminant_anti_isometry_regression() -> None:
    raw_30 = next(
        case for case in _archive()["cases"] if case["raw_id"] == "raw_30"
    )
    closure = raw_30["primitive_equivariant_K3_embedding"][
        "even_overclosure_census"
    ][1]
    anti_isometries = _all_discriminant_anti_isometries(
        sp.Matrix(closure["gram"]), sp.diag(2, 2)
    )
    assert len(anti_isometries) == 2


def test_complete_bounded_root_gate_and_non_singleton_exact_return() -> None:
    cases = {case["raw_id"]: case for case in _archive()["cases"]}
    expected = {
        "raw_22": (138, 374, 474, {"1": 20, "8": 4}),
        "raw_30": (58, 118, 36, {"1": 18, "2": 6}),
    }
    for raw_id, (square, bound, count, degrees) in expected.items():
        gate = cases[raw_id]["period_and_simultaneous_irreducibility"][
            "bounded_root_gate"
        ]
        assert gate["ample_square"] == square
        assert gate["enumeration"]["bound"] == bound
        assert gate["enumeration"]["signed_vector_count"] == count
        assert gate["proposed_root_degree_histogram"] == degrees
        assert gate["degree_zero_root_count"] == 0
        assert gate["irreducibility_obstruction_count"] == 0
        assert all(gate["checks"].values())

    payload = _archive()
    assert payload["classification_corollary"][
        "accepted_underlying_polytope_isometry_classes"
    ] == 2
    assert payload["classification_corollary"]["facet_counts"] == [22, 30]
    assert payload["scope"]["rating"] == (
        "uniform non-singleton L3 exact-return criterion and classification "
        "of the compiler-positive locus; not strict L4 under the 2026-08-09 "
        "strengthened five-part contract"
    )
    contract = payload["scope"]["contract_audit"]
    assert contract["compiler_exhaustivity"].startswith("pass:")
    assert contract["branchwise_exact_return"].startswith("pass:")
    assert contract["intrinsic_iff"].startswith("pass:")
    assert contract["class_wide_realized_branch"].startswith("not proved:")
    assert "R_1_split, not" in contract["AG_pushback_scope"]
    assert payload["promotion_theorem"]["essential_surjectivity"].startswith(
        "not proved on D_1_split:"
    )
    assert payload["all_checks_pass"] is True

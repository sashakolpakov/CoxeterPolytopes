"""Second independent hostile tests for compact narrow-facet K3 field finiteness."""

from __future__ import annotations

from scripts.audit_k3_compact_narrow_facet_field_finiteness_second_independent_hostile import (
    PINS,
    ROOT,
    build_audit,
    cartan_coboundary_certificate,
    height_and_discriminant_certificate,
    sha256,
    star_spanning_certificate,
)


def test_owner_bytes_are_terminally_frozen() -> None:
    for relative, expected in PINS.values():
        assert sha256(ROOT / relative) == expected


def test_compact_positive_dimensional_facet_star_spans() -> None:
    record = star_spanning_certificate()
    assert record["compact_Coxeter_chamber_is_simple"] is True
    assert record["positive_dimensional_facet_has_nonincident_side"] is True
    assert record["projected_star_spans_e_perp"] is True
    assert record["adding_minimal_facet_normal_spans_V"] is True
    assert len(record["exact_coordinate_rows"]) == 8
    assert "m=2" in record["dimension_one_countergate"]


def test_cycle_equalities_force_cartan_diagonal_similarity() -> None:
    record = cartan_coboundary_certificate()
    assert record["cycle_holonomy_reconstructs_vertex_coboundary"] is True
    assert record["diagonal_similarity_checked_entrywise"] is True
    assert record["poisoned_cycle_holonomy_rejected"] is True
    assert record["proper_cycle_subfield_contradiction"] is True


def test_k3_rank_and_primitive_grid_give_claimed_constant() -> None:
    record = height_and_discriminant_certificate()
    assert record["minimum_module_rank"] == 3
    assert record["nonrational_module_rank_ceiling"] == 10
    assert record["nonrational_field_degree_ceiling"] == 6
    assert record["strict_tower_generator_count_ceiling"] == 2
    assert record["embedding_pair_count_ceiling"] == 15
    assert record["conjugate_difference_height"] == 32 * 14**10
    assert record["uniform_Vandermonde_exponent"] == 30
    assert record["field_discriminant_bound"] == "(32*14^10)^30"


def test_exact_scope_and_nonclaims_are_pinned() -> None:
    audit = build_audit()
    assert audit["verdict"] == "PASS_EXACT_STATED_SCOPE"
    assert all(audit["scope"].values())

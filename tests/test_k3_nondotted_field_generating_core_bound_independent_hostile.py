"""Independent hostile tests for the nondotted field-generating-core theorem."""

from __future__ import annotations

from scripts.audit_k3_nondotted_field_generating_core_bound_independent_hostile import (
    PINS,
    ROOT,
    build_audit,
    chordless_rank_certificate,
    conjugate_height_certificate,
    primitive_grid_and_discriminant_certificate,
    sha256,
    trace_integrality_certificate,
)


def test_owner_bytes_are_frozen() -> None:
    for relative, expected in PINS.values():
        assert sha256(ROOT / relative) == expected


def test_projective_module_trace_not_directed_entry_integrality() -> None:
    record = trace_integrality_certificate()
    assert record["finite_projective_module_has_canonical_O_valued_trace"] is True
    assert record["primitive_normals_required"] is False
    assert record["individual_directed_Cartan_entries_assumed_integral"] is False
    assert record["normal_rescaling_invariance_checked"] is True


def test_chordless_minor_gives_ell_at_most_m_plus_two() -> None:
    record = chordless_rank_certificate()
    assert record["rank_inequality"] == "ell-2<=rank(G_C)<=m"
    assert record["cycle_length_bound"] == "ell<=m+2"
    assert all(row["minor_nonzero"] for row in record["exact_minor_rows"])


def test_every_embedding_gets_the_same_safe_cycle_height() -> None:
    record = conjugate_height_certificate()
    assert record["cycle_square_identity"] == "c(C)^2=product(perimeter p_ij)"
    by_rank = {row["module_rank"]: row for row in record["records"]}
    assert by_rank[3]["all_embedding_cycle_height"] == 2**5
    assert by_rank[10]["all_embedding_cycle_height"] == 2**12
    assert record["pair_generators_obey_same_final_height"] is True


def test_two_generators_grid_and_exponent_300() -> None:
    record = primitive_grid_and_discriminant_certificate()
    assert record["minimal_generating_subpacket_size_upper_bound"] == 2
    assert record["embedding_pair_count_upper_bound"] == 15
    assert record["primitive_grid"] == "x+q*y for q=0,...,15"
    assert record["Vandermonde_exponent_formula"] == "(m+7)d(d-1)"
    assert record["uniform_binary_exponent"] == 300
    assert record["field_discriminant_bound"] == "2^300"
    assert record["Hermite_Minkowski_gives_finite_field_set"] is True


def test_scope_allows_only_dots_outside_the_generating_core() -> None:
    audit = build_audit()
    assert audit["verdict"] == "PASS_CONDITIONAL_EXACT_SCOPE"
    assert all(audit["scope"].values())

"""Independent hostile tests for the chordless-cycle K3 field theorem."""

from __future__ import annotations

from pathlib import Path

from scripts.audit_k3_chordless_cycle_height_field_bound_independent_hostile import (
    OWNER_PINS,
    ROOT,
    affine_shift_audit,
    audit_payload,
    directed_chord_audit,
    generator_constant_audit,
    rank_minor_audit,
    sha256,
    trace_induction_audit,
)


REPORT = ROOT / (
    "results/k3_chordless_cycle_height_field_bound_independent_hostile_audit.md"
)
OWNER = ROOT / "results/k3_chordless_cycle_height_field_bound.md"


def test_owner_bytes_are_exactly_the_repaired_candidate() -> None:
    for item in OWNER_PINS.values():
        assert sha256(ROOT / item["path"]) == item["sha256"]

    owner = OWNER.read_text(encoding="utf-8")
    assert "Choose a totally negative codifferent scalar \\(\\alpha\\)" in owner
    assert "exact full-Picard **complex projective** K3" in owner
    assert "0<\\sigma(\\theta_C)\\le2U_C" in owner


def test_trace_induction_handles_fractional_entries_and_repeated_words() -> None:
    record = trace_induction_audit()
    assert record["integral_reflection_word_count"] == 1367
    assert record["fractional_individual_Cartan_entries_present"] > 0
    assert record["all_trace_subset_expansions_exact"] is True
    assert record["all_full_ordered_cycle_products_integral"] is True
    assert record["repeated_words_tested"] == 3


def test_cycle_reduction_keeps_directed_chord_factors() -> None:
    record = directed_chord_audit()
    assert record["directed_chord_splits_tested"] == 1092
    assert record["perimeter_directions_each_occur_once"] is True
    assert record["both_chord_directions_occur_once"] is True
    assert record["field_identity"] == "c(C1)c(C2)=c(C)c(i,j)c(j,i)"


def test_m_plus_two_bound_has_a_nonzero_universal_minor() -> None:
    record = rank_minor_audit()
    assert record["kernel_dimension_at_most"] == 2
    assert record["cycle_length_at_most_m_plus_2"] is True
    assert len(record["exact_minor_records"]) == 28
    for row in record["exact_minor_records"]:
        assert row["rank_lower_bound"] == row["length"] - 2
        assert row["witness_minor_determinant"] != 0


def test_affine_shift_and_generator_constants_are_exact() -> None:
    affine = affine_shift_audit()
    assert affine["exact_affine_place_samples"] == 240
    assert affine["large_non_distinguished_bound"] == "<=2^(ell+1)"
    assert affine["large_coordinate_dominant"] is True
    assert affine["zero_boundary_excluded_by_embedding_injectivity"] is True

    constants = generator_constant_audit()
    n2 = constants["dimension_records"][0]
    assert n2 == {"n": 2, "D": 6, "L": 5, "U": 32, "r": 2, "M_D": 16}
    assert constants["n2_primitive_bound"] == "32*max(64,C/s0)"
    assert constants["n2_discriminant_exponent"] == 30
    assert all(
        row["max_strict_tower_length"]
        == (row["degree_bound"].bit_length() - 1)
        for row in constants["tower_records"]
    )


def test_full_independent_payload_keeps_the_conditional_scope() -> None:
    payload = audit_payload()
    assert payload["scope"] == {
        "complex_projective_K3_required_for_Picard_20": True,
        "field_degree_at_least_two_required_for_cycle_length": True,
        "finite_projective_order_module_or_explicit_trace_certificate": True,
        "uniform_large_cycle_mixed_height_is_assumed": True,
        "height_bound_is_not_derived_from_reflection_maximality": True,
        "verdict": "PASS_CONDITIONAL",
    }


def test_independent_report_records_every_decisive_gate() -> None:
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    required = [
        "PASS, conditional",
        "complex projective K3 surfaces",
        r"q=\operatorname {Tr}_{k/\mathbf Q}(\alpha h)",
        "unique full subset has coefficient",
        "six genuinely fractional individual Cartan entries",
        r"c(C_1)c(C_2)=c(C)\,a_{ij}a_{ji}",
        "retains the ordered factors",
        "explicit nonzero minor",
        r"\ell\le m+2=n+3",
        r"0<\sigma(\theta_C)\le2U_C<\sigma_0(\theta_C)",
        r"H\cdot\tau_{\theta_C}(H)\le C",
        r"r=\lfloor\log_2D\rfloor",
        r"M_D={D(D-1)\over2}+1",
        r"\left(32\max\{64,C/s_0\}\right)^{30}",
        "does **not** prove",
        "height bound is not derived from reflection maximality",
    ]
    for phrase in required:
        assert phrase in text

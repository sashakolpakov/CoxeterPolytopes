from __future__ import annotations

import json

from scripts.certify_automorphic_pushback_l4_hostile_audit import (
    OUTPUT,
    analyze,
    trace_shadow_row,
)


def _payload() -> dict[str, object]:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_archive_is_exact_replay() -> None:
    assert _payload() == analyze()


def test_trace_shadow_formulas_and_reflection_eigenspace() -> None:
    for n in (4, 8, 9, 10, 12, 21, 29):
        for degree in (2, 3, 4):
            row = trace_shadow_row(n, degree)
            assert row["trace_Q_rank"] == degree * (n + 2)
            assert row["trace_signature"] == [2, degree * (n + 2) - 2]
            assert row["Hilbert_domain_codimension_in_trace_domain"] == (
                (degree - 1) * (n + 2)
            )
            assert row["facet_image_codimension_lower_bound"] == (
                (degree - 1) * (n + 2) + 1
            )
            assert row["facet_reflection_minus_eigenspace_Q_dimension"] == degree
            assert row["trace_fixed_locus_codimension"] == degree
            assert row["facet_reflection_is_Q_pseudoreflection"] is False
            assert row["trace_fixed_locus_is_divisor"] is False
            assert row["restriction_group_is_proper"] is True
            assert row["group_dimension_difference"] == (
                degree * (degree - 1) * (n + 2) ** 2 // 2
            )


def test_degree_one_is_exactly_the_classical_branch() -> None:
    row = trace_shadow_row(13, 1)
    assert row["trace_signature"] == [2, 13]
    assert row["Hilbert_domain_codimension_in_trace_domain"] == 0
    assert row["facet_reflection_minus_eigenspace_Q_dimension"] == 1
    assert row["facet_reflection_is_Q_pseudoreflection"] is True
    assert row["trace_fixed_locus_is_divisor"] is True
    assert row["restriction_group_is_proper"] is False


def test_no_candidate_class_passes_the_strict_L4_gate() -> None:
    payload = _payload()
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    candidates = payload["candidate_class_hostile_audit"]
    assert len(candidates) >= 5
    assert all(candidate["passes_L4"] is False for candidate in candidates)
    verdict = payload["verdict"]
    assert verdict["strict_noncircular_compact_L4_subclass_found"] is False
    assert verdict["new_polytope_bound_claimed"] is False
    assert verdict["new_nonexistence_theorem_for_polytopes_claimed"] is False
    assert verdict["rating"] == "L3_CARRIER_PLUS_L4_SCOPE_OBSTRUCTION_NOT_L4"


def test_two_independent_missing_lemmas_and_same_P_gate_are_explicit() -> None:
    payload = _payload()
    missing = payload["minimum_completion_package"]
    assert [item["name"] for item in missing] == [
        "automatic_bounded_purity",
        "Hilbert_reflective_finiteness",
    ]
    assert all(item["known"] is False for item in missing)

    support = payload["native_Hilbert_exact_support_gates"]
    assert support["forced_by_compact_root_reflectivity"] is False
    assert support["arbitrary_automorphic_exact_support_decidable"] is False
    assert support["compensated_divisor_counts_as_pure_support"] is False

    same_p = payload["same_polytope_return_gate"]
    assert same_p["finite_carrier_lattices_alone_imply_finite_polytopes"] is False
    assert any("complete chamber" in item for item in same_p["required"])


def test_classical_results_are_not_misapplied_to_the_Hilbert_group() -> None:
    payload = _payload()
    scope = payload["classical_theorem_scope"]
    assert scope["Ma_2018"]["usable_on_native_compact_Hilbert_carrier"] is False
    assert (
        scope["Wang_2025_2026"]["usable_on_native_compact_Hilbert_carrier"]
        is False
    )
    assert scope["Wang_2025_2026"]["trace_numerology_is_not_a_transfer"] is True
    assert (
        scope["Gritsenko_Hulek_Kodaira"][
            "strong_reflective_support_forced_by_P"
        ]
        is False
    )

from __future__ import annotations

from fractions import Fraction

from scripts.certify_compact_shimura_ci_bmy_no_go import (
    analyze,
    coefficients,
    very_ample_multiple,
)


def test_symbolic_coefficient_identity_on_independent_grid() -> None:
    for n in range(2, 25):
        for k in (1, 2, 7, 31):
            data = coefficients(n, k)
            assert data["bmy_defect"] == data["bmy_defect_expected"]
            assert data["bmy_defect"] > 0
            assert data["reduced_mirror_degree_coefficient"] >= 0


def test_demailly_multiple_and_bugaenko_specialization() -> None:
    assert very_ample_multiple(8) == 1_081_593
    row = coefficients(8, very_ample_multiple(8))
    assert row["surface_canonical"] == Fraction(6_489_566)
    assert row["surface_c2"] == Fraction(24_566_763_687_122)
    assert row["bmy_defect"] == Fraction(31_585_824_193_010)


def test_aggregate_scope_flags() -> None:
    result = analyze()
    assert result["coefficient_identity_verified"] is True
    assert result["strictly_positive_for_all_n_ge_2_k_ge_1"] is True
    assert (
        result[
            "reduced_snc_mirror_log_bmy_strictly_positive_conditional_on_snc"
        ]
        is True
    )
    assert result["unresolved_multiple_point_log_resolution_may_add_feedback"] is True
    assert result["ordinary_bmy_dimension_or_facet_bound"] is False
    assert (
        result["surface_degree_inequality_stronger_than_ambient_degree_injection"]
        is False
    )

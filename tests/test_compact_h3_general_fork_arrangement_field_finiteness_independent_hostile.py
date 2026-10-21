from scripts.audit_compact_h3_general_fork_arrangement_field_finiteness_independent_hostile import (
    ENDPOINT_EXPONENT,
    LABEL_CUTOFF,
    build_certificate,
    corrected_norm_identity,
    determinant_and_cartan_sign,
    endpoint_cutoff,
    gamma,
    h_power_bound_regression,
)


def test_cartan_cycle_sign_and_determinant_interval():
    row = determinant_and_cartan_sign()
    assert row["cycle_formula"] == "t12+t13+t23-z-4"
    assert row["cycle_formula_passed"]
    assert row["plus_z_formula_is_wrong"]
    assert row["discriminant_passed"]
    assert row["minimum_passed"]
    endpoint = row["lower_endpoint_negative_definite_example"]
    assert endpoint["passed"]
    assert endpoint["leading_principal_minors"] == ["-2", "15/4", "-9/2"]


def test_exceptional_selected_place_factor_is_present():
    row = corrected_norm_identity()
    assert row["passed"]
    assert row["selected_factor"] == "1024/(s_m0*s_n0)"
    assert row["polynomial_envelope"] == "<=4*M^4"


def test_cyclotomic_norm_factor_regression():
    assert gamma(3) == 3
    assert gamma(4) == 2
    assert gamma(5) == 5
    assert gamma(6) == 1
    assert gamma(8) == 2
    assert gamma(9) == 3
    assert gamma(12) == 1
    row = h_power_bound_regression(20_000)
    assert row["passed"]
    assert row["equality_cases"] == [3]


def test_exact_uniform_label_cutoff_endpoint():
    row = endpoint_cutoff()
    assert row["M"] == LABEL_CUTOFF == 1_000_000
    assert row["D_weakened_lower_bound"] == ENDPOINT_EXPONENT == 353
    assert row["endpoint_passed"]
    assert row["monotonicity_passed"]


def test_composite_independent_certificate():
    assert build_certificate()["all_checks_pass"]

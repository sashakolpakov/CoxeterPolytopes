from scripts.certify_h3_general_fork_field_finiteness import (
    LABEL_CUTOFF,
    build_certificate,
    capacity_ratio_bound_at_cutoff,
    euler_phi,
    monotonicity_certificate,
    three_node_determinant_identity,
    totient_lower_bound_regression,
)


def test_three_node_determinant_and_root_interval():
    row = three_node_determinant_identity()
    assert row["determinant_passed"]
    assert row["discriminant_passed"]
    assert row["beta"] == "a**2 + a*b*u + b**2 + u**2 - 4"
    assert row["cartan_cycle_formula"] == "beta=t12+t13+t23-c123-4"


def test_exact_label_cutoff_comparison():
    row = capacity_ratio_bound_at_cutoff()
    assert row["M"] == LABEL_CUTOFF == 1_000_000
    assert row["D_lower_bound"] == 353
    assert row["passed"]
    assert monotonicity_certificate()["passed"]


def test_totient_bound_regression():
    assert euler_phi(1_000_000) == 400_000
    assert totient_lower_bound_regression(25_000)


def test_composite_certificate():
    assert build_certificate()["all_checks_pass"]

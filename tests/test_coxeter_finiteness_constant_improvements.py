from scripts.certify_coxeter_finiteness_constant_improvements import (
    FORK_PAIR_BOUND,
    TAIL_BOUND,
    arithmetic_sieves,
    build_certificate,
    fork_individual_cutoff_certificate,
    fork_pair_sieve_certificate,
    gamma_from_least,
    real_cyclotomic_compositum_degree,
    tail_sieve_certificate,
)


def test_gamma_and_real_cyclotomic_degree_examples():
    phi, least = arithmetic_sieves(1000)
    assert gamma_from_least(3, least) == 3
    assert gamma_from_least(16, least) == 2
    assert gamma_from_least(12, least) == 1
    assert real_cyclotomic_compositum_degree(3, 4, phi) == 1
    assert real_cyclotomic_compositum_degree(8, 12, phi) == 4
    assert real_cyclotomic_compositum_degree(3, 420, phi) == 48


def test_individual_fork_cutoff():
    row = fork_individual_cutoff_certificate()
    assert row["passed"]
    assert row["failures"] == []
    assert row["contraction_failures"] == []
    assert row["conclusion_with_old_analytic_tail"] == "M<=840"


def test_exact_pair_sieve():
    row = fork_pair_sieve_certificate()
    assert row["passed"]
    assert row["survivor_count"] == 525
    assert row["maximum_surviving_label"] == FORK_PAIR_BOUND == 420
    assert row["endpoint_survivors"] == [{"m": 3, "n": 420, "degree": 48}]


def test_tail_sieve():
    row = tail_sieve_certificate()
    assert row["passed"]
    assert row["survivor_count"] == 66
    assert row["maximum_surviving_label"] == TAIL_BOUND == 126


def test_composite_certificate():
    assert build_certificate()["all_checks_pass"]

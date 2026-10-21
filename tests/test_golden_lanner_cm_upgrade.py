from scripts.certify_golden_lanner_cm_upgrade import analyze, golden_square_root
from coxeter_polytopes.exact.golden import Golden


def test_four_exact_quadratic_cm_upgrades() -> None:
    result = analyze()
    assert result["positive_count"] == 4
    assert result["all_positive_exact_matrix_checks_pass"] is True
    rows = {row["seed_id"]: row for row in result["rows"]}
    assert {
        seed_id
        for seed_id, row in rows.items()
        if row["status"].startswith("certified")
    } == {"triangle_255", "chain_5333", "chain_5335", "fork5_terminal5"}
    for seed_id in ("triangle_255", "chain_5333", "chain_5335", "fork5_terminal5"):
        row = rows[seed_id]
        assert row["J_square_minus_c"] is True
        assert row["J_K_skew_adjoint"] is True
        assert row["J_Q_skew_adjoint_for_trace_form"] is True
        assert row["J_commutes_with_golden_action"] is True
        assert row["determinant_CM_square_class_gate"] is True
        assert row["distinguished_positive_plane_is_J_stable"] is True
        assert row["c_totally_positive"] is True
        assert row["CM_period_domain_complex_dimension"] in {2, 3}


def test_negative_gates_are_partitioned_exactly() -> None:
    rows = {row["seed_id"]: row for row in analyze()["rows"]}
    determinant_failures = {"triangle_335", "triangle_355", "triangle_555"}
    odd_rank_failures = {
        "chain_353",
        "chain_535",
        "star_5_3_3",
        "cycle4_one5",
        "cycle4_opposite5",
    }
    assert all(
        rows[seed_id]["status"]
        == "quadratic_CM_extension_impossible_by_determinant_square_class"
        for seed_id in determinant_failures
    )
    assert all(
        rows[seed_id]["status"]
        == "quadratic_CM_extension_impossible_by_degree_divisibility"
        for seed_id in odd_rank_failures
    )


def test_exact_golden_square_gate() -> None:
    assert golden_square_root(Golden(4)) in {Golden(2), Golden(-2)}
    assert golden_square_root(Golden(4, 2)) is None
    assert golden_square_root(Golden(6) * Golden(0, 1) ** 2) is None
    assert golden_square_root(Golden(10) * Golden(0, 1) ** 2) is None

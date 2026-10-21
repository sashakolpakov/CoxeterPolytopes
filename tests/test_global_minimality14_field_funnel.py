from scripts.certify_global_minimality14_field_funnel import (
    build_certificate,
    cutoff_integer_comparison,
    euler_phi,
    gamma_value,
    h_vector_forcing_checks,
    lower_bound_totient_exact,
    small_face_forcing_checks,
    tail_capacity_strict,
    transition_routes,
)


def test_exact_cutoff_and_totient_inequality():
    assert cutoff_integer_comparison()
    assert all(lower_bound_totient_exact(n) for n in range(3, 25_000))
    assert euler_phi(20_000) == 8_000


def test_cyclotomic_norm_capacity_samples():
    assert gamma_value(7) == 7
    assert gamma_value(8) == 2
    assert gamma_value(12) == 1
    assert gamma_value(125) == 5
    assert all(tail_capacity_strict(d) for d in range(3, 20_000))


def test_every_transition_family_has_tail_or_bounded_fork():
    routes = transition_routes()
    assert set(routes) == {"Gamma1", "Gamma2", "Gamma3", "Gamma4", "Gamma5"}
    assert all(row["passed"] for row in routes.values())
    assert {name for name, row in routes.items() if row["kind"] == "fork"} == {
        "Gamma1",
        "Gamma2",
    }


def test_dimension_ten_forcing_is_uniform():
    rows = h_vector_forcing_checks(1000)
    assert rows[0]["n"] == 10
    assert rows[0]["lower"] == "9/2"
    assert rows[0]["strict_upper"] == "9/2"
    assert rows[1]["n"] == 11
    assert rows[1]["lower"] == "9/2"
    assert rows[1]["strict_upper"] == "9/2"
    assert all(row["contradiction"] for row in rows)


def test_dimension_six_forces_triangle_or_quadrilateral():
    rows = small_face_forcing_checks(1000)
    assert rows[0] == {
        "n": 6,
        "strict_upper": "5",
        "forces_triangle_or_quadrilateral": True,
    }
    assert rows[1]["n"] == 7
    assert rows[1]["strict_upper"] == "5"
    assert all(row["forces_triangle_or_quadrilateral"] for row in rows)


def test_composite_certificate():
    assert build_certificate()["all_checks_pass"]

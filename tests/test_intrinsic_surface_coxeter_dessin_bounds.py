from scripts.certify_intrinsic_surface_coxeter_dessin_bounds import (
    build_certificate,
    euler_phi,
    finite_angle_alphabet,
)


def test_totients_and_cutoff() -> None:
    assert [euler_phi(n) for n in range(1, 13)] == [
        1, 1, 2, 2, 4, 2, 6, 4, 6, 4, 10, 4
    ]
    for m in range(1, 10_001):
        assert 2 * euler_phi(m) ** 2 >= m


def test_rational_and_quadratic_angle_alphabets() -> None:
    assert finite_angle_alphabet(1) == [2, 3, 4, 6]
    assert finite_angle_alphabet(2) == [2, 3, 4, 5, 6, 8, 10, 12]


def test_k3_h9_endpoint_and_h10_obstruction() -> None:
    cert = build_certificate()["k3"]
    assert cert["H9_max_field_degree"] == 2
    assert cert["H9_compact_forces_field_degree"] == 2
    assert cert["H9_finite_angle_alphabet"] == [2, 3, 4, 5, 6, 8, 10, 12]
    assert cert["H10_max_field_degree_from_NS"] == 1
    assert cert["H10_plus_compact_full_trace_possible"] is False

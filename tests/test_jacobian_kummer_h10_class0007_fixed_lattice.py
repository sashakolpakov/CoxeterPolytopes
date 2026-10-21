from scripts.certify_jacobian_kummer_h10_class0007_fixed_lattice import (
    build_certificate,
)


def test_h10_fixed_lattice_identification() -> None:
    certificate = build_certificate()
    assert certificate["all_checks_pass"]
    assert certificate["fixed_lattice"]["determinant"] == 1024
    assert certificate["target_identification"]["lattice"] == "U + E8(-2) + <-4>"


def test_h10_primitive_mirror_types() -> None:
    certificate = build_certificate()
    assert certificate["facet_reflections"]["primitive_square_divisibility_histogram"] == {
        "-4,4": 43,
        "-4,2": 33,
        "-2,1": 8,
    }

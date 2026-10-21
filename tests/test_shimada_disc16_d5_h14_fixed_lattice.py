from __future__ import annotations

import json

from scripts.certify_shimada_disc16_d5_h14_fixed_lattice import OUTPUT


def _load() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_saturated_fixed_lattice_and_discriminant_form() -> None:
    certificate = _load()
    lattice = certificate["saturated_fixed_lattice"]

    assert lattice["definition"] == "Fix(M) = ker_Z(M^T-I) for Shimada row action M"
    assert lattice["constraint_smith_diagonal"] == [1] * 5 + [0] * 15
    assert lattice["rank"] == 15
    assert lattice["saturated"] is True
    assert lattice["basis_change_determinant"] in (-1, 1)
    assert lattice["inertia"] == [1, 14, 0]
    assert lattice["determinant"] == 128
    assert lattice["even"] is True
    assert lattice["smith_diagonal"] == [1] * 8 + [2] * 7
    assert lattice["discriminant_group"] == "(Z/2Z)^7"
    assert lattice["two_elementary_invariants"] == {"r": 15, "a": 7, "delta": 1}
    assert lattice["standard_isometry_type"] == "U(2) + E8(-1) + A1(-1)^5"

    discriminant = certificate["discriminant_form"]
    assert len(discriminant["generator_bit_vectors"]) == 7
    assert all(len(vector) == 15 for vector in discriminant["generator_bit_vectors"])
    assert discriminant["all_generators_have_order_two"] is True
    assert discriminant["delta"] == 1


def test_all_primitive_mirrors_define_integral_reflections() -> None:
    certificate = _load()
    mirrors = certificate["mirrors"]

    assert mirrors["facet_count"] == 60
    assert mirrors["primitive_type_histogram"] == {
        "square_-2_div_1": 15,
        "square_-2_div_2": 35,
        "square_-4_div_2": 10,
    }
    assert len(mirrors["records"]) == 60
    assert len(mirrors["reflection_matrices_in_reduced_basis"]) == 60
    assert mirrors["all_reflections_integral"] is True
    assert mirrors["all_reflections_preserve_fixed_gram"] is True
    assert all(record["reflection_integral"] for record in mirrors["records"])

    reflectivity = certificate["reflectivity"]
    assert reflectivity["integral_reflection_group"] is True
    assert reflectivity["finite_covolume"] is True
    assert reflectivity["fixed_lattice_reflective"] is True
    assert reflectivity["does_not_establish_literature_novelty"] is True

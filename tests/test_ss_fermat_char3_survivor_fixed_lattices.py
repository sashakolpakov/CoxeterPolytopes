from __future__ import annotations

import json

from scripts.certify_ss_fermat_char3_survivor_fixed_lattices import OUTPUT


def _payload() -> dict[str, object]:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_four_saturated_fixed_lattice_fingerprints() -> None:
    payload = _payload()
    assert payload["all_checks_pass"] is True
    cases = {case["id"]: case for case in payload["cases"]}
    assert {
        case_id: (case["hyperbolic_dimension"], case["determinant_absolute"], case["smith_diagonal_absolute"])
        for case_id, case in cases.items()
    } == {
        "char3-C2-H11": (11, 1_024, [1, 1, *([2] * 10)]),
        "char3-C2-H10": (10, 6_144, [*([2] * 10), 6]),
        "char3-V4-H11": (11, 768, [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 12]),
        "char3-V4-H10": (10, 1_024, [1, 1, *([2] * 8), 4]),
    }
    assert all(case["smith_kernel"]["saturated_smith_kernel"] for case in cases.values())


def test_v4_h10_is_the_known_vinberg_lattice() -> None:
    identification = _payload()["v4_h10_identification"]
    assert identification["target_lattice"] == "U + E8(-2) + <-4>"
    assert identification["candidate_discriminant_form"]["canonical_form"] == (
        "q_plus(F2^8) + <-1/4>"
    )
    assert identification["target_discriminant_form"]["canonical_form"] == (
        "q_plus(F2^8) + <-1/4>"
    )
    assert identification["rank"] == identification["discriminant_group_length"] + 2
    assert identification["explicit_isometry_determinant"] == 1
    assert len(identification["explicit_isometry_matrix"]) == 11
    assert identification["integrally_isometric"] is True


def test_c2_h10_half_scaled_lattice_is_mcleod_f3_h10() -> None:
    identification = _payload()["c2_h10_identification"]
    assert identification["candidate_gram_scalar"] == 2
    assert identification["target_lattice"] == "<3> + <-1>^10"
    assert identification["explicit_isometry_determinant"] == 1
    assert len(identification["explicit_isometry_matrix"]) == 11
    assert identification["integrally_isometric_after_dividing_form_by_2"] is True
    assert identification["orthogonal_group_unchanged_by_scalar_2"] is True
    assert identification["maximal_reflection_chamber_lineage"] == (
        "McLeod f_3^10, 14 walls"
    )


def test_all_196_mirrors_are_integral_roots_of_mcleod_lattice() -> None:
    mirrors = _payload()["c2_h10_196_mirrors"]
    assert mirrors["facet_count"] == 196
    assert mirrors["half_scaled_square_divisibility_histogram"] == {
        "-6,3": 80,
        "-2,1": 36,
        "-1,1": 80,
    }
    assert mirrors["original_even_square_divisibility_histogram"] == {
        "-12,6": 80,
        "-4,2": 36,
        "-2,2": 80,
    }
    assert mirrors["all_reflections_integral_involutive_isometries"] is True
    assert len(
        {tuple(row["mcleod_lattice_coordinates"]) for row in mirrors["mirrors"]}
    ) == 196


def test_all_123_mirrors_are_integral_roots_of_the_known_lattice() -> None:
    mirrors = _payload()["v4_h10_123_mirrors"]
    assert mirrors["facet_count"] == 123
    assert mirrors["primitive_square_divisibility_histogram"] == {
        "-4,2": 31,
        "-4,4": 82,
        "-2,1": 10,
    }
    assert mirrors["all_reflections_integral_involutive_isometries"] is True
    assert len({tuple(row["target_lattice_coordinates"]) for row in mirrors["mirrors"]}) == 123

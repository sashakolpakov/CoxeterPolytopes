from fractions import Fraction

from scripts.certify_rational_nef_reflection_integrality_gate import build_certificate


def certificate():
    return build_certificate()


def test_unimodular_sieve_is_exact():
    rows = certificate()["unimodular_sieve"]
    survivors = {
        (row["primitive_normal_square"], row["curve_divisibility"], row["curve_square"])
        for row in rows
        if row["survives"]
    }
    assert survivors == {(-1, 1, -1), (-2, 1, -2), (-1, 2, -4)}


def test_f4_full_nef_wall_gram_and_bad_coble_incidence():
    checks = certificate()["checks"]
    assert checks["ns_unimodular"]
    assert checks["physical_wall_gram"] == [[-4, 0, 1], [0, -1, 1], [1, 1, -1]]
    assert checks["boundary_incidence_with_minus_one_walls"] == {"E": 0, "G=F-E": 1}
    assert checks["B_not_minus_2K"]


def test_half_boundary_lattice_and_reflections():
    checks = certificate()["checks"]
    assert checks["B_reflection_not_integral_on_NS"]
    assert checks["half_lattice_determinant"] == 2
    assert checks["half_lattice_2_elementary_discriminant_order"]
    assert checks["root_span_unimodular_coordinates"]
    assert all(
        row[gate]
        for row in checks["reflection_checks"].values()
        for gate in ("integral", "involution", "isometry")
    )
    assert checks["log_class_root_pairings"] == {"h=B/2": 0, "E": -2, "G=F-E": -1}
    assert checks["boundary_reflection_fixes_log_class"]
    assert not checks["all_reflections_fix_log_class"]


def test_triangle_is_reflection_maximal_noncompact_minimum():
    checks = certificate()["checks"]
    assert checks["coxeter_labels"] == {"h,E": 2, "h,G": 3, "E,G": "infinity"}
    assert checks["triangle_area_over_pi"] == "1/6"
    assert checks["noncompact_coxeter_polygon_minimum_area_over_pi"] == "1/6"
    assert checks["reflection_maximal_by_area"]


def test_arbitrary_denominators_do_not_select_m2():
    rows = certificate()["arbitrary_denominator_family"]
    assert [row["m"] for row in rows] == list(range(2, 21))
    assert all(row["integral_isometry"] for row in rows)
    assert {tuple(tuple(x for x in line) for line in row["extension_gram"]) for row in rows} == {
        ((-2, 1), (1, 0))
    }


def test_arbitrary_commensurable_sublattices_allow_odd_squares():
    rows = certificate()["arbitrary_square_sublattice_family"]
    assert [row["physical_section_square"] for row in rows] == list(range(-1, -21, -1))
    assert all(row["integral_isometry"] for row in rows)
    assert all(row["reflection"] == [[-1, 2], [0, 1]] for row in rows)


def test_weighted_cover_noether_sieve():
    rows = certificate()["weighted_cover_samples"]
    assert certificate()["checks"]["weighted_cover_only_m2_survives_samples"]
    assert certificate()["checks"]["weighted_term_factorization_checked"]
    assert certificate()["checks"]["weighted_term_nonnegative_unique_zero_m2_checked"]
    assert next(row for row in rows if row["ramification_indices"] == [2])["c2"] == 24
    assert next(row for row in rows if row["ramification_indices"] == [2, 2])["c2"] == 24
    assert next(row for row in rows if row["ramification_indices"] == [3])["c2"] == 40

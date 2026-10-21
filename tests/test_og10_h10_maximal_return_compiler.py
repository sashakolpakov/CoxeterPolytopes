from scripts.certify_og10_h10_maximal_return_compiler import (
    binary_forms_are_gl2_isometric,
    binary_forms_are_pairwise_gl2_distinct,
    build_certificate,
    determinant_lanes,
    reduced_even_positive_binary_forms,
)


def test_reduced_even_binary_form_censuses() -> None:
    forms_12 = reduced_even_positive_binary_forms(12)
    forms_60 = reduced_even_positive_binary_forms(60)
    forms_132 = reduced_even_positive_binary_forms(132)
    forms_108 = reduced_even_positive_binary_forms(108)
    assert forms_12 == [
        [1, 0, 3],
        [2, 2, 2],
    ]
    assert forms_60 == [
        [1, 0, 15],
        [2, 2, 8],
        [3, 0, 5],
        [4, 2, 4],
    ]
    assert forms_132 == [
        [1, 0, 33],
        [2, 2, 17],
        [3, 0, 11],
        [6, 6, 7],
    ]
    assert forms_108 == [
        [1, 0, 27],
        [2, 2, 14],
        [3, 0, 9],
        [4, 2, 7],
        [6, 6, 6],
    ]
    assert all(
        binary_forms_are_pairwise_gl2_distinct(forms)
        for forms in (forms_12, forms_60, forms_108, forms_132)
    )
    assert binary_forms_are_gl2_isometric([2, -2, 8], [2, 2, 8])


def test_determinant_lane_filters() -> None:
    assert [
        (lane["glue_order"], lane["complement_determinant"])
        for lane in determinant_lanes(44)
    ] == [(44, 132)]
    assert [
        (lane["glue_order"], lane["complement_determinant"])
        for lane in determinant_lanes(20)
    ] == [(20, 60)]
    assert [
        (lane["glue_order"], lane["complement_determinant"])
        for lane in determinant_lanes(36)
    ] == [(12, 12), (36, 108)]
    assert [
        (lane["glue_order"], lane["complement_determinant"])
        for lane in determinant_lanes(4)
    ] == [(4, 12)]


def test_four_exact_rank22_replays_and_scope() -> None:
    certificate = build_certificate()
    assert all(certificate["checks"].values())

    cases = certificate["known_rank22_replays"]
    assert cases["golden_det44"]["smith_nonunits"] == [2, 22]
    assert cases["a2_coupled_det20"]["smith_nonunits"] == [2, 10]
    assert cases["norm_six_det36"]["smith_nonunits"] == [6, 6]
    assert cases["borcherds_d4_det4"]["smith_nonunits"] == [2, 2]

    assert cases["golden_det44"][
        "divisibility_three_picard_walls_automatically_absent"
    ]
    assert cases["a2_coupled_det20"][
        "divisibility_three_picard_walls_automatically_absent"
    ]
    assert not cases["norm_six_det36"][
        "divisibility_three_picard_walls_automatically_absent"
    ]

    assert "no n<=9 conclusion for NS-only carriers" in certificate["scope"][
        "not_proved"
    ]
    assert certificate["checks"][
        "coherent_full_h2_quadratic_bound_is_H9_on_OG10"
    ]

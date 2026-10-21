"""Independent regression tests for the native-host hostile audit."""

from __future__ import annotations

from fractions import Fraction

from scripts.certify_ag_field_funnel_bounded_frame_native_host_closure_independent_audit import (
    OWNER_CERTIFIER,
    OWNER_REPORT,
    OWNER_TEST,
    ROOT,
    build_certificate,
    cycle_square_from_edges,
    diagonal_rescale,
    directed_cycle_product,
    generated_algebra_dimension,
    recover_tree_scales,
    reflection_matrix,
    sha256_file,
    small_conductor_sandwich,
    core_sum_sample,
    trace_cycle_product,
)


OWNER_REPORT_SHA256 = "7ba392dd833358379ea92d1385909a52b2dc65417a9cb4ccb856170118df02a3"
OWNER_CERTIFIER_SHA256 = "44848ad5727458b1777a4989395229d847ad0e05136ae5c96d2ffd9ad12a6def"
OWNER_TEST_SHA256 = "e5b0f259ccaee50ec6c1d9d17af9d5806eed999f28c5712c330cabbc03c5ed45"
INDEPENDENT_CERTIFIER = (
    ROOT
    / "scripts/certify_ag_field_funnel_bounded_frame_native_host_closure_independent_audit.py"
)
INDEPENDENT_CERTIFIER_SHA256 = "ff96e2bd31d273ce79e9ac9ba6bf694f37564fc0e83a90f02985c368593039b9"


def lorentzian_symmetric_cartan():
    return tuple(
        tuple(Fraction(value) for value in row)
        for row in (
            (2, -1, -1, 0),
            (-1, 2, -1, -1),
            (-1, -1, 2, -1),
            (0, -1, -1, 2),
        )
    )


def test_frozen_owner_report_certifier_and_test_pins():
    assert sha256_file(OWNER_REPORT) == OWNER_REPORT_SHA256
    assert sha256_file(OWNER_CERTIFIER) == OWNER_CERTIFIER_SHA256
    assert sha256_file(OWNER_TEST) == OWNER_TEST_SHA256
    assert sha256_file(INDEPENDENT_CERTIFIER) == INDEPENDENT_CERTIFIER_SHA256


def test_projective_trace_and_symmetric_cycle_identities():
    cartan = lorentzian_symmetric_cartan()
    cycles = ((0, 1), (1, 2), (0, 2), (1, 3), (2, 3), (0, 1, 2), (1, 2, 3))
    for cycle in cycles:
        product = directed_cycle_product(cartan, cycle)
        assert trace_cycle_product(cartan, cycle) == product
        assert directed_cycle_product(cartan, tuple(reversed(cycle))) == product
        if len(cycle) >= 3:
            assert product**2 == cycle_square_from_edges(cartan, cycle)


def test_cycle_packet_recovers_every_non_tree_entry_by_diagonal_gauge():
    cartan = lorentzian_symmetric_cartan()
    scales = (Fraction(2), Fraction(3), Fraction(5), Fraction(7))
    changed = diagonal_rescale(cartan, scales)
    recovered = recover_tree_scales(cartan, changed, ((0, 1), (1, 3), (3, 2)))
    assert recovered == (
        Fraction(1),
        Fraction(3, 2),
        Fraction(5, 2),
        Fraction(7, 2),
    )
    assert diagonal_rescale(cartan, recovered) == changed


def test_connected_reflections_generate_the_full_endomorphism_algebra_sample():
    cartan = lorentzian_symmetric_cartan()
    reflections = tuple(reflection_matrix(cartan, index) for index in range(4))
    assert generated_algebra_dimension(reflections) == 4**2


def test_full_order_conductor_sandwich_and_core_sum_samples():
    sandwich = small_conductor_sandwich()
    assert sandwich["lambda_stable"] is True
    assert sandwich["two_BL_inside_L_inside_BL"] is True
    assert sandwich["BL_steps"] == [1, 1]
    assert sandwich["L_steps"] == [2, 1]
    assert sandwich["2BL_steps"] == [2, 2]

    core_sum = core_sum_sample()
    assert core_sum["translate_sum_steps"] == [1, 1]
    assert core_sum["translate_sum_is_swap_stable"] is True


def test_independent_composite_certificate_and_dimension_quantifier():
    certificate = build_certificate()
    assert certificate["all_checks_pass"] is True
    assert certificate["verdict"] == "PASS"
    assert certificate["burnside_sample"] == {
        "rank": 4,
        "generated_algebra_dimension": 16,
        "full_endomorphism_dimension": 16,
    }
    assert certificate["dimension_quantifier"] == {
        "dimensions": [2, 995],
        "ranks": [3, 996],
        "rank_count": 994,
    }

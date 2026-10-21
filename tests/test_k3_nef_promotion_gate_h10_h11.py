from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_k3_nef_promotion_gate_h10_h11 import (
    OUTPUT,
    ROOT,
    build_certificate,
)


REPORT = ROOT / "results" / "k3_nef_promotion_gate_h10_h11.md"


def archive() -> dict:
    return json.loads(Path(OUTPUT).read_text(encoding="utf-8"))


def test_archive_rebuilds_exactly() -> None:
    assert build_certificate() == archive()


def test_exact_root_hyperplane_promotion_criterion_is_recorded() -> None:
    criterion = archive()["promotion_criterion"]
    assert "if and only if" in criterion["exact_if_and_only_if"]
    assert "all >=0 or all <=0" in criterion["extreme_ray_form"]
    assert "not all zero" in criterion["extreme_ray_form"]
    assert "strict timelike points" in criterion["finite_failure_witness"]
    assert "primitive square -4 normal can never qualify" in criterion[
        "full_rank_facet_corollary"
    ]


def test_h10_16_unique_saturation_has_an_exact_strict_cut() -> None:
    row = archive()["cases"]["H10_16"]
    root = row["extra_ambient_root"]
    assert row["admissible_nontrivial_saturation_count"] == 1
    assert root["square"] == -2
    assert root["primitive"] is True
    assert root["lies_in_M"] is False
    assert root["action_orbit_size"] == 2
    assert root["inner_product_with_nontrivial_image"] == 0
    assert row["extreme_ray_test"]["pairing_histogram"] == {
        "-7": 1,
        "-5": 3,
        "-4": 2,
        "-3": 4,
        "-2": 5,
        "-1": 8,
        "0": 30,
        "1": 8,
        "2": 5,
        "3": 4,
        "4": 2,
        "5": 3,
        "7": 1,
    }
    zero = row["strict_cut_certificate"]["strict_zero_point"]
    assert zero["coordinates"] == [
        296,
        610,
        966,
        457,
        889,
        854,
        835,
        870,
        281,
        321,
        0,
    ]
    assert zero["square"] == 422572
    assert zero["minimum_old_wall_pairing"] == 16
    assert zero["extra_root_pairing"] == 0
    assert row["decision"] == (
        "not_full_invariant_nef_section_in_the_unique_rank12_saturation"
    )


def test_h10_45_all_ten_admissible_overextensions_are_strictly_cut() -> None:
    row = archive()["cases"]["H10_45"]
    extensions = row["extensions"]
    assert row["admissible_rank12_action_compatible_overextensions"] == 10
    assert row["index_histogram"] == {"2": 7, "4": 3}
    assert len(extensions) == 10
    assert [entry["subgroup_order_and_ambient_index"] for entry in extensions] == [
        *[2] * 7,
        *[4] * 3,
    ]
    assert [entry["extra_ambient_root"]["action_orbit_size"] for entry in extensions] == [
        1,
        1,
        1,
        2,
        1,
        1,
        1,
        1,
        1,
        1,
    ]
    assert [
        entry["strict_cut_certificate"]["strict_zero_point"]["square"]
        for entry in extensions
    ] == [
        2613664,
        "1463/4",
        "518399/4",
        668,
        194940,
        "440899/2",
        655486,
        38141,
        "440899/2",
        247775,
    ]
    for entry in extensions:
        assert entry["extra_ambient_root"]["square"] == -2
        assert entry["extra_ambient_root"]["primitive"] is True
        assert entry["extra_ambient_root"]["lies_in_M"] is False
        assert entry["extra_ambient_root"]["same_restriction_on_fixed_space"] is True
        zero = entry["strict_cut_certificate"]["strict_zero_point"]
        assert zero["timelike"] is True
        assert zero["fixed_by_action"] is True
        assert zero["strict_for_all_59_old_walls"] is True
        assert zero["extra_root_pairing"] == 0
        assert entry["all_checks_pass"] is True
    assert row["decision"] == (
        "no_full_invariant_nef_section_in_any_admissible_rank12_overextension"
    )


def test_h11_29_has_a_marking_independent_minus_four_facet_obstruction() -> None:
    row = archive()["cases"]["H11_29"]
    assert row["full_rank_in_ambient_positive_cone"] is True
    assert row["primitive_facet_normal_square_histogram"] == {
        "-2": 17,
        "-4": 12,
    }
    assert row["obstructing_square_minus_four_facets_one_based"] == [
        1,
        2,
        7,
        8,
        9,
        10,
        13,
        19,
        21,
        25,
        28,
        29,
    ]
    assert row["marking_independent"] is True
    assert "c^2=1/2, impossible" in row["hyperplane_argument"]
    assert row["decision"] == "not_a_full_nef_cone_for_the_rank12_ambient_NS"


def test_positive_subsystem_results_are_retained_and_scope_is_bounded() -> None:
    payload = archive()
    positive = payload["positive_context_retained"]
    assert positive["divisor_subsystem_realizations_remain_valid"] is True
    assert positive["finite_volume_certificates_remain_valid"] is True
    assert "H10/12" in positive["ambient_complete_nef_replacements"][
        "U+E8(-1)+A1(-1)^2"
    ]
    assert "H10/15" in positive["ambient_complete_nef_replacements"][
        "U+D8(-1)+A1(-1)^2"
    ]
    assert "higher-Picard" in payload["scope"]["not_proved"]
    assert payload["all_checks_pass"] is True


def test_report_states_the_result_without_overclaiming_promotion() -> None:
    report = REPORT.read_text(encoding="utf-8")
    normalized = " ".join(report.split())
    assert "None of the three new divisor/reflection carriers promotes" in normalized
    assert "if and only if no ambient" in normalized
    assert "every one of the ten" in normalized
    assert "This is a no-promotion theorem, not a retraction" in normalized
    assert "not a worldwide nonexistence theorem" in normalized

from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_r12_h10_45_bh_divisor_sublattice_carrier import (
    OUTPUT,
    ROOT,
    build_certificate,
)


REPORT = ROOT / "results" / "r12_h10_45_bh_divisor_sublattice_carrier.md"


def archive() -> dict:
    return json.loads(Path(OUTPUT).read_text(encoding="utf-8"))


def test_archive_rebuilds_exactly() -> None:
    assert build_certificate() == archive()


def test_two_genuine_order_four_k3_carriers_are_frozen() -> None:
    rows = archive()["carriers"]
    assert [row["carrier_id"] for row in rows] == [
        "BH-0.4.4.2-index2-H10-45",
        "BH-0.4.4.1-index4-H10-45",
    ]
    assert [row["ambient_NS"]["lattice"] for row in rows] == [
        "U+D8(-1)+A1(-1)^2",
        "U+E8(-1)+A1(-1)^2",
    ]
    assert [row["divisor_sublattice"]["ambient_index"] for row in rows] == [
        2,
        4,
    ]
    assert [
        row["divisor_sublattice"][
            "canonical_BH_coordinates_required_even_one_based"
        ]
        for row in rows
    ] == [[10], [10, 12]]
    for row in rows:
        assert row["brandhorst_hofmann"]["effective_group"] == "C4"
        assert row["brandhorst_hofmann"]["Picard_image"] == "C2"
        assert row["ambient_NS"]["Picard_coinvariant_lattice"] == "<-4>"
        assert row["ambient_NS"]["Picard_coinvariant_rootless"] is True
        assert row["explicit_weierstrass_carrier"][
            "global_automorphism_order"
        ] == 4
        assert row["explicit_weierstrass_carrier"]["very_general_NS_exact"] is True
        assert row["all_checks_pass"] is True


def test_integral_action_pair_conjugacy_uses_all_discriminant_and_glue_data() -> None:
    for row in archive()["carriers"]:
        pair = row["exact_H10_45_action_pair"]
        discr = pair["fixed_discriminant_isometry"]
        assert discr["source_invariant_factors"] == [2, 2, 2, 2, 4]
        assert discr["target_invariant_factors"] == [2, 2, 2, 2, 4]
        assert discr["elements_checked"] == 64
        assert discr["bijective"] is True
        assert discr["quadratic_form_preserved_mod_2Z"] is True
        assert pair["all_eigen_glue_pairs_carried"] is True
        assert pair["nikulin_uniqueness_and_surjectivity_apply"] is True
        assert pair["full_integral_lattice_action_pairs_conjugate"] is True


def test_archived_invariant_interior_and_effectivity_scope_are_literal() -> None:
    payload = archive()
    target = payload["archived_target_pair"]
    assert target["invariant_integral_interior_vector"] == [
        3,
        3,
        3,
        0,
        0,
        2,
        2,
        6,
        6,
        3,
        1,
        7,
    ]
    assert target["invariant_interior_square"] == 40
    assert target["invariant_interior_wall_pairing_minimum"] == 1
    effectivity = payload["algebro_geometric_effectivity"]
    assert effectivity["all_59_effective_signs_can_be_chosen_simultaneously"] is True
    assert effectivity["strong_Torelli_produces_effective_chamber_automorphism"] is True
    assert effectivity["formula_sigma_need_not_equal_effective_tau"] is True
    assert "divisor-subsystem" in effectivity["scope"]
    assert "not claimed to be the full nef chamber" in effectivity["scope"]


def test_finite_index_carrier_criterion_is_exhaustive_at_orders_two_and_four() -> None:
    census = archive()["finite_index_carrier_criterion"]
    assert census["fixed_discriminant_group_order"] == 64
    assert census["admissible_order_two_subgroups"] == 7
    assert census["admissible_order_four_subgroups"] == 3
    assert census["all_order_two_extensions_land_in"] == "U+D8(-1)+A1(-1)^2"
    assert census["all_order_four_extensions_land_in"] == "U+E8(-1)+A1(-1)^2"
    assert census["checks"]["all_extensions_even_hyperbolic"] is True
    assert census["checks"][
        "every_order_two_extension_has_row_0_4_4_2_form"
    ] is True
    assert census["checks"][
        "every_order_four_extension_has_row_0_4_4_1_form"
    ] is True


def test_vinberg_handoff_proves_finite_volume_and_noncompactness() -> None:
    volume = archive()["vinberg_handoff"]
    assert (volume["dimension"], volume["facets"]) == (10, 45)
    assert volume["classification"] == "finite_volume_noncompact"
    assert volume["finite_volume"] is True
    assert volume["compact"] is False
    assert (volume["elliptic_edges"], volume["finite_vertices"], volume["cusps"]) == (
        7332,
        1048,
        93,
    )
    assert volume["connected_affine_components_audited"] == 302
    assert volume["bad_edge_endpoints"] == 0
    assert volume["bad_affine_zero_neighborhoods"] == 0


def test_report_states_the_positive_result_and_its_boundary() -> None:
    report = REPORT.read_text(encoding="utf-8")
    normalized = " ".join(report.split())
    assert "TOOT" in report
    assert "two genuine algebro-geometric carriers" in normalized
    assert "effective invariant **divisor-subsystem**" in report
    assert "does not call the 59-root chamber the full nef chamber" in normalized
    assert "seven admissible subgroups of order 2" in normalized
    assert "three admissible subgroups of order 4" in normalized

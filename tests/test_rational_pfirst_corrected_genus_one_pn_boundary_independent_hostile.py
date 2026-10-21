"""Independent hostile tests for the corrected genus-one P_n boundary."""

from __future__ import annotations

import json

from scripts.audit_rational_pfirst_corrected_genus_one_pn_boundary_independent_hostile import (
    DEPENDENCY_PINS,
    OWNER_PINS,
    ROOT,
    abel_record,
    build_audit,
    mmp_orbit_record,
    positive_cone_records,
    rank_blind_record,
    sha256_file,
)


def load(path: str) -> dict[str, object]:
    return json.loads((ROOT / path).read_text())


def test_owner_and_every_cone_mmp_dependency_are_byte_pinned() -> None:
    for item in (*OWNER_PINS.values(), *DEPENDENCY_PINS.values()):
        assert sha256_file(ROOT / item["path"]) == item["sha256"]


def test_corrected_local_packet_is_rank_blind_well_beyond_the_named_family() -> None:
    records = [rank_blind_record(n) for n in range(12, 65)]
    assert {
        (record["C_square"], record["K_dot_C"], record["arithmetic_genus"], record["reflection_integral"])
        for record in records
    } == {(-2, 2, 1, True)}
    assert [records[0]["Z_square"], records[-1]["Z_square"]] == [-1, -53]


def test_positive_full_nef_range_uses_complete_cone_certificates() -> None:
    source = load(DEPENDENCY_PINS["inverse_archive"]["path"])
    records = positive_cone_records(source)
    assert [record["n"] for record in records] == [12, 13, 14, 15, 16]
    assert [record["facet_count"] for record in records] == [14, 15, 17, 18, 20]
    assert [record["primitive_candidate_ray_count"] for record in records] == [
        None,
        None,
        94,
        145,
        325,
    ]
    assert all(record["complete_nef_carrier"] for record in records)


def test_all_negative_markings_are_exact_unimodular_diagram_orbits() -> None:
    source = load(DEPENDENCY_PINS["inverse_archive"]["path"])
    p17 = load(DEPENDENCY_PINS["p17_mmp_archive"]["path"])
    p18_p19 = load(DEPENDENCY_PINS["p18_p19_mmp_archive"]["path"])
    records = [mmp_orbit_record(n, source, p17, p18_p19) for n in (17, 18, 19)]
    assert [record["canonical_class_count"] for record in records] == [8, 24, 120]
    assert [record["diagram_orbit_count"] for record in records] == [8, 24, 120]
    assert all(record["unimodular_facet_basis"] for record in records)
    assert all(record["P2_F2_reconstructions_agree"] for record in records)


def test_piecewise_Abel_quartic_formula_and_residual_degree_are_exact() -> None:
    source = load(DEPENDENCY_PINS["inverse_archive"]["path"])
    records = [abel_record(n, source) for n in (17, 18, 19)]
    assert [record["A_chain_length"] for record in records] == [17, 18, 18]
    assert records[0]["A"] == [4] + [-1] * 17
    assert records[1]["A"] == [4] + [-1] * 18
    assert records[2]["A"] == [4] + [-1] * 18 + [0]
    assert all(record["C_dot_A"] == record["residual_degree_on_C"] == 1 for record in records)


def test_full_verdict_preserves_local_family_and_universal_scopes() -> None:
    audit = build_audit()
    assert audit["verdict"] == "PASS_AFTER_TWO_SCOPE_REPAIRS"
    assert audit["scoped_carrier_iff"].endswith("iff n<=16")
    assert audit["scope"]["universal_n_le_16_claimed"] is False
    assert audit["scope"]["local_packet_gives_dimension_bound"] is False
    assert audit["scope"]["negative_no_go_characteristic_zero"] is True
    assert "Pic^1" in audit["scope"]["property_4_gate"]

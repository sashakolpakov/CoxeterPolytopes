from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_a20_pair_full_picard_k3_carriers import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def _cases(certificate: dict) -> dict[str, dict]:
    return {case["id"]: case for case in certificate["cases"]}


def test_both_generic_full_picard_carriers_are_frozen(certificate: dict) -> None:
    cases = _cases(certificate)
    assert set(cases) == {"A20_2_1_C01", "A20_2_1_C02"}
    assert cases["A20_2_1_C01"]["ledger_id"] == "A20-C2a-H12-19"
    assert cases["A20_2_1_C02"]["ledger_id"] == "A20-C2b-H11-20"
    assert [
        case["carrier_lattice"]["very_general_period_family"]["dimension"]
        for case in certificate["cases"]
    ] == [7, 8]
    assert all(
        case["carrier_lattice"]["very_general_period_family"]["full_Picard"]
        for case in certificate["cases"]
    )


def test_primitive_k3_embeddings_and_complements(certificate: dict) -> None:
    c01 = _cases(certificate)["A20_2_1_C01"]["carrier_lattice"]
    c02 = _cases(certificate)["A20_2_1_C02"]["carrier_lattice"]
    assert (c01["rank"], c01["signature"], c01["determinant"]) == (
        13,
        [1, 12, 0],
        256,
    )
    assert c01["smith_invariants"] == [1] * 6 + [2] * 6 + [4]
    assert (c02["rank"], c02["signature"], c02["determinant"]) == (
        12,
        [1, 11, 0],
        -1024,
    )
    assert c02["smith_invariants"] == [1] * 2 + [2] * 10
    for lattice, complement_signature, index in (
        (c01, [2, 7, 0], 256),
        (c02, [2, 8, 0], 1024),
    ):
        embedding = lattice["primitive_embedding_in_Lambda_K3"]
        assert embedding["orthogonal_complement_signature"] == complement_signature
        assert embedding["combined_gluing_index"] == index
        assert all(embedding["checks"].values())


def test_every_wall_is_primitive_integral_and_spans_full_picard(
    certificate: dict,
) -> None:
    cases = _cases(certificate)
    expected = {
        "A20_2_1_C01": (12, 19, {"-2,1": 7, "-4,2": 12}),
        "A20_2_1_C02": (
            11,
            20,
            {"-2,1": 3, "-2,2": 6, "-4,2": 11},
        ),
    }
    for case_id, (dimension, facets, types) in expected.items():
        chamber = cases[case_id]["oriented_chamber"]
        assert chamber["hyperbolic_dimension"] == dimension
        assert chamber["facet_count"] == len(chamber["wall_records"]) == facets
        assert chamber["wall_span_index"] == 1
        assert chamber["wall_span_smith"] == [1] * (dimension + 1)
        assert chamber["primitive_mirror_type_histogram"] == types
        assert all(
            all(record["checks"].values())
            for record in chamber["wall_records"]
        )
        assert all(chamber["checks"].values())


def test_strict_positive_classes_fix_both_orientations(certificate: dict) -> None:
    c01 = _cases(certificate)["A20_2_1_C01"]["oriented_chamber"]
    c02 = _cases(certificate)["A20_2_1_C02"]["oriented_chamber"]
    assert c01["strict_interior_square"] == 32068
    assert c01["strict_interior_pairings"] == [
        4, 8, 8, 8, 8, 8, 8, 180, 4, 180, 4, 4, 4, 4, 112, 396, 176,
        780, 780,
    ]
    assert c02["strict_interior_square"] == 5386
    assert c02["strict_interior_pairings"] == [
        4, 2, 4, 4, 2, 58, 4, 58, 4, 4, 4, 2, 2, 88, 230, 88, 218, 260,
        454, 540,
    ]
    assert c01["minimum_wall_pairing"] == 4
    assert c02["minimum_wall_pairing"] == 2


def test_exhaustive_denominator_audit_selects_unique_even_carriers(
    certificate: dict,
) -> None:
    for case in certificate["cases"]:
        audit = case["reflection_stable_overclosures"]
        exhaustive = audit["exhaustive_discriminant_module_audit"]
        maximum = audit["unique_maximal_integral_overclosure"]
        assert exhaustive["stable_integral_subgroup_orders"] == [1, 2]
        assert exhaustive["stable_even_subgroup_orders"] == [1]
        assert maximum["index"] == 2
        assert maximum["parity"] == "odd"
        assert "unique even" in audit["K3_admissible_conclusion"]
        assert all(audit["checks"].values())


def test_dual_exact_volume_proofs_have_no_hyperideal_rays(
    certificate: dict,
) -> None:
    expected = {
        "A20_2_1_C01": (527, 25, 3846, {"ff": 2589, "fi": 1146, "ii": 111}, 59),
        "A20_2_1_C02": (311, 29, 2229, {"ff": 1312, "fi": 797, "ii": 120}, 64),
    }
    for case_id, (finite, cusps, edges, endpoints, affine) in expected.items():
        volume = _cases(certificate)[case_id]["finite_volume"]
        assert volume["classification"] == "finite_volume_noncompact"
        assert (
            volume["finite_vertices"],
            volume["cusps"],
            volume["hyperideal_vertices"],
        ) == (finite, cusps, 0)
        assert volume["elliptic_rank_dimension_minus_one_edges"] == edges
        assert volume["endpoint_type_histogram"] == endpoints
        assert volume["connected_affine_components"] == affine
        assert volume["bad_edges"] == volume["bad_parabolic_neighborhoods"] == 0
        assert all(volume["checks"].values())


def test_scope_and_report_do_not_overclaim(certificate: dict) -> None:
    assert certificate["paired_gain"] == {
        "new_full_Picard_period_carriers": 2,
        "geometric_slot_count_change": 0,
        "reason": (
            "both chambers already had genuine special-period K3 "
            "automorphism carriers; the gain is deformation to full Picard"
        ),
    }
    for case in certificate["cases"]:
        scope = case["scope"]
        assert scope["not_new_abstract_polytope"] is True
        assert scope["not_equation_level"] is True
        assert "no claim" in scope["not_full_nef"]
        assert "not asserted" in scope["not_generic_surface_automorphism"]
        assert "not asserted" in scope["not_reflection_automorphisms"]
        assert case["all_checks_pass"] is True
    assert all(certificate["checks"].values())
    assert certificate["all_checks_pass"] is True

    report = (
        Path(__file__).resolve().parents[1]
        / "results/a20_pair_full_picard_k3_carriers.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "toot" in normalized
    assert "seven-dimensional" in normalized
    assert "eight-dimensional" in normalized
    assert "unique even" in normalized
    assert "no hyperideal" in normalized
    assert "not a new abstract-polytope" in normalized

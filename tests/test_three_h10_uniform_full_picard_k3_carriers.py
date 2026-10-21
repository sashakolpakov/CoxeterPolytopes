from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_three_h10_uniform_full_picard_k3_carriers import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def _chamber(certificate: dict, chamber_id: str) -> dict:
    return next(row for row in certificate["chambers"] if row["id"] == chamber_id)


def test_uniform_explicit_k3_graph_gluings(certificate: dict) -> None:
    embeddings = certificate["uniform_gluing_lemma"]["embeddings"]
    expected = {
        "D16-D5-H10-30": (2048, "T_1=<2>^2+<-2>+<-2>^8"),
        "SS2-2D-H10-47": (4096, "T_2=<2>^2+<-4>+<-2>^8"),
        "F3-C2-H10-196": (6144, "T_3=<2>^2+<-6>+<-2>^8"),
    }
    for chamber_id, (order, complement) in expected.items():
        embedding = embeddings[chamber_id]
        assert embedding["criterion"] == "explicit full-discriminant graph gluing"
        assert embedding["graph_order"] == order
        assert embedding["checked_discriminant_class_count"] == order
        assert embedding["orthogonal_complement"] == complement
        assert embedding["combined_carrier_complement_index"] == order
        assert embedding["two_primary_anti_isometry_matrix_mod_2"] == [
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 1, 1, 0, 0, 0, 0],
            [0, 0, 1, 0, 1, 1, 0, 0, 0, 0],
            [0, 0, 1, 1, 0, 1, 0, 0, 0, 0],
            [0, 0, 1, 1, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 1, 1, 1],
            [0, 0, 0, 0, 0, 0, 1, 0, 1, 1],
            [0, 0, 0, 0, 0, 0, 1, 1, 0, 1],
            [0, 0, 0, 0, 0, 0, 1, 1, 1, 0],
        ]
        assert all(embedding["checks"].values())


@pytest.mark.parametrize(
    ("chamber_id", "a", "facets", "types", "finite", "cusps", "edges"),
    [
        (
            "D16-D5-H10-30",
            1,
            30,
            {"-2,2": 18, "-4,2": 12},
            1333,
            207,
            10030,
        ),
        (
            "SS2-2D-H10-47",
            2,
            47,
            {"-2,2": 5, "-4,2": 16, "-4,4": 26},
            3136,
            216,
            19910,
        ),
        (
            "F3-C2-H10-196",
            3,
            196,
            {"-12,6": 80, "-2,2": 80, "-4,2": 36},
            212824,
            10560,
            1427720,
        ),
    ],
)
def test_full_picard_wall_and_volume_replay(
    certificate: dict,
    chamber_id: str,
    a: int,
    facets: int,
    types: dict[str, int],
    finite: int,
    cusps: int,
    edges: int,
) -> None:
    chamber = _chamber(certificate, chamber_id)
    carrier = chamber["carrier"]
    walls = chamber["oriented_wall_system"]
    assert chamber["dimension"] == 10
    assert chamber["facet_count"] == len(walls["wall_records"]) == facets
    assert carrier["name"] == f"<{2*a}>+<-2>^10"
    assert carrier["signature"] == [1, 10, 0]
    assert carrier["determinant_absolute"] == (2 * a) * 2**10
    assert carrier["discriminant"]["length"] == 11
    assert carrier["very_general_period_consequence"]["moduli_dimension"] == 9
    assert walls["primitive_mirror_type_histogram"] == types
    assert walls["wall_span"]["rank"] == 11
    assert walls["wall_span"]["index_in_carrier_lattice"] == 1
    assert all(all(record["checks"].values()) for record in walls["wall_records"])
    vinberg = chamber["finite_volume"]["literal_Vinberg"]
    cone = chamber["finite_volume"]["exact_GMP_cone"]
    assert vinberg["classification"] == cone["classification"] == (
        "finite_volume_noncompact"
    )
    assert vinberg["elliptic_rank9_edges"] == edges
    assert vinberg["finite_vertices"] == cone["finite_vertices"] == finite
    assert vinberg["cusps"] == cone["cusps"] == cusps
    assert vinberg["bad_edges"] == vinberg["bad_parabolic_neighborhoods"] == 0
    assert cone["hyperideal_rays"] == 0


def test_exact_interior_provenance(certificate: dict) -> None:
    d5 = _chamber(certificate, "D16-D5-H10-30")
    ss2 = _chamber(certificate, "SS2-2D-H10-47")
    f3 = _chamber(certificate, "F3-C2-H10-196")
    assert d5["interior_provenance"]["source_orientation_signs_match"] is True
    assert d5["oriented_wall_system"]["strict_interior_square"] == 52
    assert ss2["oriented_wall_system"]["strict_interior_square"] == 11390
    assert ss2["finite_volume"]["exact_GMP_cone"]["extreme_rays"] == 3352
    assert all(
        ss2["finite_volume"]["exact_GMP_cone"]["fresh_replay"]["checks"].values()
    )
    assert f3["interior_provenance"]["published_ample_recovered"] is True
    assert f3["oriented_wall_system"]["strict_interior_square"] == 4


def test_scope_is_machine_residue_only(certificate: dict) -> None:
    assert certificate["all_checks_pass"] is True
    assert all(certificate["checks"].values())
    scope = certificate["scope"]
    assert "presently audited repository queue" in scope["machine_ingested_residue"]
    assert scope["not_all_ledger_rows"] is True
    assert scope["not_world_literature_exhaustion"] is True
    assert scope["not_new_abstract_polytopes"] is True
    assert "not asserted" in scope["not_surface_automorphisms"]
    assert scope["worldwide_priority_claim"] == "none"


def test_report_freezes_uniform_h10_theorem() -> None:
    report = (
        Path(__file__).resolve().parents[1]
        / "results/three_h10_uniform_full_picard_k3_carriers.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "positive theorem" in normalized
    assert "d16-d5-h10-30" in normalized
    assert "ss2-2d-h10-47" in normalized
    assert "f3-c2-h10-196" in normalized
    assert "j_4-i_4" in normalized
    assert "nine-dimensional" in normalized
    assert "machine-ingested" in normalized
    assert "not full-nef" in normalized

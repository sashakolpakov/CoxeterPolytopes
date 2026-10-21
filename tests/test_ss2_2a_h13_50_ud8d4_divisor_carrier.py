from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_ss2_2a_h13_50_ud8d4_divisor_carrier import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_scaled_dual_is_actual_complex_k3_picard_lattice(certificate: dict) -> None:
    bridge = certificate["scaled_dual_bridge"]
    source = bridge["source_M"]
    dual = bridge["doubled_dual_N"]
    target = bridge["actual_K3_target"]
    isometry = bridge["explicit_N_to_K3_isometry"]
    assert source["determinant"] == -1024
    assert source["smith"] == [1] * 4 + [2] * 10
    assert dual["definition"] == "N=M^vee(2), Gram(N)=2*Gram(M)^(-1)"
    assert dual["isometry_class"] == target["isometry_class"]
    assert dual["isometry_class"] == "U+D8(-1)+D4(-1)"
    assert dual["determinant"] == -16
    assert dual["smith"] == [1] * 10 + [2] * 4
    assert len(target["genuine_BH_K3_rows"]) == 2
    assert isometry["determinant"] == -1
    assert isometry["identity"] == "X^T Gram(K3) X = 2 Gram(M)^(-1)"
    assert all(bridge["checks"].values())


def test_all_50_mirrors_generate_ns_and_reflect_integrally(certificate: dict) -> None:
    mirrors = certificate["primitive_mirrors"]
    assert mirrors["count"] == 50
    assert mirrors["source_square_divisibility_histogram"] == {
        "-4,2": 40,
        "-2,1": 10,
    }
    assert mirrors["doubled_dual_coordinate_content_histogram"] == {
        "2": 40,
        "1": 10,
    }
    assert mirrors["K3_square_divisibility_histogram"] == {
        "-2,1": 40,
        "-4,2": 10,
    }
    assert mirrors["wall_lattice_index_in_K3_NS"] == 1
    assert mirrors["wall_lattice_hnf"] == [
        [1 if row == column else 0 for column in range(14)] for row in range(14)
    ]
    assert mirrors["all_50_reflections_integral"] is True
    assert all(mirrors["checks"].values())
    assert all(all(row["checks"].values()) for row in mirrors["records"])


def test_strict_interior_and_independent_exact_cone(certificate: dict) -> None:
    interior = certificate["strict_carrier_interior"]
    cone = certificate["independent_exact_cone_cross_check"]
    assert interior["square"] == 2012
    assert interior["minimum_wall_pairing"] == 6
    assert all(interior["checks"].values())
    assert cone["canonical_input_sha256"] == (
        "96e2abe698f9f80a41f1cc983c35fb5b094f31a6de5b7ec8c76351e2de5dab3c"
    )
    assert cone["ray_count"] == 48552
    assert (cone["timelike_rays"], cone["isotropic_rays"], cone["hyperideal_rays"]) == (
        45954,
        2598,
        0,
    )
    assert cone["ideal_ray_active_facet_histogram"] == {
        "14": 1656,
        "15": 872,
        "16": 68,
        "17": 2,
    }
    assert all(cone["checks"].values())


def test_literal_vinberg_gate_and_scope(certificate: dict) -> None:
    vinberg = certificate["vinberg_finite_volume"]
    scope = certificate["scope"]
    assert vinberg["classification"] == "finite_volume_noncompact_Coxeter_polytope"
    assert (vinberg["dimension"], vinberg["facets"]) == (13, 50)
    assert (vinberg["finite_vertices"], vinberg["cusps"]) == (45954, 2598)
    assert vinberg["elliptic_edges"] == 391614
    assert vinberg["connected_affine_components"] == 6210
    assert vinberg["bad_edges"] == vinberg["bad_parabolic_neighborhoods"] == 0
    assert all(vinberg["checks"].values())
    assert "not a same-form" in scope["scaled_dual_not_overlattice"]
    assert "square-minus-four" in scope["not_full_nef"]
    assert "no Brandhorst" in scope["not_automorphism_fixed"]
    assert "no abstract-polytope priority" in scope["novelty"]
    assert certificate["hostile_audit"]["status"] == "green"
    assert certificate["all_checks_pass"] is True


def test_report_freezes_gain_and_boundaries() -> None:
    report = (
        Path(__file__).resolve().parents[1]
        / "results"
        / "ss2_2a_h13_50_ud8d4_divisor_carrier.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "complex-k3 carrier" in normalized
    assert "scaled dual" in normalized
    assert "40 primitive `(-2)`" in normalized
    assert "10 primitive `(-4)`" in normalized
    assert "48,552" in normalized
    assert "not a same-form overlattice" in normalized
    assert "not the full nef cone" in normalized
    assert "not a new abstract polytope" in normalized


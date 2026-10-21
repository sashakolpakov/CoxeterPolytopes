"""Regression tests for exact-trace cusp ground-field collapse."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

import pytest

from scripts.certify_exact_trace_noncompact_ground_field_collapse import (
    bilinear,
    build_certificate,
    cube_cusp_record,
    digest,
    encoded_certificate,
    identity,
    matmul,
    matrix_rank,
    reflection_matrix,
    transpose,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/exact_trace_noncompact_ground_field_collapse.md"
CERTIFIER = ROOT / "scripts/certify_exact_trace_noncompact_ground_field_collapse.py"
ARCHIVE = ROOT / "data/enumerations/exact_trace_noncompact_ground_field_collapse.json"


def standard_frame(dimension: int):
    size = dimension + 1
    form = identity(size)
    form[-1][-1] = Fraction(-1)
    cusp = [Fraction(0) for _ in range(size)]
    cusp[-2] = cusp[-1] = Fraction(1)
    normals = []
    for coordinate in range(dimension - 1):
        for sign in (1, -1):
            normal = cusp[:]
            normal[coordinate] = Fraction(sign)
            normals.append(normal)
    return form, cusp, normals


def test_cube_cusp_frames_cover_simple_and_nonsimple_pressure_points():
    for dimension in range(2, 21):
        record = cube_cusp_record(dimension)
        assert record["ambient_rank"] == dimension + 1
        assert record["incident_facet_count"] == 2 * (dimension - 1)
        assert record["incident_normal_rank"] == dimension
        assert record["common_fixed_space_dimension"] == 1
        assert record["cusp_square"] == 0
        assert record["nonsimple_when_incident_count_exceeds_dimension"] == (
            dimension > 2
        )


def test_incident_normals_span_u_perp_and_have_positive_dependence():
    for dimension in range(2, 11):
        form, cusp, normals = standard_frame(dimension)
        assert bilinear(cusp, form, cusp) == 0
        assert all(bilinear(cusp, form, normal) == 0 for normal in normals)
        assert all(bilinear(normal, form, normal) == 1 for normal in normals)
        assert matrix_rank(normals) == dimension
        summed = [sum(normal[j] for normal in normals) for j in range(dimension + 1)]
        assert summed == [2 * (dimension - 1) * entry for entry in cusp]


def test_each_model_wall_is_an_exact_integral_reflection():
    for dimension in range(2, 9):
        form, _, normals = standard_frame(dimension)
        unit = identity(dimension + 1)
        for normal in normals:
            reflection = reflection_matrix(form, normal)
            assert matmul(reflection, reflection) == unit
            assert matmul(transpose(reflection), matmul(form, reflection)) == form
            assert all(entry.denominator == 1 for row in reflection for entry in row)


def test_h3_square_cusp_is_explicitly_nonsimple():
    record = cube_cusp_record(3)
    assert record["horosphere_section"] == "[-1,1]^2"
    assert record["incident_facet_count"] == 4
    assert record["simple_cusp_facet_count"] == 3
    assert record["nonsimple_when_incident_count_exceeds_dimension"] is True
    assert record["positive_dependence"] == "sum(e_i)=4*u"


def test_certificate_archive_replays_exactly():
    envelope = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = build_certificate()
    assert envelope["payload"] == payload
    assert envelope["payload_sha256"] == digest(payload)
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()


def test_conclusions_and_nonclaims_are_literal():
    payload = build_certificate()
    conclusions = payload["conclusions"]
    assert conclusions["returned_coefficient_field_is_Q"] is True
    assert conclusions["full_trace_K3_noncompact_field_set_is_singleton_Q"] is True
    assert conclusions["full_trace_K3_rank_identity_after_collapse"] == "n+1=rho(X)<=20"
    assert conclusions["full_trace_K3_hyperbolic_dimension_at_most_19"] is True
    assert conclusions["does_not_supersede_rational_weak_log_n_at_most_25"] is True
    assert conclusions["requires_reflection_maximality"] is False
    scope = payload["scope"]
    assert scope["finite_projective_order_module_allowed"] is True
    assert scope["free_module_required"] is False
    assert scope["primitive_wall_normals_required"] is False
    assert scope["simple_cusp_required"] is False
    assert scope["exact_full_incident_cusp_packet_required"] is True
    assert scope["upper_half_space_argument_includes_nonsimple_cusps"] is True
    assert scope["partial_affine_subdiagram_suffices"] is False
    assert scope["compact_chambers_closed"] is False
    assert scope["Property_4_proved"] is False
    assert scope["raw_chamber_dimension_bound_proved"] is False


def test_report_contains_every_field_and_cusp_pressure_point():
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    for phrase in (
        "proved in every hyperbolic dimension",
        "finite projective",
        "positive definite for every",
        "This includes nonsimple cusps",
        r"C\times(T,\infty)",
        "An unbounded full-dimensional Euclidean polyhedron has infinite",
        r"\sum_i\lambda_i a_i\ne0",
        r"\bigcap_{i\in I}\ker(s_i-1)",
        "Extension of scalars",
        r"h(u,u)=0\quad\hbox{in }k",
        r"u^\sigma\ne0",
        r"\boxed{k=\mathbf Q}",
        r"n+1=\rho(X)\le20",
        "does not supersede the distinct rational weak-log branch",
        "partial wall subsystem",
    ):
        assert phrase in text


def test_transport_and_no_classification_controls_are_clean():
    combined = REPORT.read_bytes() + CERTIFIER.read_bytes()
    assert b"\x08" not in combined
    text = combined.decode("utf-8")
    assert ",quad" not in text
    assert ",qquad" not in text
    for forbidden in ("Nikulin", "Vinberg", "ABSW"):
        assert forbidden not in text


def test_dimension_one_is_rejected_by_the_model_certifier():
    with pytest.raises(ValueError):
        cube_cusp_record(1)


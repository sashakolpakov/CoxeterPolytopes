"""Independent hostile tests for compact narrow-facet K3 field finiteness."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

from scripts.audit_k3_compact_narrow_facet_field_finiteness_independent_hostile import (
    build_certificate,
    coboundary_similarity_record,
    constant_record,
    primitive_grid_rows,
    proper_extension_rows,
    rank_degree_rows,
    star_span_rows,
    trace_cycle_record,
)


ROOT = Path(__file__).resolve().parents[1]
OWNER_REPORT = ROOT / "results/k3_compact_narrow_facet_field_finiteness.md"
OWNER_CERT = ROOT / "scripts/certify_k3_compact_narrow_facet_field_finiteness.py"
OWNER_ARCHIVE = (
    ROOT / "data/enumerations/k3_compact_narrow_facet_field_finiteness.json"
)
OWNER_TEST = ROOT / "tests/test_k3_compact_narrow_facet_field_finiteness.py"
INDEPENDENT_REPORT = (
    ROOT
    / "results/k3_compact_narrow_facet_field_finiteness_independent_hostile_audit.md"
)
INDEPENDENT_SCRIPT = (
    ROOT
    / "scripts/audit_k3_compact_narrow_facet_field_finiteness_independent_hostile.py"
)

OWNER_PINS = {
    OWNER_REPORT: "39c9a3a8422110bcdcf751ed7bae90ac86cdf3a720d220ee164ff241bd9ac9b5",
    OWNER_CERT: "d754d69fad85366ad6a2774be2034c24ebf50a9620ec79618df4b18a7f42d9e7",
    OWNER_ARCHIVE: "3d0a2dff63bae86421a286131026f75fc4dbfaae9dec72be10cc4f38f5cd7e09",
    OWNER_TEST: "10ddc82f12c88d7e718f62ee859812ac22dbd6f37d17450e2f5c59280dce8a45",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_frozen_owner_bytes_are_exact() -> None:
    assert {path: sha256(path) for path in OWNER_PINS} == OWNER_PINS


def test_star_dimension_count_closes_every_rank_three_through_twenty() -> None:
    rows = star_span_rows()
    assert [row["module_rank"] for row in rows] == list(range(3, 21))
    for row in rows:
        assert row["incident_tangent_rank"] + 1 == row["projected_star_rank"]
        assert row["projected_star_rank"] + 1 == row["normal_star_rank"]
        assert row["normal_star_rank"] == row["module_rank"]
        assert row["nonincident_side_adds_one"] is True


def test_rank_one_products_have_exact_cyclic_traces_with_repetitions() -> None:
    rows = trace_cycle_record()
    assert len(rows) == 4
    assert all(row["equal"] for row in rows)
    assert any(len(set(row["word"])) < len(row["word"]) for row in rows)
    assert rows[0]["trace"] == "3"


def test_coboundary_similarity_preserves_lorentz_characteristic_data() -> None:
    record = coboundary_similarity_record()
    assert record["connected_support"] is True
    assert record["two_cycles_preserved"] is True
    assert record["characteristic_polynomial_preserved"] is True
    assert record["leading_principal_minors"] == ["2", "-5", "-12"]
    assert record["ldl_pivots"] == ["2", "-5/2", "12/5"]
    assert record["negative_pivots"] == 1


def test_proper_subfield_has_at_least_two_real_extensions() -> None:
    rows = proper_extension_rows()
    assert rows
    assert all(row["real_extensions"] >= 2 for row in rows)
    assert {row["total_degree"] for row in rows} == {2, 3, 4, 5, 6}


def test_k3_trace_rank_forces_d_six_m_ten_and_two_generators() -> None:
    rows = rank_degree_rows()
    assert max(row["field_degree"] for row in rows) == 6
    assert max(row["module_rank"] for row in rows) == 10
    assert max(row["strict_adjunction_count_bound"] for row in rows) == 2
    assert max(row["discriminant_exponent"] for row in rows) == 30
    assert all(row["trace_rank"] <= 20 for row in rows)
    assert all(row["hyperbolic_dimension"] >= 2 for row in rows)


def test_sixteen_point_primitive_grid_beats_all_embedding_pairs() -> None:
    rows = primitive_grid_rows()
    assert rows[-1] == {
        "degree": 6,
        "embedding_pairs": 15,
        "candidate_count_needed": 16,
        "uniform_grid_size_16_suffices": True,
    }
    assert all(row["uniform_grid_size_16_suffices"] for row in rows)


def test_uniform_height_and_discriminant_constants_are_exact() -> None:
    record = constant_record()
    assert record["cycle_height"] == 14**10 == 289_254_654_976
    assert record["discriminant_base"] == 32 * 14**10 == 9_256_148_959_232
    assert record["discriminant_exponent"] == 30
    assert record["discriminant_expression"] == "(32*14^10)^30"
    assert record["discriminant_decimal_digits"] == len(
        str((32 * 14**10) ** 30)
    )


def test_independent_script_replays_without_owner_imports() -> None:
    completed = subprocess.run(
        [sys.executable, str(INDEPENDENT_SCRIPT)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(completed.stdout)
    assert payload["verdict"] == "PASS_EXACT_COMPACT_FULL_TRACE_K3_N_GE_2"
    assert payload["field_descent"]["coefficient_field_equals_cycle_field"]
    assert not payload["field_descent"]["minimal_field_hypothesis_needed"]
    assert payload["constants"]["discriminant_expression"] == "(32*14^10)^30"
    source = INDEPENDENT_SCRIPT.read_text(encoding="utf-8")
    assert "certify_k3_compact_narrow_facet_field_finiteness" not in source


def test_certificate_scope_keeps_the_n_one_poison_and_nonclaims() -> None:
    payload = build_certificate()
    scope = payload["scope"]
    assert scope["hyperbolic_dimension_at_least_two_required"]
    assert scope["dimension_one_infinite_field_control_retained"]
    assert scope["admissible_totally_real_k_form_required"]
    assert scope["k_linear_integral_reflections_required"]
    assert not scope["reflection_maximality_used"]
    assert not scope["saturation_used"]
    assert not scope["carrier_construction_proved"]
    assert not scope["property4_proved"]
    assert not scope["finite_hosts_or_chambers_proved"]


def test_report_states_every_pressure_point_and_scope_boundary() -> None:
    text = " ".join(INDEPENDENT_REPORT.read_text(encoding="utf-8").split())
    for phrase in (
        "The narrow-facet star contains an actual",
        "some side of",
        "The intrinsic facet",
        "tree block",
        "trivial holonomy",
        "cannot be a gratuitous proper scalar extension",
        "finite projective",
        r"|\tau(c)|<14^m",
        "sixteen integers",
        "genuine discriminant bound",
        r"gate \(n\ge2\) is essential",
        "does not construct that carrier",
    ):
        assert phrase in text


def test_tex_transport_and_dependency_controls_are_clean() -> None:
    combined = INDEPENDENT_REPORT.read_bytes() + INDEPENDENT_SCRIPT.read_bytes()
    assert not any(byte < 32 and byte not in (9, 10, 13) for byte in combined)
    text = combined.decode("utf-8")
    assert ",quad" not in text
    assert ",qquad" not in text
    assert "|le" not in text
    assert "TODO" not in text
    assert "FIXME" not in text

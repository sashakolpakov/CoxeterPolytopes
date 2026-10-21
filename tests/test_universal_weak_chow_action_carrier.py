import hashlib
import json
from pathlib import Path

from scripts.certify_universal_weak_chow_action_carrier import build_certificate


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/universal_weak_chow_action_carrier.json"
REPORT = ROOT / "results/universal_weak_chow_action_carrier_and_leverage_no_go.md"


def test_certificate_matches_frozen_archive() -> None:
    live = build_certificate()
    frozen = json.loads(ARCHIVE.read_text())
    assert live == frozen
    assert live["all_checks_pass"] is True


def test_triangle_chow_action_replay() -> None:
    result = build_certificate()
    triangle = result["triangle_calibration"]
    checks = triangle["checks"]
    assert triangle["name"] == "arithmetic (2,3,infinity) triangle"
    assert triangle["surface"] == "Bl_9(P^2)"
    assert triangle["projector_rank"] == 3
    assert triangle["all_checks_pass"] is True
    assert all(checks.values())
    assert checks["projector_cycle_idempotent"] is True
    assert checks["lifted_cycles_self_transpose"] is True
    assert checks["all_pair_cycle_composition_identities"] is True
    assert checks["lifted_m2_relation"] is True
    assert checks["lifted_m3_relation"] is True
    assert checks["lifted_parallel_product_nontrivial_on_summand"] is True


def test_uniform_published_bound_is_imported_not_derived() -> None:
    bound = build_certificate()["unconditional_published_bound_corollary"]
    assert bound["published_degree_bound"] == 25
    assert bound["published_n_plus_one_bound"] == 30
    assert bound["trace_rank_bound"] == 750
    assert bound["uniform_blowup_bound"] == 2997
    assert bound["uniform_picard_rank_bound"] == 2998
    assert bound["classical_bounds_imported_not_derived"] is True
    assert bound["weak_carrier_alone_bounds_dimension_or_field"] is False


def test_dependency_pin_and_scope() -> None:
    result = build_certificate()
    for relative, expected in result["dependency_sha256"].items():
        actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        assert actual == expected
    scope = result["scope"]
    assert scope["finite_volume_needed_for_carrier"] is False
    assert scope["reflection_maximality_needed_for_carrier"] is False
    assert scope["effectivity_or_nef_return_claimed"] is False
    assert scope["regular_or_birational_surface_action_claimed"] is False


def test_report_keeps_the_outer_and_leverage_quantifiers_separate() -> None:
    text = REPORT.read_text()
    required = (
        "unconditional in the marked rational-Chow category",
        "W(P)\\longrightarrow",
        "b\\le2997",
        "imports",
        "not a new proof",
        "Finite volume and reflection maximality are not needed",
        "Existence of a carrier satisfying Theorem 1",
        "implies neither a dimension bound nor a finite",
        "No universal effective/nef or mixed positive-characteristic theorem is",
    )
    for phrase in required:
        assert phrase in text

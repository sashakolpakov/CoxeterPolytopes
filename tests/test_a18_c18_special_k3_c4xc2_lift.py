from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_a18_c18_special_k3_c4xc2_lift import build_certificate


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT / "data" / "enumerations" / "a18_c18_special_k3_c4xc2_lift.json"
)
REPORT = ROOT / "results" / "a18_c18_special_k3_c4xc2_lift.md"
LEDGER = ROOT / "results" / "higher_dimensional_polytope_lineage_ledger.md"


def _archived() -> dict:
    return json.loads(CERTIFICATE.read_text(encoding="utf-8"))


def test_simultaneous_c4xc2_glue_replays_exactly() -> None:
    archived = _archived()
    assert build_certificate() == archived
    assert archived["all_checks_pass"]
    assert archived["certificate_payload_sha256"] == (
        "8b36468302ee6bc3b9cd3bd199578201cea8e26091ef0416cdd8937fe7da7cc3"
    )
    group = archived["k3_group_action"]
    assert group["abstract_group"] == "C4 x C2"
    assert group["order"] == 8
    assert group["restriction_to_N"] == "V4"
    assert "common square" in group["central_kernel_on_N"]


def test_special_period_and_common_fixed_picard_lattice_are_exact() -> None:
    archived = _archived()
    period = archived["period_and_picard_lattice"]
    assert period["common_period_in_T_tensor_C"] == ["-i", "-i", "1", "1"]
    assert period["three_eigenvalues"] == ["1", "i", "-i"]
    assert period["period_square"] == "0"
    assert period["period_positive_norm"] == "8"
    assert period["picard_rank"] == 20
    assert period["picard_signature"] == [1, 19]
    assert period["picard_determinant"] == -16
    assert period["picard_smith_diagonal_absolute"] == [*[1] * 18, 4, 4]
    assert period["common_fixed_NS_is_exactly_archived_C18_fixed_lattice"]
    assert period["common_fixed_NS_rank"] == 11
    assert period["common_fixed_NS_determinant"] == 128


def test_rootless_coinvariant_closes_the_group_kahler_gate() -> None:
    archived = _archived()
    audit = archived["rootless_coinvariant_and_torelli"]
    assert audit["common_algebraic_coinvariant_rank"] == 9
    assert audit["common_algebraic_coinvariant_signature"] == [0, 9]
    assert audit["common_algebraic_coinvariant_determinant"] == -512
    assert audit["common_algebraic_coinvariant_smith_diagonal_absolute"] == [
        1,
        1,
        2,
        2,
        2,
        2,
        2,
        2,
        8,
    ]
    roots = audit["minus_two_vector_audit_on_negative_coinvariant"]
    assert roots["target"] == 2
    assert roots["vector_count"] == 0
    assert roots["vectors"] == []
    assert audit["rootless"]
    assert audit["source_48_wall_chamber_is_preserved"]
    assert "genuine C4 x C2 automorphism group" in audit[
        "global_torelli_conclusion"
    ]


def test_h10_21_vinberg_handoff_and_scope_are_frozen() -> None:
    archived = _archived()
    chamber = archived["hyperbolic_fixed_section"]
    assert (chamber["dimension"], chamber["facet_count"]) == (10, 21)
    assert chamber["finite_volume"] and not chamber["compact"]
    assert (chamber["finite_vertex_count"], chamber["ideal_vertex_count"]) == (
        329,
        22,
    )
    assert chamber["elliptic_edge_count"] == 2088
    assert chamber["bad_edge_count"] == 0
    assert chamber["bad_parabolic_neighborhood_count"] == 0
    assert "not a new abstract polytope" in archived["novelty_scope"]
    assert "not asserted to be the literal nef cone" in archived["nef_scope_caution"]


def test_report_and_ledger_state_the_actual_group_and_special_scope() -> None:
    report = REPORT.read_text(encoding="utf-8")
    ledger = LEDGER.read_text(encoding="utf-8")
    assert "C4 x C2" in report
    assert "restriction to (N) is the archived (V_4)" in report
    assert "Picard rank 20" in report
    assert "rank-9 algebraic coinvariant" in report
    assert "zero vectors of square (-2)" in report
    assert "not the literal nef cone" in report
    assert "A18-V4-H10-21" in ledger
    assert "special Picard-rank-20 K3" in ledger
    assert "C4 x C2" in ledger


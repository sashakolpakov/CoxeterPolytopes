import hashlib
import json
from pathlib import Path

from scripts.certify_bounded_universal_mixed_ag_carrier_no_go_and_next_species import (
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/bounded_universal_mixed_ag_carrier_no_go_and_next_species.json"
REPORT = ROOT / "results/bounded_universal_mixed_ag_carrier_no_go_and_next_species.md"


def test_live_certificate_matches_archive() -> None:
    live = build_certificate()
    assert live == json.loads(ARCHIVE.read_text())
    assert live["all_replay_checks_pass"] is True


def test_dependency_pins() -> None:
    for relative, expected in build_certificate()["dependency_sha256"].items():
        assert hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == expected


def test_both_controls_pass_only_the_mixed_species() -> None:
    controls = build_certificate()["controls"]
    assert controls["X3"]["mixed_k3_fm_pass"] is True
    assert controls["X3"]["full_nef_k3_pass"] is False
    assert controls["X3"]["absolute_maximality_claimed"] is False
    assert controls["H21"]["mixed_k3_fm_pass"] is True
    assert controls["H21"]["og10_full_cone_pass"] is False
    assert controls["H21"]["absolute_maximality_proved_in_dependency"] is True


def test_no_global_completion_is_claimed() -> None:
    status = build_certificate()["terminal_theorem_status"]
    assert status["full_trace_mukai_rank_ceiling"] == 24
    assert status["full_trace_necessary_inequality"] == "[k:Q](n+1) <= 24"
    assert status["mixed_k3_fm_species_currently_total"] is False
    assert status["unconditional_global_AG_finiteness_proved"] is False
    assert status["classical_reflection_rank_bound_imported"] is False
    assert not any(build_certificate()["unproved_universal_gates"].values())


def test_report_contains_rigorous_fm_lift_and_missing_theorem() -> None:
    text = REPORT.read_text()
    for phrase in (
        "Proposition 3 (Fourier--Mukai lift of every long `X3` reflection)",
        "r_v(v/3)=-v/3",
        "orientation-preserving Mukai Hodge isometries",
        "mixed K3 curve/Fourier--Mukai-action",
        "Uniform Mukai essential-surjectivity conjecture",
        "d(n+1)=\\operatorname {rk}_{\\mathbf Z}T_P\\le24",
        "the literal user acceptance criterion is not met",
    ):
        assert phrase in text

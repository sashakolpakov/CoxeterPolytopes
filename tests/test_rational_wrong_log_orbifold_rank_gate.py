import json
from fractions import Fraction
from pathlib import Path

from scripts.certify_rational_wrong_log_orbifold_rank_gate import (
    ARCHIVE,
    admissible_rank_pairs,
    build_certificate,
    canonical_degree,
    certify,
    d_square,
    exact_picard_ceiling,
    natural_wrong_log_pairing,
    orbifold_c2_times_four,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_wrong_log_orbifold_rank_gate.md"


def test_adjunction_and_natural_pairing_table():
    assert canonical_degree(-4) == 2
    assert canonical_degree(-2) == 0
    assert canonical_degree(-1) == -1
    assert natural_wrong_log_pairing(-4, -4, 0) == 0
    assert natural_wrong_log_pairing(-4, -4, 2) == -1
    assert natural_wrong_log_pairing(-2, 0) == 0
    assert natural_wrong_log_pairing(-2, 2) == -1
    assert natural_wrong_log_pairing(-1, 0) == 1
    assert natural_wrong_log_pairing(-1, 1) == Fraction(1, 2)
    assert natural_wrong_log_pairing(-1, 2) == 0
    assert natural_wrong_log_pairing(-1, 3) == Fraction(-1, 2)


def test_symbolic_rank_elimination_and_bounded_poison_checks():
    assert exact_picard_ceiling() == 26
    pairs = admissible_rank_pairs(limit=26)
    assert [(rho, b) for rho, b in pairs if rho == 26] == [(26, 16)]
    assert d_square(26, 16) == 0
    assert orbifold_c2_times_four(26, 16) == 0
    for limit in (27, 100, 500):
        poison_pairs = admissible_rank_pairs(limit=limit)
        assert max(rho for rho, _ in poison_pairs) == 26


def test_certificate_controls_and_archive():
    payload = build_certificate()
    certify(payload)
    assert payload["all_checks_pass"]
    assert payload["F4_control"]["passes"]
    assert payload["F4_control"]["K_plus_half_B_pairings"] == {
        "B": "0",
        "E": "-1",
        "G": "-1/2",
    }
    assert payload["F4_control"]["D_square"] == "8"
    assert payload["F4_control"]["orbifold_c2"] == "13/4"
    assert all(row["passes"] for row in payload["geometric_hostile_controls"])
    distributed = payload["geometric_hostile_controls"][-1]
    assert distributed["component_boundary_degrees"] == [1, 1, 1]
    assert distributed["aggregate_boundary_degree"] == 3
    assert distributed["forced_exceptional_coefficient_lower_bound"] == "2"
    assert json.loads(ARCHIVE.read_text(encoding="utf-8")) == payload


def test_report_scope_and_clean_bytes():
    raw = REPORT.read_bytes()
    assert b"\r" not in raw
    assert not any(byte < 32 and byte not in (9, 10) for byte in raw)
    text = raw.decode("utf-8")
    assert "rho(S)\\leq26" in text
    assert "\\frac14(1,1)" in text
    assert "rho(Y)=\\rho(S)-b" in text
    assert "klt surface base-point-free theorem" in text
    assert "one-string Hirzebruch--Jung contraction theorem" in text
    assert "e(Y_{\\rm reg})" in text
    assert "thm-nef-anti-canonical" in text
    assert "Nothing in the argument below" in text
    assert "not essential surjectivity" in text
    assert "Four geometric hostile controls" in text
    assert "aggregate degree is \\(B.E=1+1+1=3\\)" in text
    assert "The maximal triangle is therefore a positive calibration" in text

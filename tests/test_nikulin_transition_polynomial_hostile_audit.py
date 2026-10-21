from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/independent_nikulin_transition_polynomial_audit.py"
ARCHIVE = ROOT / "data/enumerations/nikulin_transition_polynomial_audit.json"
REPORT = ROOT / "results/nikulin_transition_polynomial_hostile_audit.md"


def load_module():
    spec = importlib.util.spec_from_file_location(
        "nikulin_transition_polynomial_hostile_audit", SCRIPT
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_frozen_archive_is_exact_replay() -> None:
    module = load_module()
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert archived == module.build_certificate()
    assert archived["all_checks_pass"] is True


def test_arxiv_failure_and_exact_true_maximum() -> None:
    cert = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    printed = cert["printed_certificate"]
    assert printed["degree"] == 18
    assert printed["rational_counterexample"] == {
        "x": "3/4",
        "absolute_value": "6471355275/68719476736",
        "strictly_exceeds_claimed_supremum": True,
    }
    assert printed["true_supremum"]["minimal_polynomial_coefficients"] == [
        387420489,
        -84059036,
        3538944,
    ]
    bound = printed["lemma_7_at_t_14"]
    assert bound["integer_degree_bound_from_this_certificate"] == 27
    assert bound["printed_claim_is_established_by_printed_certificate"] is False


def test_journal_repair_is_exact() -> None:
    cert = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    repair = cert["minimal_repair"]
    assert repair["degree"] == 20
    assert repair["exact_supremum"] == "884736/9765625"
    assert repair["identity"] == "P_repaired(x)=P_equation_13(x+2)"
    assert repair["lemma_7_at_t_14"]["integer_degree_bound"] == 23
    assert cert["source"]["journal_version"]["corrected_exponent"] == 3
    assert cert["propagation_scope"]["global_N14_theorem_disproved"] is False


def test_report_has_version_and_scope_guards() -> None:
    raw = REPORT.read_bytes()
    assert b"\r" not in raw and b"\x08" not in raw and b"\x0b" not in raw
    text = raw.decode("utf-8")
    assert "arXiv-v4 certificate is false as printed" in text
    assert "journal silently repairs it" in text
    assert "does **not** create a gap in the final journal theorem" in text
    assert "hostile regression corpus" in text
    assert "no Nikulin/Vinberg classification result is allowed" in text


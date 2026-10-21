"""Lean independent pin and scope checks for the compact K3 compiler audit."""

from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_compact_full_picard_raw_compiler_criterion_independent_hostile_audit.md"

PINS = {
    "results/k3_compact_full_picard_raw_compiler_criterion.md":
        "3d8c1fc26d6d8242bf77095640a54d3690672ae76367b0508d89a10ca59a745f",
    "scripts/certify_k3_compact_full_picard_raw_compiler_criterion.py":
        "2598775af699d6bfbbe395e4923a37f924706cdf0224e990fb5b1ecdbf095816",
    "data/enumerations/k3_compact_full_picard_raw_compiler_criterion.json":
        "f272928df529a923d01ce5e46b2e9fe95f215ed26fb0ddd33524a199a896b173",
    "tests/test_k3_compact_full_picard_raw_compiler_criterion.py":
        "6ebc9d56fea4fa5515c9386bd48639450424bc372691f71821829d643167231a",
    "results/exact_trace_noncompact_ground_field_collapse.md":
        "d738018bacc3803f78fa725da375bb8045ecce92a1c879b03b30a9998e980255",
    "results/noncompact_k3_bounded_spanning_root_lemma.md":
        "de3b78a9070dc7035a8ecf51c5aa697b0736ad1559f0db8867b6629eb3cc4446",
    "results/noncompact_k3_bounded_spanning_root_lemma_independent_hostile_audit.md":
        "70ae2c1e09ddd0fe8cd1999087172ed14499dbd47709ebbc2236aa77ad985bde",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_all_audited_inputs_are_exactly_pinned() -> None:
    assert {relative: sha256(ROOT / relative) for relative in PINS} == PINS


def test_owner_contract_contains_the_repaired_ample_eigenframe() -> None:
    owner = (ROOT / "results/k3_compact_full_picard_raw_compiler_criterion.md").read_text(
        encoding="utf-8"
    )
    for needle in (
        "p\\in E\\cap\\operatorname {Amp}(S)",
        "[p]\\in P^\\circ",
        "\\overline{\\mathbf P\\bigl(\\operatorname {Amp}(S)\\cap E\\bigr)}",
        "The compiler can return `EMPTY`",
        "not Property 4 by itself",
    ):
        assert needle in owner


def test_independent_audit_states_the_exact_scope_and_pressure_points() -> None:
    audit = REPORT.read_text(encoding="utf-8")
    for needle in (
        "PASS in the exact pinned-realization scope",
        "common-square-class conditions",
        "F_0\\subseteq M\\subseteq F_0^\\#",
        "q_T\\simeq-q_N",
        "Q_p(\\delta)\\le2\\cosh(2R)",
        "codimension greater than one",
        "0\\le E_i.E_j\\le13",
        "The output list can be empty",
        "not universal Property 4",
    ):
        assert needle in audit

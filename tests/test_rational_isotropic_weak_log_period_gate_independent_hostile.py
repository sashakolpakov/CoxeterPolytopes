"""Independent hostile tests for the rational isotropic period gate.

The lattice and cohomology recurrence below are reconstructed without using
the producer's implementations.  Importing the producer is reserved for the
separate byte-exact archive replay.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

from scripts.certify_rational_isotropic_weak_log_period_gate import (
    encoded_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_isotropic_weak_log_period_gate.md"
CERTIFIER = ROOT / "scripts/certify_rational_isotropic_weak_log_period_gate.py"
ARCHIVE = ROOT / "data/enumerations/rational_isotropic_weak_log_period_gate.json"
ORIGINAL_TEST = ROOT / "tests/test_rational_isotropic_weak_log_period_gate.py"
AUDIT = ROOT / "results/rational_isotropic_weak_log_period_gate_independent_hostile_audit.md"
THIS_TEST = Path(__file__).resolve()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def intersection(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    assert len(left) == len(right) == 10
    return left[0] * right[0] - sum(
        a * b for a, b in zip(left[1:], right[1:])
    )


def independent_profile(order: int | None, horizon: int) -> list[tuple[int, int, int, int]]:
    """Return (r,target h0,restriction rank,h0(rF)) from the proof logic."""

    assert horizon >= 0
    assert order is None or order >= 1
    h0 = 1
    rows: list[tuple[int, int, int, int]] = [(0, 0, 0, h0)]
    for r in range(1, horizon + 1):
        target_h0 = int(order is not None and r % order == 0)
        # At a multiple, a power of the first lift restricts nontrivially;
        # otherwise the degree-zero target has no section.
        restriction_rank = target_h0
        h0 += restriction_rank
        rows.append((r, target_h0, restriction_rank, h0))
    return rows


def test_independent_I_1_9_reconstruction_and_rr_identity():
    canonical = (-3,) + (1,) * 9
    fibre = (3,) + (-1,) * 9
    assert fibre == tuple(-value for value in canonical)
    assert intersection(fibre, fibre) == 0
    assert intersection(canonical, fibre) == 0
    assert intersection(canonical, canonical) == 0
    assert __import__("math").gcd(*map(abs, fibre)) == 1
    for index in range(1, 10):
        exceptional = tuple(int(position == index) for position in range(10))
        assert intersection(exceptional, exceptional) == -1
        assert intersection(fibre, exceptional) == 1
    for r in range(0, 4097):
        rr_numerator = r * r * intersection(fibre, fibre) - r * intersection(
            fibre, canonical
        )
        assert rr_numerator == 0
        assert 1 + rr_numerator // 2 == 1


def test_independent_torsion_and_nontorsion_recursions_include_edge_cases():
    for order in range(1, 65):
        rows = independent_profile(order, 512)
        for r, target_h0, rank, h0 in rows:
            if r == 0:
                assert (target_h0, rank, h0) == (0, 0, 1)
                continue
            assert target_h0 == int(r % order == 0)
            assert rank == target_h0
            assert h0 == r // order + 1
            # chi=1 and h2=0 independently force h1=h0-1.
            assert h0 - (r // order) == 1
        assert rows[order][3] == 2
        assert all(row[3] == 1 for row in rows[1:order])

    non_torsion = independent_profile(None, 2048)
    assert all(row[3] == 1 for row in non_torsion)
    assert all(row[1] == row[2] == 0 for row in non_torsion)


def test_archive_is_byte_exact_and_labels_finite_replay_as_nonproof():
    expected = encoded_certificate()
    assert ARCHIVE.read_text(encoding="utf-8") == expected
    envelope = json.loads(expected)
    payload = envelope["payload"]
    canonical = json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode()
    assert envelope["payload_sha256"] == hashlib.sha256(canonical).hexdigest()
    assert payload["scope"]["computer_replay_replaces_geometric_proof"] is False


def test_report_has_every_hypothesis_proof_bridge_and_scope_barrier():
    normalized = " ".join(REPORT.read_text(encoding="utf-8").split())
    proof_bridges = (
        "smooth connected genus-one curve",
        r"N:=\mathcal O_F(F)\in\operatorname {Pic}^0(F)",
        r"\mathcal O_S((r-1)F)",
        r"\longrightarrow N^r",
        r"H^0(S,\mathcal O_S(K_S-rF))^\vee=0",
        r"H^1(S,(m-1)F)=0",
        r"the section \(u^q\) restricts nontrivially",
        "base-point-free pencil",
        r"\mathcal O_C(3)(-p_1-\cdots-p_9)",
        "translate of the Abel--Jacobi isomorphism",
        "not a claim that the two complete nef cones are marked-isomorphic",
    )
    for phrase in proof_bridges:
        assert phrase in normalized

    exclusions = (
        "every isotropic weak-log class has a smooth genus-one representative",
        "displayed blowups realize a prescribed Coxeter polytope",
        "complete nef chamber or all effective facet curves",
        "finite hyperbolic volume",
        "reflection group is saturated or maximal",
        "exhaustive for raw arithmetic polytopes",
        "finite replay does not replace the cohomological or moduli arguments",
    )
    for phrase in exclusions:
        assert phrase in normalized


def test_independent_audit_is_not_the_manual_audit_and_freezes_current_pins():
    text = AUDIT.read_text(encoding="utf-8")
    normalized = " ".join(text.split())
    required = (
        "PASS, in exactly the stated smooth complex genus-one scope",
        "does not use the conclusions of the earlier manual audit as premises",
        "There is no exponent shift in the last term",
        "does not assume connected fibres",
        "translate of the negative Abel--Jacobi isomorphism",
        "not itself a polytope-to-surface existence theorem",
        "No Nikulin, Vinberg, or ABSW classification is imported",
    )
    for phrase in required:
        assert phrase in normalized

    rows = dict(
        re.findall(r"\| ([^|]+?) \| `([0-9a-f]{64})` \|", text)
    )
    expected = {
        "theorem report": sha256(REPORT),
        "exact certifier": sha256(CERTIFIER),
        "deterministic archive": sha256(ARCHIVE),
        "original regression tests": sha256(ORIGINAL_TEST),
        "independent hostile tests": sha256(THIS_TEST),
    }
    assert rows == expected
    for placeholder in (
        "REPORT_SHA256",
        "CERTIFIER_SHA256",
        "ARCHIVE_SHA256",
        "ORIGINAL_TEST_SHA256",
        "INDEPENDENT_TEST_SHA256",
    ):
        assert placeholder not in text


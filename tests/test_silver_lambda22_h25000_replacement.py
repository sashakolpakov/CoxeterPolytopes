"""Frozen checks for the silver Lambda22 H25000 persistence certificate."""

from __future__ import annotations

import json
from math import isqrt
from pathlib import Path
import runpy


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/silver_lambda22_h25000_replacement.json"
SCRIPT = ROOT / "scripts/audit_silver_lambda22_h25000_replacement.py"


def payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def exact_sign(rational: int, radical: int) -> int:
    sr = (rational > 0) - (rational < 0)
    sb = (radical > 0) - (radical < 0)
    if not sr:
        return sb
    if not sb or sr == sb:
        return sr
    return sr if rational * rational - 2 * radical * radical > 0 else sb


def test_frozen_frontier_counts_hash_and_scope() -> None:
    frozen = payload()
    frontier = frozen["frontier"]
    assert frontier["control_count"] == 301
    assert frontier["controls_by_norm"] == {
        "1": 41,
        "2": 83,
        "3": 15,
        "4": 41,
        "4 + 2*sqrt(2)": 58,
        "6": 28,
        "12": 15,
        "12 + 6*sqrt(2)": 20,
    }
    assert frontier["retained_gate_compatible_type_count"] == 5_575
    assert frontier["retained_types_by_norm"] == {
        "2": 2_046,
        "4": 2_069,
        "6": 730,
        "12": 730,
    }
    assert frontier["locally_admissible_extension_count_with_multiplicity"] == 2_663_705
    assert frontier["locally_admissible_extensions_by_norm"] == {
        "2": 1_592_152,
        "4": 797_565,
        "6": 182_676,
        "12": 91_312,
    }
    assert frontier["retained_type_stream_sha256"] == (
        "21240a46b23221f1c167ffd5c6267c71eb328a57b0a067fa70e488b3955f2e5f"
    )
    assert frontier["nonnegative_type_count"] == 2
    assert frozen["logical_scope"]["bounded_prefix_only"] is True
    assert frozen["logical_scope"]["vinberg_termination"] is False
    assert frozen["logical_scope"]["global_nonreflectivity"] is False


def test_two_threats_are_exactly_frozen_and_resolved() -> None:
    first, second = payload()["nonnegative_threat_resolutions"]
    assert first["height"]["exact"] == {"a": "6392", "b": "4520"}
    assert first["norm"]["display"] == "4"
    assert first["gram"] == [12772, 6400, 9040]
    assert first["locally_admissible_extension_count"] == 437
    assert first["nonnegative_extension_count"] == 1
    assert first["unique_nonnegative_extension"]["d_pairings"] == [226, 160]
    assert first["forced_relation"] == "v=40d"
    assert first["short_vector"] == "e=u-56d"
    assert first["short_vector_norm"] == 4
    assert first["gate_sublattice_minimum"] == 12

    assert second["height"]["exact"] == {"a": "9588", "b": "6780"}
    assert second["norm"]["display"] == "6"
    assert second["gram"] == [28734, 14400, 20340]
    assert second["locally_admissible_extension_count"] == 327
    assert second["nonnegative_extension_count"] == 1
    assert second["unique_nonnegative_extension"]["d_pairings"] == [339, 240]
    assert second["forced_relation"] == "v=60d"
    assert second["short_vector"] == "e=u-84d"
    assert second["short_vector_norm"] == 6
    assert second["gate_sublattice_minimum"] == 24


def test_critical_extension_theorems_replay_without_the_full_scan() -> None:
    module = runpy.run_path(str(SCRIPT))
    resolutions = [
        module["threat_resolution"](row) for row in module["EXPECTED_THREATS"]
    ]
    assert [row["locally_admissible_extension_count"] for row in resolutions] == [
        437,
        327,
    ]
    assert [row["locally_admissible_stream_sha256"] for row in resolutions] == [
        "238503f78fae77ed01d5453c01ab22ec47691fbe633af77d78245e06c6a0d412",
        "446c1f99685dea1a8b48639b4b4b217002d368d86f876c9b01a2412f5c784716",
    ]
    assert [row["zero_determinant_extension_count"] for row in resolutions] == [6, 4]


def test_rectangular_brute_force_confirms_each_unique_nonnegative_extension() -> None:
    cases = (
        ((12772, 6400, 9040), 2, (176, 124), (226, 160, 0)),
        ((28734, 14400, 20340), 3, (264, 186), (339, 240, 0)),
    )
    expected_counts = (437, 327)
    for (gram, modulus, control, expected), expected_count in zip(
        cases, expected_counts, strict=True
    ):
        u_square, v_square, dot = gram
        rows = []
        nonnegative = []
        for p in range(-isqrt(4 * u_square), isqrt(4 * u_square) + 1):
            if p % modulus:
                continue
            for q in range(-isqrt(4 * v_square), isqrt(4 * v_square) + 1):
                if q % modulus:
                    continue
                determinant = (
                    4 * (u_square * v_square - dot * dot)
                    - v_square * p * p
                    - u_square * q * q
                    + 2 * dot * p * q
                )
                if determinant < 0:
                    continue
                rows.append((p, q, determinant))
                # 87 times the exact margin.  Since
                # 87*(sqrt(2)-1)*(a+b*sqrt(2)) is integral here, compare
                # the resulting coefficient pair without floating point.
                ca, cb = control
                threshold_a = 87 * (-ca + 2 * cb)
                threshold_b = 87 * (ca - cb)
                margin_a = 28 * p - threshold_a
                margin_b = 28 * q - threshold_b
                if exact_sign(margin_a, margin_b) >= 0:
                    nonnegative.append((p, q, determinant))
        assert len(rows) == expected_count
        assert nonnegative == [expected]


def test_h10000_dependency_is_independently_pinned() -> None:
    dependency = payload()["dependency"]
    assert dependency["independent_H10000_verdict"] == "PASS"
    assert dependency["production_H10000_sha256"] == (
        "776a9a58f7d4b323e8167b087cd420def4c7a620d90973b2e99b8ac2c3bd4ad5"
    )
    assert dependency["independent_H10000_audit_sha256"] == (
        "6c2e5782d84180ee53cac58742388e50f86bdd6c98b863b1c48adab889e693e4"
    )

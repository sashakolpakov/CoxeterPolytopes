"""Frozen checks for the clean-room silver Lambda22 H10000 audit."""

from __future__ import annotations

import json
from math import isqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "data/enumerations/silver_lambda22_h10000_replacement_independent_hostile_audit.json"
)
SOURCE = (
    ROOT
    / "scripts/independent_hostile_audit_silver_lambda22_h10000_replacement.py"
)


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


def test_frontier_counts_and_bounded_scope() -> None:
    frozen = payload()
    interval = frozen["new_interval"]
    assert frozen["verdict"] == "PASS"
    assert interval["control_count"] == 129
    assert interval["retained_type_count"] == 1_555
    assert interval["relaxed_nonnegative_type_count"] == 1
    assert interval["controls_by_norm"] == {
        "1": 18,
        "2": 35,
        "3": 6,
        "4": 18,
        "4 + 2*sqrt(2)": 26,
        "6": 12,
        "12": 6,
        "12 + 6*sqrt(2)": 8,
    }
    assert interval["retained_types_by_norm"] == {
        "2": 569,
        "4": 574,
        "6": 206,
        "12": 206,
    }
    assert frozen["cutoff_semantics"]["lower"] == "strict at the identity embedding"
    assert "inclusive" in frozen["cutoff_semantics"]["upper"]
    assert frozen["logical_scope"] == {
        "why_local_gate_is_safe": (
            "Any PSD extension rejected by the explicit [-3,3]^3 gate already "
            "displays an actual nonzero integral Lambda22 combination of norm "
            "1--3, impossible because min(Lambda22)=4."
        ),
        "why_unique_survivor_is_empty": (
            "Its zero v,d Cauchy minor forces v=20d in the positive-definite "
            "lattice; then u-28d is a nonzero norm-two Lambda22 vector."
        ),
        "bounded_prefix_only": True,
        "vinberg_termination": False,
        "global_nonreflectivity": False,
    }


def test_critical_665_theorem_is_frozen_exactly() -> None:
    theorem = payload()["critical_extension_theorem"]
    assert theorem["gram"] == [3194, 1600, 2260]
    assert theorem["integral_PSD_extension_count"] == 665
    assert theorem["local_gate_rejection_count"] == 0
    assert theorem["locally_admissible_extension_count"] == 665
    assert theorem["locally_admissible_stream_sha256"] == (
        "3bc3dc04c50765b15ba1be158bdefc6d6bbbbae0e8e90b6d721a659b85556594"
    )
    assert theorem["zero_determinant_extension_count"] == 4
    assert theorem["nonnegative_extension_count"] == 1
    extension = theorem["unique_nonnegative_extension"]
    assert extension == {
        "d_pairings": [113, 80],
        "extension_determinant": 0,
        "margin": {"a": "32/87", "b": "-22/87"},
        "margin_numerator": [32, -22],
        "margin_numerator_norm": 56,
    }
    assert theorem["forced_relation_check"] == {
        "norm_v_minus_20d": 0,
        "u_dot_v_minus_20d": 0,
        "positive_definiteness_forces": "v=20d",
    }
    assert theorem["short_vector"] == {
        "vector": "u-28d",
        "norm": 2,
        "contradicts_Lambda22_minimum": 4,
    }


def test_rectangular_brute_replay_has_665_and_one_nonnegative() -> None:
    """Use a second, deliberately simple loop rather than the audit ellipse."""

    u_square, v_square, dot = 3194, 1600, 2260
    psd = []
    nonnegative = []
    for p in range(-isqrt(4 * u_square), isqrt(4 * u_square) + 1):
        for q in range(-isqrt(4 * v_square), isqrt(4 * v_square) + 1):
            determinant = (
                4 * (u_square * v_square - dot * dot)
                - v_square * p * p
                - u_square * q * q
                + 2 * dot * p * q
            )
            if determinant < 0:
                continue
            psd.append((p, q, determinant))
            margin_a = 28 * p - 3132
            margin_b = 28 * q - 2262
            if exact_sign(margin_a, margin_b) >= 0:
                nonnegative.append((p, q, determinant, margin_a, margin_b))
    assert len(psd) == 665
    assert nonnegative == [(113, 80, 0, 32, -22)]


def test_source_is_clean_room_and_continuity_is_pinned() -> None:
    source = SOURCE.read_text()
    assert "from coxeter_polytopes" not in source
    assert "import coxeter_polytopes" not in source
    assert "def enumerate_controls" in source
    assert "def candidate_grams" in source
    assert "def extensions" in source
    assert "exact_ellipse_inequality" in source
    frozen = payload()
    assert frozen["continuity"]["base_ray_strict_on_all_four_accepted_cutters"] is True
    assert frozen["continuity"]["Lambda22_minimum"] == 4
    assert frozen["replacement_ray"]["identity_sign"] == 1

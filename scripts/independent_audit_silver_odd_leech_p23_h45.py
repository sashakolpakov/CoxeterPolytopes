#!/usr/bin/env python3
"""Independent exact regression for the silver odd-Leech H45 certificate."""

from __future__ import annotations

import json
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_e8_cubes import _bounded_for_caps


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/silver_odd_leech_p23_first_packet.json"
REPORT = ROOT / "results/silver_odd_leech_p23_h45_independent_audit.md"

A = Silver(-3, 4)
EPS = Silver(1, 1)
NORMS = (Silver(1), Silver(2, 1), Silver(2), Silver(5, 1), Silver(12, 7), Silver(10, 2))
P23 = frozenset(NORMS[3:])
ENDPOINT = Silver(45, 32) / 2


def triples(target: Silver) -> tuple[tuple[int, int, int], ...]:
    if not target.is_integral or int(target.b) % 2:
        return ()
    aa = int(target.a)
    bb = int(target.b) // 2
    output = []
    for vv in range(aa // 2 + 1):
        uu = aa - 2 * vv
        if (uu not in (0,) and uu < 3) or (vv not in (0,) and vv < 3):
            continue
        if bb * bb > uu * vv:
            continue
        determinant = uu * vv - bb * bb
        if determinant == 0:
            if gcd(uu, gcd(vv, abs(bb))) < 3:
                continue
        elif any(
            0 < i * i * uu + j * j * vv + 2 * i * j * bb < 3
            for i in range(-3, 4)
            for j in range(-3, 4)
            if i or j
        ):
            continue
        output.append((uu, vv, bb))
    return tuple(output)


def controls() -> tuple[tuple[Silver, Silver, Silver, Silver, tuple[tuple[int, int, int], ...]], ...]:
    output = set()
    for norm in NORMS:
        for control in _bounded_for_caps((norm * 46 / A, norm / (-A))):
            if control.sign(0) <= 0 or not (2 * A * control / norm).is_integral:
                continue
            target = norm + A * control * control
            if target.sign(1) < 0:
                continue
            height = A * control * control / norm
            if (ENDPOINT - height).sign(0) < 0:
                continue
            reduced = target / (A * A) if norm in P23 else target
            output.add((height, norm, control, reduced, triples(reduced)))
    return tuple(sorted(output, key=lambda row: row[:3]))


def joint_maximum(
    aa: int, cc: int, bb: int, coefficient: Silver
) -> tuple[int, int, Silver]:
    best = None
    for mm in range(-isqrt(4 * aa), isqrt(4 * aa) + 1):
        for nn in range(-isqrt(4 * cc), isqrt(4 * cc) + 1):
            determinant = 4 * aa * cc + 2 * bb * mm * nn - aa * nn * nn - cc * mm * mm - 4 * bb * bb
            if determinant < 0:
                continue
            value = Silver(mm) + coefficient * nn
            if best is None or (value - best[2]).sign(0) > 0:
                best = (mm, nn, value)
    assert best is not None
    return best


def main() -> None:
    rows = controls()
    nonempty = tuple(row for row in rows if row[-1])
    assert len(rows) == 22
    assert nonempty == (
        (Silver(7, 6) / 2, Silver(2), EPS, Silver(9, 6), ((3, 3, 3),)),
        (
            ENDPOINT,
            Silver(2),
            Silver(3, 2),
            Silver(47, 32),
            ((27, 10, 16), (25, 11, 16), (23, 12, 16)),
        ),
    )
    maxima = (
        joint_maximum(3, 10, -4, Silver(2, 1)),
        joint_maximum(4, 11, 5, EPS),
        joint_maximum(3, 12, 4, EPS),
    )
    assert maxima == (
        (-2, 6, Silver(10, 6)),
        (3, 6, Silver(9, 6)),
        (3, 6, Silver(9, 6)),
    )
    scale = Fraction(5, 6)
    wall = A * Silver(3, 2)
    margins = tuple(wall - scale * row[2] for row in maxima)
    assert margins == (
        Silver(Fraction(-4, 3), 1),
        Silver(Fraction(-1, 2), 1),
        Silver(Fraction(-1, 2), 1),
    )
    assert all(margin.sign(0) > 0 for margin in margins)
    first_margin = A - scale
    excess = Silver(4) * scale * scale - A
    assert first_margin.sign(0) > 0
    assert excess == Silver(Fraction(52, 9), -4) and excess.sign(0) > 0

    payload = json.loads(ARTIFACT.read_text())
    frozen = payload["h45_necessary_wall_prefix"]
    assert frozen["complete_control_count"] == 22
    assert frozen["replacement_witness_scale"] == "5/6"
    assert frozen["global_minimum_wall_margin"] == {"a": "-4/3", "b": "1"}
    assert frozen["hyperideal_excess"] == {"a": "52/9", "b": "-4"}

    REPORT.write_text(
        """# Independent audit: silver odd-Leech p23 prefix through H45

PASS.

This script does not import the odd-Leech backend.  It independently repeats
the six mirror-norm control boxes, the p23 target division, the binary Gram
census, and the finite integral `3 x 3` positive-semidefinite extension scan.

- complete normalized controls through `45/2+16sqrt(2)`: **22**;
- nonempty necessary levels: the 4096-root first packet and the H45 endpoint;
- H45 necessary Gram types: `(27,10,16)`, `(25,11,16)`, `(23,12,16)`;
- joint unscaled PSD maxima: `10+6sqrt(2)`, `9+6sqrt(2)`, `9+6sqrt(2)`;
- the unscaled direction is not used for the H45 conclusion;
- replacement: `z=(5/6)d`;
- exact H45 margins: `-4/3+sqrt(2)`, `-1/2+sqrt(2)`, `-1/2+sqrt(2)`;
- exact hyperideal excess: `52/9-4sqrt(2)>0`.

The independent values match the frozen JSON byte-for-value at the critical
replacement scale, minimum margin, and hyperideal excess.
"""
    )
    print(REPORT.relative_to(ROOT))


if __name__ == "__main__":
    main()

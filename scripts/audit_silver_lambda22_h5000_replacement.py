#!/usr/bin/env python3
"""Certify that the H1598 replacement ray survives through real height 5000.

This continues the exact bounded-prefix certificate for
``<-(sqrt(2)-1)> + Lambda22``.  The control bound is the rational field
element 5000, hence the identity-embedding height cutoff is exactly 5000.
It is not Vinberg termination or a global nonreflectivity theorem.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import lcm
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_lambda22 import (
    LAMBDA22_COEFFICIENT,
    enumerate_lambda22_controls,
    gated_relaxed_witness_bound,
)
from coxeter_polytopes.search.silver_leech import candidate_leech_gram_triples
from coxeter_polytopes.search.silver_leech_manifest import (
    has_short_binary_combination,
)


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/silver_lambda22_h5000_replacement.json"
REPORT_PATH = ROOT / "results/silver_lambda22_h5000_replacement.md"
START = Silver(1598, 1130)
BOUND = Silver(5000)
SCALE = Silver(Fraction(10, 31))
NORM_GATES = {
    Silver(1): (1, 1),
    Silver(2): (1, 1),
    Silver(4, 2): (2, 1),
    Silver(4): (2, 2),
    Silver(3): (3, 3),
    Silver(6): (3, 3),
    Silver(12, 6): (6, 3),
    Silver(12): (6, 6),
}
MINIMUM = {1: 4, 2: 12, 3: 24, 6: 96}
EXPECTED = (
    (Silver(1687, 1193), Silver(2), Silver(64, 45), (1688, 844, 1193), (82, 58), Silver(Fraction(14, 31), Fraction(-9, 31))),
    (Silver(1687, 1193), Silver(4), Silver(90, 64), (3376, 1688, 2386), (116, 82), Silver(Fraction(-18, 31), Fraction(14, 31))),
    (Silver(1912, 1352), Silver(2), Silver(68, 48), (1898, 964, 1352), (87, 62), Silver(Fraction(2, 31))),
    (Silver(2009, 1421), Silver(4), Silver(98, 70), (3976, 2032, 2842), (126, 90), Silver(Fraction(-42, 31), Fraction(32, 31))),
    (Silver(2151, 1521), Silver(2), Silver(72, 51), (2120, 1092, 1521), (92, 66), Silver(Fraction(-10, 31), Fraction(9, 31))),
    (Silver(2151, 1521), Silver(4), Silver(102, 72), (4368, 2120, 3042), (132, 92), Silver(Fraction(18, 31), Fraction(-10, 31))),
    (Silver(2254, 1594), Silver(2), Silver(74, 52), (2262, 1124, 1594), (95, 67), Silver(Fraction(20, 31), Fraction(-12, 31))),
)


def norm_gate(square: int, modulus: int) -> bool:
    return (
        (not square or square >= MINIMUM[modulus])
        and square % (2 * modulus) == 0
        and (modulus % 2 or square % 8 in (0, 4))
    )


def value(value: Silver) -> dict[str, object]:
    return {"exact": value.to_json(), "display": str(value), "signs": [value.sign(0), value.sign(1)]}


def scan() -> tuple[int, int, tuple[tuple[object, ...], ...]]:
    controls = retained = 0
    rows = []
    for level in enumerate_lambda22_controls(BOUND):
        if level.height <= START:
            continue
        controls += 1
        um, vm = NORM_GATES[level.norm]
        for item in candidate_leech_gram_triples(level.target):
            if has_short_binary_combination(item):
                continue
            gram = (item.u_square, item.v_square, item.u_dot_v)
            if (
                not norm_gate(gram[0], um)
                or not norm_gate(gram[1], vm)
                or gram[2] % lcm(um, vm)
            ):
                continue
            try:
                du, dv, maximum = gated_relaxed_witness_bound(*gram, um, vm)
            except AssertionError:
                continue
            retained += 1
            margin = SCALE * maximum - LAMBDA22_COEFFICIENT * level.control
            if margin.sign(0) >= 0:
                rows.append((level.height, level.norm, level.control, gram, (du, dv), margin))
    result = tuple(rows)
    if (controls, retained, result) != (69, 635, EXPECTED):
        raise AssertionError(f"H5000 frontier changed: {controls}, {retained}, {result!r}")
    return controls, retained, result


def resolution(index: int) -> str:
    if index in (0, 2):
        return "u-v-6d has norm 2, contradicting min(Lambda22)=4."
    if index == 4:
        return "u-v-7d has norm 2, contradicting min(Lambda22)=4."
    if index == 6:
        return "u-v-7d has norm 2 (equivalently another short combination does), contradicting min(Lambda22)=4."
    if index == 1:
        return "u-v-8d lies in S2 and has norm 4, contradicting min(S2)=12."
    if index == 3:
        return "u-v-8d lies in S2 and has norm 4, contradicting min(S2)=12."
    if index == 5:
        return "u-v-10d lies in S2 and has norm 4, contradicting min(S2)=12."
    raise AssertionError(index)


def build_payload() -> dict[str, object]:
    controls, retained, rows = scan()
    norm = Silver(4) * SCALE * SCALE - LAMBDA22_COEFFICIENT
    if norm != Silver(Fraction(1361, 961), -1) or norm.sign(0) <= 0:
        raise AssertionError("replacement norm changed")
    return {
        "schema_version": "0.1.0",
        "seed_id": "silver-H22-Lambda22-a-sqrt2-1",
        "scope": "Exact necessary-superset obstruction through identity-embedding height 5000; bounded prefix only.",
        "height_bound": value(BOUND),
        "replacement_ray": "(1,(10/31)d), d^2=4",
        "replacement_norm": value(norm),
        "new_control_count": controls,
        "retained_gate_compatible_type_count": retained,
        "nonnegative_relaxed_row_count": len(rows),
        "rows": [
            {
                "height": value(row[0]),
                "norm": value(row[1]),
                "control": value(row[2]),
                "gram": list(row[3]),
                "maximizing_d_pairings": list(row[4]),
                "margin": value(row[5]),
                "resolution": resolution(index),
            }
            for index, row in enumerate(rows)
        ],
        "conclusion": "All seven relaxed threats are arithmetically empty; the strict hyperideal ray survives every possible wall through real height 5000.",
        "discovery": None,
    }


def report() -> str:
    return """# Silver Lambda22 replacement survives through real height 5000

For `<-(sqrt(2)-1)> + Lambda22` in `H22`, the exact ray

`z=(1,(10/31)d)`,  `d^2=4`,

has norm `1361/961-sqrt(2)>0`.  Starting immediately after the accepted
wall at `1598+1130sqrt(2)`, the complete exact control scan through the
rational field-element bound `H=5000` contains 69 new controls and 635
gate-compatible necessary Gram types.

Only seven relaxed types can pair nonnegatively with `z`.  Four force a
nonzero Lambda22 vector of norm two.  The other three force a vector in the
dyadic radical sublattice `S2` of norm four, below `min(S2)=12`.  Therefore
all seven are empty and `z` is strict against every possible wall through
the cutoff.

The bounded prefix still has infinite volume.  This is not Vinberg
termination and does not prove global nonreflectivity.

```sh
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h5000_replacement.py
```
"""


def main() -> None:
    JSON_PATH.write_text(json.dumps(build_payload(), indent=2, sort_keys=True) + "\n")
    REPORT_PATH.write_text(report())
    print(JSON_PATH)
    print(REPORT_PATH)


if __name__ == "__main__":
    main()

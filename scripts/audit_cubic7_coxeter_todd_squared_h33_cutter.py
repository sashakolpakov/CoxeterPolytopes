#!/usr/bin/env python3
"""Freeze the first cutter and replacement for cubic-seven ``K12^2``."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_coxeter_todd_squared import (
    CUTTING_HEIGHT,
    OLD_WITNESS_SCALE,
    REPLACEMENT_WITNESS_SCALE,
    TWO_THETA_K12_SEED,
    cubic7_k12_controls,
    cubic7_k12_necessary_gram_types,
    cutting_prefix_directional_maxima,
    cutting_root,
    cutting_root_audit,
    dual_witness,
    dual_witness_norm,
    replacement_prefix_pairing_bounds,
)


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/cubic7_coxeter_todd_squared_h33_cutter.json"
REPORT_PATH = ROOT / "results/cubic7_coxeter_todd_squared_h33_cutter.md"
DEPENDENCY = ROOT / "data/enumerations/cubic7_coxeter_todd_squared_h100.json"


def _value(value: Cubic7) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(place) for place in range(3)],
    }


def build_payload() -> dict[str, object]:
    if not DEPENDENCY.exists():
        raise AssertionError("the H100 coefficient-screen dependency is missing")
    audit = cutting_root_audit()
    controls = tuple(
        level
        for level in cubic7_k12_controls(TWO_THETA_K12_SEED, 1000)
        if level.height <= CUTTING_HEIGHT
    )
    nonempty = tuple(
        (level, types)
        for level in controls
        if (types := cubic7_k12_necessary_gram_types(TWO_THETA_K12_SEED, level))
    )
    observed = tuple(
        (level.height, level.norm, len(types)) for level, types in nonempty
    )
    expected = (
        (Cubic7(5, 14, 6), Cubic7(2), 1),
        (Cubic7(16, 45, 20), Cubic7(2), 1),
        (Cubic7(25, 70, 31), Cubic7(2), 1),
        (CUTTING_HEIGHT, Cubic7(2), 2),
    )
    if observed != expected:
        raise AssertionError(f"the cutter chronology changed: {observed!r}")
    old_norm = dual_witness_norm(OLD_WITNESS_SCALE)
    replacement_norm = dual_witness_norm(REPLACEMENT_WITNESS_SCALE)
    directions = cutting_prefix_directional_maxima()
    replacement_bounds = replacement_prefix_pairing_bounds()
    if old_norm.sign(0) <= 0 or replacement_norm.sign(0) <= 0:
        raise AssertionError("a displayed witness ceased to be hyperideal")
    if audit.old_witness_pairing != Cubic7(0, 1) or audit.old_witness_pairing.sign(0) <= 0:
        raise AssertionError("the explicit wall no longer cuts the old witness")
    if audit.replacement_witness_pairing.sign(0) >= 0:
        raise AssertionError("the replacement is not strict at the cutter")
    if any(
        value.sign(0) > 0
        for value in (
            audit.first_packet_maximum,
            audit.second_packet_maximum,
            audit.third_packet_maximum,
            audit.other_cutting_type_maximum,
            audit.same_height_maximum,
        )
    ):
        raise AssertionError("the explicit cutter lost chronological acceptance")
    return {
        "schema": "coxeter.cubic7_k12_squared_h33_cutter.v1",
        "host": "<-2*theta> + K12 + K12",
        "hyperbolic_dimension": 24,
        "dependency": {
            "path": str(DEPENDENCY.relative_to(ROOT)),
            "sha256": hashlib.sha256(DEPENDENCY.read_bytes()).hexdigest(),
        },
        "scope": (
            "Complete primitive control and necessary coefficient-Gram chronology "
            "through H=33+94theta+42theta^2; explicit accepted first cutter of "
            "the 3/2 dual witness and a strict 7/5 replacement through the same "
            "height. This is not Vinberg termination."
        ),
        "chronology": {
            "control_count": len(controls),
            "nonempty_necessary_level_count": len(nonempty),
            "necessary_gram_type_count": sum(len(types) for _level, types in nonempty),
            "nonempty_levels": [
                {
                    "height": _value(level.height),
                    "norm": _value(level.norm),
                    "control": _value(level.control),
                    "target": _value(level.target),
                    "gram_types": [
                        {
                            "matrix": gram.matrix,
                            "determinant": gram.determinant,
                            "rank": gram.rank,
                        }
                        for gram in types
                    ],
                }
                for level, types in nonempty
            ],
        },
        "first_cutting_wall": {
            "height": _value(audit.height),
            "control": _value(audit.control),
            "norm": _value(audit.norm),
            "coefficient_gram": audit.coefficient_gram,
            "coefficient_relation": "2U-2V+3W=0",
            "h_pairings_U_V_W": audit.h_pairings,
            "root_vector": [value.to_json() for value in cutting_root()],
            "primitive": True,
            "primitive_proof": (
                "The root has norm 2. A nonunit divisor would have square ideal "
                "dividing inert (2), impossible in O_K; equivalently (2) is not "
                "the square of an integral nonunit ideal."
            ),
            "reflection_integral": True,
            "reflection_integrality_proof": (
                "For norm 2 the reflection coefficient is the host pairing. "
                "The time pairing 2theta*x and every K12 Gram pairing are integral."
            ),
            "old_witness_pairing": _value(audit.old_witness_pairing),
            "earlier_packet_maxima": [
                _value(audit.first_packet_maximum),
                _value(audit.second_packet_maximum),
                _value(audit.third_packet_maximum),
            ],
            "same_height": {
                "other_necessary_type_schur_candidate_count": (
                    audit.other_cutting_type_schur_count
                ),
                "maximum_pairing_with_other_necessary_type": _value(
                    audit.other_cutting_type_maximum
                ),
                "schur_divisibility_candidate_count": audit.same_height_schur_count,
                "minimum_four_candidate_count": audit.same_height_minimum_gate_count,
                "maximum_pairing_with_distinct_root": _value(audit.same_height_maximum),
            },
            "accepted": True,
        },
        "replacement_hyperideal_ray": {
            "construction": "z'=(1,(7/5)h,0), h=G_K12^-1*e1",
            "vector": [
                value.to_json() for value in dual_witness(REPLACEMENT_WITNESS_SCALE)
            ],
            "scale": str(REPLACEMENT_WITNESS_SCALE),
            "norm": _value(replacement_norm),
            "pairing_with_cutter": _value(audit.replacement_witness_pairing),
            "necessary_directional_maxima_by_gram_type": [
                _value(value) for value in directions
            ],
            "strict_pairing_bounds_by_gram_type": [
                _value(value) for value in replacement_bounds
            ],
            "strict_against_every_possible_wall_through_bound": True,
            "finite_volume_prefix": False,
            "compact_prefix": False,
        },
        "conclusion": {
            "old_witness_is_cut": True,
            "replacement_is_strict_hyperideal": True,
            "prefix_through_bound_has_infinite_volume": True,
            "vinberg_terminated": False,
            "compact_example": False,
            "finite_volume_example": False,
        },
    }


def render_report(payload: dict[str, object]) -> str:
    cut = payload["first_cutting_wall"]
    replacement = payload["replacement_hyperideal_ray"]
    chronology = payload["chronology"]
    return f"""# Cubic-seven Coxeter--Todd-squared host through its first cutter

## Outcome

For

```text
<-2*theta> + K12 + K12  in H24,
```

the point `z=(1,(3/2)h,0)`, `h=G_K12^-1 e1`, is first cut at

```text
H = {cut['height']['display']}.
```

An explicit primitive accepted norm-two wall pairs with `z` by exactly
`theta>0`.  The replacement `z'=(1,(7/5)h,0)` has positive norm
`{replacement['norm']['display']}` and is strict against every possible wall
through the same height.  Therefore this complete prefix still has infinite
volume.  This is a bounded obstruction, not a compact example and not a
Vinberg-termination claim.

## Complete chronology

There are {chronology['control_count']} locally admissible controls through the
cutting height.  Exactly {chronology['nonempty_necessary_level_count']} levels,
containing {chronology['necessary_gram_type_count']} necessary coefficient
Gram types, survive.  The first three are the previously accepted packets at

```text
5+14theta+6theta^2,
16+45theta+20theta^2,
25+70theta+31theta^2.
```

The last height has two necessary types.  The cutting type is realized by

```text
U = {cut['root_vector'][1:13]}
```

in archived coordinate form, with coefficient Gram

```text
[[12,30,12],[30,84,36],[12,36,16]],
2U-2V+3W=0,
(h.U,h.V,h.W)=(4,10,4).
```

Consequently

```text
r=(1+6theta+3theta^2; U+theta*V+theta^2*W, 0),  r^2=2,
```

and `(z,r)=theta`.

## Chronological acceptance

Complete exact cross-Gram maximization over the minimum-four necessary
supersets gives maximum root pairings

```text
0,
-3theta-2theta^2,
-4-8theta-3theta^2
```

against the three earlier packets.  Thus every earlier wall is nonpositive.
Against the other necessary type at the cutting height, the full relaxed
Schur/divisibility superset has maximum
`{cut['same_height']['maximum_pairing_with_other_necessary_type']['display']}<0`.
For its own type, a rank-two Schur/divisibility enumeration retains
{cut['same_height']['schur_divisibility_candidate_count']} cross blocks; the
exact integral minimum-four gate retains
{cut['same_height']['minimum_four_candidate_count']}.  Excluding the identical
root, their maximum is `{cut['same_height']['maximum_pairing_with_distinct_root']['display']}<0`.
The cutter is therefore accepted, not merely locally realizable.

## Replacement

The replacement has

```text
(z',z') = 196/75 - 2theta > 0,
(z',r) = -2/5 - (2/5)theta^2 < 0.
```

The exact necessary directional maxima at the three earlier types and the
two cutting-height types are respectively `(0+4theta+2theta^2)`,
`(2+6theta+3theta^2)`, `(2+8theta+4theta^2)`,
`(3+10theta+5theta^2)`, and `(4+10theta+4theta^2)`.  The fourth entry is
the other necessary H33 type, not the explicit cutter type.  Inserting
scale `7/5` gives a strict negative pairing at all five types.  Empty
necessary levels contain no wall.  Hence `z'` survives the complete prefix.

## Reproduction

```bash
PYTHONPATH=src python3 scripts/audit_cubic7_coxeter_todd_squared_h33_cutter.py --verify
PYTHONPATH=src pytest -q tests/test_cubic7_coxeter_todd_squared_h33_cutter.py
```
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    report = render_report(payload)
    if args.verify:
        if JSON_PATH.read_text() != encoded or REPORT_PATH.read_text() != report:
            raise SystemExit("frozen cubic-seven K12^2 cutter artifact drifted")
        print("PASS")
        return
    JSON_PATH.write_text(encoded)
    REPORT_PATH.write_text(report)
    print(JSON_PATH)
    print(REPORT_PATH)


if __name__ == "__main__":
    main()

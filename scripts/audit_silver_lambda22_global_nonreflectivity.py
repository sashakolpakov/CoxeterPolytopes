#!/usr/bin/env python3
"""Freeze the global nonreflectivity theorem for the silver Lambda22 host."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from coxeter_polytopes.classification.silver_lambda22_global_nonreflective import (
    BASEPOINT,
    DIRECTION,
    FRONTIER_ARTIFACT,
    FRONTIER_HEIGHT_BOUND,
    GLOBAL_NONREFLECTIVITY_PROOF,
    ISOMETRY_A,
    ISOMETRY_B,
    ISOMETRY_DOT,
    ISOMETRY_SHIFT,
    ORBIT_ONE,
    SEPARATING_CONTROL_BOUND,
    SEPARATING_HEIGHT_BOUND,
    TIME_ZERO_NORM_CASES,
    infinite_order_audit,
    isometry_audit,
    no_time_zero_mirror_audit,
    orbit_recurrence,
    separator_implies_positive_frontier_pairing,
    separating_mirror_bound,
)
from coxeter_polytopes.exact.silver import Silver


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/silver_lambda22_global_nonreflectivity.json"
REPORT_PATH = ROOT / "results/silver_lambda22_global_nonreflectivity.md"
H14_ARTIFACT = ROOT / "data/enumerations/silver_lambda22_h14.json"


def value(item: Silver) -> dict[str, object]:
    return {
        "exact": item.to_json(),
        "display": str(item),
        "signs": [item.sign(0), item.sign(1)],
    }


def vector(item: tuple[Silver, ...]) -> list[dict[str, str]]:
    return [entry.to_json() for entry in item]


def build_payload() -> dict[str, object]:
    if not isometry_audit():
        raise AssertionError("the full integral-isometry audit failed")
    if not infinite_order_audit():
        raise AssertionError("the infinite-order audit failed")
    if not no_time_zero_mirror_audit():
        raise AssertionError("the regular-basepoint audit failed")
    if not separator_implies_positive_frontier_pairing():
        raise AssertionError("the separating-mirror implication failed")
    frontier = json.loads(FRONTIER_ARTIFACT.read_text())
    h14 = json.loads(H14_ARTIFACT.read_text())
    if (
        frontier["height_bound"]["display"] != "150"
        or frontier["necessary_superset_scan"]["post_H14_nonnegative_type_count"] != 7
        or "Every possible post-H14 root through H<=150" not in frontier["conclusion"]
        or not h14["scope"].startswith("Complete exact Vinberg prefix through H=")
        or h14["strict_hyperideal_witness"]["ordinary_max_pairing"]["signs"][0] >= 0
        or h14["strict_hyperideal_witness"]["dyadic_max_pairing"]["signs"][0] >= 0
    ):
        raise AssertionError("the frozen H150 frontier no longer proves its scope")
    bound = separating_mirror_bound()
    if bound.height_bound != SEPARATING_HEIGHT_BOUND or not (
        bound.height_bound < bound.frontier_bound
    ):
        raise AssertionError("the separator bound no longer lies in the frontier")
    recurrence = []
    for exponent in range(7):
        x, coefficient = orbit_recurrence(exponent)
        recurrence.append(
            {
                "exponent": exponent,
                "time": value(x),
                "direction_coefficient": value(coefficient),
            }
        )
    return {
        "schema_version": "0.1.0",
        "seed_id": "silver-H22-Lambda22-global-nonreflectivity",
        "field": "Q(sqrt(2))",
        "ring": "Z[sqrt(2)]",
        "hyperbolic_dimension": 22,
        "lattice": "<-(sqrt(2)-1)> + Lambda22",
        "arithmetic_host": {
            "identity_signature": [22, 1],
            "conjugate_positive_definite": True,
            "integral_orthogonal_group_discrete": True,
            "arithmetic_quotient_cocompact": True,
        },
        "scope": (
            "Global theorem for the full reflection subgroup of the arithmetic "
            "host, using the complete H14 base packet and H14<H<=150 frontier "
            "only as a finite separating-mirror lemma."
        ),
        "direction": {
            "basis_coordinates": list(DIRECTION),
            "norm": 4,
        },
        "isometry": {
            "formula": (
                "g(x,y)=((5+4sqrt2)x+(4+3sqrt2)(d,y), "
                "y+((2+sqrt2)x+(1+sqrt2)(d,y))d)"
            ),
            "A": value(ISOMETRY_A),
            "B": value(ISOMETRY_B),
            "dot_coefficient": value(ISOMETRY_DOT),
            "shift_coefficient": value(ISOMETRY_SHIFT),
            "basepoint": vector(BASEPOINT),
            "basepoint_image": vector(ORBIT_ONE),
            "integral": True,
            "integral_inverse": True,
            "full_basis_form_preservation": True,
            "plane_determinant": value(Silver(1)),
            "plane_trace": value(2 * ISOMETRY_A),
            "identity_embedding_trace_greater_than_two": True,
            "infinite_order": True,
            "orbit_recurrence": recurrence,
        },
        "regular_basepoint": {
            "complete_norm_class_count": len(TIME_ZERO_NORM_CASES),
            "time_zero_cases": [
                {
                    "norm": value(q),
                    "trace_half_target": target,
                    "u_minimum": u_minimum,
                    "v_minimum": v_minimum,
                    "impossible": True,
                }
                for q, target, u_minimum, v_minimum in TIME_ZERO_NORM_CASES
            ],
            "no_reflective_root_with_time_coordinate_zero": True,
        },
        "separating_mirror_lemma": {
            "endpoint_inequality": (
                "(d,y) >= ((4-sqrt(2))/2)x for an oriented root whose mirror "
                "separates v0 and gv0"
            ),
            "cauchy_control_square_bound": value(
                SEPARATING_CONTROL_BOUND**2
            ),
            "control_bound": value(SEPARATING_CONTROL_BOUND),
            "height_bound": value(SEPARATING_HEIGHT_BOUND),
            "positive_z_pairing_lower_coefficient": value(Silver(10, -7) / 6),
            "frontier_height_bound": value(FRONTIER_HEIGHT_BOUND),
            "separator_is_inside_frontier": True,
        },
        "frontier_dependency": {
            "path": str(FRONTIER_ARTIFACT.relative_to(ROOT)),
            "sha256": hashlib.sha256(FRONTIER_ARTIFACT.read_bytes()).hexdigest(),
            "base_H14_path": str(H14_ARTIFACT.relative_to(ROOT)),
            "base_H14_sha256": hashlib.sha256(H14_ARTIFACT.read_bytes()).hexdigest(),
            "base_H14_and_post_H14_complete_through_150": True,
        },
        "chamber_argument": {
            "v0_and_gv0_same_regular_reflection_chamber": True,
            "g_permutes_complete_primitive_reflective_root_set": True,
            "g_normalizes_full_reflection_group": True,
            "g_stabilizes_chamber": True,
            "ambient_integral_isometry_group_discrete": True,
            "infinite_cyclic_chamber_stabilizer": True,
            "finite_volume_chamber_has_finite_stabilizer": True,
        },
        "proof": GLOBAL_NONREFLECTIVITY_PROOF,
        "conclusion": {
            "reflection_chamber_has_infinite_volume": True,
            "reflection_subgroup_has_infinite_index": True,
            "host_is_reflective": False,
            "global_nonreflectivity_theorem": True,
        },
    }


def report(payload: dict[str, object]) -> str:
    dep = payload["frontier_dependency"]
    return f"""# Global nonreflectivity of the silver Lambda22 host

Let `a=sqrt(2)-1`, `L=Lambda22`, and

`M=<-a> + L` over `Z[sqrt(2)]`.

At the identity embedding this form has signature `(22,1)`.  At the
conjugate embedding the time coefficient becomes `1+sqrt(2)>0`, while
`Lambda22` stays positive definite, so the conjugate form is positive
definite.  Thus `O+(M)` is a discrete cocompact arithmetic hyperbolic
lattice.

This is a **global theorem**, not a bounded-prefix observation: the full
reflection subgroup of `O+(M)` has infinite index, and its fundamental
reflection chamber has infinite hyperbolic volume.

## Infinite-order chamber stabilizer

Put `d=e1`, so `d^2=4`.  The following integral map is an isometry of `M`:

`g(x,y)=((5+4sqrt(2))x+(4+3sqrt(2))(d,y),`

`         y+((2+sqrt(2))x+(1+sqrt(2))(d,y))d)`.

Its inverse is given by reversing the two off-diagonal signs.  Exact
calculation on all 23 lattice basis vectors proves `g^T B g=B`.  On the
Lorentz plane generated by `(1,0)` and `(0,d)`, its determinant is one and
its trace is `10+8sqrt(2)>2`; therefore `g` has infinite order.

## Why no mirror separates the endpoints

Let `v0=(1,0)`.  Then

`g v0=(5+4sqrt(2),(2+sqrt(2))d`.

Orient a reflective root `r=(x,y)` by `x>0`.  If its mirror separated `v0`
from `g v0`, then

`(d,y) >= ((4-sqrt(2))/2)x`.

Cauchy and `y^2=q+(sqrt(2)-1)x^2` imply the exact bounds

`x^2/q <= 136+96sqrt(2)=(8+6sqrt(2))^2`,

`H=(sqrt(2)-1)x^2/q <= 56+40sqrt(2) < 150`.

Moreover, for the frozen hyperideal test direction `z=(1,d/3)`,

`(z,r) >= ((10-7sqrt(2))/6)x > 0`.

The complete H14 packet together with the complete post-H14 necessary-wall
frontier through H=150 excludes every such reflective root.  The frozen
frontier dependency is `{dep['path']}` with SHA-256
`{dep['sha256']}`.

The complete eight-class norm alphabet also proves that no reflective root
has `x=0`; hence `v0` is regular.  Thus `v0` and `g v0` lie in the same
reflection chamber.  Since an integral isometry permutes the complete set of
primitive reflective roots and conjugates each reflection to the reflection
in the image root, `g` normalizes the full reflection group and maps chambers
to chambers.  Therefore `g` stabilizes this chamber.

## Global conclusion

The integral orthogonal group is discrete and `g` is a hyperbolic
translation.  Proper discontinuity gives a sufficiently small open ball
`U` whose distinct translates `g^n U` are pairwise disjoint.  Because
`gC=C`, choosing `U` inside the chamber puts every translate inside `C`.
They all have the same positive volume.  Therefore the chamber has infinite
volume and the reflection subgroup has infinite index.

Equivalently, for the complete locally finite mirror arrangement, the
standard reflection-chamber theorem makes the reflection group act simply
transitively on chambers.  Its quotient in the ambient arithmetic group is
the chamber stabilizer, which contains the infinite cyclic group generated
by `g`; hence the reflection subgroup has infinite index.  Since the ambient
arithmetic quotient is cocompact and has positive finite volume, the
infinite-degree chamber cover also has infinite volume.

```sh
PYTHONPATH=src python3 scripts/audit_silver_lambda22_global_nonreflectivity.py --freeze
PYTHONPATH=src python3 scripts/audit_silver_lambda22_global_nonreflectivity.py --verify
PYTHONPATH=src pytest -q tests/test_silver_lambda22_global_nonreflectivity.py
```
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--freeze", action="store_true")
    mode.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    rendered = report(payload)
    if args.verify:
        if not JSON_PATH.exists() or not REPORT_PATH.exists():
            raise SystemExit("FAIL: frozen production files are missing")
        if JSON_PATH.read_text() != serialized or REPORT_PATH.read_text() != rendered:
            raise SystemExit("FAIL: frozen global nonreflectivity files changed")
        print("PASS")
        return
    JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    JSON_PATH.write_text(serialized)
    REPORT_PATH.write_text(rendered)
    print(JSON_PATH)
    print(REPORT_PATH)


if __name__ == "__main__":
    main()

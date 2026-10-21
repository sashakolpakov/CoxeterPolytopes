#!/usr/bin/env python3
"""Freeze the universal odd-unimodular Pell-boost theorem."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from coxeter_polytopes.classification.pell_boost_nonreflectivity import (
    CRITERION_THEOREM,
)
from coxeter_polytopes.classification.silver_odd_unimodular_global_nonreflective import (
    CATALOGUE_HOSTS,
    CONTROL_SQUARE_BOUND,
    DIRECTION_FROM_DIRECTION,
    DIRECTION_NORM,
    GLOBAL_CHAMBER_ARGUMENT,
    GLOBAL_FAMILY_THEOREM,
    PELL_A,
    PELL_U,
    SEPARATOR_HEIGHT_BOUND,
    SEPARATOR_SLOPE,
    SPATIAL_FROM_TIME,
    TIME_FROM_DIRECTION,
    UNIVERSAL_NORM_CLASS_PROOF,
    WITNESS_SCALE,
    WITNESS_SEPARATOR_COEFFICIENT,
    boost_audit,
    boost_plane_matrix,
    plane_isometry_audit,
)
from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_s26 import (
    NORM_CLASSES,
    RHO,
    height_zero_gated_gram_triples,
    minimal_direction_audit_through,
)


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/silver_odd_unimodular_global_nonreflectivity.json"
REPORT_PATH = ROOT / "results/silver_odd_unimodular_global_nonreflectivity.md"
NAMED_FAMILY_ARTIFACT = (
    ROOT / "data/enumerations/silver_odd_unimodular_h27_h29_h1000.json"
)
S26_ARTIFACT = ROOT / "data/enumerations/silver_s26_h26_first_packet.json"


def value(item: Silver) -> dict[str, object]:
    return {
        "exact": item.to_json(),
        "display": str(item),
        "signs": [item.sign(0), item.sign(1)],
    }


def build_payload() -> dict[str, object]:
    if not boost_audit() or not plane_isometry_audit():
        raise AssertionError("the universal Pell boost failed")
    if any(height_zero_gated_gram_triples(row) for row in NORM_CLASSES):
        raise AssertionError("a time-zero mirror class survived")
    frontier = minimal_direction_audit_through(SEPARATOR_HEIGHT_BOUND)
    if not frontier.strict_hyperideal:
        raise AssertionError("the witness failed at the separator cutoff")
    expected_counts = (139, 108, 32, 297, 26)
    actual_counts = (
        frontier.control_count,
        frontier.distinct_height_count,
        frontier.nonempty_control_count,
        frontier.gated_binary_type_count,
        frontier.nonempty_height_count,
    )
    if actual_counts != expected_counts:
        raise AssertionError(f"the separator frontier changed: {actual_counts}")
    worst = frontier.worst_row
    return {
        "schema_version": "0.1.0",
        "seed_id": "silver-odd-unimodular-min3-global-nonreflectivity",
        "field": "Q(sqrt(2))",
        "ring": "Z[sqrt(2)]",
        "coefficient": value(RHO),
        "scope": (
            "Universal global theorem for every positive-definite integral "
            "unimodular spatial lattice L of minimum three, with any chosen "
            "d in L of square three. Named catalogue hosts are audited instances, "
            "not an asserted classification of every such L."
        ),
        "hypotheses": {
            "spatial_lattice": "positive-definite integral unimodular",
            "minimum": 3,
            "chosen_direction_square": 3,
            "rank_restriction": None,
            "shell_or_orbit_assumption": False,
        },
        "reusable_criterion": {
            "statement": CRITERION_THEOREM,
            "pell_identity": "A^2-a*m*U^2=1",
            "integrality_conditions": ["a*U integral", "(A-1)/m integral"],
            "separator_height_formula": "H_separator=A^2-1",
        },
        "boost": {
            "formula": (
                "g(x,y)=(A*x+U*(d,y), "
                "y+((a*U)*x+((A-1)/m)*(d,y))*d)"
            ),
            "A": value(PELL_A),
            "U": value(PELL_U),
            "time_from_direction_U": value(TIME_FROM_DIRECTION),
            "spatial_from_time_aU": value(SPATIAL_FROM_TIME),
            "A_minus_1_over_m": value(DIRECTION_FROM_DIRECTION),
            "plane_matrix": [
                [value(entry) for entry in row] for row in boost_plane_matrix()
            ],
            "pell_identity": True,
            "integral": True,
            "integral_inverse": True,
            "form_preserving": True,
            "identity_embedding_trace_greater_than_two": True,
            "infinite_order": True,
        },
        "separator": {
            "direction_slope": value(SEPARATOR_SLOPE),
            "control_square_bound": value(CONTROL_SQUARE_BOUND),
            "height_bound": value(SEPARATOR_HEIGHT_BOUND),
            "equals_A_square_minus_one": True,
            "witness": f"z=(1,{WITNESS_SCALE}*d)",
            "positive_witness_lower_coefficient": value(
                WITNESS_SEPARATOR_COEFFICIENT
            ),
        },
        "complete_reflective_arithmetic": {
            "norm_class_count": len(NORM_CLASSES),
            "norm_classes": [value(row.norm) for row in NORM_CLASSES],
            "proof": UNIVERSAL_NORM_CLASS_PROOF,
            "time_zero_surviving_class_count": 0,
            "basepoint_regular": True,
        },
        "universal_frontier": {
            "height_bound": value(frontier.complete_height_bound),
            "control_count": frontier.control_count,
            "distinct_height_count": frontier.distinct_height_count,
            "nonempty_control_count": frontier.nonempty_control_count,
            "necessary_binary_type_count": frontier.gated_binary_type_count,
            "nonempty_height_count": frontier.nonempty_height_count,
            "witness_scale": value(frontier.witness_scale),
            "witness_lorentz_norm": value(frontier.witness_lorentz_norm),
            "worst_row": {
                "height": value(worst.control.height),
                "norm": value(worst.control.norm_class.norm),
                "control": value(worst.control.control),
                "gram": [
                    worst.triple.u_square,
                    worst.triple.v_square,
                    worst.triple.u_dot_v,
                ],
                "maximum_direction_pairing": value(
                    worst.pairing.maximum_pairing
                ),
                "margin": value(worst.witness_margin),
            },
            "all_margins_strict": frontier.all_margins_strict,
            "rank_independent_reason": (
                "The coefficient/control box, binary Gram gate and augmented "
                "pairing maximum use only integrality, unimodularity, minimum "
                "three, d^2=3 and positive semidefiniteness. They retain a "
                "superset of actual roots in every rank."
            ),
        },
        "named_audited_instances": {
            "hosts": [
                {"name": name, "hyperbolic_dimension": dimension}
                for name, dimension in CATALOGUE_HOSTS
            ],
            "S26_artifact": {
                "path": str(S26_ARTIFACT.relative_to(ROOT)),
                "sha256": hashlib.sha256(S26_ARTIFACT.read_bytes()).hexdigest(),
            },
            "T27_U28_dim29odd_artifact": {
                "path": str(NAMED_FAMILY_ARTIFACT.relative_to(ROOT)),
                "sha256": hashlib.sha256(
                    NAMED_FAMILY_ARTIFACT.read_bytes()
                ).hexdigest(),
            },
        },
        "published_classification_corollary": {
            "scope": (
                "The rootless odd-unimodular orthogonal silver slice only; "
                "not all arithmetic hosts in any dimension."
            ),
            "rank_counts": {
                "23": 1,
                "24": 1,
                "25": 0,
                "26": 1,
                "27": 3,
                "28": 38,
                "29": 10092,
            },
            "total": 10136,
            "eliminated_by_universal_theorem": True,
            "source_roles": [
                {
                    "url": "https://arxiv.org/abs/math/9912236",
                    "role": "Borcherds: unique rootless rank-26 unimodular lattice",
                },
                {
                    "url": "https://arxiv.org/abs/math/0012231",
                    "role": (
                        "King mass formula verifies Bacher--Venkov rootless "
                        "enumerations in ranks 27 and 28"
                    ),
                },
                {
                    "url": "https://arxiv.org/abs/2410.19569",
                    "role": (
                        "Allombert--Chenevier complete rank <=29 context, "
                        "38 rootless rank-28 classes and 10092 rank-29 classes"
                    ),
                },
            ],
            "dependency_note": (
                "The universal Pell-boost theorem is independent of these "
                "classification counts; only this finite-slice corollary uses them."
            ),
        },
        "chamber_logic": {
            "proof": GLOBAL_CHAMBER_ARGUMENT,
            "ambient_group": "Gamma=O^+(M)",
            "reflection_subgroup": (
                "W=<all integral reflections in Gamma>"
            ),
            "no_separating_reflective_mirror": True,
            "boost_normalizes_complete_reflection_group": True,
            "boost_stabilizes_complete_reflection_chamber": True,
            "reflection_subgroup_normal": True,
            "reflection_subgroup_chamber_stabilizer_trivial": True,
            "boost_powers_inject_into_Gamma_mod_W": True,
            "infinite_order_chamber_stabilizer": True,
            "conjugate_form_positive_definite": True,
            "host_anisotropic": True,
            "ambient_arithmetic_quotient_cocompact": True,
            "reflection_chamber_infinite_volume": True,
            "reflection_subgroup_infinite_index": True,
        },
        "theorem": GLOBAL_FAMILY_THEOREM,
        "conclusion": {
            "universal_minimum_three_unimodular_family_nonreflective": True,
            "named_catalogue_instances_nonreflective": [
                name for name, _ in CATALOGUE_HOSTS
            ],
            "field_wide_nonexistence": False,
        },
    }


def report(payload: dict[str, object]) -> str:
    frontier = payload["universal_frontier"]
    production = payload["named_audited_instances"]
    return f"""# A Pell-boost global nonreflectivity criterion and silver family

## General criterion

Let `M=<-a>+L` over a real quadratic integer ring and let `d in L` have
square `m`.  Suppose integral `A,U` satisfy

`A^2-a*m*U^2=1`, `a*U` is integral, and `(A-1)/m` is integral.

Then

`g(x,y)=(A*x+U(d,y), y+((a*U)x+((A-1)/m)(d,y))d)`

is an integral Lorentz isometry with integral inverse.  Any oriented mirror
separating `(1,0)` from its image has height at most

`H_sep=A^2-1`.

For a spacelike test ray `z=(1,t*d)`, the same endpoint inequality supplies
an explicit positive lower bound for `(z,r)`.  Consequently a complete
strict-ray certificate through `H_sep`, together with a complete mirror-norm
alphabet excluding time-zero mirrors, proves that `g` stabilizes the full
reflection chamber.  Since `A>1`, it is an infinite-order hyperbolic
translation; the chamber has infinite volume and the reflection subgroup
has infinite index.  More explicitly, if `Gamma=O^+(M)` and `W` is generated
by all integral reflections, then `W` is normal because `Gamma` permutes all
reflective mirrors.  The `W`-stabilizer of a chamber is trivial, so the powers
of the infinite-order `g` inject into `Gamma/W`.  Thus `[Gamma:W]=infinity`.
For an admissible host `Gamma` is an arithmetic lattice, and infinite index
forces the `W`-chamber to have infinite volume.

## Universal silver theorem

Take

`a=2sqrt(2)-1`, `m=3`, `A=16+12sqrt(2)`, `U=7+5sqrt(2)`.

Then

- `aU=13+9sqrt(2)`;
- `(A-1)/3=5+4sqrt(2)`;
- separating slope `(d,y)/x >= 8-4sqrt(2)`;
- `x^2/q <= 297+210sqrt(2)`;
- `H_sep=543+384sqrt(2)=A^2-1`.

For `z=(1,(43/55)d)`, every separator would satisfy

`(z,r) >= ((399-282sqrt(2))/55)x > 0`.

The complete six-class reflection arithmetic has no time-zero mirror for
any positive-definite integral unimodular `L` of minimum three.  The exact
rank-independent frontier through the separator cutoff has
{frontier['control_count']} controls, {frontier['necessary_binary_type_count']}
necessary binary Gram types and {frontier['nonempty_height_count']} nonempty
heights.  Every augmented Gram pairing is strict.  The worst row remains

`H=455+322sqrt(2), q=2, x=16+11sqrt(2), Gram=(460,226,322)`,

with maximum `37+26sqrt(2)` and margin
`(-51+37sqrt(2))/55>0`.

Therefore:

> For every positive-definite integral unimodular lattice `L` with
> `min(L)=3`, the arithmetic host
> `<-(2sqrt(2)-1)>+L` has an infinite-volume complete reflection chamber;
> its full integral reflection subgroup has infinite index.

This statement is independent of rank and makes no shell-transitivity,
kissing-number, or automorphism-orbit assumption.

At the nonidentity embedding, `-(2sqrt(2)-1)` and the spatial form are both
positive definite.  Hence the host is anisotropic, `Gamma` is cocompact, and
the group-theoretic infinite-index conclusion above applies without a cusp
or finite-volume caveat.

## Published-classification corollary

Combining the universal theorem with the published complete classifications
gives a finite-slice consequence.  The numbers of rootless odd unimodular
lattices in ranks `23,24,25,26,27,28,29` are respectively

`1, 1, 0, 1, 3, 38, 10092`,

for a total of **10,136**.  Every associated orthogonal silver host
`<-(2sqrt(2)-1)>+L` is nonreflective by the theorem above.  This eliminates
the complete rootless odd-unimodular orthogonal slice in hyperbolic
dimensions `23,24,26,27,28,29` (rank 25 is vacuous).  It is emphatically
not whole-dimension or field-wide nonexistence.

The count corollary uses [Borcherds, rank 26](https://arxiv.org/abs/math/9912236),
[King/Bacher--Venkov, ranks 27--28](https://arxiv.org/abs/math/0012231), and
[Allombert--Chenevier, ranks 28--29](https://arxiv.org/abs/2410.19569).
The core Pell-boost theorem does not depend on those classification counts.

The explicitly audited catalogue instances are `S26`, `T27`, `U28`, and
`dim29odd`, pinned by `{production['S26_artifact']['path']}` and
`{production['T27_U28_dim29odd_artifact']['path']}`.  These names are examples,
not an asserted classification of every minimum-three unimodular lattice.
The theorem is not field-wide nonexistence over `Q(sqrt(2))`.

```sh
PYTHONPATH=src python3 scripts/audit_silver_odd_unimodular_global_nonreflectivity.py --freeze
PYTHONPATH=src python3 scripts/audit_silver_odd_unimodular_global_nonreflectivity.py --verify
PYTHONPATH=src pytest -q tests/test_pell_boost_nonreflectivity.py tests/test_silver_odd_unimodular_global_nonreflectivity.py
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
            raise SystemExit("FAIL: frozen universal theorem is missing")
        if JSON_PATH.read_text() != serialized or REPORT_PATH.read_text() != rendered:
            raise SystemExit("FAIL: frozen universal theorem changed")
        print("PASS_GLOBAL_ODD_UNIMODULAR_NONREFLECTIVITY")
        return
    JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    JSON_PATH.write_text(serialized)
    REPORT_PATH.write_text(rendered)
    print(JSON_PATH)
    print(REPORT_PATH)


if __name__ == "__main__":
    main()

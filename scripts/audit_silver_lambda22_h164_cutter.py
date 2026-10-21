#!/usr/bin/env python3
"""Certify the first explicit cutter and a replacement ray for silver Lambda22.

The fixed H14 ray ``(1,d/3)`` is first cut by the explicit H164 root built
from a short pair ``d,e``.  An exact 4-by-4 Schur necessary-superset scan
shows that it is acute to every earlier possible root: nine earlier rows are
only tangent and no earlier row is positive.  A smaller replacement ray
``(1,13d/40)`` stays strict through the whole H164 prefix.

This remains a bounded Vinberg-prefix certificate, not a termination or
global nonreflectivity theorem.
"""

from __future__ import annotations

import json
from fractions import Fraction
from hashlib import sha256
from math import isqrt, lcm
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_lambda22 import (
    FIRST_PACKET_HEIGHT,
    LAMBDA22_COEFFICIENT,
    cross_packet_pairing_from_dot,
    dyadic_packet_pairing_from_dot,
    enumerate_lambda22_controls,
    gated_relaxed_witness_bound,
    lambda22_gram,
    lambda22_pair_silver,
)
from coxeter_polytopes.search.cubic7_lambda22 import (
    lambda22_integer_gram_minimum_at_least,
)
from coxeter_polytopes.search.silver_leech import candidate_leech_gram_triples
from coxeter_polytopes.search.silver_leech_manifest import (
    has_short_binary_combination,
)


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/silver_lambda22_h164_cutter.json"
REPORT_PATH = ROOT / "results/silver_lambda22_h164_cutter.md"
H150_JSON_PATH = ROOT / "data/enumerations/silver_lambda22_h150_frontier.json"
HEIGHT = Silver(164, 116)
CONTROL = Silver(28, 20)
ROOT_NORM = Silver(4)
TARGET = Silver(660, 464)
D = (1,) + (0,) * 21
E = (-6, 31, 21, 15, 10, 20, 14, 7, 4, 2, 2, 1) + (0,) * 10
U = tuple(8 * D[i] - E[i] for i in range(22))
V = tuple(6 * D[i] - E[i] for i in range(22))
SPATIAL = tuple(Silver(U[i], V[i]) for i in range(22))
REPLACEMENT_SCALE = Fraction(13, 40)
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


def value_json(value: Silver) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(0), value.sign(1)],
    }


def norm_gate(square: int, modulus: int) -> bool:
    return (
        (not square or square >= MINIMUM[modulus])
        and square % (2 * modulus) == 0
        and (modulus % 2 or square % 8 in (0, 4))
    )


def cross_extensions(
    u_square: int,
    v_square: int,
    u_dot_v: int,
    u_modulus: int,
    v_modulus: int,
    other_control: Silver,
) -> list[tuple[Silver, int, int, int, int, int, int, int]]:
    """Enumerate nonnegative pairings allowed by the exact 4x4 Schur test.

    The fixed block is ``Gram(d,e)=[[4,-2],[-2,12]]``.  Multiplying the
    Schur complement by its determinant 44 gives diagonal numerators A,B
    and cross numerator N.  PSD is exactly ``A>=0,B>=0,AB>=N^2``.
    """

    threshold = LAMBDA22_COEFFICIENT * CONTROL * other_control
    u_rows = []
    for d_dot_u in range(-isqrt(4 * u_square), isqrt(4 * u_square) + 1):
        if d_dot_u % u_modulus:
            continue
        for e_dot_u in range(-isqrt(12 * u_square), isqrt(12 * u_square) + 1):
            if e_dot_u % lcm(2, u_modulus):
                continue
            a = (
                44 * u_square
                - 12 * d_dot_u * d_dot_u
                - 4 * d_dot_u * e_dot_u
                - 4 * e_dot_u * e_dot_u
            )
            if a < 0:
                continue
            value = Silver(8 * d_dot_u - e_dot_u, 6 * d_dot_u - e_dot_u)
            u_rows.append((value, d_dot_u, e_dot_u, a))
    v_rows = []
    for d_dot_v in range(-isqrt(4 * v_square), isqrt(4 * v_square) + 1):
        if d_dot_v % v_modulus:
            continue
        for e_dot_v in range(-isqrt(12 * v_square), isqrt(12 * v_square) + 1):
            if e_dot_v % lcm(2, v_modulus):
                continue
            b = (
                44 * v_square
                - 12 * d_dot_v * d_dot_v
                - 4 * d_dot_v * e_dot_v
                - 4 * e_dot_v * e_dot_v
            )
            if b < 0:
                continue
            value = Silver(12 * d_dot_v - 2 * e_dot_v, 8 * d_dot_v - e_dot_v)
            v_rows.append((value, d_dot_v, e_dot_v, b))
    u_rows.sort(reverse=True)
    v_rows.sort(reverse=True)
    output = []
    if not u_rows or not v_rows:
        return output
    for u_value, d_dot_u, e_dot_u, a in u_rows:
        if (u_value + v_rows[0][0] - threshold).sign(0) < 0:
            break
        for v_value, d_dot_v, e_dot_v, b in v_rows:
            pairing = u_value + v_value - threshold
            if pairing.sign(0) < 0:
                break
            n = 44 * u_dot_v - (
                12 * d_dot_u * d_dot_v
                + 2 * (d_dot_u * e_dot_v + e_dot_u * d_dot_v)
                + 4 * e_dot_u * e_dot_v
            )
            if a * b < n * n:
                continue
            output.append(
                (pairing, d_dot_u, d_dot_v, e_dot_u, e_dot_v, a, b, n)
            )
    return output


def chronological_cross_scan() -> list[dict[str, object]]:
    rows = []
    for level in enumerate_lambda22_controls(HEIGHT):
        if level.height <= FIRST_PACKET_HEIGHT:
            continue
        u_modulus, v_modulus = NORM_GATES[level.norm]
        for item in candidate_leech_gram_triples(level.target):
            if has_short_binary_combination(item):
                continue
            gram = (item.u_square, item.v_square, item.u_dot_v)
            if (
                not norm_gate(gram[0], u_modulus)
                or not norm_gate(gram[1], v_modulus)
                or gram[2] % lcm(u_modulus, v_modulus)
            ):
                continue
            extensions = cross_extensions(
                *gram, u_modulus, v_modulus, level.control
            )
            if not extensions:
                continue
            rows.append(
                {
                    "height": value_json(level.height),
                    "norm": value_json(level.norm),
                    "control": value_json(level.control),
                    "gram": list(gram),
                    "extensions": [
                        {
                            "pairing": value_json(extension[0]),
                            "d_dot_u": extension[1],
                            "d_dot_v": extension[2],
                            "e_dot_u": extension[3],
                            "e_dot_v": extension[4],
                            "Schur_A": extension[5],
                            "Schur_B": extension[6],
                            "Schur_N": extension[7],
                        }
                        for extension in extensions
                    ],
                }
            )
    if len(rows) != 12:
        raise AssertionError(f"the H164 chronological frontier changed: {len(rows)}")
    if any(
        extension["pairing"]["signs"][0] != 0
        for row in rows[:-1]
        for extension in row["extensions"]
    ):
        raise AssertionError("an earlier possible root became positive")
    final = rows[-1]
    if (
        final["height"]["display"] != str(HEIGHT)
        or final["norm"]["display"] != "4"
        or final["gram"] != [300, 180, 232]
        or final["extensions"]
        != [
            {
                "pairing": value_json(Silver(4)),
                "d_dot_u": 34,
                "d_dot_v": 26,
                "e_dot_u": -28,
                "e_dot_v": -24,
                "Schur_A": 0,
                "Schur_B": 0,
                "Schur_N": 0,
            }
        ]
    ):
        raise AssertionError("the unique H164 positive extension changed")
    return rows


def explicit_root_audit() -> dict[str, object]:
    gram = lambda22_gram()
    d_square = lambda22_pair_silver(D, D)
    e_square = lambda22_pair_silver(E, E)
    d_dot_e = lambda22_pair_silver(D, E)
    spatial_norm = lambda22_pair_silver(SPATIAL, SPATIAL)
    root_norm = -LAMBDA22_COEFFICIENT * CONTROL**2 + spatial_norm
    spatial_pairings = tuple(
        sum((Silver(gram[i][j]) * SPATIAL[j] for j in range(22)), start=Silver())
        for i in range(22)
    )
    reflection_integral = (
        (2 * LAMBDA22_COEFFICIENT * CONTROL / ROOT_NORM).is_integral
        and all((2 * value / ROOT_NORM).is_integral for value in spatial_pairings)
    )
    # The displayed spatial vector actually contains the unit -(1+sqrt(2)),
    # so no nonunit of Z[sqrt(2)] divides the root.
    primitive_unit_coordinate = Silver(-1, -1) in SPATIAL
    primitive = primitive_unit_coordinate
    fixed_ray_pairing = (
        -LAMBDA22_COEFFICIENT * CONTROL
        + lambda22_pair_silver(tuple(Silver(value) / 3 for value in D), SPATIAL)
    )
    if (
        (d_square, e_square, d_dot_e) != (4, 12, -2)
        or (lambda22_pair_silver(U, U), lambda22_pair_silver(V, V), lambda22_pair_silver(U, V))
        != (300, 180, 232)
        or (lambda22_pair_silver(D, U), lambda22_pair_silver(D, V)) != (34, 26)
        or root_norm != ROOT_NORM
        or LAMBDA22_COEFFICIENT * CONTROL**2 / ROOT_NORM != HEIGHT
        or not reflection_integral
        or not primitive
        or fixed_ray_pairing != Silver(-2, 2) / 3
        or fixed_ray_pairing.sign(0) <= 0
    ):
        raise AssertionError("the explicit H164 cutter arithmetic changed")
    return {
        "d": list(D),
        "e": list(E),
        "Gram_d_e": [[4, -2], [-2, 12]],
        "u": list(U),
        "v": list(V),
        "coefficient_gram": [[300, 232], [232, 180]],
        "d_pairings": [34, 26],
        "control": value_json(CONTROL),
        "spatial": "y=u+sqrt(2)v=(8+6sqrt(2))d-(1+sqrt(2))e",
        "norm": value_json(root_norm),
        "height": value_json(HEIGHT),
        "primitive": primitive,
        "primitive_unit_coordinate": str(Silver(-1, -1)),
        "reflection_integral": reflection_integral,
        "fixed_ray_pairing": value_json(fixed_ray_pairing),
    }


def first_packet_acuteness_audit() -> dict[str, object]:
    """Maximize the cutter pairing over both complete H14 shell packets.

    If ``w`` is a shell vector, write ``p=d.w`` and ``q=e.w``.  Cauchy
    bounds make the scan finite.  For every integral extension Gram we use
    the exact quotient-lattice minimum decider, so the rejection of the two
    apparent PSD-only maxima is complete rather than a small-box heuristic.
    """

    output = []
    for name, shell_square, control, pairing_modulus in (
        ("ordinary_q2", 10, Silver(6, 4), 1),
        ("dyadic_q4_S2", 20, Silver(8, 6), 2),
    ):
        best: tuple[Silver, int, int] | None = None
        tested = 0
        minimum_survivors = 0
        for d_dot_w in range(-isqrt(4 * shell_square), isqrt(4 * shell_square) + 1):
            if d_dot_w % pairing_modulus:
                continue
            for e_dot_w in range(-isqrt(12 * shell_square), isqrt(12 * shell_square) + 1):
                if e_dot_w % pairing_modulus:
                    continue
                gram = (
                    (4, -2, d_dot_w),
                    (-2, 12, e_dot_w),
                    (d_dot_w, e_dot_w, shell_square),
                )
                # Principal minors plus the exact minimum routine together
                # cover positive-definite and rank-deficient PSD extensions.
                if (
                    d_dot_w * d_dot_w > 4 * shell_square
                    or e_dot_w * e_dot_w > 12 * shell_square
                ):
                    continue
                determinant = (
                    4 * 12 * shell_square
                    + 2 * (-2) * d_dot_w * e_dot_w
                    - 4 * e_dot_w * e_dot_w
                    - 12 * d_dot_w * d_dot_w
                    - shell_square * 4
                )
                if determinant < 0:
                    continue
                tested += 1
                if not lambda22_integer_gram_minimum_at_least(gram, 4):
                    continue
                minimum_survivors += 1
                spatial_pairing = (
                    Silver(20 * d_dot_w - 3 * e_dot_w, 14 * d_dot_w - 2 * e_dot_w)
                )
                if best is None or spatial_pairing > best[0]:
                    best = spatial_pairing, d_dot_w, e_dot_w
        if best is None:
            raise AssertionError("an H14 shell lost every extension")
        threshold = LAMBDA22_COEFFICIENT * CONTROL * control
        margin = best[0] - threshold
        output.append(
            {
                "packet": name,
                "shell_square": shell_square,
                "pairing_modulus": pairing_modulus,
                "PSD_extension_count": tested,
                "minimum_four_extension_count": minimum_survivors,
                "max_spatial_pairing": value_json(best[0]),
                "maximizer": {"d_dot_w": best[1], "e_dot_w": best[2]},
                "timelike_threshold": value_json(threshold),
                "max_root_pairing": value_json(margin),
            }
        )
    if [
        (row["max_spatial_pairing"]["display"], row["maximizer"], row["max_root_pairing"]["display"])
        for row in output
    ] != [
        (
            "124 + 86*sqrt(2)",
            {"d_dot_w": 5, "e_dot_w": -8},
            "-12 - 10*sqrt(2)",
        ),
        (
            "190 + 132*sqrt(2)",
            {"d_dot_w": 8, "e_dot_w": -10},
            "-2 - 4*sqrt(2)",
        ),
    ]:
        raise AssertionError("the complete H14 packet acuteness maxima changed")
    if any(row["max_root_pairing"]["signs"][0] >= 0 for row in output):
        raise AssertionError("the H164 cutter is no longer acute to H14")
    return {"packets": output, "both_strict": True}


def replacement_ray_audit() -> dict[str, object]:
    nonnegative = []
    retained = 0
    for level in enumerate_lambda22_controls(HEIGHT):
        u_modulus, v_modulus = NORM_GATES[level.norm]
        for item in candidate_leech_gram_triples(level.target):
            if has_short_binary_combination(item):
                continue
            u_square, v_square, u_dot_v = (
                item.u_square,
                item.v_square,
                item.u_dot_v,
            )
            if (
                not norm_gate(u_square, u_modulus)
                or not norm_gate(v_square, v_modulus)
                or u_dot_v % lcm(u_modulus, v_modulus)
            ):
                continue
            try:
                d_dot_u, d_dot_v, maximum = gated_relaxed_witness_bound(
                    u_square,
                    v_square,
                    u_dot_v,
                    u_modulus,
                    v_modulus,
                )
            except AssertionError:
                continue
            retained += 1
            margin = REPLACEMENT_SCALE * maximum - LAMBDA22_COEFFICIENT * level.control
            if margin.sign(0) >= 0:
                nonnegative.append(
                    {
                        "height": value_json(level.height),
                        "norm": value_json(level.norm),
                        "control": value_json(level.control),
                        "gram": [u_square, v_square, u_dot_v],
                        "relaxed_pairings": [d_dot_u, d_dot_v],
                        "margin": value_json(margin),
                    }
                )
    expected = [
        ("14 + 10*sqrt(2)", [28, 16, 20]),
        ("21 + 15*sqrt(2)", [60, 36, 45]),
        ("42 + 30*sqrt(2)", [228, 144, 180]),
        ("56 + 40*sqrt(2)", [100, 64, 80]),
    ]
    if [(row["height"]["display"], row["gram"]) for row in nonnegative] != expected:
        raise AssertionError("the replacement-ray relaxed frontier changed")
    replacement_norm = Silver(4) * Silver(REPLACEMENT_SCALE) ** 2 - LAMBDA22_COEFFICIENT
    cutter_pairing = (
        -LAMBDA22_COEFFICIENT * CONTROL
        + Silver(REPLACEMENT_SCALE) * lambda22_pair_silver(D, SPATIAL)
    )
    alternate_h164_margin = (
        Silver(REPLACEMENT_SCALE) * Silver(36, 24)
        - LAMBDA22_COEFFICIENT * CONTROL
    )
    if (
        retained != 112
        or replacement_norm != Silver(569) / 400 - Silver(0, 1)
        or replacement_norm.sign(0) <= 0
        or cutter_pairing != Silver(-19, 9) / 20
        or cutter_pairing.sign(0) >= 0
        or alternate_h164_margin != Silver(-3, -2) / 10
        or alternate_h164_margin.sign(0) >= 0
    ):
        raise AssertionError(
            "the replacement-ray certificate changed: "
            f"retained={retained}, norm={replacement_norm}, "
            f"cutter={cutter_pairing}, alternate={alternate_h164_margin}"
        )
    return {
        "vector": "z'=(1,(13/40)d)",
        "scale": "13/40",
        "norm": value_json(replacement_norm),
        "control_count": len(enumerate_lambda22_controls(HEIGHT)),
        "retained_necessary_type_count": retained,
        "relaxed_nonnegative_rows": nonnegative,
        "resolution": (
            "The H14 q4 row is controlled by the complete frozen H14 shell; "
            "the H21, H42, and H56 rows are unrealizable by the S3/S6/S2 "
            "difference-minimum exclusions frozen in the H150 audit."
        ),
        "explicit_cutter_pairing": value_json(cutter_pairing),
        "alternate_H164_relaxed_margin": value_json(alternate_h164_margin),
    }


def build_payload() -> dict[str, object]:
    h150_sha = sha256(H150_JSON_PATH.read_bytes()).hexdigest()
    chronological = chronological_cross_scan()
    explicit = explicit_root_audit()
    first_packet = first_packet_acuteness_audit()
    replacement = replacement_ray_audit()
    return {
        "schema_version": "0.1.0",
        "seed_id": "silver-H22-Lambda22-a-sqrt2-1",
        "scope": (
            "Exact explicit accepted cutter and strict replacement hyperideal ray "
            "through H=164+116sqrt(2). This is a bounded prefix result, not "
            "Vinberg termination or global nonreflectivity."
        ),
        "base_H150_artifact": {
            "path": str(H150_JSON_PATH.relative_to(ROOT)),
            "sha256": h150_sha,
        },
        "explicit_cutter": explicit,
        "complete_H14_packet_acuteness": first_packet,
        "chronological_acceptance": {
            "method": (
                "Complete necessary binary Gram census with exact 4x4 PSD Schur "
                "extensions by d,e, full pairing gates and S2 mod-8 residues."
            ),
            "nonnegative_row_count": len(chronological),
            "pre_H164_tangent_row_count": sum(
                row["height"]["display"] != str(HEIGHT) for row in chronological
            ),
            "H164_q2_tangent_row_count": sum(
                row["height"]["display"] == str(HEIGHT)
                and row["norm"]["display"] == "2"
                for row in chronological
            ),
            "H164_positive_row_count": 1,
            "rows": chronological,
            "conclusion": (
                "No earlier possible root pairs positively with the cutter. The "
                "unique positive same-height extension is the cutter itself."
            ),
        },
        "replacement_hyperideal_ray": replacement,
        "conclusion": (
            "The H164 root genuinely cuts the old d/3 ray, but the replacement "
            "13d/40 ray proves the complete bounded prefix still has infinite volume."
        ),
        "discovery": None,
    }


def report(payload: dict[str, object]) -> str:
    return """# Silver Lambda22: first accepted cutter at H164, but no closure

For the arithmetic host

`<-(sqrt(2)-1)> + Lambda22`

in hyperbolic dimension 22, the first explicit root found that cuts the
frozen ray `z=(1,d/3)` occurs at

`H=164+116sqrt(2)`.

Take `d=e1` and

`e=(-6,31,21,15,10,20,14,7,4,2,2,1,0^10)`.

Then `d^2=4`, `e^2=12`, `d.e=-2`, and

`u=8d-e`, `v=6d-e`, `r=(28+20sqrt(2),u+sqrt(2)v)`

has `r^2=4`.  It is primitive and reflection-integral.  Its pairing with
the old ray is `(-2+2sqrt(2))/3>0`, so it is a genuine cutter.

## Chronological acceptance

An exact necessary-superset scan uses every reflective control through H164,
every binary coefficient Gram type surviving the local gates and minima,
and every positive-semidefinite 4-by-4 extension by the fixed pair `(d,e)`.
There are nine pre-H164 nonnegative rows, all tangent, and no positive row.
At H164 there are two q=2 tangent rows and one positive q=4 extension.  The
positive extension has zero Schur complement and is exactly the displayed
root.  Hence no earlier wall rejects it.

## Replacement ray

The cutter still does not close the chamber.  The exact ray

`z'=(1,(13/40)d)`

has norm `569/400-sqrt(2)>0`.  A complete gated necessary-superset scan of
91 controls and 112 gate-compatible binary types through H164 leaves only four relaxed
nonnegative rows; all four are already excluded by the frozen H14/H150
shell and minimum arguments.  The explicit cutter pairs with `z'` by
`-19/20+9sqrt(2)/20<0`.

Thus the exact prefix through H164 still has infinite volume.  This is a
bounded-prefix result, not a global nonreflectivity theorem.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h164_cutter.py
PYTHONPATH=src pytest -q tests/test_silver_lambda22_h164_cutter.py
```
"""


def main() -> None:
    payload = build_payload()
    JSON_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    REPORT_PATH.write_text(report(payload))
    print(JSON_PATH)
    print(REPORT_PATH)


if __name__ == "__main__":
    main()

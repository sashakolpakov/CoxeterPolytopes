#!/usr/bin/env python3
"""Independent exact audit of the first ray-cutting silver Lambda22 wall.

The host is ``<-(sqrt(2)-1)> + Lambda22``.  This script deliberately does
not import an H164 owner artifact.  Starting from the catalogue Gram matrix,
it checks the explicit root and replays a complete necessary-superset
chronology through its height.  It also certifies a replacement hyperideal
ray against an even larger necessary superset.
"""

from __future__ import annotations

import json
from functools import lru_cache
from hashlib import sha256
from math import gcd, isqrt, lcm
from pathlib import Path

from coxeter_polytopes.exact.silver import SILVER_ZERO, SQRT2, Silver
from coxeter_polytopes.search.silver_lambda22 import (
    FIRST_PACKET_HEIGHT,
    GRAM_FILE_SHA256,
    LAMBDA22_COEFFICIENT,
    SOURCE_URL,
    enumerate_lambda22_controls,
    lambda22_gram,
    lambda22_lorentz_pair_silver,
    lambda22_pair_silver,
)
from coxeter_polytopes.search.silver_leech import candidate_leech_gram_triples
from coxeter_polytopes.search.silver_leech_manifest import (
    has_short_binary_combination,
)


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = (
    ROOT / "data/enumerations/silver_lambda22_h164_independent_audit.json"
)
REPORT_PATH = ROOT / "results/silver_lambda22_h164_independent_audit.md"
BASE_H14_PATH = (
    ROOT / "data/enumerations/silver_lambda22_h14_independent_hostile_audit.json"
)

D = (1,) + (0,) * 21
E = (-6, 31, 21, 15, 10, 20, 14, 7, 4, 2, 2, 1) + (0,) * 10
U0 = tuple(8 * D[index] - E[index] for index in range(22))
V0 = tuple(6 * D[index] - E[index] for index in range(22))
X0 = Silver(28, 20)
Q0 = Silver(4)
H0 = Silver(164, 116)
TARGET0 = Silver(660, 464)

MINIMUM_BY_MODULUS = {1: 4, 2: 12, 3: 24, 6: 96}


def value_json(value: Silver) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(0), value.sign(1)],
    }


def gate_data(level: object) -> tuple[int, int]:
    """Pairing-radical moduli for ``u,v`` in ``y=u+sqrt(2)v``."""

    u_modulus = v_modulus = 1
    if level.dyadic_valuation == 3:
        u_modulus = 2
    elif level.dyadic_valuation == 4:
        u_modulus = v_modulus = 2
    if level.ternary_gate:
        u_modulus = lcm(u_modulus, 3)
        v_modulus = lcm(v_modulus, 3)
    return u_modulus, v_modulus


def basic_necessary_types(level: object) -> tuple[tuple[int, int, int], ...]:
    """Every binary Gram type allowed by exact elementary obstructions."""

    u_modulus, v_modulus = gate_data(level)
    cross_modulus = lcm(u_modulus, v_modulus)
    output = []
    for item in candidate_leech_gram_triples(level.target):
        if has_short_binary_combination(item):
            continue
        u_square, v_square, u_dot_v = (
            item.u_square,
            item.v_square,
            item.u_dot_v,
        )
        if (
            (u_square and u_square < MINIMUM_BY_MODULUS[u_modulus])
            or (v_square and v_square < MINIMUM_BY_MODULUS[v_modulus])
            or u_square % u_modulus
            or v_square % v_modulus
            or u_dot_v % cross_modulus
        ):
            continue
        # The four classes of S2/2Lambda22 have norm residues 0,4,4,4.
        if u_modulus % 2 == 0 and u_square % 8 not in (0, 4):
            continue
        if v_modulus % 2 == 0 and v_square % 8 not in (0, 4):
            continue
        output.append((u_square, v_square, u_dot_v))
    return tuple(output)


@lru_cache(maxsize=None)
def extension_pairs(
    vector_square: int, gate_modulus: int
) -> tuple[tuple[int, int, int], ...]:
    """All exact PSD extensions of Gram(d,e) by one gated vector.

    The fixed Gram is ``[[4,-2],[-2,12]]``.  Rows are
    ``(d.vector,e.vector,22*SchurResidual)``.  Since ``e in S2``, its
    pairing also has the lcm of 2 with the vector's gate modulus.
    """

    output = []
    e_modulus = lcm(2, gate_modulus)
    for d_dot in range(-isqrt(4 * vector_square), isqrt(4 * vector_square) + 1):
        if d_dot % gate_modulus:
            continue
        for e_dot in range(
            -isqrt(12 * vector_square), isqrt(12 * vector_square) + 1
        ):
            if e_dot % e_modulus:
                continue
            residual = 22 * vector_square - 2 * (
                3 * d_dot * d_dot + d_dot * e_dot + e_dot * e_dot
            )
            if residual >= 0:
                output.append((d_dot, e_dot, residual))
    return tuple(output)


def nonnegative_cross_extensions(
    level: object, u_square: int, v_square: int, u_dot_v: int
) -> tuple[dict[str, object], ...]:
    """Exact 4x4 PSD extensions pairing nonnegatively with the H164 root."""

    u_modulus, v_modulus = gate_data(level)
    left = extension_pairs(u_square, u_modulus)
    right = tuple(
        sorted(
            extension_pairs(v_square, v_modulus),
            key=lambda row: Silver(
                12 * row[0] - 2 * row[1], 8 * row[0] - row[1]
            ),
            reverse=True,
        )
    )
    threshold = LAMBDA22_COEFFICIENT * X0 * level.control
    output = []
    for d_dot_u, e_dot_u, residual_u in left:
        first = Silver(8 * d_dot_u - e_dot_u, 6 * d_dot_u - e_dot_u)
        if right:
            best_right = Silver(
                12 * right[0][0] - 2 * right[0][1],
                8 * right[0][0] - right[0][1],
            )
            if first + best_right < threshold:
                continue
        for d_dot_v, e_dot_v, residual_v in right:
            spatial = first + Silver(
                12 * d_dot_v - 2 * e_dot_v,
                8 * d_dot_v - e_dot_v,
            )
            margin = spatial - threshold
            if margin.sign(0) < 0:
                break
            cross_residual = 22 * u_dot_v - (
                6 * d_dot_u * d_dot_v
                + d_dot_u * e_dot_v
                + e_dot_u * d_dot_v
                + 2 * e_dot_u * e_dot_v
            )
            if residual_u * residual_v < cross_residual * cross_residual:
                continue
            output.append(
                {
                    "cross": [d_dot_u, d_dot_v, e_dot_u, e_dot_v],
                    "schur_scaled22": [
                        residual_u,
                        residual_v,
                        cross_residual,
                    ],
                    "pairing": value_json(margin),
                }
            )
    return tuple(output)


def h14_packet_maximum(shell_norm: int, d_modulus: int) -> dict[str, object]:
    """Maximize the candidate pairing over a complete H14 shell.

    If ``w`` is the underlying shell vector, put ``p=d.w,r=e.w``.  The
    candidate's spatial pairing is
    ``(20p-3r)+(14p-2r)sqrt(2)``.  Cauchy, ``e in S2``, and exact PSD of
    Gram(d,e,w) give a finite superset.  For norm ten, minimum four removes
    ``|p|=6`` because one of ``w+d,w-d`` would have norm two.
    """

    if (shell_norm, d_modulus) not in ((10, 1), (20, 2)):
        raise ValueError("this audit has exactly the two H14 shells")
    best: tuple[Silver, int, int, int] | None = None
    d_bound = isqrt(4 * shell_norm)
    if shell_norm == 10:
        d_bound = 5
    for d_dot in range(-d_bound, d_bound + 1):
        if d_dot % d_modulus:
            continue
        for e_dot in range(-isqrt(12 * shell_norm), isqrt(12 * shell_norm) + 1):
            if e_dot % 2:
                continue
            determinant = (
                44 * shell_norm
                - 12 * d_dot * d_dot
                - 4 * d_dot * e_dot
                - 4 * e_dot * e_dot
            )
            if determinant < 0:
                continue
            spatial = Silver(
                20 * d_dot - 3 * e_dot,
                14 * d_dot - 2 * e_dot,
            )
            row = (spatial, d_dot, e_dot, determinant)
            if best is None or spatial > best[0]:
                best = row
    if best is None:
        raise AssertionError("the H14 shell bound was empty")
    control = Silver(6, 4) if shell_norm == 10 else Silver(8, 6)
    threshold = LAMBDA22_COEFFICIENT * X0 * control
    return {
        "shell_norm": shell_norm,
        "d_pairing_modulus": d_modulus,
        "maximizer": {"d_dot": best[1], "e_dot": best[2]},
        "gram_determinant": best[3],
        "spatial_maximum": value_json(best[0]),
        "lorentz_threshold": value_json(threshold),
        "pairing_maximum": value_json(best[0] - threshold),
    }


def chronology_audit() -> dict[str, object]:
    controls = enumerate_lambda22_controls(H0)
    type_count = 0
    earlier_zero = []
    earlier_positive = []
    simultaneous_zero = []
    simultaneous_positive = []
    for level in controls:
        if level.height <= FIRST_PACKET_HEIGHT:
            continue
        for gram in basic_necessary_types(level):
            type_count += 1
            extensions = nonnegative_cross_extensions(level, *gram)
            for extension in extensions:
                row = {
                    "height": value_json(level.height),
                    "norm": value_json(level.norm),
                    "control": value_json(level.control),
                    "gram": list(gram),
                    **extension,
                }
                positive = extension["pairing"]["signs"][0] > 0
                simultaneous = level.height == H0
                if simultaneous and positive:
                    simultaneous_positive.append(row)
                elif simultaneous:
                    simultaneous_zero.append(row)
                elif positive:
                    earlier_positive.append(row)
                else:
                    earlier_zero.append(row)

    if len(controls) != 91 or type_count != 108:
        raise AssertionError(
            f"the H164 census changed: {len(controls)} controls, {type_count} types"
        )
    if earlier_positive or len(earlier_zero) != 9:
        raise AssertionError("the strict pre-H164 cross frontier changed")
    if len(simultaneous_zero) != 2 or len(simultaneous_positive) != 1:
        raise AssertionError("the simultaneous H164 cross frontier changed")
    unique = simultaneous_positive[0]
    if (
        unique["norm"]["display"] != "4"
        or unique["gram"] != [300, 180, 232]
        or unique["cross"] != [34, 26, -28, -24]
        or unique["schur_scaled22"] != [0, 0, 0]
        or unique["pairing"]["display"] != "4"
    ):
        raise AssertionError("the unique positive H164 extension changed")
    return {
        "control_count": len(controls),
        "post_H14_necessary_binary_type_count": type_count,
        "strictly_positive_pre_H164_count": len(earlier_positive),
        "pre_H164_tangency_count": len(earlier_zero),
        "pre_H164_tangencies": earlier_zero,
        "H164_tangent_q2_count": len(simultaneous_zero),
        "H164_tangent_q2_rows": simultaneous_zero,
        "H164_positive_count": len(simultaneous_positive),
        "H164_positive_rows": simultaneous_positive,
        "unique_positive_resolution": (
            "The zero Schur residuals force u and v into span(d,e). Solving "
            "the two exact 2x2 systems gives u=8d-e and v=6d-e, so the "
            "unique positive extension is the candidate paired with itself."
        ),
    }


def modulus_class(value: int) -> int:
    if value % 6 == 0:
        return 6
    if value % 3 == 0:
        return 3
    if value % 2 == 0:
        return 2
    return 1


def rank_one_gate_ok(
    u_square: int,
    v_square: int,
    u_dot_v: int,
    u_modulus: int,
    v_modulus: int,
) -> bool:
    """Apply the exact primitive-direction gate to a rank-one binary Gram."""

    primitive_norm = gcd(gcd(u_square, v_square), abs(u_dot_v))
    if primitive_norm == 0:
        return True
    u_scale = isqrt(u_square // primitive_norm)
    v_scale = isqrt(v_square // primitive_norm)
    if u_dot_v < 0:
        v_scale = -v_scale
    if (
        primitive_norm * u_scale * u_scale != u_square
        or primitive_norm * v_scale * v_scale != v_square
        or primitive_norm * u_scale * v_scale != u_dot_v
        or gcd(abs(u_scale), abs(v_scale)) != 1
    ):
        raise AssertionError("rank-one primitive-direction recovery failed")
    forced = lcm(
        u_modulus // gcd(u_modulus, abs(u_scale)),
        v_modulus // gcd(v_modulus, abs(v_scale)),
    )
    category = modulus_class(forced)
    if primitive_norm < MINIMUM_BY_MODULUS[category]:
        return False
    return category % 2 or primitive_norm % 8 in (0, 4)


def complete_binary_gate_ok(
    u_square: int,
    v_square: int,
    u_dot_v: int,
    u_modulus: int,
    v_modulus: int,
) -> bool:
    """Check all binary combinations against S2/S3/S6 minima and residues."""

    determinant = u_square * v_square - u_dot_v * u_dot_v
    if determinant == 0:
        return rank_one_gate_ok(
            u_square, v_square, u_dot_v, u_modulus, v_modulus
        )
    # A forbidden norm is at most 95.  Completing squares gives the exact
    # finite coefficient bounds below; the extra two endpoints are slack.
    u_bound = isqrt((95 * v_square) // determinant) + 2
    v_bound = isqrt((95 * u_square) // determinant) + 2
    for u_coefficient in range(-u_bound, u_bound + 1):
        for v_coefficient in range(-v_bound, v_bound + 1):
            if not (u_coefficient or v_coefficient):
                continue
            square = (
                u_coefficient * u_coefficient * u_square
                + v_coefficient * v_coefficient * v_square
                + 2 * u_coefficient * v_coefficient * u_dot_v
            )
            forced = gcd(
                abs(u_coefficient) * u_modulus,
                abs(v_coefficient) * v_modulus,
            )
            category = modulus_class(forced)
            if 0 < square < MINIMUM_BY_MODULUS[category]:
                return False
            if category % 2 == 0 and square % 8 not in (0, 4):
                return False
    return True


def replacement_ray_audit() -> dict[str, object]:
    """Bound the replacement ray against every necessary type through H164."""

    retained = []
    for level in enumerate_lambda22_controls(H0):
        u_modulus, v_modulus = gate_data(level)
        for gram in basic_necessary_types(level):
            if not complete_binary_gate_ok(
                *gram, u_modulus, v_modulus
            ):
                continue
            u_square, v_square, u_dot_v = gram
            best: tuple[Silver, int, int, int] | None = None
            for d_dot_u in range(
                -isqrt(4 * u_square), isqrt(4 * u_square) + 1
            ):
                if d_dot_u % u_modulus:
                    continue
                for d_dot_v in range(
                    -isqrt(4 * v_square), isqrt(4 * v_square) + 1
                ):
                    if d_dot_v % v_modulus:
                        continue
                    determinant = (
                        4 * (u_square * v_square - u_dot_v * u_dot_v)
                        - v_square * d_dot_u * d_dot_u
                        - u_square * d_dot_v * d_dot_v
                        + 2 * u_dot_v * d_dot_u * d_dot_v
                    )
                    if determinant < 0:
                        continue
                    # This is 40 times <(1,13d/40),root>.
                    margin40 = Silver(13 * d_dot_u, 13 * d_dot_v) - (
                        40 * LAMBDA22_COEFFICIENT * level.control
                    )
                    row = (margin40, d_dot_u, d_dot_v, determinant)
                    if best is None or margin40 > best[0]:
                        best = row
            if best is not None:
                retained.append(
                    {
                        "height": value_json(level.height),
                        "norm": value_json(level.norm),
                        "control": value_json(level.control),
                        "gram": list(gram),
                        "d_cross": [best[1], best[2]],
                        "gram_determinant": best[3],
                        "pairing_times_40": value_json(best[0]),
                    }
                )
    retained.sort(
        key=lambda row: Silver.from_json(row["pairing_times_40"]["exact"]),
        reverse=True,
    )
    if len(retained) != 62:
        raise AssertionError(f"the replacement-ray type count changed: {len(retained)}")
    maximum = retained[0]
    if (
        maximum["pairing_times_40"]["display"] != "5 - 6*sqrt(2)"
        or maximum["pairing_times_40"]["signs"][0] >= 0
        or maximum["height"]["display"] != "164 + 116*sqrt(2)"
        or maximum["norm"]["display"] != "2"
        or maximum["gram"] != [162, 84, 116]
        or maximum["d_cross"] != [25, 18]
    ):
        raise AssertionError("the replacement-ray maximum changed")

    replacement = (Silver(1), *((Silver(13 * entry) / 40) for entry in D))
    replacement_norm = lambda22_lorentz_pair_silver(replacement, replacement)
    candidate = (X0, *(Silver(U0[index], V0[index]) for index in range(22)))
    candidate_pairing = lambda22_lorentz_pair_silver(replacement, candidate)
    if replacement_norm != Silver(569, 0) / 400 - SQRT2:
        raise AssertionError("the replacement-ray norm changed")
    if replacement_norm.sign(0) <= 0 or candidate_pairing.sign(0) >= 0:
        raise AssertionError("the replacement ray is not strict and hyperideal")

    q2_h14 = -LAMBDA22_COEFFICIENT * Silver(6, 4) + (
        Silver(65, 65) / 40
    )
    q4_h14 = -LAMBDA22_COEFFICIENT * Silver(8, 6) + (
        Silver(104, 104) / 40
    )
    return {
        "vector": "z'=(1,13d/40)",
        "norm": value_json(replacement_norm),
        "pairing_with_H164_candidate": value_json(candidate_pairing),
        "complete_binary_gate_type_count": len(retained),
        "global_necessary_superset_maximum": maximum,
        "global_pairing_maximum": value_json(
            Silver.from_json(maximum["pairing_times_40"]["exact"]) / 40
        ),
        "strict_against_every_necessary_type": True,
        "actual_H14_packet_maxima": {
            "q2_norm10_shell": value_json(q2_h14),
            "q4_S2_norm20_shell": value_json(q4_h14),
        },
    }


def candidate_audit() -> dict[str, object]:
    gram = lambda22_gram()
    if (
        lambda22_pair_silver(D, D),
        lambda22_pair_silver(E, E),
        lambda22_pair_silver(D, E),
    ) != (4, 12, -2):
        raise AssertionError("the fixed d,e Gram changed")
    if any(
        sum(gram[row][column] * E[column] for column in range(22)) % 2
        for row in range(22)
    ):
        raise AssertionError("e is no longer in S2")
    ordinary_gram = (
        int(lambda22_pair_silver(U0, U0).a),
        int(lambda22_pair_silver(V0, V0).a),
        int(lambda22_pair_silver(U0, V0).a),
    )
    if ordinary_gram != (300, 180, 232):
        raise AssertionError("the candidate spatial Gram changed")
    spatial = tuple(Silver(U0[index], V0[index]) for index in range(22))
    root = (X0, *spatial)
    root_norm = lambda22_lorentz_pair_silver(root, root)
    height = LAMBDA22_COEFFICIENT * X0 * X0 / root_norm
    if root_norm != Q0 or height != H0:
        raise AssertionError("the candidate norm or height changed")

    time_reflection_coefficient = 2 * LAMBDA22_COEFFICIENT * X0 / Q0
    spatial_reflection_coefficients = []
    for row in range(22):
        pairing = sum(
            (
                Silver(gram[row][column]) * spatial[column]
                for column in range(22)
                if gram[row][column]
            ),
            start=SILVER_ZERO,
        )
        coefficient = 2 * pairing / Q0
        if not coefficient.is_integral:
            raise AssertionError("the candidate reflection is not integral")
        spatial_reflection_coefficients.append(coefficient)
    if not time_reflection_coefficient.is_integral:
        raise AssertionError("the candidate time reflection is not integral")
    if spatial[11] != Silver(-1, -1) or abs(spatial[11].norm()) != 1:
        raise AssertionError("the advertised primitive unit coordinate changed")

    old_witness = (Silver(1), *((Silver(entry) / 3) for entry in D))
    old_pairing = lambda22_lorentz_pair_silver(old_witness, root)
    if old_pairing != Silver(-2, 2) / 3 or old_pairing.sign(0) <= 0:
        raise AssertionError("the H164 root no longer cuts the old ray")
    return {
        "d": list(D),
        "e": list(E),
        "u": list(U0),
        "v": list(V0),
        "control": value_json(X0),
        "spatial_ordinary_gram": list(ordinary_gram),
        "spatial_square": value_json(TARGET0),
        "root_norm": value_json(root_norm),
        "height": value_json(height),
        "primitive_certificate": (
            "The twelfth catalogue coordinate of y is -(1+sqrt(2)), "
            "a unit of Z[sqrt(2)]."
        ),
        "reflection_integral": True,
        "time_reflection_coefficient": value_json(time_reflection_coefficient),
        "spatial_reflection_coefficients": [
            value_json(value) for value in spatial_reflection_coefficients
        ],
        "old_ray": {
            "vector": "z=(1,d/3)",
            "pairing": value_json(old_pairing),
            "cut": True,
        },
    }


def build_payload() -> dict[str, object]:
    h14_maxima = [h14_packet_maximum(10, 1), h14_packet_maximum(20, 2)]
    expected_h14 = ["-12 - 10*sqrt(2)", "-2 - 4*sqrt(2)"]
    if [row["pairing_maximum"]["display"] for row in h14_maxima] != expected_h14:
        raise AssertionError("the independent H14 cross maxima changed")
    if not BASE_H14_PATH.exists():
        raise FileNotFoundError("the frozen independent H14 audit is missing")
    return {
        "schema_version": "0.1.0",
        "status": "EXACT_INDEPENDENT_PASS",
        "scope": (
            "Independent exact acceptance audit of the explicit H164 silver "
            "Lambda22 root, plus a complete necessary-superset replacement-ray "
            "certificate through the same height. This is a bounded Vinberg "
            "prefix and makes no global termination claim."
        ),
        "field": "Q(sqrt(2))",
        "hyperbolic_dimension": 22,
        "lattice": "<-(sqrt(2)-1)> + Lambda22",
        "lambda22_source": SOURCE_URL,
        "gram_file_sha256": GRAM_FILE_SHA256,
        "base_H14_independent_audit": {
            "path": str(BASE_H14_PATH.relative_to(ROOT)),
            "sha256": sha256(BASE_H14_PATH.read_bytes()).hexdigest(),
        },
        "candidate": candidate_audit(),
        "H14_packet_cross_maxima": h14_maxima,
        "chronology": chronology_audit(),
        "replacement_hyperideal_ray": replacement_ray_audit(),
        "conclusion": (
            "The H164 q=4 root is a genuine chronologically accepted wall and "
            "cuts z=(1,d/3). The replacement z'=(1,13d/40) has positive norm "
            "and pairs strictly negatively with every necessary wall type "
            "through H164, so this enlarged prefix still has infinite volume."
        ),
    }


def report(payload: dict[str, object]) -> str:
    chronology = payload["chronology"]
    replacement = payload["replacement_hyperideal_ray"]
    return f"""# Independent silver Lambda22 H164 wall audit

The explicit vector

`r=(28+20sqrt(2), u+sqrt(2)v)`,

with `d=e1`, the frozen catalogue vector `e`, `u=8d-e`, and `v=6d-e`,
is a primitive reflective root of norm 4 at height `164+116sqrt(2)`.
Its pairing with the old hyperideal ray `z=(1,d/3)` is
`-2/3+2sqrt(2)/3>0`, so it genuinely cuts that ray.

## Chronological acceptance

The independent replay enumerates {chronology['control_count']} exact controls and
{chronology['post_H14_necessary_binary_type_count']} gated necessary binary Gram
types after the complete H14 packet.  Exact 4-by-4 Schur-complement enumeration
finds no strictly positive pairing before H164.  There are
{chronology['pre_H164_tangency_count']} possible earlier tangencies, all harmless.
At H164 there are two possible q=2 tangencies and one positive q=4 extension.
The latter has Gram `(300,180,232)`, cross data `(34,26,-28,-24)`, and zero
Schur residual.  Hence it forces `u=8d-e,v=6d-e`: it is the root paired with
itself, with pairing 4.  Therefore the displayed root is accepted.

Against the two complete H14 shells, the exact pairing maxima are respectively
`-12-10sqrt(2)` and `-2-4sqrt(2)`, both strict.

## Replacement hyperideal ray

The ray `z'=(1,13d/40)` has norm `569/400-sqrt(2)>0`.  After all binary
S2/S3/S6 minimum and residue gates, {replacement['complete_binary_gate_type_count']}
necessary types remain through H164.  Its largest possible pairing on that
superset is `(5-6sqrt(2))/40<0`.  In particular, its pairing with the new
H164 wall is `-19/20+9sqrt(2)/20<0`.  Thus the enlarged bounded prefix still
has infinite volume; this is not a global nonreflectivity result.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/independent_audit_silver_lambda22_h164.py
PYTHONPATH=src pytest -q tests/test_silver_lambda22_h164_independent_audit.py
```

Source: {SOURCE_URL}
"""


def main() -> None:
    payload = build_payload()
    JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    JSON_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    REPORT_PATH.write_text(report(payload))
    print(JSON_PATH)
    print(REPORT_PATH)


if __name__ == "__main__":
    main()

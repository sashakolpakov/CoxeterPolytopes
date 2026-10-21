#!/usr/bin/env python3
"""Certify the silver Lambda22 replacement from H=10000 through H=25000.

The inherited ray is ``z=(1,(28/87)d)``, where ``d`` is a fixed minimal
vector of ``Lambda22``.  This script exhausts the exact reflective-control
alphabet, the binary coefficient-Gram necessary superset, and every integral
PSD extension by ``d`` surviving the explicit ``[-3,3]^3`` minimum-four
gate.  It proves a bounded persistence statement only; it does not assert
Vinberg termination or global nonreflectivity.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import isqrt, lcm
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_lambda22 import (
    LAMBDA22_COEFFICIENT,
    enumerate_lambda22_controls,
)
from coxeter_polytopes.search.silver_leech import candidate_leech_gram_triples
from coxeter_polytopes.search.silver_leech_manifest import has_short_binary_combination


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "data/enumerations/silver_lambda22_h10000_replacement.json"
BASE_AUDIT_PATH = (
    ROOT
    / "data/enumerations/silver_lambda22_h10000_replacement_independent_hostile_audit.json"
)
JSON_PATH = ROOT / "data/enumerations/silver_lambda22_h25000_replacement.json"
REPORT_PATH = ROOT / "results/silver_lambda22_h25000_replacement.md"

START = Silver(10_000)
BOUND = Silver(25_000)
SCALE = Silver(Fraction(28, 87))
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
LOCAL_COEFFICIENTS = tuple(
    (first, second, third)
    for first in range(-3, 4)
    for second in range(-3, 4)
    for third in range(-3, 4)
    if first or second or third
)

EXPECTED_COUNTS = (301, 5_575, 2_663_705, 2)
EXPECTED_THREATS = (
    {
        "height": Silver(6392, 4520),
        "norm": Silver(4),
        "control": Silver(176, 124),
        "target": Silver(25572, 18080),
        "gram": (12772, 6400, 9040),
        "modulus": 2,
        "extension_count": 437,
        "pairings": (226, 160),
        "margin": Silver(Fraction(64, 87), Fraction(-44, 87)),
        "multiple": 40,
        "shift": 56,
        "short_norm": 4,
        "sublattice_minimum": 12,
    },
    {
        "height": Silver(9588, 6780),
        "norm": Silver(6),
        "control": Silver(264, 186),
        "target": Silver(57534, 40680),
        "gram": (28734, 14400, 20340),
        "modulus": 3,
        "extension_count": 327,
        "pairings": (339, 240),
        "margin": Silver(Fraction(32, 29), Fraction(-22, 29)),
        "multiple": 60,
        "shift": 84,
        "short_norm": 6,
        "sublattice_minimum": 24,
    },
)


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


def extension_determinant(
    gram: tuple[int, int, int], d_dot_u: int, d_dot_v: int
) -> int:
    u_square, v_square, dot = gram
    return (
        4 * (u_square * v_square - dot * dot)
        - v_square * d_dot_u * d_dot_u
        - u_square * d_dot_v * d_dot_v
        + 2 * dot * d_dot_u * d_dot_v
    )


def locally_admissible_extension(
    gram: tuple[int, int, int], d_dot_u: int, d_dot_v: int
) -> bool:
    """Reject any displayed norm-one, -two, or -three lattice vector."""

    u_square, v_square, dot = gram
    return not any(
        0
        < (
            first * first * u_square
            + second * second * v_square
            + 4 * third * third
            + 2 * first * second * dot
            + 2 * first * third * d_dot_u
            + 2 * second * third * d_dot_v
        )
        < 4
        for first, second, third in LOCAL_COEFFICIENTS
    )


def ceil_div(numerator: int, denominator: int) -> int:
    return -((-numerator) // denominator)


def extensions(
    gram: tuple[int, int, int],
    u_modulus: int,
    v_modulus: int,
    *,
    apply_local_gate: bool = True,
) -> tuple[tuple[int, int, int], ...]:
    """Enumerate every congruence-gated integral PSD extension by ``d^2=4``.

    For ``U>0``, the determinant condition is the exact ellipse

    ``(U*q-M*p)^2 <= (U*V-M^2)*(4*U-p^2)``.
    """

    u_square, v_square, dot = gram
    base_determinant = u_square * v_square - dot * dot
    if base_determinant < 0:
        return ()
    p_bound = isqrt(4 * u_square)
    q_bound = isqrt(4 * v_square)
    output: list[tuple[int, int, int]] = []
    if not u_square:
        if dot:
            return ()
        q_start = ceil_div(-q_bound, v_modulus) * v_modulus
        for q in range(q_start, q_bound + 1, v_modulus):
            if not apply_local_gate or locally_admissible_extension(gram, 0, q):
                output.append((0, q, 0))
        return tuple(output)
    p_start = ceil_div(-p_bound, u_modulus) * u_modulus
    for p in range(p_start, p_bound + 1, u_modulus):
        radius_square = base_determinant * (4 * u_square - p * p)
        if radius_square < 0:
            continue
        radius = isqrt(radius_square)
        lower = max(-q_bound, ceil_div(dot * p - radius, u_square))
        upper = min(q_bound, (dot * p + radius) // u_square)
        q_start = ceil_div(lower, v_modulus) * v_modulus
        for q in range(q_start, upper + 1, v_modulus):
            determinant = extension_determinant(gram, p, q)
            if determinant < 0:
                raise AssertionError("the exact ellipse admitted a negative determinant")
            if apply_local_gate and not locally_admissible_extension(gram, p, q):
                continue
            output.append((p, q, determinant))
    return tuple(output)


def ray_margin(control: Silver, p: int, q: int) -> Silver:
    return SCALE * Silver(p, q) - LAMBDA22_COEFFICIENT * control


def threat_resolution(expected: dict[str, object]) -> dict[str, object]:
    gram = expected["gram"]
    modulus = expected["modulus"]
    if not isinstance(gram, tuple) or not isinstance(modulus, int):
        raise TypeError("invalid frozen threat")
    psd = extensions(gram, modulus, modulus, apply_local_gate=False)
    local = extensions(gram, modulus, modulus)
    control = expected["control"]
    if not isinstance(control, Silver):
        raise TypeError("invalid frozen threat control")
    nonnegative = tuple(
        (p, q, determinant, ray_margin(control, p, q))
        for p, q, determinant in local
        if ray_margin(control, p, q).sign(0) >= 0
    )
    p, q = expected["pairings"]  # type: ignore[misc]
    margin = expected["margin"]
    if not isinstance(margin, Silver):
        raise TypeError("invalid frozen threat margin")
    if (
        len(psd) != expected["extension_count"]
        or len(local) != expected["extension_count"]
        or nonnegative != ((p, q, 0, margin),)
    ):
        raise AssertionError("a critical post-H10000 extension census changed")
    u_square, v_square, dot = gram
    multiple = int(expected["multiple"])
    shift = int(expected["shift"])
    relation_square = v_square + multiple * multiple * 4 - 2 * multiple * q
    relation_cross = dot - multiple * p
    short_square = u_square - 2 * shift * p + shift * shift * 4
    short_d_dot = p - 4 * shift
    if (
        relation_square,
        relation_cross,
        short_square,
        margin.sign(0),
    ) != (0, 0, expected["short_norm"], 1):
        raise AssertionError("a post-H10000 short-vector reduction changed")
    stream = sha256()
    for first, second, determinant in local:
        stream.update(f"{first},{second},{determinant}\n".encode())
    return {
        "height": value_json(expected["height"]),  # type: ignore[arg-type]
        "norm": value_json(expected["norm"]),  # type: ignore[arg-type]
        "control": value_json(control),
        "target": value_json(expected["target"]),  # type: ignore[arg-type]
        "gram": list(gram),
        "pairing_moduli": [modulus, modulus],
        "integral_PSD_extension_count": len(psd),
        "local_gate_rejection_count": len(psd) - len(local),
        "locally_admissible_extension_count": len(local),
        "locally_admissible_stream_sha256": stream.hexdigest(),
        "zero_determinant_extension_count": sum(row[2] == 0 for row in local),
        "nonnegative_extension_count": 1,
        "unique_nonnegative_extension": {
            "d_pairings": [p, q],
            "extension_determinant": 0,
            "margin": value_json(margin),
        },
        "forced_relation": f"v={multiple}d",
        "forced_relation_norm": relation_square,
        "forced_relation_cross": relation_cross,
        "short_vector": f"e=u-{shift}d",
        "short_vector_norm": short_square,
        "d_dot_short_vector": short_d_dot,
        "gate_sublattice": f"S{modulus}",
        "gate_sublattice_minimum": expected["sublattice_minimum"],
        "resolution": (
            f"The zero Cauchy minor forces v={multiple}d. The vector "
            f"e=u-{shift}d is nonzero of norm {short_square}; because u is "
            f"in S{modulus} and {shift}d is visibly in S{modulus}, e is in "
            f"S{modulus}, contradicting min(S{modulus})="
            f"{expected['sublattice_minimum']}."
        ),
    }


def scan() -> dict[str, object]:
    levels = tuple(
        level
        for level in enumerate_lambda22_controls(BOUND)
        if START < level.height <= BOUND
    )
    controls_by_norm: dict[str, int] = {}
    retained_by_norm: dict[str, int] = {}
    extensions_by_norm: dict[str, int] = {}
    retained = extension_count = 0
    threats: list[dict[str, object]] = []
    type_stream = sha256()
    for index, level in enumerate(levels, start=1):
        label = str(level.norm)
        controls_by_norm[label] = controls_by_norm.get(label, 0) + 1
        u_modulus, v_modulus = NORM_GATES[level.norm]
        for item in candidate_leech_gram_triples(level.target):
            if has_short_binary_combination(item):
                continue
            gram = item.u_square, item.v_square, item.u_dot_v
            if (
                not norm_gate(gram[0], u_modulus)
                or not norm_gate(gram[1], v_modulus)
                or gram[2] % lcm(u_modulus, v_modulus)
            ):
                continue
            rows = extensions(gram, u_modulus, v_modulus)
            if not rows:
                continue
            retained += 1
            extension_count += len(rows)
            retained_by_norm[label] = retained_by_norm.get(label, 0) + 1
            extensions_by_norm[label] = extensions_by_norm.get(label, 0) + len(rows)
            type_stream.update(
                (
                    f"{level.norm.a},{level.norm.b};"
                    f"{level.height.a},{level.height.b};"
                    f"{level.control.a},{level.control.b};"
                    f"{gram[0]},{gram[1]},{gram[2]};{len(rows)}\n"
                ).encode()
            )
            nonnegative = tuple(
                (p, q, determinant, ray_margin(level.control, p, q))
                for p, q, determinant in rows
                if ray_margin(level.control, p, q).sign(0) >= 0
            )
            if nonnegative:
                threats.append(
                    {
                        "height": level.height,
                        "norm": level.norm,
                        "control": level.control,
                        "target": level.target,
                        "gram": gram,
                        "extension_count": len(rows),
                        "nonnegative": nonnegative,
                    }
                )
        if index % 50 == 0:
            print(
                "PROGRESS",
                index,
                len(levels),
                retained,
                extension_count,
                len(threats),
                flush=True,
            )
    if (len(levels), retained, extension_count, len(threats)) != EXPECTED_COUNTS:
        raise AssertionError(
            "the post-H10000 H25000 frontier changed: "
            f"{len(levels)}, {retained}, {extension_count}, {len(threats)}"
        )
    for row, expected in zip(threats, EXPECTED_THREATS, strict=True):
        if any(
            row[key] != expected[key]
            for key in ("height", "norm", "control", "target", "gram", "extension_count")
        ):
            raise AssertionError("the chronological post-H10000 threat changed")
        p, q = expected["pairings"]  # type: ignore[misc]
        if row["nonnegative"] != ((p, q, 0, expected["margin"]),):
            raise AssertionError("the unique nonnegative extension changed")
    return {
        "control_count": len(levels),
        "controls_by_norm": controls_by_norm,
        "first_control_height": value_json(levels[0].height),
        "last_control_height": value_json(levels[-1].height),
        "retained_gate_compatible_type_count": retained,
        "retained_types_by_norm": retained_by_norm,
        "locally_admissible_extension_count_with_multiplicity": extension_count,
        "locally_admissible_extensions_by_norm": extensions_by_norm,
        "retained_type_stream_sha256": type_stream.hexdigest(),
        "nonnegative_type_count": len(threats),
    }


def dependency_audit() -> dict[str, object]:
    base = json.loads(BASE_PATH.read_text())
    hostile = json.loads(BASE_AUDIT_PATH.read_text())
    if (
        base["replacement_ray"] != "(1,(28/87)d), d^2=4"
        or base["replacement_norm"]["signs"][0] != 1
        or hostile["verdict"] != "PASS"
        or hostile["replacement_ray"]["identity_sign"] != 1
    ):
        raise AssertionError("the independently audited H10000 dependency changed")
    return {
        "production_H10000_artifact": str(BASE_PATH.relative_to(ROOT)),
        "production_H10000_sha256": sha256(BASE_PATH.read_bytes()).hexdigest(),
        "independent_H10000_audit": str(BASE_AUDIT_PATH.relative_to(ROOT)),
        "independent_H10000_audit_sha256": sha256(
            BASE_AUDIT_PATH.read_bytes()
        ).hexdigest(),
        "independent_H10000_verdict": hostile["verdict"],
    }


def build_payload() -> dict[str, object]:
    frontier = scan()
    resolutions = [threat_resolution(expected) for expected in EXPECTED_THREATS]
    ray_norm = 4 * SCALE * SCALE - LAMBDA22_COEFFICIENT
    if ray_norm != Silver(Fraction(10705, 7569), -1) or ray_norm.sign(0) <= 0:
        raise AssertionError("the inherited ray ceased to be hyperideal")
    return {
        "schema_version": "0.1.0",
        "seed_id": "silver-H22-Lambda22-a-sqrt2-1",
        "scope": (
            "Complete gated necessary-wall obstruction strictly after rational "
            "identity height 10000 through inclusive rational identity height "
            "25000; bounded prefix only."
        ),
        "dependency": dependency_audit(),
        "height_start_exclusive": value_json(START),
        "height_bound_inclusive": value_json(BOUND),
        "replacement_ray": "(1,(28/87)d), d^2=4",
        "replacement_norm": value_json(ray_norm),
        "frontier": frontier,
        "nonnegative_threat_resolutions": resolutions,
        "logical_scope": {
            "complete_control_alphabet": True,
            "complete_binary_Gram_necessary_superset": True,
            "all_integral_PSD_d_extensions_enumerated": True,
            "explicit_local_minimum_gate": "[-3,3]^3",
            "bounded_prefix_only": True,
            "vinberg_termination": False,
            "global_nonreflectivity": False,
        },
        "conclusion": (
            "Exactly two necessary types can pair nonnegatively with the ray. "
            "Both are arithmetically empty by the exact S2/S3 minimum gates, "
            "so the strict hyperideal ray survives every possible wall through "
            "rational identity height 25000."
        ),
        "discovery": None,
    }


def report() -> str:
    return """# Silver Lambda22 replacement survives through real height 25000

For `<-(sqrt(2)-1)> + Lambda22` in `H22`, inherit the independently audited
strict hyperideal ray

`z=(1,(28/87)d)`, with `d^2=4`.

The complete interval strictly after rational identity height `10000` and
through inclusive height `25000` contains 301 reflective controls, 5,575
gate-compatible binary coefficient-Gram types, and 2,663,705 locally
admissible integral PSD extensions by `d` (counted with multiplicity).

Exactly two types have any nonnegative extension.

1. At `H=6392+4520sqrt(2)`, norm `q=4`, the Gram type
   `(12772,6400,9040)` has 437 extensions and exactly one nonnegative one:
   `(d.u,d.v)=(226,160)`, of determinant zero and margin
   `(64-44sqrt(2))/87>0`.  Cauchy equality forces `v=40d`; then
   `e=u-56d` has norm four.  The norm-four reflection gate gives `u in S2`,
   while `56d in S2`, so `e in S2`, contradicting `min(S2)=12`.

2. At `H=9588+6780sqrt(2)`, norm `q=6`, the Gram type
   `(28734,14400,20340)` has 327 extensions and exactly one nonnegative one:
   `(339,240)`, again of determinant zero, with margin
   `(32-22sqrt(2))/29>0`.  It forces `v=60d`; now `e=u-84d` has norm six
   and lies in `S3`, contradicting `min(S3)=24`.

There are no further nonnegative necessary types through the cutoff.  Thus
the exact ray remains strict against every possible wall through height
25000, so this bounded chamber prefix still contains an open hyperideal
patch.  This is not Vinberg termination and does not prove global
nonreflectivity.

```sh
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h25000_replacement.py --freeze
PYTHONPATH=src pytest -q tests/test_silver_lambda22_h25000_replacement.py
```
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", action="store_true")
    arguments = parser.parse_args()
    payload = build_payload()
    if arguments.freeze:
        JSON_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        REPORT_PATH.write_text(report())
        print(JSON_PATH)
        print(REPORT_PATH)
    else:
        print(
            json.dumps(
                {
                    "controls": payload["frontier"]["control_count"],
                    "types": payload["frontier"]["retained_gate_compatible_type_count"],
                    "extensions": payload["frontier"][
                        "locally_admissible_extension_count_with_multiplicity"
                    ],
                    "threats": payload["frontier"]["nonnegative_type_count"],
                },
                indent=2,
            )
        )


if __name__ == "__main__":
    main()

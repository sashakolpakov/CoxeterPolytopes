#!/usr/bin/env python3
"""Independent hostile audit of the norm-223 H=114998 Leech root.

The frozen H=250000 necessary-prefix artifact is used only as the list of
arithmetic components to challenge.  This script independently reconstructs
the binary Leech Gram census, its norm-four pairing maxima, the explicit MOG
root, and the two global inequalities proving chronological acceptance.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
from math import gcd, isqrt
from pathlib import Path
import sys
from typing import Sequence


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from coxeter_polytopes.exact.silver import Silver  # noqa: E402


BASE_PATH = ROOT / "data/enumerations/silver_leech_norm223_h250000.json"
PRODUCTION_PATH = (
    ROOT / "data/enumerations/silver_leech_norm223_h114998_limiter.json"
)
EXPECTED_PRODUCTION_SHA256 = (
    "6250394c4c24b129bb3fff6839ba5a020708985b5c96756c903f71a158285c4f"
)
JSON_PATH = (
    ROOT
    / "data/enumerations/"
    "silver_leech_norm223_h114998_independent_hostile_audit.json"
)
REPORT_PATH = (
    ROOT / "results/silver_leech_norm223_h114998_independent_hostile_audit.md"
)

COEFFICIENT = Silver(-13, 14)
LIMITING_HEIGHT = Silver(114_998, 81_316)
LIMITING_CONTROL = Silver(130, 92)
LIMITING_TARGET = Silver(229_998, 162_632)
LIMITING_MAXIMUM = Silver(679, 479)
REALIZED_TRIPLE = (115_266, 57_366, 81_316)
IMPOSSIBLE_TRIPLE = (115_270, 57_364, 81_316)
LOCAL_CLASSES = (
    (Silver(1), Silver(1)),
    (Silver(2), Silver(1)),
    (Silver(15, 1), Silver(15, 1)),
    (Silver(30, 2), Silver(15, 1)),
)

# Raw standard MOG coordinates are divided by sqrt(8).  P is type (4^2),
# Q is all-positive on a standard Golay octad, E=P+Q, and D is type (4^2).
MOG_D = (4, -4) + (0,) * 22
MOG_P = (-4, -4) + (0,) * 22
MOG_Q = (
    0,
    2,
    2,
    0,
    0,
    0,
    2,
    2,
    2,
    0,
    2,
    0,
    2,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2,
)
MOG_E = tuple(left + right for left, right in zip(MOG_P, MOG_Q))
MOG_U = tuple(170 * d + e for d, e in zip(MOG_D, MOG_E))
MOG_V = tuple(120 * d + e for d, e in zip(MOG_D, MOG_E))

_GOLAY23_GENERATOR = sum(
    1 << exponent for exponent in (11, 9, 7, 6, 5, 1, 0)
)


def _rational_string(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def encode(value: object) -> object:
    if isinstance(value, Silver):
        return value.to_json()
    if isinstance(value, Fraction):
        return _rational_string(value)
    if isinstance(value, dict):
        return {str(key): encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(item) for item in value]
    return value


def decode_silver(value: dict[str, str]) -> Silver:
    return Silver(Fraction(value["a"]), Fraction(value["b"]))


def scaled_dot(left: Sequence[int], right: Sequence[int]) -> int:
    numerator = sum(a * b for a, b in zip(left, right))
    if len(left) != 24 or len(right) != 24 or numerator % 8:
        raise AssertionError("invalid scaled MOG pairing")
    return numerator // 8


def golay_words() -> frozenset[int]:
    words = set()
    for message in range(1 << 12):
        word = 0
        for bit in range(12):
            if message & (1 << bit):
                word ^= _GOLAY23_GENERATOR << bit
        word |= (word.bit_count() & 1) << 23
        words.add(word)
    distribution = Counter(word.bit_count() for word in words)
    if len(words) != 4096 or distribution != Counter(
        {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}
    ):
        raise AssertionError("the independent Golay code reconstruction failed")
    return frozenset(words)


def candidate_triples(target: Silver) -> tuple[tuple[int, int, int], ...]:
    """All binary Grams for y=u+sqrt(2)v before the minimum-four gate."""

    if not target.is_integral or target.b.denominator != 1:
        return ()
    rational = int(target.a)
    if rational < 0 or rational % 2 or target.b.denominator != 1:
        return ()
    dot_fraction = target.b / 2
    if dot_fraction.denominator != 1:
        return ()
    dot = int(dot_fraction)
    discriminant = rational * rational - 8 * dot * dot
    if discriminant < 0:
        return ()
    radical_bound = isqrt(discriminant)
    lower = max(0, (rational - radical_bound) // 4 - 2)
    upper = min(rational // 2, (rational + radical_bound) // 4 + 2)
    output = []
    for v_square in range(lower, upper + 1):
        if v_square % 2:
            continue
        u_square = rational - 2 * v_square
        if u_square < 0 or u_square % 2:
            continue
        if (u_square and u_square < 4) or (v_square and v_square < 4):
            continue
        if dot * dot <= u_square * v_square:
            output.append((u_square, v_square, dot))
    return tuple(output)


def combination_square(triple: tuple[int, int, int], a: int, b: int) -> int:
    u_square, v_square, dot = triple
    return a * a * u_square + b * b * v_square + 2 * a * b * dot


def has_short_combination(triple: tuple[int, int, int]) -> bool:
    u_square, v_square, dot = triple
    determinant = u_square * v_square - dot * dot
    if determinant < 0:
        return True
    if determinant == 0:
        return gcd(gcd(u_square, v_square), abs(dot)) < 4
    u_bound = isqrt((3 * v_square) // determinant) + 1
    v_bound = isqrt((3 * u_square) // determinant) + 1
    return any(
        0 < combination_square(triple, a, b) < 4
        for a in range(-u_bound, u_bound + 1)
        for b in range(-v_bound, v_bound + 1)
        if a or b
    )


def reduced_basis(
    triple: tuple[int, int, int],
) -> tuple[int, int, int, tuple[int, int], tuple[int, int]]:
    """Gauss-reduce a positive binary Gram while tracking a GL(2,Z) basis."""

    u_square, v_square, dot = triple
    first = (1, 0)
    second = (0, 1)

    def inner(left: tuple[int, int], right: tuple[int, int]) -> int:
        return (
            left[0] * right[0] * u_square
            + left[1] * right[1] * v_square
            + (left[0] * right[1] + left[1] * right[0]) * dot
        )

    for _ in range(100):
        a = inner(first, first)
        b = inner(second, second)
        c = inner(first, second)
        if a <= 0:
            raise AssertionError("a retained Gram is not positive definite")
        if b < a:
            first, second = second, first
            continue
        if 2 * abs(c) <= a:
            if abs(first[0] * second[1] - first[1] * second[0]) != 1:
                raise AssertionError("binary reduction lost unimodularity")
            return a, b, c, first, second
        quotient = (
            (c + a // 2) // a if c > 0 else -((-c + a // 2) // a)
        )
        second = (
            second[0] - quotient * first[0],
            second[1] - quotient * first[1],
        )
    raise AssertionError("binary reduction did not terminate")


def silver_integer_sign(rational: int, radical: int) -> int:
    """Sign of rational+radical*sqrt(2), with no field-object allocation."""

    if rational == 0:
        return (radical > 0) - (radical < 0)
    if radical == 0:
        return (rational > 0) - (rational < 0)
    if (rational > 0) == (radical > 0):
        return (rational > 0) - (rational < 0)
    norm = rational * rational - 2 * radical * radical
    norm_sign = (norm > 0) - (norm < 0)
    return norm_sign if rational > 0 else -norm_sign


def d4_pairing_data(
    triple: tuple[int, int, int],
) -> tuple[Silver, tuple[tuple[int, int, int], ...]]:
    """Independent exact norm-four pairing bound and its maximizing rows."""

    a, b, c, first, second = reduced_basis(triple)
    determinant = a * b - c * c
    basis_determinant = first[0] * second[1] - first[1] * second[0]
    best: tuple[int, int] | None = None
    maximizers: list[tuple[int, int, int]] = []
    for first_dot in range(-isqrt(4 * a), isqrt(4 * a) + 1):
        if any(
            value != 0 and value < 4
            for value in (4 + a + 2 * first_dot, 4 + a - 2 * first_dot)
        ):
            continue
        for second_dot in range(-isqrt(4 * b), isqrt(4 * b) + 1):
            if any(
                value != 0 and value < 4
                for value in (
                    4 + b + 2 * second_dot,
                    4 + b - 2 * second_dot,
                )
            ):
                continue
            augmented = (
                4 * determinant
                - first_dot * (first_dot * b - second_dot * c)
                + second_dot * (first_dot * c - second_dot * a)
            )
            if augmented < 0:
                continue
            d_dot_u = (
                second[1] * first_dot - first[1] * second_dot
            ) // basis_determinant
            d_dot_v = (
                -second[0] * first_dot + first[0] * second_dot
            ) // basis_determinant
            row = (d_dot_u, d_dot_v, augmented)
            if best is None or silver_integer_sign(
                d_dot_u - best[0], d_dot_v - best[1]
            ) > 0:
                best = (d_dot_u, d_dot_v)
                maximizers = [row]
            elif (d_dot_u, d_dot_v) == best:
                maximizers.append(row)
    if best is None:
        raise AssertionError("no norm-four extension survived")
    return Silver(*best), tuple(sorted(set(maximizers)))


def square_roots_in_ring(value: Silver) -> tuple[Silver, ...]:
    if not value.is_integral or value.a < 0:
        return ()
    rational = int(value.a)
    radical = int(value.b)
    output = set()
    bound = isqrt(rational // 2) + 1
    for n in range(-bound, bound + 1):
        if n == 0:
            if radical:
                continue
            m = isqrt(rational)
            for candidate in {m, -m}:
                if candidate * candidate == rational:
                    output.add(Silver(candidate))
            continue
        if radical % (2 * n):
            continue
        m = radical // (2 * n)
        if m * m + 2 * n * n == rational:
            output.add(Silver(m, n))
    return tuple(sorted(output))


def explicit_root_audit() -> dict[str, object]:
    support = sum(1 << index for index, value in enumerate(MOG_Q) if value)
    if support.bit_count() != 8 or support not in golay_words():
        raise AssertionError("Q is not a standard Golay-octad Leech vector")
    if (
        scaled_dot(MOG_D, MOG_D),
        scaled_dot(MOG_P, MOG_P),
        scaled_dot(MOG_Q, MOG_Q),
        scaled_dot(MOG_P, MOG_Q),
    ) != (4, 4, 4, -1):
        raise AssertionError("the explicit minimal-vector construction failed")
    if (
        scaled_dot(MOG_E, MOG_E),
        scaled_dot(MOG_D, MOG_E),
    ) != (6, -1):
        raise AssertionError("the fixed (d,e) Gram changed")
    gram = (
        scaled_dot(MOG_U, MOG_U),
        scaled_dot(MOG_V, MOG_V),
        scaled_dot(MOG_U, MOG_V),
    )
    direction_pairing = Silver(
        scaled_dot(MOG_D, MOG_U), scaled_dot(MOG_D, MOG_V)
    )
    spatial_square = Silver(gram[0] + 2 * gram[1], 2 * gram[2])
    root_square = spatial_square - COEFFICIENT * LIMITING_CONTROL**2
    height = COEFFICIENT * LIMITING_CONTROL**2 / root_square
    primitive = scaled_dot(MOG_D, MOG_U) % 2 == 1
    reflection_integral = root_square == 2
    if (
        gram != REALIZED_TRIPLE
        or direction_pairing != LIMITING_MAXIMUM
        or spatial_square != LIMITING_TARGET
        or root_square != 2
        or height != LIMITING_HEIGHT
        or not primitive
        or not reflection_integral
    ):
        raise AssertionError("the explicit root audit failed")
    return {
        "d_scaled_sqrt8": MOG_D,
        "e_scaled_sqrt8": MOG_E,
        "u_scaled_sqrt8": MOG_U,
        "v_scaled_sqrt8": MOG_V,
        "d_e_gram": ((4, -1), (-1, 6)),
        "u_v_gram": gram,
        "direction_pairing": direction_pairing,
        "spatial_square": spatial_square,
        "root_square": root_square,
        "height": height,
        "primitive": primitive,
        "primitivity_reason": (
            "A nonunit divisor of a norm-two root could only be sqrt(2). "
            "That would force u into 2*Lambda, but d.u=679 is odd."
        ),
        "reflection_integral": reflection_integral,
    }


def independent_component_scan(
    components: list[dict[str, object]],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    output = []
    total_raw = 0
    total_surviving = 0
    limiting_hits = []
    for component in components:
        target = decode_silver(component["target_after_gate"])
        raw = candidate_triples(target)
        survivors = tuple(triple for triple in raw if not has_short_combination(triple))
        bounds = tuple((triple, *d4_pairing_data(triple)) for triple in survivors)
        maximum = max(row[1] for row in bounds)
        declared = decode_silver(component["maximum_pairing"])
        if (
            len(raw) != component["raw_triple_count"]
            or len(survivors) != component["surviving_triple_count"]
            or maximum > declared
        ):
            raise AssertionError(
                f"independent component replay failed at {component['control_ordinal']}"
            )
        total_raw += len(raw)
        total_surviving += len(survivors)
        output.append(
            {
                "control_ordinal": component["control_ordinal"],
                "height": decode_silver(component["height"]),
                "raw_type_count": len(raw),
                "surviving_type_count": len(survivors),
                "independent_maximum": maximum,
                "frozen_conservative_maximum": declared,
                "independent_tightening": maximum != declared,
            }
        )
        if decode_silver(component["height"]) == LIMITING_HEIGHT:
            limiting_hits = [row for row in bounds if row[1] == LIMITING_MAXIMUM]
    # One further raw type lies in a control whose entire component is deleted
    # by the exhaustive minimum-four gate, so it is absent from the frozen
    # list of 35 nonempty components.
    if total_raw != 2874 or total_surviving != 2855:
        raise AssertionError("the independent binary census totals changed")
    if tuple(row[0] for row in limiting_hits) != (
        IMPOSSIBLE_TRIPLE,
        REALIZED_TRIPLE,
    ):
        raise AssertionError("the limiting maximizing locus changed")
    maximizers = []
    for triple, _maximum, extensions in limiting_hits:
        if len(extensions) != 1 or extensions[0][:2] != (679, 479):
            raise AssertionError("a limiting extension ceased to be unique")
        d_dot_u, d_dot_v, augmented = extensions[0]
        e_square = triple[0] - 2 * 170 * d_dot_u + 170**2 * 4
        f_square = triple[1] - 2 * 120 * d_dot_v + 120**2 * 4
        e_dot_f = (
            triple[2] - 170 * d_dot_v - 120 * d_dot_u + 170 * 120 * 4
        )
        difference_square = e_square + f_square - 2 * e_dot_f
        maximizers.append(
            {
                "triple": triple,
                "d_dot_u": d_dot_u,
                "d_dot_v": d_dot_v,
                "augmented_determinant": augmented,
                "shifted_gram": (e_square, f_square, e_dot_f),
                "shifted_difference_square": difference_square,
                "resolution": (
                    "impossible: a nonzero Leech vector would have norm two"
                    if difference_square == 2
                    else "realized: the zero difference forces the common vector e"
                ),
            }
        )
    if (
        maximizers[0]["shifted_gram"],
        maximizers[0]["shifted_difference_square"],
        maximizers[1]["shifted_gram"],
        maximizers[1]["shifted_difference_square"],
    ) != ((10, 4, 6), 2, (6, 6, 6), 0):
        raise AssertionError("the limiting maximizer resolutions changed")
    return output, {
        "raw_type_count_in_nonempty_components": total_raw,
        "raw_type_count_in_complete_prefix": 2875,
        "short_combination_exclusion_count_in_nonempty_components": (
            total_raw - total_surviving
        ),
        "short_combination_exclusion_count_in_complete_prefix": 20,
        "fully_deleted_raw_type_count": 1,
        "surviving_type_count": total_surviving,
        "independently_tightened_component_count": sum(
            row["independent_tightening"] for row in output
        ),
        "limiting_surviving_type_count": 216,
        "limiting_maximum": LIMITING_MAXIMUM,
        "maximizers": maximizers,
    }


def earlier_envelope(
    components: list[dict[str, object]],
    component_rows: list[dict[str, object]],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    c = LIMITING_MAXIMUM / 4
    beta = Silver(1, 1)
    perpendicular_square = Fraction(23, 4)
    rows = []
    for component in components:
        height = decode_silver(component["height"])
        if height >= LIMITING_HEIGHT:
            continue
        divisor = decode_silver(component["odd_spatial_divisor"])
        spatial_square = decode_silver(component["target_after_gate"]) * divisor**2
        direction_bound = decode_silver(component["maximum_pairing"]) * divisor
        residual = spatial_square - direction_bound**2 / 4
        gap = (
            COEFFICIENT
            * LIMITING_CONTROL
            * decode_silver(component["control"])
            - c * direction_bound
        )
        monotonicity = (
            16 * c**2 * spatial_square
            - direction_bound**2
            * (4 * c**2 + beta**2 * perpendicular_square)
        )
        rejection = gap**2 - beta**2 * perpendicular_square * residual
        if any(
            value.sign(0) <= 0
            for value in (residual, gap, monotonicity, rejection)
        ):
            raise AssertionError(
                f"an earlier envelope failed at {component['control_ordinal']}"
            )
        rows.append(
            {
                "control_ordinal": component["control_ordinal"],
                "height": height,
                "actual_spatial_square": spatial_square,
                "actual_direction_bound": direction_bound,
                "perpendicular_residual": residual,
                "threshold_gap": gap,
                "monotonicity_discriminant": monotonicity,
                "rejection_discriminant": rejection,
            }
        )
    tight = min(rows, key=lambda row: row["rejection_discriminant"])
    if (
        len(rows) != 34
        or sum(row["surviving_type_count"] for row in component_rows[:-1])
        != 2639
        or tight["control_ordinal"] != 565
        or tight["threshold_gap"]
        != Silver(Fraction(63, 4), Fraction(45, 4))
        or tight["perpendicular_residual"] != Silver(Fraction(55, 4), 9)
        or tight["rejection_discriminant"] != Silver(57, 41)
    ):
        raise AssertionError("the tight earlier-wall certificate changed")
    return rows, {
        "decomposition": (
            "Y0=((679+479*sqrt(2))/4)d+(1+sqrt(2))(e+d/4), "
            "with d^2=4 and (e+d/4)^2=23/4"
        ),
        "component_count": len(rows),
        "binary_type_count": 2639,
        "all_strictly_negative": True,
        "tight_control_ordinal": tight["control_ordinal"],
        "tight_threshold_gap": tight["threshold_gap"],
        "tight_perpendicular_residual": tight["perpendicular_residual"],
        "tight_rejection_discriminant": tight["rejection_discriminant"],
    }


def same_height_and_height_zero() -> dict[str, object]:
    height_zero_targets = tuple(norm / divisor**2 for norm, divisor in LOCAL_CLASSES)
    if any(candidate_triples(target) for target in height_zero_targets):
        raise AssertionError("a height-zero necessary Gram type appeared")

    squared_controls = []
    for norm, _divisor in LOCAL_CLASSES:
        control_square = LIMITING_HEIGHT * norm / COEFFICIENT
        positive_roots = tuple(
            root for root in square_roots_in_ring(control_square) if root.sign(0) > 0
        )
        squared_controls.append((norm, control_square, positive_roots))
    if tuple((row[0], row[2]) for row in squared_controls) != (
        (Silver(1), (Silver(92, 65),)),
        (Silver(2), (LIMITING_CONTROL,)),
        (Silver(15, 1), ()),
        (Silver(30, 2), ()),
    ):
        raise AssertionError("the complete same-height control packet changed")
    q1_control = Silver(92, 65)
    q1_target = Silver(1) + COEFFICIENT * q1_control**2
    if candidate_triples(q1_target):
        raise AssertionError("the simultaneous norm-one component is nonempty")

    conjugate_target = LIMITING_TARGET.conjugate()
    conjugate_cap_clearance = Silver(7) - 4 * conjugate_target
    if conjugate_target.sign(0) <= 0 or conjugate_cap_clearance.sign(0) <= 0:
        raise AssertionError("the exact same-height conjugate cap failed")
    # If a distinct same-height root paired positively with the displayed
    # root, Delta^2<4 at the selected place.  At the conjugate place Cauchy
    # gives Delta'^2 <= 4*N' < 7, hence Tr(Delta^2)<11.  Writing
    # Delta=p+sqrt(2)q gives Tr=2(p^2+2q^2).  Leech minimum four forces this
    # to be at least 16 when q!=0.  For q=0 it is 2p^2; below 11 the only
    # nonzero possibility is p^2=4, contradicting Delta^2<4.
    return {
        "height_zero_targets": height_zero_targets,
        "height_zero_empty": True,
        "same_height_control_squares": squared_controls,
        "same_height_q1_control": q1_control,
        "same_height_q1_target": q1_target,
        "same_height_q1_empty": True,
        "same_height_q2_complete_trace_lemma": {
            "conjugate_spatial_target": conjugate_target,
            "exact_clearance_for_4N_conjugate_below_7": conjugate_cap_clearance,
            "trace_upper_bound": "Tr(Delta^2)<11",
            "q_nonzero_trace_lower_bound": 16,
            "q_zero_only_sub_11_case": "p^2=4 and Delta^2=4",
            "positive_distinct_pairing_impossible": True,
        },
        "same_height_acute": True,
    }


@lru_cache(maxsize=1)
def build_payload() -> dict[str, object]:
    base_bytes = BASE_PATH.read_bytes()
    base = json.loads(base_bytes)["audit"]
    components = base["components"]
    if len(components) != 35:
        raise AssertionError("the frozen dependency lost its 35 components")
    independent_rows, census = independent_component_scan(components)
    envelope_rows, envelope = earlier_envelope(components, independent_rows)
    root = explicit_root_audit()
    endpoint = same_height_and_height_zero()
    accepted = (
        root["primitive"]
        and root["reflection_integral"]
        and envelope["all_strictly_negative"]
        and endpoint["height_zero_empty"]
        and endpoint["same_height_q1_empty"]
        and endpoint["same_height_acute"]
    )
    if not accepted:
        raise AssertionError("the independent acceptance verdict failed")
    # This comparison is deliberately last: no production field participates
    # in the hostile replay above.
    production_bytes = PRODUCTION_PATH.read_bytes()
    production_sha = sha256(production_bytes).hexdigest()
    production = json.loads(production_bytes)["audit"]
    production_agrees = (
        production_sha == EXPECTED_PRODUCTION_SHA256
        and production["accepted"] is True
        and production["limiting_height"] == LIMITING_HEIGHT.to_json()
        and production["limiting_control"] == LIMITING_CONTROL.to_json()
        and production["realized_u_v_gram"]
        == {"u_square": 115266, "v_square": 57366, "u_dot_v": 81316}
        and production["earlier_component_count"] == 34
        and production["tight_earlier_control_ordinal"] == 565
        and production["tight_earlier_rejection_discriminant"]
        == Silver(57, 41).to_json()
        and production["same_height_acute"] is True
    )
    if not production_agrees:
        raise AssertionError("the independent verdict disagrees with production")
    return {
        "schema": "coxeter.silver_leech_norm223_h114998.independent_hostile.v1",
        "dependency": {
            "path": str(BASE_PATH.relative_to(ROOT)),
            "sha256": sha256(base_bytes).hexdigest(),
            "use": "frozen list of necessary arithmetic components only",
        },
        "source_sha256": {
            "src/coxeter_polytopes/exact/silver.py": sha256(
                (ROOT / "src/coxeter_polytopes/exact/silver.py").read_bytes()
            ).hexdigest(),
            "scripts/independent_hostile_audit_silver_leech_norm223_h114998.py": (
                sha256(Path(__file__).read_bytes()).hexdigest()
            ),
        },
        "production_cross_check": {
            "path": str(PRODUCTION_PATH.relative_to(ROOT)),
            "sha256": production_sha,
            "expected_sha256": EXPECTED_PRODUCTION_SHA256,
            "consulted_only_after_independent_replay": True,
            "verdict_and_key_exact_fields_agree": production_agrees,
        },
        "independent_binary_replay": {
            **census,
            "component_rows": independent_rows,
        },
        "limiting_root": root,
        "earlier_chronology": {
            **envelope,
            "certificate_rows": envelope_rows,
        },
        "endpoint_chronology": endpoint,
        "verdict": {
            "primitive": True,
            "reflection_integral": True,
            "chronologically_vinberg_accepted": True,
            "new_finite_volume_example": False,
            "new_compact_example": False,
        },
        "scope": (
            "Complete acceptance audit through H=114998+81316*sqrt(2), "
            "relative to the frozen complete necessary arithmetic prefix. "
            "This does not assert termination, finite volume, compactness, "
            "or global nonreflectivity."
        ),
    }


def render(payload: dict[str, object]) -> str:
    chronology = payload["earlier_chronology"]
    census = payload["independent_binary_replay"]
    return f"""# Independent hostile audit: norm-223 H114998 root

The explicit norm-two root at
`H=114998+81316*sqrt(2)` is **chronologically Vinberg-accepted**.

This replay independently rebuilt all `{census['surviving_type_count']}`
surviving binary Leech Gram types in the 35 frozen nonempty components and
conservative norm-four pairing maxima (tightening two frozen bounds). The
endpoint has 216 types. Exactly two
attain `679+479*sqrt(2)`: one is impossible because it creates a Leech
norm-two difference; the other is the explicit MOG realization
`u=170d+e`, `v=120d+e`, with `d^2=4`, `e^2=6`, and `d.e=-1`.

For the `{chronology['binary_type_count']}` types in all
`{chronology['component_count']}` earlier components, a single exact
Schur--Cauchy envelope is strict. Its tight row is control ordinal
`{chronology['tight_control_ordinal']}`:

```text
gap = (63+45*sqrt(2))/4
gap^2 - radicand = 57+41*sqrt(2) > 0.
```

The simultaneous norm-one component is empty. For distinct norm-two roots at
the endpoint, a positive pairing would give `Delta^2<4` at the selected place
and `Delta'^2<7` at the conjugate place, hence `Tr(Delta^2)<11`. But
`Delta=p+sqrt(2)q` has trace `2(p^2+2q^2)`: it is at least 16 if `q!=0`,
while the only nonzero sub-11 case for `q=0` is `p^2=4`, contradicting
`Delta^2<4`. Thus every distinct same-height wall is nonacute.

The root is primitive (`d.u=679` is odd) and its norm two makes reflection
integrality automatic. This upgrades the necessary-prefix limiter to a
genuine accepted wall, but it is not a new finite-volume or compact example.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    payload = encode(build_payload())
    artifact = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
    report = render(payload).encode()
    if args.verify:
        if JSON_PATH.read_bytes() != artifact or REPORT_PATH.read_bytes() != report:
            raise SystemExit("frozen independent norm-223 audit differs")
        print("PASS", sha256(artifact).hexdigest())
        return
    JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    JSON_PATH.write_bytes(artifact)
    REPORT_PATH.write_bytes(report)
    print(sha256(artifact).hexdigest())


if __name__ == "__main__":
    main()

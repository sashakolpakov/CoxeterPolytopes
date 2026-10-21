#!/usr/bin/env python3
"""Certify the accepted H1598 cutter and its replacement hyperideal ray.

The host is ``<-(sqrt(2)-1)> + Lambda22`` in H22.  This continues the
frozen H738 certificate.  The default audit checks a compact, exact record
of the exhaustive scans.  ``--full`` replays all 284 reflective controls,
all necessary coefficient Gram types, and every relevant Schur extension.

This is only a bounded Vinberg-prefix certificate.  It is not a finite-
volume example, Vinberg termination, or a global nonreflectivity theorem.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from math import isqrt, lcm, sqrt
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.cubic7_lambda22 import (
    lambda22_integer_gram_minimum_at_least,
)
from coxeter_polytopes.search.silver_lambda22 import (
    LAMBDA22_COEFFICIENT,
    enumerate_lambda22_controls,
    gated_relaxed_witness_bound,
    lambda22_gram,
    lambda22_pair_silver,
)
from coxeter_polytopes.search.silver_leech import candidate_leech_gram_triples
from coxeter_polytopes.search.silver_leech_manifest import (
    has_short_binary_combination,
)


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/silver_lambda22_h1598_cutter.json"
REPORT_PATH = ROOT / "results/silver_lambda22_h1598_cutter.md"
H738_PATH = ROOT / "data/enumerations/silver_lambda22_h738_cutter.json"

HEIGHT = Silver(1598, 1130)
CONTROL = Silver(62, 44)
ROOT_NORM = Silver(2)
D = (1,) + (0,) * 21
E = (0, 0, -1) + (0,) * 13 + (-1,) + (0,) * 5
U = tuple(20 * D[i] + E[i] for i in range(22))
V = tuple(14 * D[i] + E[i] for i in range(22))
SPATIAL = tuple(Silver(U[i], V[i]) for i in range(22))
OLD_SCALE = Fraction(32, 99)
NEW_SCALE = Fraction(10, 31)

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


# Each extension is (pairing, d.u, d.v, e.u, e.v, Schur_A, Schur_B,
# Schur_N).  For Gram(d,e)=diag(4,10), the Schur numerators are multiplied
# by det=40.  These 32 rows (34 extensions) are the complete nonnegative
# output of the exhaustive --full replay.
PINNED_CHRONOLOGICAL_ROWS = (
    ("21 + 15*sqrt(2)", "6", "12 + 9*sqrt(2)", (60, 36, 45),
     ((Silver(6, 6), 15, 12, 6, 0, 6, 0, 0),
      (Silver(3, 3), 15, 12, 3, 0, 114, 0, 0),
      (Silver(0), 15, 12, 0, 0, 150, 0, 0))),
    ("42 + 30*sqrt(2)", "12", "24 + 18*sqrt(2)", (228, 144, 180),
     ((Silver(0), 30, 24, 0, 0, 120, 0, 0),)),
    ("56 + 40*sqrt(2)", "4", "16 + 12*sqrt(2)", (100, 64, 80),
     ((Silver(0), 20, 16, 0, 0, 0, 0, 0),)),
    ("274 + 194*sqrt(2)", "4", "36 + 26*sqrt(2)", (580, 260, 388),
     ((Silver(2), 48, 32, 6, 6, 16, 16, 16),)),
    ("369 + 261*sqrt(2)", "2", "30 + 21*sqrt(2)", (368, 186, 261),
     ((Silver(0), 38, 27, 8, 6, 24, 6, -12),)),
    ("369 + 261*sqrt(2)", "4", "42 + 30*sqrt(2)", (744, 368, 522),
     ((Silver(0), 54, 38, 12, 8, 24, 24, -24),)),
    ("478 + 338*sqrt(2)", "2", "34 + 24*sqrt(2)", (470, 244, 338),
     ((Silver(0), 43, 31, 8, 6, 54, 6, -2),)),
    ("527 + 373*sqrt(2)", "2", "36 + 25*sqrt(2)", (536, 260, 373),
     ((Silver(0), 46, 32, 8, 6, 24, 16, 8),)),
    ("527 + 373*sqrt(2)", "2", "36 + 25*sqrt(2)", (508, 274, 373),
     ((Silver(0), 45, 33, 4, 4, 6, 6, 6),)),
    ("527 + 373*sqrt(2)", "4", "50 + 36*sqrt(2)", (1040, 536, 746),
     ((Silver(0), 64, 46, 12, 8, 64, 24, 16),)),
    ("601 + 425*sqrt(2)", "2", "38 + 27*sqrt(2)", (584, 310, 425),
     ((Silver(0), 48, 35, 8, 6, 64, 6, 8),)),
    ("601 + 425*sqrt(2)", "4", "54 + 38*sqrt(2)", (1240, 584, 850),
     ((Silver(0), 70, 48, 12, 8, 24, 64, 16),)),
    ("656 + 464*sqrt(2)", "2", "40 + 28*sqrt(2)", (658, 328, 464),
     ((Silver(0), 51, 36, 8, 6, 54, 16, 8),)),
    ("738 + 522*sqrt(2)", "2", "42 + 30*sqrt(2)", (710, 384, 522),
     ((Silver(0), 53, 39, 8, 6, 54, 6, 18),)),
    ("799 + 565*sqrt(2)", "2", "44 + 31*sqrt(2)", (792, 404, 565),
     ((Silver(0), 56, 40, 8, 6, 64, 16, 8),)),
    ("799 + 565*sqrt(2)", "4", "62 + 44*sqrt(2)", (1616, 792, 1130),
     ((Silver(0), 80, 56, 12, 8, 64, 64, 16),)),
    ("862 + 610*sqrt(2)", "2", "46 + 32*sqrt(2)", (878, 424, 610),
     ((Silver(0), 59, 41, 8, 6, 54, 6, 18),)),
    ("956 + 676*sqrt(2)", "2", "48 + 34*sqrt(2)", (938, 488, 676),
     ((Silver(0), 61, 44, 8, 6, 54, 16, 8),)),
    ("956 + 676*sqrt(2)", "4", "68 + 48*sqrt(2)", (1876, 976, 1352),
     ((Silver(0), 86, 62, 16, 12, 56, 24, -8),)),
    ("1025 + 725*sqrt(2)", "2", "50 + 35*sqrt(2)", (1032, 510, 725),
     ((Silver(0), 64, 45, 8, 6, 64, 6, 8),)),
    ("1025 + 725*sqrt(2)", "4", "70 + 50*sqrt(2)", (2040, 1032, 1450),
     ((Silver(0), 90, 64, 12, 8, 24, 64, 16),)),
    ("1127 + 797*sqrt(2)", "2", "52 + 37*sqrt(2)", (1096, 580, 797),
     ((Silver(0), 66, 48, 8, 6, 24, 16, 8),)),
    ("1127 + 797*sqrt(2)", "4", "74 + 52*sqrt(2)", (2320, 1096, 1594),
     ((Silver(0), 96, 66, 12, 8, 64, 24, 16),)),
    ("1202 + 850*sqrt(2)", "2", "54 + 38*sqrt(2)", (1198, 604, 850),
     ((Silver(0), 69, 49, 8, 6, 54, 6, -2),)),
    ("1202 + 850*sqrt(2)", "4", "76 + 54*sqrt(2)", (2332, 1240, 1700),
     ((Silver(0), 96, 70, 16, 12, 96, 24, 32),)),
    ("1312 + 928*sqrt(2)", "4", "80 + 56*sqrt(2)", (2628, 1312, 1856),
     ((Silver(0), 102, 72, 16, 12, 56, 64, 32),)),
    ("1393 + 985*sqrt(2)", "2", "58 + 41*sqrt(2)", (1376, 706, 985),
     ((Silver(0), 74, 53, 8, 6, 24, 6, -12),)),
    ("1393 + 985*sqrt(2)", "4", "82 + 58*sqrt(2)", (2824, 1376, 1970),
     ((Silver(0), 106, 74, 12, 8, 24, 24, -24),)),
    ("1434 + 1014*sqrt(2)", "6", "102 + 72*sqrt(2)", (4218, 2196, 3042),
     ((Silver(0), 129, 93, 24, 18, 6, 54, -18),)),
    ("1598 + 1130*sqrt(2)", "2", "62 + 44*sqrt(2)", (1610, 794, 1130),
     ((Silver(2), 80, 56, 10, 10, 0, 0, 0),)),
    ("1598 + 1130*sqrt(2)", "4", "88 + 62*sqrt(2)", (3164, 1616, 2260),
     ((Silver(0), 112, 80, 16, 12, 96, 64, 32),)),
    ("1598 + 1130*sqrt(2)", "4", "88 + 62*sqrt(2)", (3156, 1620, 2260),
     ((Silver(2), 112, 80, 14, 14, 16, 16, 16),)),
)


# (height, norm, control, Gram, maximizing (d.u,d.v), relaxed margin,
# exact exclusion).  This is the complete nonnegative output for 10d/31.
PINNED_REPLACEMENT_ROWS = (
    ("14 + 10*sqrt(2)", "4", "8 + 6*sqrt(2)", (28, 16, 20), (10, 8),
     Silver(Fraction(-24, 31), Fraction(18, 31)),
     "The complete H14 S2 norm-20 shell has |d.w|<=8; its actual maximum pairing is strict."),
    ("21 + 15*sqrt(2)", "6", "12 + 9*sqrt(2)", (60, 36, 45), (15, 12),
     Silver(Fraction(-36, 31), Fraction(27, 31)),
     "u-v is a nonzero S3 vector of norm 6 below min(S3)=24."),
    ("42 + 30*sqrt(2)", "12", "24 + 18*sqrt(2)", (228, 144, 180), (30, 24),
     Silver(Fraction(-72, 31), Fraction(54, 31)),
     "u-v is a nonzero S6 vector of norm 12 below min(S6)=96."),
    ("56 + 40*sqrt(2)", "4", "16 + 12*sqrt(2)", (100, 64, 80), (20, 16),
     Silver(Fraction(-48, 31), Fraction(36, 31)),
     "u-v is a nonzero S2 vector of norm 4 below min(S2)=12."),
    ("656 + 464*sqrt(2)", "4", "56 + 40*sqrt(2)", (1372, 628, 928), (74, 50),
     Silver(Fraction(-4, 31), Fraction(4, 31)),
     "2u-3v lies in S2 and has norm 4 below min(S2)=12."),
    ("1202 + 850*sqrt(2)", "4", "76 + 54*sqrt(2)", (2500, 1156, 1700), (100, 68),
     Silver(Fraction(8, 31), Fraction(-2, 31)),
     "Schur zero forces u=25d,v=17d; u violates the q4 gate because Gd is not even."),
    ("1393 + 985*sqrt(2)", "2", "58 + 41*sqrt(2)", (1376, 706, 985), (74, 53),
     Silver(Fraction(-4, 31), Fraction(3, 31)),
     "v-u+5d has norm 2 below the ambient minimum 4."),
    ("1393 + 985*sqrt(2)", "4", "82 + 58*sqrt(2)", (2824, 1376, 1970), (106, 74),
     Silver(Fraction(6, 31), Fraction(-4, 31)),
     "u-v+8d lies in S2 and has norm 4 below min(S2)=12."),
    ("1476 + 1044*sqrt(2)", "2", "60 + 42*sqrt(2)", (1490, 732, 1044), (77, 54),
     Silver(Fraction(26, 31), Fraction(-18, 31)),
     "v-u+6d has norm 2 below the ambient minimum 4."),
    ("1598 + 1130*sqrt(2)", "2", "62 + 44*sqrt(2)", (1566, 816, 1130), (79, 57),
     Silver(Fraction(-16, 31), Fraction(12, 31)),
     "v-u+5d has norm 2 below the ambient minimum 4."),
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


def compact_chronological_rows(rows: object = PINNED_CHRONOLOGICAL_ROWS) -> list[dict[str, object]]:
    return [
        {
            "height": height,
            "norm": norm,
            "control": control,
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
        for height, norm, control, gram, extensions in rows
    ]


@lru_cache(maxsize=None)
def one_vector_rows(square: int, modulus: int, kind: int) -> tuple[tuple[float, Silver, int, int, int], ...]:
    rows = []
    for d_dot in range(-isqrt(4 * square), isqrt(4 * square) + 1):
        if d_dot % modulus:
            continue
        for e_dot in range(-isqrt(10 * square), isqrt(10 * square) + 1):
            if e_dot % modulus:
                continue
            schur = 40 * square - 10 * d_dot * d_dot - 4 * e_dot * e_dot
            if schur < 0:
                continue
            value = (
                Silver(20 * d_dot + e_dot, 14 * d_dot + e_dot)
                if kind == 0
                else Silver(28 * d_dot + 2 * e_dot, 20 * d_dot + e_dot)
            )
            numeric = float(value.a) + sqrt(2) * float(value.b)
            rows.append((numeric, value, d_dot, e_dot, schur))
    return tuple(sorted(rows, reverse=True))


def cross_extensions(
    u_square: int,
    v_square: int,
    u_dot_v: int,
    u_modulus: int,
    v_modulus: int,
    other_control: Silver,
) -> tuple[tuple[Silver, int, int, int, int, int, int, int], ...]:
    """All nonnegative 4x4 Schur extensions by the fixed pair (d,e)."""

    threshold = LAMBDA22_COEFFICIENT * CONTROL * other_control
    threshold_numeric = float(threshold.a) + sqrt(2) * float(threshold.b)
    u_rows = one_vector_rows(u_square, u_modulus, 0)
    v_rows = one_vector_rows(v_square, v_modulus, 1)
    output = []
    if not u_rows or not v_rows:
        return ()
    for u_numeric, u_value, d_u, e_u, schur_a in u_rows:
        if u_numeric + v_rows[0][0] < threshold_numeric - 1e-8:
            break
        for v_numeric, v_value, d_v, e_v, schur_b in v_rows:
            if u_numeric + v_numeric < threshold_numeric - 1e-8:
                break
            pairing = u_value + v_value - threshold
            if pairing.sign(0) < 0:
                break
            schur_n = 40 * u_dot_v - 10 * d_u * d_v - 4 * e_u * e_v
            if schur_a * schur_b < schur_n * schur_n:
                continue
            output.append(
                (pairing, d_u, d_v, e_u, e_v, schur_a, schur_b, schur_n)
            )
    return tuple(output)


def full_chronological_scan() -> tuple[int, tuple[object, ...]]:
    """Expensive clean replay of all types through H1598."""

    type_count = 0
    rows = []
    for level in enumerate_lambda22_controls(HEIGHT):
        if level.height <= Silver(14, 10):
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
            type_count += 1
            extensions = cross_extensions(
                *gram, u_modulus, v_modulus, level.control
            )
            if extensions:
                rows.append(
                    (str(level.height), str(level.norm), str(level.control), gram, extensions)
                )
    return type_count, tuple(rows)


def chronological_acceptance(full: bool) -> dict[str, object]:
    if full:
        type_count, rows = full_chronological_scan()
        if type_count != 1134 or rows != PINNED_CHRONOLOGICAL_ROWS:
            raise AssertionError(
                f"the full H1598 chronological replay changed: {type_count}, {rows!r}"
            )
    rows = compact_chronological_rows()
    extensions = [extension for row in rows for extension in row["extensions"]]
    positives = [
        (row, extension)
        for row in rows
        for extension in row["extensions"]
        if extension["pairing"]["signs"][0] > 0
    ]
    expected_positives = [
        ("21 + 15*sqrt(2)", "6", [60, 36, 45]),
        ("21 + 15*sqrt(2)", "6", [60, 36, 45]),
        ("274 + 194*sqrt(2)", "4", [580, 260, 388]),
        ("1598 + 1130*sqrt(2)", "2", [1610, 794, 1130]),
        ("1598 + 1130*sqrt(2)", "4", [3156, 1620, 2260]),
    ]
    if (
        len(rows) != 32
        or len(extensions) != 34
        or [(r["height"], r["norm"], r["gram"]) for r, _ in positives]
        != expected_positives
    ):
        raise AssertionError("the pinned H1598 chronological census is malformed")

    # Exact resolution of the only new same-height positive impostor.
    # From its Gram and pairings, (u-v-8d)^2=0, hence u-v=8d.
    # With h=v-20d, h^2=20 and e.h=14, so (h-e)^2=2.
    q4 = positives[-1][0]
    U0, V0, C0 = q4["gram"]
    if (
        U0 + V0 - 2 * C0 != 256
        or 112 - 80 != 32
        or V0 + 400 * 4 - 40 * 80 != 20
        or 20 + 10 - 2 * 14 != 2
    ):
        raise AssertionError("the H1598 q4 norm-two contradiction changed")
    return {
        "control_count": 284,
        "gate_compatible_post_H14_type_count": 1134,
        "nonnegative_row_count": 32,
        "nonnegative_extension_count": 34,
        "positive_extension_count": 5,
        "method": (
            "Every reflective control, necessary binary Gram type, and exact 4x4 "
            "Schur extension by Gram(d,e)=diag(4,10)."
        ),
        "rows": rows,
        "positive_resolutions": [
            "The two H21 extensions are empty because u-v is in S3 with norm 6<24.",
            "The H274 extension is empty because u-12d is in S2 with norm 4<12.",
            "The H1598 q2 extension is the displayed root itself.",
            "The H1598 q4 extension forces (h-e)^2=2 in Lambda22 and is empty.",
        ],
        "conclusion": (
            "The unique realizable positive extension is the displayed norm-two "
            "root, so it is chronologically accepted."
        ),
    }


def explicit_root() -> dict[str, object]:
    gram = lambda22_gram()
    if len(E) != 22:
        raise AssertionError("the explicit e vector changed rank")
    coefficient_gram = (
        lambda22_pair_silver(U, U),
        lambda22_pair_silver(V, V),
        lambda22_pair_silver(U, V),
    )
    root_norm = (
        -LAMBDA22_COEFFICIENT * CONTROL * CONTROL
        + lambda22_pair_silver(SPATIAL, SPATIAL)
    )
    spatial_pairings = tuple(
        sum(
            (Silver(gram[i][j]) * SPATIAL[j] for j in range(22)),
            start=Silver(),
        )
        for i in range(22)
    )
    reflection_integral = (
        (2 * LAMBDA22_COEFFICIENT * CONTROL / root_norm).is_integral
        and all((2 * pairing / root_norm).is_integral for pairing in spatial_pairings)
    )
    primitive = Silver(-1, -1) in SPATIAL
    old_pairing = (
        -LAMBDA22_COEFFICIENT * CONTROL
        + Silver(OLD_SCALE) * lambda22_pair_silver(D, SPATIAL)
    )
    if (
        (lambda22_pair_silver(D, D), lambda22_pair_silver(E, E), lambda22_pair_silver(D, E))
        != (4, 10, 0)
        or coefficient_gram != (1610, 794, 1130)
        or lambda22_pair_silver(D, SPATIAL) != Silver(80, 56)
        or root_norm != ROOT_NORM
        or LAMBDA22_COEFFICIENT * CONTROL * CONTROL / root_norm != HEIGHT
        or not reflection_integral
        or not primitive
        or old_pairing != Silver(Fraction(-14, 99), Fraction(10, 99))
        or old_pairing.sign(0) <= 0
    ):
        raise AssertionError("the explicit H1598 root arithmetic changed")
    return {
        "d": list(D),
        "e": list(E),
        "Gram_d_e": [[4, 0], [0, 10]],
        "u": list(U),
        "v": list(V),
        "coefficient_gram": [[1610, 1130], [1130, 794]],
        "d_pairing": value_json(Silver(80, 56)),
        "control": value_json(CONTROL),
        "spatial": "y=(20+14sqrt(2))d+(1+sqrt(2))e",
        "norm": value_json(root_norm),
        "height": value_json(HEIGHT),
        "primitive": primitive,
        "reflection_integral": reflection_integral,
        "old_32d_over_99_ray_pairing": value_json(old_pairing),
    }


def packet_and_prior_wall_acuteness() -> dict[str, object]:
    """Complete H14 shell maxima and direct H164/H738 pairings."""

    packets = []
    for name, shell_square, control, modulus in (
        ("ordinary_q2", 10, Silver(6, 4), 1),
        ("dyadic_q4_S2", 20, Silver(8, 6), 2),
    ):
        best = None
        psd = minimum_four = 0
        for d_dot in range(-isqrt(4 * shell_square), isqrt(4 * shell_square) + 1):
            if d_dot % modulus:
                continue
            for e_dot in range(-isqrt(10 * shell_square), isqrt(10 * shell_square) + 1):
                if e_dot % modulus or 40 * shell_square - 10 * d_dot * d_dot - 4 * e_dot * e_dot < 0:
                    continue
                psd += 1
                gram = ((4, 0, d_dot), (0, 10, e_dot), (d_dot, e_dot, shell_square))
                if not lambda22_integer_gram_minimum_at_least(gram, 4):
                    continue
                minimum_four += 1
                spatial = Silver(1, 1) * (
                    Silver(20, 14) * d_dot + Silver(1, 1) * e_dot
                )
                if best is None or spatial > best[0]:
                    best = (spatial, d_dot, e_dot)
        if best is None:
            raise AssertionError("an H14 packet lost every extension")
        margin = best[0] - LAMBDA22_COEFFICIENT * CONTROL * control
        packets.append(
            {
                "packet": name,
                "PSD_extension_count": psd,
                "minimum_four_extension_count": minimum_four,
                "maximizer": {"d_dot_w": best[1], "e_dot_w": best[2]},
                "maximum_pairing": value_json(margin),
            }
        )
    if [
        (row["PSD_extension_count"], row["minimum_four_extension_count"], row["maximizer"], row["maximum_pairing"]["display"])
        for row in packets
    ] != [
        (201, 165, {"d_dot_w": 5, "e_dot_w": 5}, "-45 - 32*sqrt(2)"),
        (103, 89, {"d_dot_w": 8, "e_dot_w": 4}, "-28 - 20*sqrt(2)"),
    ]:
        raise AssertionError("the complete H14 acuteness maxima changed")

    # Direct pairings with the two accepted post-H14 walls.
    e164 = (-6, 31, 21, 15, 10, 20, 14, 7, 4, 2, 2, 1) + (0,) * 10
    g738 = (0, -1) + (0,) * 19 + (-1,)
    h164_spatial = tuple(
        Silver(8 * D[i] - e164[i], 6 * D[i] - e164[i]) for i in range(22)
    )
    h738_spatial = tuple(
        Silver(14 * D[i] + g738[i], 10 * D[i] + g738[i]) for i in range(22)
    )
    h164_pair = -LAMBDA22_COEFFICIENT * CONTROL * Silver(28, 20) + lambda22_pair_silver(SPATIAL, h164_spatial)
    h738_pair = -LAMBDA22_COEFFICIENT * CONTROL * Silver(42, 30) + lambda22_pair_silver(SPATIAL, h738_spatial)
    if h164_pair != Silver(-34, -24) or h738_pair != Silver(-37, -26):
        raise AssertionError("prior accepted-wall acuteness changed")
    return {
        "complete_H14_packets": packets,
        "H164_cutter_pairing": value_json(h164_pair),
        "H738_cutter_pairing": value_json(h738_pair),
        "all_strict": True,
    }


def full_replacement_scan() -> tuple[int, tuple[object, ...]]:
    retained = 0
    rows = []
    scale = Silver(NEW_SCALE)
    for level in enumerate_lambda22_controls(HEIGHT):
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
            try:
                d_u, d_v, maximum = gated_relaxed_witness_bound(
                    *gram, u_modulus, v_modulus
                )
            except AssertionError:
                continue
            retained += 1
            margin = scale * maximum - LAMBDA22_COEFFICIENT * level.control
            if margin.sign(0) >= 0:
                rows.append(
                    (str(level.height), str(level.norm), str(level.control), gram, (d_u, d_v), margin)
                )
    return retained, tuple(rows)


def replacement_ray(full: bool) -> dict[str, object]:
    compact = tuple(row[:6] for row in PINNED_REPLACEMENT_ROWS)
    if full:
        retained, rows = full_replacement_scan()
        if retained != 1138 or rows != compact:
            raise AssertionError(
                f"the full H1598 replacement replay changed: {retained}, {rows!r}"
            )
    scale = Silver(NEW_SCALE)
    norm = Silver(4) * scale * scale - LAMBDA22_COEFFICIENT
    prior = {
        "H164": -LAMBDA22_COEFFICIENT * Silver(28, 20) + scale * Silver(34, 26),
        "H738": -LAMBDA22_COEFFICIENT * Silver(42, 30) + scale * Silver(54, 38),
        "H1598": -LAMBDA22_COEFFICIENT * CONTROL + scale * Silver(80, 56),
    }
    if (
        norm != Silver(Fraction(1361, 961), -1)
        or norm.sign(0) <= 0
        or prior != {
            "H164": Silver(Fraction(-32, 31), Fraction(12, 31)),
            "H738": Silver(Fraction(-18, 31), Fraction(8, 31)),
            "H1598": Silver(Fraction(-6, 31), Fraction(2, 31)),
        }
        or any(value.sign(0) >= 0 for value in prior.values())
    ):
        raise AssertionError("the 10d/31 replacement arithmetic changed")
    rows = [
        {
            "height": height,
            "norm": root_norm,
            "control": control,
            "gram": list(gram),
            "relaxed_pairings": list(pairings),
            "margin": value_json(margin),
            "resolution": resolution,
        }
        for height, root_norm, control, gram, pairings, margin, resolution in PINNED_REPLACEMENT_ROWS
    ]
    if len(rows) != 10 or any(not row["resolution"] for row in rows):
        raise AssertionError("the pinned replacement frontier is malformed")
    return {
        "vector": "z'''=(1,(10/31)d)",
        "scale": "10/31",
        "norm": value_json(norm),
        "control_count": 284,
        "retained_necessary_type_count": 1138,
        "relaxed_nonnegative_row_count": 10,
        "relaxed_nonnegative_rows": rows,
        "accepted_wall_pairings": {name: value_json(value) for name, value in prior.items()},
        "conclusion": "All ten relaxed threats are exactly excluded, so the ray is strict through H1598.",
    }


def build_payload(full: bool = False) -> dict[str, object]:
    return {
        "schema_version": "0.1.0",
        "seed_id": "silver-H22-Lambda22-a-sqrt2-1",
        "scope": (
            "Exact accepted H1598 cutter and strict replacement hyperideal ray "
            "through H1598; bounded prefix only, not finite volume, Vinberg "
            "termination, or global nonreflectivity."
        ),
        "base_H738_artifact": {
            "path": str(H738_PATH.relative_to(ROOT)),
            "sha256": sha256(H738_PATH.read_bytes()).hexdigest(),
        },
        "full_replay_requested": full,
        "explicit_cutter": explicit_root(),
        "complete_acuteness": packet_and_prior_wall_acuteness(),
        "chronological_acceptance": chronological_acceptance(full),
        "replacement_hyperideal_ray": replacement_ray(full),
        "conclusion": (
            "A third genuine accepted wall cuts the 32d/99 ray, but the strict "
            "hyperideal replacement 10d/31 proves that the bounded prefix still "
            "has infinite volume."
        ),
        "discovery": None,
    }


def report() -> str:
    return """# Silver Lambda22: third accepted cutter at H1598

For `<-(sqrt(2)-1)> + Lambda22` in H22, let `d=e1` and
`e=-e3-e17`.  Then `d^2=4`, `e^2=10`, and `d.e=0`.  The root

`r=(62+44sqrt(2),(20d+e)+sqrt(2)(14d+e))`

is primitive, reflection-integral, has norm 2, and occurs at
`H=1598+1130sqrt(2)`.  It pairs positively with the previous `32d/99`
hyperideal ray, so it cuts that ray.

## Chronological acceptance

The exact 4-by-4 scan covers 284 controls and 1,134 gated post-H14
necessary Gram types.  It has 32 nonnegative rows, 34 extensions, and five
positive extensions.  The two H21 positives and the H274 positive are the
previously certified empty S3/S2 rows.  At H1598, the norm-two positive is
the displayed root.  The only other positive, a norm-four row, forces
`u-v=8d`; putting `h=v-20d` gives `h^2=20`, `e.h=14`, and therefore
`(h-e)^2=2`, impossible in Lambda22.  Thus the displayed wall is accepted.

Complete minimum-four extension scans show it is strictly acute to both H14
packets, and direct calculation gives pairings `-34-24sqrt(2)` and
`-37-26sqrt(2)` with the accepted H164 and H738 walls.

## Replacement ray

The ray `z'''=(1,(10/31)d)` has norm `1361/961-sqrt(2)>0` and pairs
strictly negatively with all three accepted cutters.  A complete gated
necessary scan has 1,138 retained types and only ten relaxed nonnegative
rows.  All ten are excluded by exact ambient/S2/S3/S6 minimum arguments,
the complete H14 shell bound, or the q4 reflection gate.  Hence this ray is
strict through H1598.

The bounded prefix therefore still has infinite volume.  This is not a
finite-volume example, Vinberg termination, or a global nonreflectivity
claim.

```sh
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h1598_cutter.py
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h1598_cutter.py --full
PYTHONPATH=src pytest -q tests/test_silver_lambda22_h1598_cutter.py
```
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full", action="store_true", help="replay the expensive exhaustive scans")
    arguments = parser.parse_args()
    payload = build_payload(arguments.full)
    JSON_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    REPORT_PATH.write_text(report())
    print(JSON_PATH)
    print(REPORT_PATH)


if __name__ == "__main__":
    main()

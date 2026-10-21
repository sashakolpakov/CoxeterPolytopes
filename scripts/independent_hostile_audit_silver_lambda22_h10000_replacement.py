#!/usr/bin/env python3
"""Clean-room audit of the silver Lambda22 replacement through H=10000.

No project module is imported.  The replay implements its own exact
``Q(sqrt(2))`` arithmetic, proves a complete control coefficient box,
reconstructs the binary Gram necessary superset and congruence gates, and
enumerates every gated integral positive-semidefinite extension by a fixed
norm-four direction.  Its critical subcertificate separately exhausts all
extensions of ``(3194,1600,2260)`` and the exact ray margins.

This is a bounded necessary-superset audit.  It inherits the complete mirror
norm alphabet and lattice minima from earlier frozen certificates; it is not
Vinberg termination or a global nonreflectivity proof.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from functools import cmp_to_key
from hashlib import sha256
import json
from math import gcd, isqrt, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRODUCTION = ROOT / "data/enumerations/silver_lambda22_h10000_replacement.json"
BASE = ROOT / "data/enumerations/silver_lambda22_h2786_cutter.json"
MINIMA = ROOT / "data/enumerations/silver_lambda22_h14_independent_hostile_audit.json"
OUTPUT = (
    ROOT
    / "data/enumerations/silver_lambda22_h10000_replacement_independent_hostile_audit.json"
)


@dataclass(frozen=True, slots=True)
class K2:
    """Small standalone model of ``a+b*sqrt(2)``."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    @classmethod
    def make(
        cls, a: int | str | Fraction = 0, b: int | str | Fraction = 0
    ) -> "K2":
        return cls(Fraction(a), Fraction(b))

    @classmethod
    def from_exact(cls, value: dict[str, str]) -> "K2":
        return cls.make(value["a"], value["b"])

    def __add__(self, other: "K2") -> "K2":
        return K2(self.a + other.a, self.b + other.b)

    def __sub__(self, other: "K2") -> "K2":
        return K2(self.a - other.a, self.b - other.b)

    def __neg__(self) -> "K2":
        return K2(-self.a, -self.b)

    def __mul__(self, other: "K2") -> "K2":
        return K2(
            self.a * other.a + 2 * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    def __truediv__(self, other: "K2") -> "K2":
        norm = other.a * other.a - 2 * other.b * other.b
        if not norm:
            raise ZeroDivisionError
        return self * K2(other.a / norm, -other.b / norm)

    def sign(self, embedding: int = 0) -> int:
        if embedding not in (0, 1):
            raise ValueError("embedding must be zero or one")
        a, b = (self.a, self.b) if embedding == 0 else (self.a, -self.b)
        sa = (a > 0) - (a < 0)
        sb = (b > 0) - (b < 0)
        if not sa:
            return sb
        if not sb or sa == sb:
            return sa
        return sa if a * a - 2 * b * b > 0 else sb

    @property
    def integral(self) -> bool:
        return self.a.denominator == self.b.denominator == 1

    def exact(self) -> dict[str, str]:
        return {"a": str(self.a), "b": str(self.b)}

    def display(self) -> str:
        if not self.b:
            return str(self.a)
        if not self.a:
            return f"{self.b}*sqrt(2)"
        operation = "+" if self.b > 0 else "-"
        return f"{self.a} {operation} {abs(self.b)}*sqrt(2)"


ZERO = K2.make()
A = K2.make(-1, 1)
START = K2.make(2786, 1970)
BOUND = K2.make(10_000)
SCALE = K2.make(Fraction(28, 87))
COEFFICIENT_BOX = 354
THREAT_GRAM = (3194, 1600, 2260)
THREAT_CONTROL = K2.make(88, 62)
THREAT_HEIGHT = K2.make(3196, 2260)
THREAT_MARGIN = K2.make(Fraction(32, 87), Fraction(-22, 87))

NORMS = (
    K2.make(1),
    K2.make(2),
    K2.make(4, 2),
    K2.make(4),
    K2.make(3),
    K2.make(6),
    K2.make(12, 6),
    K2.make(12),
)
NORM_GATES = {
    K2.make(1): (1, 1),
    K2.make(2): (1, 1),
    K2.make(4, 2): (2, 1),
    K2.make(4): (2, 2),
    K2.make(3): (3, 3),
    K2.make(6): (3, 3),
    K2.make(12, 6): (6, 3),
    K2.make(12): (6, 6),
}
SUBLATTICE_MINIMUM = {1: 4, 2: 12, 3: 24, 6: 96}
LOCAL_COEFFICIENTS = tuple(
    (first, second, third)
    for first in range(-3, 4)
    for second in range(-3, 4)
    for third in range(-3, 4)
    if first or second or third
)


@dataclass(frozen=True, slots=True)
class Control:
    norm: K2
    x: K2
    target: K2
    height: K2


def compare_k2(left: K2, right: K2) -> int:
    return (left - right).sign(0)


def compare_controls(left: Control, right: Control) -> int:
    for first, second in (
        (left.height, right.height),
        (left.norm, right.norm),
        (left.x, right.x),
    ):
        comparison = compare_k2(first, second)
        if comparison:
            return comparison
    return 0


def sign_integer_pair(rational: int, radical: int, embedding: int = 0) -> int:
    """Exact sign of ``rational+radical*sqrt(2)`` without floats."""

    if embedding == 1:
        radical = -radical
    elif embedding != 0:
        raise ValueError("embedding must be zero or one")
    sa = (rational > 0) - (rational < 0)
    sb = (radical > 0) - (radical < 0)
    if not sa:
        return sb
    if not sb or sa == sb:
        return sa
    return sa if rational * rational - 2 * radical * radical > 0 else sb


def coefficient_box_proof() -> dict[str, object]:
    """Prove the integral ``x=m+n*sqrt(2)`` search box is exhaustive."""

    # H<=10000 implies A*x_0^2<=10000*q_0.  The uniform 705 bound is
    # strict for every canonical norm.  Target conjugate positivity gives
    # q_1+A_1*x_1^2>=0, and q+9A<0 at the conjugate place gives |x_1|<3.
    for norm in NORMS:
        if norm.sign(0) <= 0 or norm.sign(1) <= 0:
            raise AssertionError("a mirror norm ceased to be totally positive")
        if (K2.make(705 * 705) * A - K2.make(10_000) * norm).sign(0) <= 0:
            raise AssertionError("identity control bound 705 is insufficient")
        if (norm + K2.make(9) * A).sign(1) >= 0:
            raise AssertionError("conjugate control bound 3 is insufficient")
    # m=(x0+x1)/2 gives |m|<354.  Also |n|sqrt(2)<354 and
    # 251sqrt(2)>354, so |n|<=250.  The symmetric box 354 has slack.
    if (K2.make(0, 251) - K2.make(354)).sign(0) <= 0:
        raise AssertionError("radical coefficient bound was not proved")
    return {
        "identity_embedding_absolute_bound": "<705",
        "conjugate_embedding_absolute_bound": "<3",
        "derived_rational_coefficient_bound": "<=353",
        "derived_radical_coefficient_bound": "<=250",
        "enumerated_symmetric_coefficient_box": COEFFICIENT_BOX,
    }


def enumerate_controls() -> tuple[Control, ...]:
    """Reconstruct every reflective control in the open/closed interval."""

    coefficient_box_proof()
    rows: set[Control] = set()
    for norm in NORMS:
        qa, qb = int(norm.a), int(norm.b)
        qnorm = qa * qa - 2 * qb * qb
        for rational in range(-COEFFICIENT_BOX, COEFFICIENT_BOX + 1):
            for radical in range(-COEFFICIENT_BOX, COEFFICIENT_BOX + 1):
                if sign_integer_pair(rational, radical) <= 0:
                    continue
                # 2*A*x=c+d*sqrt(2); exact division by q is integral iff
                # both quotient numerators are divisible by N(q).
                c = -2 * rational + 4 * radical
                d = 2 * rational - 2 * radical
                if (c * qa - 2 * d * qb) % qnorm:
                    continue
                if (d * qa - c * qb) % qnorm:
                    continue
                # target=q+A*x^2, expanded over ordinary integers.
                target_a = (
                    qa
                    - rational * rational
                    - 2 * radical * radical
                    + 4 * rational * radical
                )
                target_b = (
                    qb
                    + rational * rational
                    + 2 * radical * radical
                    - 2 * rational * radical
                )
                if sign_integer_pair(target_a, target_b, 1) < 0:
                    continue
                height_a = target_a - qa
                height_b = target_b - qb
                start_q_a = 2786 * qa + 3940 * qb
                start_q_b = 2786 * qb + 1970 * qa
                if sign_integer_pair(
                    height_a - start_q_a, height_b - start_q_b
                ) <= 0:
                    continue
                if sign_integer_pair(height_a - 10_000 * qa, height_b - 10_000 * qb) > 0:
                    continue
                rows.add(
                    Control(
                        norm,
                        K2.make(rational, radical),
                        K2.make(target_a, target_b),
                        K2.make(height_a, height_b) / norm,
                    )
                )
    return tuple(sorted(rows, key=cmp_to_key(compare_controls)))


def combination_square(
    gram: tuple[int, int, int], first: int, second: int
) -> int:
    u_square, v_square, dot = gram
    return (
        first * first * u_square
        + second * second * v_square
        + 2 * first * second * dot
    )


def has_forbidden_binary_combination(gram: tuple[int, int, int]) -> bool:
    """Exactly detect every binary integral combination of norm 1--3."""

    u_square, v_square, dot = gram
    determinant = u_square * v_square - dot * dot
    if determinant < 0:
        return True
    if determinant == 0:
        # A rank-one integral form has primitive direction square
        # gcd(U,V,M), attained by a Bezout combination.
        return gcd(gcd(u_square, v_square), abs(dot)) < 4
    first_bound = isqrt((3 * v_square) // determinant) + 1
    second_bound = isqrt((3 * u_square) // determinant) + 1
    return any(
        0 < combination_square(gram, first, second) < 4
        for first in range(-first_bound, first_bound + 1)
        for second in range(-second_bound, second_bound + 1)
        if first or second
    )


def candidate_grams(target: K2) -> tuple[tuple[int, int, int], ...]:
    """Rebuild the even/minimum-four/Cauchy binary necessary census."""

    if not target.integral or target.b.numerator % 2:
        return ()
    rational_target = int(target.a)
    dot = int(target.b / 2)
    rows = []
    for v_square in range(0, rational_target // 2 + 1, 2):
        u_square = rational_target - 2 * v_square
        if u_square < 0 or u_square % 2:
            continue
        if u_square and u_square < 4:
            continue
        if v_square and v_square < 4:
            continue
        gram = (u_square, v_square, dot)
        if dot * dot > u_square * v_square:
            continue
        if has_forbidden_binary_combination(gram):
            continue
        rows.append(gram)
    return tuple(rows)


def norm_gate(square: int, modulus: int) -> bool:
    return (
        (not square or square >= SUBLATTICE_MINIMUM[modulus])
        and square % (2 * modulus) == 0
        and (modulus % 2 or square % 8 in (0, 4))
    )


def extension_determinant(
    gram: tuple[int, int, int], p: int, q: int
) -> int:
    u_square, v_square, dot = gram
    return (
        4 * (u_square * v_square - dot * dot)
        - v_square * p * p
        - u_square * q * q
        + 2 * dot * p * q
    )


def locally_admissible_extension(
    gram: tuple[int, int, int], p: int, q: int
) -> bool:
    """Apply exactly the declared finite ``[-3,3]^3`` short-vector gate."""

    u_square, v_square, dot = gram
    for first, second, third in LOCAL_COEFFICIENTS:
        square = (
            first * first * u_square
            + second * second * v_square
            + 4 * third * third
            + 2 * first * second * dot
            + 2 * first * third * p
            + 2 * second * third * q
        )
        if 0 < square < 4:
            return False
    return True


def ceil_div(numerator: int, denominator: int) -> int:
    return -((-numerator) // denominator)


def extensions(
    gram: tuple[int, int, int],
    p_modulus: int,
    q_modulus: int,
    *,
    apply_local_gate: bool = True,
) -> tuple[tuple[int, int, int], ...]:
    """Enumerate every gated integral PSD extension by ``d^2=4``.

    For ``U>0``, the determinant inequality is evaluated through the exact
    equivalent

    ``(U*q-M*p)^2 <= (U*V-M^2)*(4*U-p^2)``.

    This derives the complete integer q interval for each p without a loose
    rectangular scan.
    """

    if p_modulus <= 0 or q_modulus <= 0:
        raise ValueError("pairing moduli must be positive")
    u_square, v_square, dot = gram
    base_determinant = u_square * v_square - dot * dot
    if base_determinant < 0:
        return ()
    output = []
    p_bound = isqrt(4 * u_square)
    q_bound = isqrt(4 * v_square)
    if not u_square:
        if dot:
            return ()
        for q in range(ceil_div(-q_bound, q_modulus) * q_modulus, q_bound + 1, q_modulus):
            row = (0, q, 0)
            if not apply_local_gate or locally_admissible_extension(gram, 0, q):
                output.append(row)
        return tuple(output)
    p_start = ceil_div(-p_bound, p_modulus) * p_modulus
    for p in range(p_start, p_bound + 1, p_modulus):
        radius_square = base_determinant * (4 * u_square - p * p)
        if radius_square < 0:
            continue
        radius = isqrt(radius_square)
        lower = max(-q_bound, ceil_div(dot * p - radius, u_square))
        upper = min(q_bound, (dot * p + radius) // u_square)
        q_start = ceil_div(lower, q_modulus) * q_modulus
        for q in range(q_start, upper + 1, q_modulus):
            determinant = extension_determinant(gram, p, q)
            if determinant < 0:
                raise AssertionError("exact ellipse interval admitted a negative determinant")
            if apply_local_gate and not locally_admissible_extension(gram, p, q):
                continue
            output.append((p, q, determinant))
    return tuple(output)


def ray_margin(control: K2, p: int, q: int) -> K2:
    return SCALE * K2.make(p, q) - A * control


def threat_extension_theorem() -> dict[str, object]:
    """Exhaust and resolve every integral extension of the sole threat."""

    psd = extensions(THREAT_GRAM, 1, 1, apply_local_gate=False)
    local = extensions(THREAT_GRAM, 1, 1, apply_local_gate=True)
    nonnegative = tuple(
        (p, q, determinant, ray_margin(THREAT_CONTROL, p, q))
        for p, q, determinant in local
        if ray_margin(THREAT_CONTROL, p, q).sign(0) >= 0
    )
    expected = ((113, 80, 0, THREAT_MARGIN),)
    if len(psd) != 665 or len(local) != 665 or nonnegative != expected:
        raise AssertionError(
            "critical extension census changed: "
            f"PSD={len(psd)}, local={len(local)}, nonnegative={nonnegative!r}"
        )
    digest = sha256()
    for p, q, determinant in local:
        digest.update(f"{p},{q},{determinant}\n".encode())
    p, q, determinant, margin = nonnegative[0]
    u_square, v_square, dot = THREAT_GRAM
    v_minus_20d_square = v_square + 20 * 20 * 4 - 40 * q
    u_dot_v_minus_20d = dot - 20 * p
    short_square = u_square - 56 * p + 28 * 28 * 4
    if (
        v_minus_20d_square,
        u_dot_v_minus_20d,
        short_square,
        32 * 32 - 2 * 22 * 22,
    ) != (0, 0, 2, 56):
        raise AssertionError("the exact threat reduction changed")
    return {
        "gram": list(THREAT_GRAM),
        "two_by_two_pairing_bounds": {"d_dot_u": [-113, 113], "d_dot_v": [-80, 80]},
        "exact_ellipse_inequality": (
            "(3194*q-2260*p)^2 <= 2800*(12776-p^2)"
        ),
        "integral_PSD_extension_count": len(psd),
        "explicit_local_gate": "all (a,b,c) in [-3,3]^3 except (0,0,0)",
        "local_gate_rejection_count": len(psd) - len(local),
        "locally_admissible_extension_count": len(local),
        "locally_admissible_stream_sha256": digest.hexdigest(),
        "zero_determinant_extension_count": sum(row[2] == 0 for row in local),
        "nonnegative_extension_count": len(nonnegative),
        "unique_nonnegative_extension": {
            "d_pairings": [p, q],
            "extension_determinant": determinant,
            "margin": margin.exact(),
            "margin_numerator": [32, -22],
            "margin_numerator_norm": 56,
        },
        "forced_relation_check": {
            "norm_v_minus_20d": v_minus_20d_square,
            "u_dot_v_minus_20d": u_dot_v_minus_20d,
            "positive_definiteness_forces": "v=20d",
        },
        "short_vector": {
            "vector": "u-28d",
            "norm": short_square,
            "contradicts_Lambda22_minimum": 4,
        },
    }


def scan_frontier() -> dict[str, object]:
    """Clean-room reconstruction of all H2786--H10000 necessary types."""

    controls = enumerate_controls()
    controls_by_norm: dict[str, int] = {}
    retained_by_norm: dict[str, int] = {}
    retained = 0
    locally_admissible_extensions = 0
    threats = []
    digest = sha256()
    for control in controls:
        label = control.norm.display()
        controls_by_norm[label] = controls_by_norm.get(label, 0) + 1
        u_modulus, v_modulus = NORM_GATES[control.norm]
        for gram in candidate_grams(control.target):
            u_square, v_square, dot = gram
            if (
                not norm_gate(u_square, u_modulus)
                or not norm_gate(v_square, v_modulus)
                or dot % lcm(u_modulus, v_modulus)
            ):
                continue
            rows = extensions(gram, u_modulus, v_modulus)
            if not rows:
                continue
            retained += 1
            retained_by_norm[label] = retained_by_norm.get(label, 0) + 1
            locally_admissible_extensions += len(rows)
            best_p, best_q, best_det = rows[0]
            best_value = K2.make(best_p, best_q)
            for p, q, determinant in rows[1:]:
                value = K2.make(p, q)
                if compare_k2(value, best_value) > 0:
                    best_p, best_q, best_det = p, q, determinant
                    best_value = value
            best_margin = ray_margin(control.x, best_p, best_q)
            digest.update(
                (
                    f"{control.norm.a},{control.norm.b};"
                    f"{control.x.a},{control.x.b};"
                    f"{u_square},{v_square},{dot};"
                    f"{best_p},{best_q},{best_det};{len(rows)}\n"
                ).encode()
            )
            if best_margin.sign(0) < 0:
                continue
            nonnegative = []
            for p, q, determinant in rows:
                margin = ray_margin(control.x, p, q)
                if margin.sign(0) >= 0:
                    nonnegative.append(
                        {
                            "d_pairings": [p, q],
                            "extension_determinant": determinant,
                            "margin": margin.exact(),
                        }
                    )
            threats.append(
                {
                    "height": control.height.exact(),
                    "norm": control.norm.exact(),
                    "control": control.x.exact(),
                    "gram": list(gram),
                    "maximizing_d_pairings": [best_p, best_q],
                    "maximizing_extension_determinant": best_det,
                    "maximizing_margin": best_margin.exact(),
                    "nonnegative_extension_count": len(nonnegative),
                    "nonnegative_extensions": nonnegative,
                }
            )
    if (len(controls), retained, len(threats)) != (129, 1555, 1):
        raise AssertionError(
            f"independent frontier drifted: {len(controls)}, {retained}, {len(threats)}"
        )
    expected_threat = threats[0]
    if (
        K2.from_exact(expected_threat["height"]) != THREAT_HEIGHT
        or K2.from_exact(expected_threat["norm"]) != K2.make(2)
        or K2.from_exact(expected_threat["control"]) != THREAT_CONTROL
        or expected_threat["gram"] != list(THREAT_GRAM)
        or expected_threat["maximizing_d_pairings"] != [113, 80]
        or expected_threat["nonnegative_extension_count"] != 1
    ):
        raise AssertionError("the independently replayed unique threat changed")
    return {
        "control_count": len(controls),
        "controls_by_norm": controls_by_norm,
        "first_control_height": controls[0].height.exact(),
        "last_control_height": controls[-1].height.exact(),
        "retained_type_count": retained,
        "retained_types_by_norm": retained_by_norm,
        "locally_admissible_extension_count_with_multiplicity": locally_admissible_extensions,
        "retained_type_stream_sha256": digest.hexdigest(),
        "relaxed_nonnegative_type_count": len(threats),
        "threats": threats,
    }


def compare_production(frontier: dict[str, object], theorem: dict[str, object]) -> None:
    production = json.loads(PRODUCTION.read_text())
    if production["new_control_count"] != frontier["control_count"]:
        raise AssertionError("production control count disagrees")
    if production["retained_gate_compatible_type_count"] != frontier["retained_type_count"]:
        raise AssertionError("production retained-type count disagrees")
    if production["nonnegative_relaxed_row_count"] != frontier["relaxed_nonnegative_type_count"]:
        raise AssertionError("production relaxed-threat count disagrees")
    threat = production["unique_threat"]
    if threat["gram"] != theorem["gram"]:
        raise AssertionError("production threat Gram disagrees")
    if threat["maximizing_d_pairings"] != [113, 80]:
        raise AssertionError("production threat pairings disagree")
    if threat["locally_admissible_extension_count"] != theorem["locally_admissible_extension_count"]:
        raise AssertionError("production extension count disagrees")
    if threat["nonnegative_extension_count"] != theorem["nonnegative_extension_count"]:
        raise AssertionError("production nonnegative extension count disagrees")
    frozen_extension = threat["all_nonnegative_extensions"][0]
    replay_extension = theorem["unique_nonnegative_extension"]
    if frozen_extension["extension_determinant"] != replay_extension["extension_determinant"]:
        raise AssertionError("production extension determinant disagrees")
    if frozen_extension["margin"]["exact"] != replay_extension["margin"]:
        raise AssertionError("production exact margin disagrees")


def continuity_audit(ray_norm: K2) -> dict[str, object]:
    base = json.loads(BASE.read_text())
    replacement = base["replacement_hyperideal_ray"]
    if replacement["scale"] != "28/87":
        raise AssertionError("the inherited ray scale changed")
    if K2.from_exact(replacement["norm"]["exact"]) != ray_norm:
        raise AssertionError("the inherited ray norm changed")
    if any(
        K2.from_exact(value["exact"]).sign(0) >= 0
        for value in replacement["accepted_wall_pairings"].values()
    ):
        raise AssertionError("the inherited ray is not strict on every accepted cutter")
    minima = json.loads(MINIMA.read_text())
    if minima["verdict"] != "PASS":
        raise AssertionError("the independent lattice-minimum audit is not passing")
    if minima["pari_count_only"]["ordinary_cumulative"]["4"][1] != 4:
        raise AssertionError("the independently audited Lambda22 minimum changed")
    return {
        "base_H2786_sha256": sha256(BASE.read_bytes()).hexdigest(),
        "base_ray_strict_on_all_four_accepted_cutters": True,
        "independent_minima_sha256": sha256(MINIMA.read_bytes()).hexdigest(),
        "Lambda22_minimum": 4,
    }


def build_payload() -> dict[str, object]:
    frontier = scan_frontier()
    theorem = threat_extension_theorem()
    compare_production(frontier, theorem)
    ray_norm = K2.make(1, -1) + K2.make(4) * SCALE * SCALE
    expected_ray_norm = K2.make(Fraction(10705, 7569), -1)
    if ray_norm != expected_ray_norm or ray_norm.sign(0) <= 0 or ray_norm.sign(1) <= 0:
        raise AssertionError("the replacement ray is no longer spacelike")
    continuity = continuity_audit(ray_norm)
    return {
        "schema_version": "0.1.0",
        "verdict": "PASS",
        "scope": (
            "Independent bounded necessary-superset replay strictly after "
            "2786+1970sqrt(2) through inclusive rational identity height 10000."
        ),
        "independence": (
            "No project imports: standalone Q(sqrt(2)), exhaustive control box, "
            "binary Gram census, congruence/minimum gates, exact PSD ellipse "
            "enumeration, finite local short-vector gate, and ray signs."
        ),
        "sources": {
            "production_artifact": str(PRODUCTION.relative_to(ROOT)),
            "production_sha256": sha256(PRODUCTION.read_bytes()).hexdigest(),
        },
        "cutoff_semantics": {
            "lower": "strict at the identity embedding",
            "upper": "inclusive rational field element 10000 at the identity embedding",
            "coefficient_box_proof": coefficient_box_proof(),
        },
        "continuity": continuity,
        "replacement_ray": {
            "vector": "(1,(28/87)d), d^2=4",
            "norm": ray_norm.exact(),
            "identity_sign": ray_norm.sign(0),
        },
        "new_interval": frontier,
        "critical_extension_theorem": theorem,
        "logical_scope": {
            "why_local_gate_is_safe": (
                "Any PSD extension rejected by the explicit [-3,3]^3 gate "
                "already displays an actual nonzero integral Lambda22 combination "
                "of norm 1--3, impossible because min(Lambda22)=4."
            ),
            "why_unique_survivor_is_empty": (
                "Its zero v,d Cauchy minor forces v=20d in the positive-definite "
                "lattice; then u-28d is a nonzero norm-two Lambda22 vector."
            ),
            "bounded_prefix_only": True,
            "vinberg_termination": False,
            "global_nonreflectivity": False,
        },
        "conclusion": (
            "The production H10000 bounded-prefix obstruction passes a clean-room "
            "replay. The sole nonnegative necessary extension is impossible, so the "
            "inherited strict hyperideal ray survives the necessary wall superset "
            "through H=10000; no claim beyond that cutoff is made."
        ),
    }


def encoded_payload() -> bytes:
    return (json.dumps(build_payload(), indent=2, sort_keys=True) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    arguments = parser.parse_args()
    encoded = encoded_payload()
    if arguments.verify:
        if not OUTPUT.exists() or OUTPUT.read_bytes() != encoded:
            raise SystemExit("FAIL: silver H10000 hostile replay differs from frozen JSON")
        status = "PASS"
    else:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_bytes(encoded)
        status = "WROTE"
    payload = json.loads(encoded)
    print(
        json.dumps(
            {
                "status": status,
                "output": str(OUTPUT),
                "sha256": sha256(encoded).hexdigest(),
                "controls": payload["new_interval"]["control_count"],
                "retained_types": payload["new_interval"]["retained_type_count"],
                "critical_extensions": payload["critical_extension_theorem"][
                    "locally_admissible_extension_count"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

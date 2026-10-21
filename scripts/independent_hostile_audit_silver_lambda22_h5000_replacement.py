#!/usr/bin/env python3
"""Clean-room audit of the silver Lambda22 replacement frontier to H=5000.

The production certificate delegates its control and Gram searches to shared
search modules.  This replay deliberately imports no project code.  It uses a
small standalone ``Q(sqrt(2))`` implementation and reconstructs:

* every reflective control strictly after ``1598+1130*sqrt(2)`` and through
  identity-embedding height 5000;
* every binary coefficient Gram type surviving the stated necessary gates;
* every locally admissible integral PSD extension by a norm-four direction;
* all extensions, rather than only the maximizing extension, whose pairing
  with ``(1,(10/31)d)`` is nonnegative.

This is a bounded necessary-superset audit.  It inherits the complete mirror
norm alphabet and the Lambda22/S2 minima from the earlier independently
audited H14 lattice certificate; it is not Vinberg termination.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from functools import cmp_to_key
from hashlib import sha256
from math import gcd, isqrt, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRODUCTION_PATH = ROOT / "data/enumerations/silver_lambda22_h5000_replacement.json"
PREVIOUS_PATH = ROOT / "data/enumerations/silver_lambda22_h1598_cutter.json"
MINIMA_AUDIT_PATH = (
    ROOT / "data/enumerations/silver_lambda22_h14_independent_hostile_audit.json"
)
OUTPUT_PATH = (
    ROOT
    / "data/enumerations/silver_lambda22_h5000_replacement_independent_hostile_audit.json"
)


@dataclass(frozen=True, slots=True)
class K2:
    """A deliberately local exact model of ``a+b*sqrt(2)``."""

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
        operator = "+" if self.b > 0 else "-"
        return f"{self.a} {operator} {abs(self.b)}*sqrt(2)"


ZERO = K2.make()
ONE = K2.make(1)
TWO = K2.make(2)
SQRT2 = K2.make(0, 1)
A = K2.make(-1, 1)
START = K2.make(1598, 1130)
BOUND = K2.make(5000)
SCALE = K2.make(Fraction(10, 31))
COEFFICIENT_BOX = 252
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
    """Exact sign of an integral coefficient pair, without Fraction churn."""

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


def control_coefficient_box_proof() -> dict[str, object]:
    """Prove that the fixed coefficient box contains every control.

    At the identity place, height at most 5000 gives ``|x_0|<500`` for
    every norm class.  Positivity of the spatial target at the conjugate
    place gives ``|x_1|<3``.  Hence ``|m|<=251`` and ``|n|<=177`` for
    ``x=m+n*sqrt(2)``.  The symmetric box 252 is therefore strict slack.
    """

    for norm in NORMS:
        if norm.sign(0) <= 0 or norm.sign(1) <= 0:
            raise AssertionError("mirror norm ceased to be totally positive")
        if (K2.make(250_000) * A - K2.make(5000) * norm).sign(0) <= 0:
            raise AssertionError("identity control bound 500 is too small")
        if (norm + K2.make(9) * A).sign(1) >= 0:
            raise AssertionError("conjugate control bound 3 is too small")
    if (K2.make(0, 178) - K2.make(Fraction(503, 2))).sign(0) <= 0:
        raise AssertionError("radical coefficient bound was not proved")
    return {
        "identity_embedding_absolute_bound": "<500",
        "conjugate_embedding_absolute_bound": "<3",
        "derived_rational_coefficient_bound": "<=251",
        "derived_radical_coefficient_bound": "<=177",
        "enumerated_symmetric_coefficient_box": COEFFICIENT_BOX,
    }


def enumerate_controls() -> tuple[Control, ...]:
    control_coefficient_box_proof()
    rows: set[Control] = set()
    for norm in NORMS:
        qa, qb = int(norm.a), int(norm.b)
        qnorm = qa * qa - 2 * qb * qb
        for rational in range(-COEFFICIENT_BOX, COEFFICIENT_BOX + 1):
            for radical in range(-COEFFICIENT_BOX, COEFFICIENT_BOX + 1):
                if sign_integer_pair(rational, radical) <= 0:
                    continue
                # 2*A*x = c+d*sqrt(2).  Division by qa+qb*sqrt(2)
                # is integral precisely when the following two numerators
                # are divisible by the rational field norm of q.
                c = -2 * rational + 4 * radical
                d = 2 * rational - 2 * radical
                if (c * qa - 2 * d * qb) % qnorm or (d * qa - c * qb) % qnorm:
                    continue
                # target=q+A*x^2, expanded entirely over the integers.
                target_a = qa - rational * rational - 2 * radical * radical + 4 * rational * radical
                target_b = qb + rational * rational + 2 * radical * radical - 2 * rational * radical
                if sign_integer_pair(target_a, target_b, 1) < 0:
                    continue
                height_numerator_a = target_a - qa
                height_numerator_b = target_b - qb
                start_times_q_a = 1598 * qa + 2260 * qb
                start_times_q_b = 1598 * qb + 1130 * qa
                if sign_integer_pair(
                    height_numerator_a - start_times_q_a,
                    height_numerator_b - start_times_q_b,
                ) <= 0:
                    continue
                if sign_integer_pair(
                    height_numerator_a - 5000 * qa,
                    height_numerator_b - 5000 * qb,
                ) > 0:
                    continue
                x = K2.make(rational, radical)
                target = K2.make(target_a, target_b)
                height = K2.make(height_numerator_a, height_numerator_b) / norm
                rows.add(Control(norm, x, target, height))
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
    """Independently test every binary combination that could have norm <4."""

    u_square, v_square, dot = gram
    determinant = u_square * v_square - dot * dot
    if determinant < 0:
        return True
    if determinant == 0:
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
    if not target.integral or target.b.denominator != 1 or target.b.numerator % 2:
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


def locally_admissible_extension(
    gram: tuple[int, int, int], p: int, q: int
) -> bool:
    """Replay the production certificate's explicit [-3,3]^3 minimum gate."""

    u_square, v_square, dot = gram
    return not any(
        0
        < (
            first * first * u_square
            + second * second * v_square
            + 4 * third * third
            + 2 * first * second * dot
            + 2 * first * third * p
            + 2 * second * third * q
        )
        < 4
        for first in range(-3, 4)
        for second in range(-3, 4)
        for third in range(-3, 4)
        if first or second or third
    )


def extensions(
    gram: tuple[int, int, int], u_modulus: int, v_modulus: int
) -> tuple[tuple[int, int, int], ...]:
    """Enumerate all gated integral PSD extensions by ``d^2=4``."""

    u_square, v_square, dot = gram
    output = []
    for p in range(-isqrt(4 * u_square), isqrt(4 * u_square) + 1):
        if p % u_modulus:
            continue
        for q in range(-isqrt(4 * v_square), isqrt(4 * v_square) + 1):
            if q % v_modulus:
                continue
            determinant = (
                4 * (u_square * v_square - dot * dot)
                - v_square * p * p
                - u_square * q * q
                + 2 * dot * p * q
            )
            if determinant < 0 or not locally_admissible_extension(gram, p, q):
                continue
            output.append((p, q, determinant))
    return tuple(output)


def short_resolution(
    norm: K2, gram: tuple[int, int, int], p: int, q: int
) -> dict[str, object]:
    """Find the full-lattice contradiction for a nonnegative extension."""

    u_square, v_square, dot = gram
    candidates = []
    for k in range(-50, 51):
        square = (
            u_square
            + v_square
            - 2 * dot
            + 4 * k * k
            - 2 * k * (p - q)
        )
        if norm == K2.make(2) and 0 < square < 4:
            candidates.append((k, square, "Lambda22", 4))
        if norm == K2.make(4) and k % 2 == 0 and 0 < square < 12:
            candidates.append((k, square, "S2", 12))
    if not candidates:
        raise AssertionError("a nonnegative extension has no certified short vector")
    k, square, lattice, minimum = min(candidates, key=lambda row: (row[1], abs(row[0])))
    return {
        "combination": f"u-v-({k})d",
        "k": k,
        "square": square,
        "lattice": lattice,
        "certified_minimum": minimum,
    }


def scan() -> dict[str, object]:
    controls = enumerate_controls()
    retained = 0
    threats = []
    retained_by_norm: dict[str, int] = {}
    controls_by_norm: dict[str, int] = {}
    for control in controls:
        norm_label = control.norm.display()
        controls_by_norm[norm_label] = controls_by_norm.get(norm_label, 0) + 1
        u_modulus, v_modulus = NORM_GATES[control.norm]
        for gram in candidate_grams(control.target):
            u_square, v_square, dot = gram
            if (
                not norm_gate(u_square, u_modulus)
                or not norm_gate(v_square, v_modulus)
                or dot % lcm(u_modulus, v_modulus)
            ):
                continue
            surviving_extensions = extensions(gram, u_modulus, v_modulus)
            if not surviving_extensions:
                continue
            retained += 1
            retained_by_norm[norm_label] = retained_by_norm.get(norm_label, 0) + 1
            nonnegative = []
            x_a, x_b = int(control.x.a), int(control.x.b)
            for p, q, determinant in surviving_extensions:
                # 31*margin = 10*(p+q*sqrt(2))-31*A*x.
                margin_a = 10 * p + 31 * x_a - 62 * x_b
                margin_b = 10 * q - 31 * x_a + 31 * x_b
                if sign_integer_pair(margin_a, margin_b) < 0:
                    continue
                margin = K2.make(Fraction(margin_a, 31), Fraction(margin_b, 31))
                nonnegative.append(
                    {
                        "d_dot_u": p,
                        "d_dot_v": q,
                        "extension_determinant": determinant,
                        "margin": margin.exact(),
                        "resolution": short_resolution(control.norm, gram, p, q),
                    }
                )
            if nonnegative:
                threats.append(
                    {
                        "height": control.height.exact(),
                        "norm": control.norm.exact(),
                        "control": control.x.exact(),
                        "gram": list(gram),
                        "nonnegative_extension_count": len(nonnegative),
                        "nonnegative_extensions": nonnegative,
                    }
                )
    if (len(controls), retained, len(threats)) != (69, 635, 7):
        raise AssertionError(
            f"independent frontier drifted: {len(controls)}, {retained}, {len(threats)}"
        )
    if [row["nonnegative_extension_count"] for row in threats] != [1] * 7:
        raise AssertionError("a threat has an unresolved alternative nonnegative extension")
    return {
        "controls": controls,
        "control_count": len(controls),
        "controls_by_norm": controls_by_norm,
        "retained_type_count": retained,
        "retained_types_by_norm": retained_by_norm,
        "threats": threats,
    }


def compare_production(scan_result: dict[str, object]) -> None:
    production = json.loads(PRODUCTION_PATH.read_text())
    if production["new_control_count"] != scan_result["control_count"]:
        raise AssertionError("production control count disagrees")
    if production["retained_gate_compatible_type_count"] != scan_result["retained_type_count"]:
        raise AssertionError("production retained-type count disagrees")
    threats = scan_result["threats"]
    if production["nonnegative_relaxed_row_count"] != len(threats):
        raise AssertionError("production threat count disagrees")
    if len(production["rows"]) != len(threats):
        raise AssertionError("production threat rows disagree")
    for frozen, replay in zip(production["rows"], threats, strict=True):
        extension = replay["nonnegative_extensions"][0]
        comparisons = (
            (frozen["height"]["exact"], replay["height"]),
            (frozen["norm"]["exact"], replay["norm"]),
            (frozen["control"]["exact"], replay["control"]),
            (frozen["margin"]["exact"], extension["margin"]),
        )
        if any(left != right for left, right in comparisons):
            raise AssertionError("production exact threat value disagrees")
        if frozen["gram"] != replay["gram"]:
            raise AssertionError("production threat Gram disagrees")
        if frozen["maximizing_d_pairings"] != [
            extension["d_dot_u"],
            extension["d_dot_v"],
        ]:
            raise AssertionError("production maximizing extension disagrees")


def continuity_audit() -> dict[str, object]:
    previous = json.loads(PREVIOUS_PATH.read_text())
    prior_ray = previous["replacement_hyperideal_ray"]
    if prior_ray["scale"] != "10/31" or "strict through H1598" not in prior_ray["conclusion"]:
        raise AssertionError("the H1598 continuity certificate drifted")
    if K2.from_exact(prior_ray["norm"]["exact"]).sign(0) <= 0:
        raise AssertionError("the inherited ray is no longer hyperideal")
    if any(
        K2.from_exact(value["exact"]).sign(0) >= 0
        for value in prior_ray["accepted_wall_pairings"].values()
    ):
        raise AssertionError("the inherited ray is not strict on an accepted cutter")
    minima = json.loads(MINIMA_AUDIT_PATH.read_text())
    if minima["verdict"] != "PASS":
        raise AssertionError("the independent Lambda22 minima audit does not pass")
    if minima["pari_count_only"]["S2_cumulative"]["10"] != [0, 0]:
        raise AssertionError("the inherited S2 minimum certificate drifted")
    if minima["pari_count_only"]["S2_cumulative"]["12"][1] != 12:
        raise AssertionError("the inherited S2 minimum is no longer 12")
    return {
        "prior_H1598_artifact_sha256": sha256(PREVIOUS_PATH.read_bytes()).hexdigest(),
        "independent_minima_artifact_sha256": sha256(
            MINIMA_AUDIT_PATH.read_bytes()
        ).hexdigest(),
        "prior_ray_strict_through_H1598": True,
        "ambient_minimum": 4,
        "S2_minimum": 12,
    }


def audit() -> dict[str, object]:
    scan_result = scan()
    compare_production(scan_result)
    continuity = continuity_audit()
    controls = scan_result.pop("controls")
    return {
        "schema_version": "0.1.0",
        "verdict": "PASS",
        "scope": (
            "Clean-room exact necessary-superset replay strictly after "
            "H=1598+1130sqrt(2) through inclusive identity height 5000."
        ),
        "independence": (
            "No project imports: standalone Q(sqrt(2)), coefficient-box proof, "
            "control census, binary Gram census, gate checks, and complete "
            "integral PSD extension enumeration."
        ),
        "sources": {
            "production_artifact": str(PRODUCTION_PATH.relative_to(ROOT)),
            "production_sha256": sha256(PRODUCTION_PATH.read_bytes()).hexdigest(),
        },
        "cutoff_semantics": {
            "lower": "strictly greater than 1598+1130sqrt(2) at the identity embedding",
            "upper": "less than or equal to the rational field element 5000 at the identity embedding",
            "first_replayed_height": controls[0].height.exact(),
            "last_replayed_height": controls[-1].height.exact(),
            "coefficient_box_proof": control_coefficient_box_proof(),
        },
        "continuity": continuity,
        "new_interval": scan_result,
        "logical_check": {
            "all_nonnegative_extensions_enumerated": True,
            "nonnegative_extensions_per_threat": [
                row["nonnegative_extension_count"] for row in scan_result["threats"]
            ],
            "all_nonnegative_extensions_short_vector_impossible": True,
        },
        "conclusion": (
            "The production H5000 bounded-prefix claim passes.  Combined with "
            "the inherited H1598 certificate, the ray is strict against every "
            "wall in the necessary superset through identity height 5000.  "
            "This does not prove global nonreflectivity or Vinberg termination."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    arguments = parser.parse_args()
    payload = audit()
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "verdict": payload["verdict"],
                "controls": payload["new_interval"]["control_count"],
                "retained_types": payload["new_interval"]["retained_type_count"],
                "threats": len(payload["new_interval"]["threats"]),
                "nonnegative_extensions_per_threat": payload["logical_check"][
                    "nonnegative_extensions_per_threat"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

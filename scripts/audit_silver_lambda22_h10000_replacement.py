#!/usr/bin/env python3
"""Certify the post-H2786 silver Lambda22 replacement through real H=10000.

The ray ``(1,(28/87)d)`` was proved strict through the fourth accepted wall
at ``2786+1970*sqrt(2)``.  This script exhausts the complete gated necessary
wall superset after that wall and through the rational identity-embedding
height bound 10000.  The conclusion is bounded; it is neither Vinberg
termination nor global nonreflectivity.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import isqrt, lcm
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_lambda22 import (
    LAMBDA22_COEFFICIENT,
    enumerate_lambda22_controls,
    gated_relaxed_witness_bound,
)
from coxeter_polytopes.search.silver_leech import candidate_leech_gram_triples
from coxeter_polytopes.search.silver_leech_manifest import (
    has_short_binary_combination,
)


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "data/enumerations/silver_lambda22_h2786_cutter.json"
JSON_PATH = ROOT / "data/enumerations/silver_lambda22_h10000_replacement.json"
REPORT_PATH = ROOT / "results/silver_lambda22_h10000_replacement.md"
START = Silver(2786, 1970)
BOUND = Silver(10000)
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
EXPECTED_THREAT = (
    Silver(3196, 2260),
    Silver(2),
    Silver(88, 62),
    (3194, 1600, 2260),
    (113, 80),
    Silver(Fraction(32, 87), Fraction(-22, 87)),
)


def norm_gate(square: int, modulus: int) -> bool:
    return (
        (not square or square >= MINIMUM[modulus])
        and square % (2 * modulus) == 0
        and (modulus % 2 or square % 8 in (0, 4))
    )


def value_json(value: Silver) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(0), value.sign(1)],
    }


def locally_admissible_extension(
    gram: tuple[int, int, int], d_dot_u: int, d_dot_v: int
) -> bool:
    """Apply the explicit ``[-3,3]^3`` minimum-four gate to ``<u,v,d>``.

    Failure is already a rigorous exclusion: it exhibits a nonzero integral
    combination of ``u,v,d`` having norm one, two, or three, impossible in
    ``Lambda22``.  Passing is only a necessary local condition, so every
    genuine lattice extension is retained for the exhaustive threat check.
    """

    u_square, v_square, u_dot_v = gram
    return not any(
        0
        < (
            first * first * u_square
            + second * second * v_square
            + 4 * third * third
            + 2 * first * second * u_dot_v
            + 2 * first * third * d_dot_u
            + 2 * second * third * d_dot_v
        )
        < 4
        for first in range(-3, 4)
        for second in range(-3, 4)
        for third in range(-3, 4)
        if first or second or third
    )


def nonnegative_threat_extensions() -> tuple[
    int, tuple[tuple[int, int, int, Silver], ...]
]:
    """Enumerate *all* extensions surviving the explicit local gate.

    Merely excluding the relaxed maximizing extension would not exclude a
    second nonnegative lattice extension.  This finite enumeration closes
    that logical gap.
    """

    gram = EXPECTED_THREAT[3]
    u_square, v_square, u_dot_v = gram
    locally_admissible = 0
    nonnegative: list[tuple[int, int, int, Silver]] = []
    for d_dot_u in range(-isqrt(4 * u_square), isqrt(4 * u_square) + 1):
        for d_dot_v in range(-isqrt(4 * v_square), isqrt(4 * v_square) + 1):
            determinant = (
                4 * (u_square * v_square - u_dot_v * u_dot_v)
                - v_square * d_dot_u * d_dot_u
                - u_square * d_dot_v * d_dot_v
                + 2 * u_dot_v * d_dot_u * d_dot_v
            )
            if determinant < 0 or not locally_admissible_extension(
                gram, d_dot_u, d_dot_v
            ):
                continue
            locally_admissible += 1
            maximum = Silver(d_dot_u, d_dot_v)
            margin = SCALE * maximum - LAMBDA22_COEFFICIENT * EXPECTED_THREAT[2]
            if margin.sign(0) >= 0:
                nonnegative.append((d_dot_u, d_dot_v, determinant, margin))
    result = tuple(nonnegative)
    if locally_admissible != 665 or result != (
        (113, 80, 0, Silver(Fraction(32, 87), Fraction(-22, 87))),
    ):
        raise AssertionError(
            "the exhaustive H10000 threat-extension enumeration changed: "
            f"{locally_admissible}, {result!r}"
        )
    return locally_admissible, result


def scan() -> tuple[int, int, tuple[tuple[object, ...], ...]]:
    controls = retained = 0
    threats: list[tuple[object, ...]] = []
    for level in enumerate_lambda22_controls(BOUND):
        if level.height <= START:
            continue
        controls += 1
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
            try:
                d_dot_u, d_dot_v, maximum = gated_relaxed_witness_bound(
                    *gram, u_modulus, v_modulus
                )
            except AssertionError:
                continue
            retained += 1
            margin = SCALE * maximum - LAMBDA22_COEFFICIENT * level.control
            if margin.sign(0) >= 0:
                threats.append(
                    (
                        level.height,
                        level.norm,
                        level.control,
                        gram,
                        (d_dot_u, d_dot_v),
                        margin,
                    )
                )
    result = tuple(threats)
    if (controls, retained, result) != (129, 1555, (EXPECTED_THREAT,)):
        raise AssertionError(
            f"the post-H2786 H10000 frontier changed: {controls}, {retained}, {result!r}"
        )
    return controls, retained, result


def threat_resolution() -> dict[str, object]:
    height, norm, control, gram, pairings, margin = EXPECTED_THREAT
    u_square, v_square, u_dot_v = gram
    d_dot_u, d_dot_v = pairings
    # The relaxed maximum has v^2=1600 and d.v=80 with d^2=4, hence
    # equality in Cauchy and v=20d.  Put e=u-28d.  Its exact norm is two.
    v_schur = 4 * v_square - d_dot_v * d_dot_v
    cross_schur = 4 * u_dot_v - d_dot_u * d_dot_v
    short_square = u_square - 56 * d_dot_u + 28 * 28 * 4
    if (v_schur, cross_schur, short_square) != (0, 0, 2):
        raise AssertionError("the unique H10000 threat reduction changed")
    locally_admissible, nonnegative = nonnegative_threat_extensions()
    d_dot_u_all, d_dot_v_all, determinant, margin_all = nonnegative[0]
    if (d_dot_u_all, d_dot_v_all, margin_all) != (
        d_dot_u,
        d_dot_v,
        margin,
    ):
        raise AssertionError("the exhaustive threat extension disagrees with the maximum")
    return {
        "height": value_json(height),
        "norm": value_json(norm),
        "control": value_json(control),
        "gram": list(gram),
        "maximizing_d_pairings": list(pairings),
        "relaxed_margin": value_json(margin),
        "locally_admissible_extension_count": locally_admissible,
        "local_extension_gate": (
            "PSD plus absence of a norm 1--3 combination in the explicit "
            "coefficient box [-3,3]^3. Failure supplies a rigorous min-4 "
            "contradiction; passing is only necessary."
        ),
        "nonnegative_extension_count": len(nonnegative),
        "all_nonnegative_extensions": [
            {
                "d_pairings": [d_dot_u_all, d_dot_v_all],
                "extension_determinant": determinant,
                "margin": value_json(margin_all),
            }
        ],
        "v_Schur_numerator": v_schur,
        "cross_Schur_numerator": cross_schur,
        "forced_relation": "v=20d",
        "short_vector": "e=u-28d",
        "short_vector_norm": short_square,
        "resolution": (
            "Exhaustive enumeration finds exactly one nonnegative locally "
            "admissible extension.  It forces a nonzero "
            "Lambda22 vector of norm 2, contradicting min(Lambda22)=4."
        ),
    }


def build_payload() -> dict[str, object]:
    if not BASE_PATH.exists():
        raise FileNotFoundError(BASE_PATH)
    base = json.loads(BASE_PATH.read_text())
    replacement = base["replacement_hyperideal_ray"]
    if replacement["scale"] != "28/87" or replacement["norm"]["signs"][0] != 1:
        raise AssertionError("the frozen H2786 replacement dependency changed")
    controls, retained, threats = scan()
    ray_norm = 4 * SCALE * SCALE - LAMBDA22_COEFFICIENT
    if ray_norm != Silver(Fraction(10705, 7569), -1) or ray_norm.sign(0) <= 0:
        raise AssertionError("the H10000 replacement norm changed")
    return {
        "schema_version": "0.1.0",
        "seed_id": "silver-H22-Lambda22-a-sqrt2-1",
        "scope": (
            "Complete gated necessary-wall obstruction strictly after H2786 "
            "and through rational identity height 10000; bounded prefix only."
        ),
        "base_H2786_artifact": str(BASE_PATH.relative_to(ROOT)),
        "height_start_exclusive": value_json(START),
        "height_bound_inclusive": value_json(BOUND),
        "replacement_ray": "(1,(28/87)d), d^2=4",
        "replacement_norm": value_json(ray_norm),
        "new_control_count": controls,
        "retained_gate_compatible_type_count": retained,
        "nonnegative_relaxed_row_count": len(threats),
        "unique_threat": threat_resolution(),
        "conclusion": (
            "The unique relaxed threat is arithmetically empty, so the strict "
            "hyperideal ray survives every possible wall through real H=10000."
        ),
        "discovery": None,
    }


def report() -> str:
    return """# Silver Lambda22 replacement survives through real height 10000

For `<-(sqrt(2)-1)> + Lambda22` in `H22`, the fourth accepted cutter at
`2786+1970sqrt(2)` leaves the exact ray

`z=(1,(28/87)d)`,  `d^2=4`,

with norm `10705/7569-sqrt(2)>0`.  Starting strictly after that cutter, the
complete control scan through the rational identity-height bound `H=10000`
contains 129 controls and 1,555 gate-compatible necessary Gram types.

Only one relaxed type can pair nonnegatively with `z`.  It occurs at
`3196+2260sqrt(2)`, has coefficient Gram `(3194,1600,2260)`, and maximizing
pairings `(d.u,d.v)=(113,80)`.  An exhaustive integral PSD enumeration,
followed by the explicit `[-3,3]^3` short-combination gate, checks all 665
extensions not already contradicted by `min(Lambda22)=4` and
finds exactly one nonnegative extension, namely that maximizer.  Equality
in the `v` and cross Schur bounds
forces `v=20d`; then `e=u-28d` has norm two, contradicting
`min(Lambda22)=4`.  Thus the type is empty.

The complete bounded prefix through real height 10000 therefore retains a
strict hyperideal patch and is not finite-volume.  This is not Vinberg
termination and does not prove global nonreflectivity.

```sh
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h10000_replacement.py
PYTHONPATH=src pytest -q tests/test_silver_lambda22_h10000_replacement.py
```
"""


def main() -> None:
    JSON_PATH.write_text(json.dumps(build_payload(), indent=2, sort_keys=True) + "\n")
    REPORT_PATH.write_text(report())
    print(JSON_PATH)
    print(REPORT_PATH)


if __name__ == "__main__":
    main()

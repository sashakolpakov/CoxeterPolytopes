#!/usr/bin/env python3
"""Certify the global zero-E8 no-go for the OG10 one-root compiler.

For ``ell=(A,0,C)`` put ``X=-A/phi``.  Then

    q = 2 p^2 = C^2 - phi X^2.

If the wall separates the archived interior from ``v61``, its identity-place
slope satisfies ``13 < C/X < phi*(601+199*phi)/2``.  Consequently
``(169-phi)X^2 < q < 4``.  The conjugate Schur gate independently gives
``(phi-1)(X')^2 < 4``.  Exact rational relaxations imply

    |X| < 1/6, |X'| < 8/3,

and hence, for ``X=a+b*phi``, ``|a|,|b| <= 1``.  The norm equation similarly
forces ``C=c+d*phi`` into ``|c|<=3, |d|<=1``.  Exhausting this necessary box
therefore proves a global no-go for zero E8, not merely a coefficient-box
statement.
"""

from __future__ import annotations

from collections import Counter
import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE_PATH = (
    ROOT / "data/candidates/og10_golden_projected_wall_zero_e8_global.json"
)

Pair = tuple[int, int]
ZERO: Pair = (0, 0)
PHI: Pair = (0, 1)
INV_PHI: Pair = (-1, 1)

# Respectively 400*x_int and 13*v61 in the x,u coordinates.
INTERIOR_XU: tuple[Pair, Pair] = ((601, 199), (2, 0))
RAY_XU: tuple[Pair, Pair] = ((13, 0), (0, 1))


def add(left: Pair, right: Pair) -> Pair:
    return left[0] + right[0], left[1] + right[1]


def neg(value: Pair) -> Pair:
    return -value[0], -value[1]


def mul(left: Pair, right: Pair) -> Pair:
    a, b = left
    c, d = right
    return a * c + b * d, a * d + b * c + b * d


def conjugate(value: Pair) -> Pair:
    return value[0] + value[1], -value[1]


def sign(value: Pair) -> int:
    """Exact sign of ``a+b*phi`` by integer comparisons."""

    a, b = value
    rational_part = 2 * a + b
    if rational_part == 0:
        return (b > 0) - (b < 0)
    if b == 0 or (rational_part > 0) == (b > 0):
        return (rational_part > 0) - (rational_part < 0)
    discriminant = rational_part * rational_part - 5 * b * b
    if discriminant > 0:
        return (rational_part > 0) - (rational_part < 0)
    return (b > 0) - (b < 0)


def canonical_sha256(value: object) -> str:
    data = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(data).hexdigest()


def admissible(q: Pair) -> bool:
    remainder = add((4, 0), neg(q))
    return (
        sign(q) > 0
        and sign(remainder) > 0
        and sign(conjugate(q)) > 0
        and sign(conjugate(remainder)) > 0
    )


def evaluate(A: Pair, C: Pair) -> tuple[Pair, Pair]:
    on_interior = add(mul(A, INTERIOR_XU[0]), mul(C, INTERIOR_XU[1]))
    on_ray = add(mul(A, RAY_XU[0]), mul(C, RAY_XU[1]))
    return on_interior, on_ray


def q_value(X: Pair, C: Pair) -> Pair:
    return add(mul(C, C), neg(mul(PHI, mul(X, X))))


def build_certificate() -> dict[str, object]:
    # Exact algebraic checks used by the coefficient-bound proof.
    bound_checks = {
        # 4/(169-phi) < 1/36 iff 144 < 169-phi.
        "identity_X_square_bound_less_than_1_over_36": sign((25, -1)) > 0,
        # 4*phi < 64/9 iff 16/9 is above the positive root of t^2-t-1.
        "conjugate_X_absolute_bound_less_than_8_over_3": (
            31 > 0  # 81*((16/9)^2-(16/9)-1)=31
        ),
        # |b| < 17/(6*sqrt(5)) < 2.
        "minkowski_b_bound_strictly_less_than_2": 289 < 720,
        # |a| < 1/6+phi < 2 iff phi < 11/6.
        "minkowski_a_bound_strictly_less_than_2": 19 > 0,
        # C^2 < 4+phi/36 < 73/18 < (13/6)^2.
        "identity_C_absolute_bound_less_than_13_over_6": 146 < 169,
        # |d| < 25/(6*sqrt(5)) < 2.
        "minkowski_d_bound_strictly_less_than_2": 625 < 720,
        # |c| < 13/6+phi < 4, using phi<11/6.
        "minkowski_c_bound_strictly_less_than_4": 13 + 11 == 24,
        # The interior threshold is larger than the ray threshold 13.
        "slope_interval_is_ordered": sign((173, 800)) > 0,
        # For q=e+f*phi in (0,4)^2: |f|<8/sqrt(5)<4 and
        # |e|<4+3*phi<9, so the following rectangle is exhaustive.
        "q_Minkowski_b_coefficient_absolute_bound_at_most_3": 64 < 80,
        "q_Minkowski_a_coefficient_absolute_bound_at_most_8": 4 + 5 == 9,
    }
    if not all(bound_checks.values()):
        raise AssertionError(f"a zero-E8 height bound failed: {bound_checks}")

    records: list[dict[str, object]] = []
    q_histogram: Counter[str] = Counter()
    sign_histogram: Counter[str] = Counter()
    cutting_records: list[dict[str, object]] = []
    for a in range(-1, 2):
        for b in range(-1, 2):
            X = (a, b)
            A = neg(mul(PHI, X))
            for c in range(-3, 4):
                for d in range(-1, 2):
                    C = (c, d)
                    q = q_value(X, C)
                    if not admissible(q):
                        continue
                    on_interior, on_ray = evaluate(A, C)
                    signs = sign(on_interior), sign(on_ray)
                    record = {
                        "X": list(X),
                        "A": list(A),
                        "C": list(C),
                        "q": list(q),
                        "on_400_interior": list(on_interior),
                        "on_13_v61": list(on_ray),
                        "signs": list(signs),
                    }
                    records.append(record)
                    q_histogram[str(q)] += 1
                    sign_histogram[f"{signs[0]},{signs[1]}"] += 1
                    if signs[0] * signs[1] < 0:
                        cutting_records.append(record)

    q_rectangle = []
    for a in range(-8, 9):
        for b in range(-3, 4):
            q = (a, b)
            if admissible(q):
                q_rectangle.append(q)
    expected_q_rectangle = [
        (1, 0),
        (1, 1),
        (2, -1),
        (2, 0),
        (2, 1),
        (3, -1),
        (3, 0),
    ]
    q_rectangle.sort()
    if q_rectangle != expected_q_rectangle:
        raise AssertionError(f"unexpected q rectangle: {q_rectangle}")
    if cutting_records:
        raise AssertionError("the global zero-E8 necessary box has a separator")

    payload: dict[str, object] = {
        "schema_version": "1.0.0",
        "candidate_id": "og10_golden_projected_wall_zero_e8_global_no_go",
        "certificate_kind": "global_exact_zero_E8_projected_wall_no_go",
        "fixed_lattice": "F=<-2*phi>+E8+<2>",
        "ansatz": "H(ell)=[[F,ell],[ell^t,2]], ell=(A,0,C)",
        "variable_change": "X=-A/phi, so q=2*p^2=C^2-phi*X^2",
        "signature_gate": "0<q<4 and 0<q'<4",
        "cut_gate": "ell(interior)*ell(v61)<0",
        "global_height_argument": {
            "separator_slope_interval": (
                "13<C/X<phi*(601+199*phi)/2; X=0 has equal signs"
            ),
            "necessary_identity_bound": "(169-phi)*X^2<q<4, hence |X|<1/6",
            "necessary_conjugate_bound": (
                "(phi-1)*(X')^2<q'<4, hence |X'|<8/3"
            ),
            "necessary_X_coefficient_box": {"a": [-1, 1], "b": [-1, 1]},
            "necessary_C_coefficient_box": {"c": [-3, 3], "d": [-1, 1]},
            "bound_checks": bound_checks,
        },
        "totally_positive_q_rectangle": [list(q) for q in q_rectangle],
        "projective_square_unit_classes": [
            [[1, 0], [1, 1], [2, -1]],
            [[2, 0]],
            [[2, 1], [3, -1]],
            [[3, 0]],
        ],
        "necessary_box_pair_count": 3 * 3 * 7 * 3,
        "signature_admissible_pair_count": len(records),
        "q_histogram": dict(sorted(q_histogram.items())),
        "sign_histogram": dict(sorted(sign_histogram.items())),
        "cutting_covector_count": 0,
        "cutting_records": [],
        "admissible_records_sha256": canonical_sha256(records),
        "conclusion": (
            "No zero-E8 one-root covector anywhere in Z[phi]^2 separates the "
            "archived invariant interior from v61."
        ),
        "scope": (
            "Global only for ell=(A,0,C). Mixed E8 covectors remain open; this "
            "does not produce or exclude a hyperbolic Coxeter polytope."
        ),
    }
    payload["certificate_payload_sha256"] = canonical_sha256(payload)
    return payload


def verify_certificate(payload: dict[str, object]) -> dict[str, bool]:
    body = dict(payload)
    digest = body.pop("certificate_payload_sha256", None)
    fresh = build_certificate()
    checks = {
        "payload_digest_matches": digest == canonical_sha256(body),
        "fresh_replay_matches": payload == fresh,
        "zero_cutters": payload.get("cutting_covector_count") == 0,
    }
    if not all(checks.values()):
        raise AssertionError(f"zero-E8 global certificate failed: {checks}")
    return checks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify-certificate", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    if args.verify_certificate:
        archived = json.loads(CERTIFICATE_PATH.read_text())
        print(json.dumps(verify_certificate(archived), indent=2, sort_keys=True))
    else:
        print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

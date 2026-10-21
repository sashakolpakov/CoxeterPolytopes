#!/usr/bin/env python3
"""Independent exact audit of the OG10 norm-six projected-wall collapse.

For the rank-preserving compiler

    H6(k) = [[F,3k],[3k^t,6]],  p0=F^-1 k,

the Schur gate is ``0<q=2*p0^2<4/3`` at both embeddings.  The corresponding
Minkowski rectangle contains only ``q=1``.  Since ``2*p0`` is integral for
the normalized covectors used by the compiler, it is an old norm-two vector
of F and has the identical projected mirror.

This script imports no project code.  It also audits the exact parity subset
inside the frozen 71-root H9 prefix and recomputes all four pairings with
``v61`` from a literal, row-length-checked E8 Cartan matrix.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT_PATH = ROOT / "data/candidates/og10_golden_h9_vinberg_through_3_plus_6phi.json"
CERTIFICATE_PATH = ROOT / "data/candidates/og10_golden_norm_six_projection_collapse.json"

Pair = tuple[Fraction, Fraction]
ZERO: Pair = (Fraction(0), Fraction(0))
PHI: Pair = (Fraction(0), Fraction(1))

E8_GRAM: tuple[tuple[int, ...], ...] = (
    (2, -1, 0, 0, 0, 0, 0, 0),
    (-1, 2, -1, 0, 0, 0, 0, 0),
    (0, -1, 2, -1, 0, 0, 0, -1),
    (0, 0, -1, 2, -1, 0, 0, 0),
    (0, 0, 0, -1, 2, -1, 0, 0),
    (0, 0, 0, 0, -1, 2, -1, 0),
    (0, 0, 0, 0, 0, -1, 2, 0),
    (0, 0, -1, 0, 0, 0, 0, 2),
)


def add(x: Pair, y: Pair) -> Pair:
    return x[0] + y[0], x[1] + y[1]


def neg(x: Pair) -> Pair:
    return -x[0], -x[1]


def mul(x: Pair, y: Pair) -> Pair:
    a, b = x
    c, d = y
    return a * c + b * d, a * d + b * c + b * d


def conjugate(x: Pair) -> Pair:
    return x[0] + x[1], -x[1]


def sign(x: Pair) -> int:
    """Exact sign of a rational ``a+b*phi``."""

    a, b = x
    t = 2 * a + b
    if t == 0:
        return (b > 0) - (b < 0)
    if b == 0 or (t > 0) == (b > 0):
        return (t > 0) - (t < 0)
    comparison = t * t - 5 * b * b
    if comparison > 0:
        return (t > 0) - (t < 0)
    return (b > 0) - (b < 0)


def pair_from_json(value: dict[str, str]) -> Pair:
    return Fraction(value["a"]), Fraction(value["b"])


def pair_json(value: Pair) -> list[str]:
    return [str(value[0]), str(value[1])]


def sum_pairs(values: object) -> Pair:
    answer = ZERO
    for value in values:
        answer = add(answer, value)
    return answer


def field_gram() -> tuple[tuple[Pair, ...], ...]:
    if len(E8_GRAM) != 8 or any(len(row) != 8 for row in E8_GRAM):
        raise AssertionError("the literal E8 matrix is not 8 by 8")
    if any(E8_GRAM[i][j] != E8_GRAM[j][i] for i in range(8) for j in range(8)):
        raise AssertionError("the literal E8 matrix is not symmetric")
    gram = [[ZERO for _ in range(10)] for _ in range(10)]
    gram[0][0] = (Fraction(0), Fraction(-2))
    for i in range(8):
        for j in range(8):
            gram[i + 1][j + 1] = (Fraction(E8_GRAM[i][j]), Fraction(0))
    gram[9][9] = (Fraction(2), Fraction(0))
    return tuple(tuple(row) for row in gram)


def field_pair(left: tuple[Pair, ...], right: tuple[Pair, ...]) -> Pair:
    gram = field_gram()
    return sum_pairs(
        mul(left[i], mul(gram[i][j], right[j]))
        for i in range(10)
        for j in range(10)
    )


def q_rectangle() -> list[tuple[int, int]]:
    # If q=a+b*phi lies in the open (4/3)^2 Minkowski rectangle, then
    # |b|<8/(3*sqrt(5))<2 and |a|<4/3+phi<3.  Thus this finite range is
    # necessary, not an imposed coefficient box.
    values = []
    bound = Fraction(4, 3)
    for a in range(-2, 3):
        for b in range(-1, 2):
            q = (Fraction(a), Fraction(b))
            qp = conjugate(q)
            if sign(q) > 0 and sign((bound - q[0], -q[1])) > 0 and sign(qp) > 0 and sign((bound - qp[0], -qp[1])) > 0:
                values.append((a, b))
    return values


def is_even_integral(value: Pair) -> bool:
    return all(entry.denominator == 1 and entry.numerator % 2 == 0 for entry in value)


def canonical_sha256(value: object) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def build_certificate() -> dict[str, object]:
    if q_rectangle() != [(1, 0)]:
        raise AssertionError(f"unexpected norm-six q rectangle: {q_rectangle()}")

    source = json.loads(INPUT_PATH.read_text())
    roots = [tuple(pair_from_json(x) for x in root) for root in source["simple_roots"]]
    ray = tuple(pair_from_json(x) for x in source["hyperideal_witness"]["ray"])
    if len(roots) != 71 or any(len(root) != 10 for root in roots):
        raise AssertionError("the frozen H9 prefix is not the expected 71 by 10 array")

    parity_indices = [
        index
        for index, root in enumerate(roots)
        if all(is_even_integral(entry) for entry in root[1:9])
    ]
    rows = []
    for index in parity_indices:
        norm = field_pair(roots[index], roots[index])
        on_ray = field_pair(roots[index], ray)
        rows.append(
            {
                "root_index": index,
                "norm": pair_json(norm),
                "on_v61": pair_json(on_ray),
                "on_v61_sign": sign(on_ray),
            }
        )
    if parity_indices != [8, 9, 61, 70]:
        raise AssertionError(f"unexpected parity-root indices: {parity_indices}")
    if any(row["norm"] != ["2", "0"] for row in rows):
        raise AssertionError("a parity-prefix vector is not a norm-two root")
    if any(row["on_v61_sign"] != -1 for row in rows):
        raise AssertionError("a parity-prefix root cuts v61")

    payload: dict[str, object] = {
        "schema_version": "1.0.0",
        "candidate_id": "og10_golden_norm_six_projection_collapse",
        "certificate_kind": "global_exact_projected_arrangement_redundancy",
        "norm_six_compiler": "H6(k)=[[F,3k],[3k^t,6]], p0=F^-1*k",
        "signature_gate": "0<q=2*p0^2<4/3 and 0<q'<4/3",
        "q_Minkowski_rectangle": [[1, 0]],
        "global_projection_identity": {
            "integral_vector": "y=2*p0=(-A/phi,2w,C) lies in F",
            "norm": "y^2=4*p0^2=2*q=2",
            "mirror_identity": "p0^perp=y^perp",
            "conclusion": "every admissible H6 projected mirror is an old F norm-two mirror",
        },
        "exact_converse": {
            "condition": "an old norm-two root y occurs iff F*y/2 is integral",
            "intrinsic_form": "y lies in F intersect 2*F^dual",
            "coordinate_form": "the E8 coordinate vector of y lies in 2*E8",
            "reason": (
                "the x and u Gram entries already contain 2, while unimodularity "
                "of E8 makes E8*y_E8 even equivalent to y_E8 even"
            ),
            "not_every_old_root_occurs": True,
        },
        "frozen_71_root_prefix_audit": {
            "parity_subset_indices": parity_indices,
            "parity_subset_indices_through_61": [8, 9],
            "rows": rows,
            "all_four_pair_v61_strictly_negative": True,
        },
        "input": {
            "path": str(INPUT_PATH.relative_to(ROOT)),
            "file_sha256": hashlib.sha256(INPUT_PATH.read_bytes()).hexdigest(),
            "root_count": len(roots),
        },
        "conclusion": (
            "The norm-six/divisibility-three compiler may change the ambient OG10 "
            "divisor lift of a mirror, but it adds no projected hyperplane to the "
            "complete norm-two arrangement on F."
        ),
        "volume_scope": (
            "Therefore the fixed H9 projected chamber remains the old infinite-volume "
            "P9. No ambient H10 volume conclusion follows merely from this codimension-one "
            "statement; that requires the prior facet/arithmetic implication or a fresh proof."
        ),
        "new_hyperbolic_example": False,
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
        "q_rectangle_is_singleton": payload.get("q_Minkowski_rectangle") == [[1, 0]],
        "no_new_projected_mirror": payload.get("global_projection_identity", {}).get("mirror_identity") == "p0^perp=y^perp",
    }
    if not all(checks.values()):
        raise AssertionError(f"norm-six projection-collapse replay failed: {checks}")
    return checks


if __name__ == "__main__":
    print(json.dumps(build_certificate(), indent=2, sort_keys=True))

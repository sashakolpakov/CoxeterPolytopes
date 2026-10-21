#!/usr/bin/env python3
"""Independent weak-log coloring audit for the standard-host P17--P19 packets.

The frozen producer correctly avoids reusing canonical MMP orbits when a
primitive norm-one ray is physically scaled from a (-1)-curve to a boundary
B=2r.  This independent replay pins that producer quartet and enumerates
every scaling directly.  Adjunction requires K.r=-1 or +1 on norm-one rays
and K.r=0 on norm-two rays.  Exact linear relations among all facet vectors
reject every sign assignment.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any, Sequence

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/rational_i1n_inverse_design.json"
SOURCE_SHA256 = "351a0a8525209313de13f9dd3752320bd074daf0985da704fc634684c75efc18"
PRODUCER_ARTIFACTS = {
    "report": (
        ROOT / "results/p17_p19_three_branch_property4_escape.md",
        "6bc82285d0ca42c4e7eeaed40d2b6302182083557e88db9c57e20e842d67b4b0",
    ),
    "certificate_generator": (
        ROOT / "scripts/certify_p17_p19_three_branch_property4_escape.py",
        "0755edff0bf042970ba19c87b01a79b2a44f5e21fd36662e570b2282953ad5e7",
    ),
    "archive": (
        ROOT / "data/enumerations/p17_p19_three_branch_property4_escape.json",
        "005b24f6ab15ee4114f25194fc38f605e456ecc693faaca327ad1119a03810af",
    ),
    "tests": (
        ROOT / "tests/test_p17_p19_three_branch_property4_escape.py",
        "3c57d95d4ec1ac858277716c8d3d44000d4d582afd62cb588619980ac829d746",
    ),
}
OUTPUT = ROOT / "data/enumerations/p17_p19_weak_log_coloring_independent_audit.json"
SCHEMA = "coxeter.p17_p19.weak_log_coloring_independent_audit.v2"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _dot(left: Sequence[int], right: Sequence[int]) -> int:
    return left[0] * right[0] - sum(a * b for a, b in zip(left[1:], right[1:]))


def _primitive_integer_row(row: Sequence[sp.Rational]) -> tuple[int, ...] | None:
    if all(value == 0 for value in row):
        return None
    denominator_lcm = math.lcm(*(int(value.q) for value in row))
    values = [int(value * denominator_lcm) for value in row]
    divisor = math.gcd(*(abs(value) for value in values))
    values = [value // divisor for value in values]
    first = next(value for value in values if value)
    if first < 0:
        values = [-value for value in values]
    return tuple(values)


def _case(n: int, roots: list[list[int]]) -> dict[str, Any]:
    rank = n + 1
    norm_one = [index for index, root in enumerate(roots) if _dot(root, root) == -1]
    norm_two = [index for index, root in enumerate(roots) if _dot(root, root) == -2]
    if len(norm_one) + len(norm_two) != len(roots):
        raise ArithmeticError(f"P{n} has a wall outside the norm-one/two alphabet")
    norm_two_span_rank = sp.Matrix([roots[index] for index in norm_two]).rank()

    # A row of A evaluates the I_(1,n) pairing K.r.
    evaluation = sp.Matrix(
        [[root[0], *[-entry for entry in root[1:]]] for root in roots]
    )
    if evaluation.rank() != rank:
        raise ArithmeticError(f"P{n} walls do not span the declared host")
    independent_rows = list(evaluation.T.rref()[1])
    basis_evaluation = evaluation.extract(independent_rows, range(rank))
    transport = evaluation * basis_evaluation.inv()

    # E sends the sign vector on the norm-one facets to the desired full
    # pairing vector b.  Consistency is b=A A_basis^{-1} b_basis.
    selector = sp.zeros(len(roots), len(norm_one))
    for column, wall_index in enumerate(norm_one):
        selector[wall_index, column] = 1
    basis_selector = selector.extract(independent_rows, range(len(norm_one)))
    obstruction = transport * basis_selector - selector
    constraint_set: set[tuple[int, ...]] = set()
    for row in range(obstruction.rows):
        normalized = _primitive_integer_row(list(obstruction.row(row)))
        if normalized is not None:
            constraint_set.add(normalized)
    constraints = sorted(constraint_set)

    survivors: list[list[int]] = []
    first_failed_constraint_counts: dict[str, int] = {}
    for signs in itertools.product((-1, 1), repeat=len(norm_one)):
        failed = next(
            (
                constraint
                for constraint in constraints
                if sum(coefficient * sign for coefficient, sign in zip(constraint, signs))
                != 0
            ),
            None,
        )
        if failed is None:
            survivors.append(list(signs))
        else:
            key = ",".join(str(value) for value in failed)
            first_failed_constraint_counts[key] = first_failed_constraint_counts.get(key, 0) + 1

    return {
        "n": n,
        "rank": rank,
        "facet_count": len(roots),
        "norm_one_facet_labels": [index + 1 for index in norm_one],
        "norm_two_facet_count": len(norm_two),
        "norm_two_span_rank": norm_two_span_rank,
        "sign_assignment_count": 2 ** len(norm_one),
        "independent_wall_labels": [index + 1 for index in independent_rows],
        "primitive_consistency_constraints": [list(row) for row in constraints],
        "constraint_count": len(constraints),
        "surviving_sign_assignments": survivors,
        "survivor_count": len(survivors),
        "first_failed_constraint_counts": first_failed_constraint_counts,
        "conclusion": (
            "No linear functional K has K.r in {-1,+1} on every primitive "
            "norm-one wall and K.r=0 on every norm-two wall.  Thus no exact "
            "standard-host smooth-rational (-1)/(-2)/(-4) physical coloring "
            "exists, before MMP, primality, or semiampleness."
        ),
    }


def build_certificate() -> dict[str, Any]:
    if _sha256(SOURCE) != SOURCE_SHA256:
        raise ArithmeticError("the frozen P17--P19 wall source changed")
    producer_digests = {
        name: {
            "path": str(path.relative_to(ROOT)),
            "sha256": _sha256(path),
            "expected_sha256": expected,
        }
        for name, (path, expected) in PRODUCER_ARTIFACTS.items()
    }
    if any(
        item["sha256"] != item["expected_sha256"]
        for item in producer_digests.values()
    ):
        raise ArithmeticError("the frozen P17--P19 producer quartet changed")
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    rows = []
    for n in (17, 18, 19):
        roots = [item["vector"] for item in source["dimensions"][str(n)]["roots"]]
        rows.append(_case(n, roots))

    checks = {
        "source_digest_pinned": _sha256(SOURCE) == SOURCE_SHA256,
        "producer_quartet_digests_pinned": all(
            item["sha256"] == item["expected_sha256"]
            for item in producer_digests.values()
        ),
        "dimensions_exact": [row["n"] for row in rows] == [17, 18, 19],
        "assignment_counts_exact": [row["sign_assignment_count"] for row in rows]
        == [4, 32768, 32],
        "all_sign_assignments_rejected": all(row["survivor_count"] == 0 for row in rows),
        "norm_two_walls_already_span_full_hosts": all(
            row["norm_two_span_rank"] == row["rank"] for row in rows
        ),
        "nonempty_exact_constraint_systems": all(row["constraint_count"] > 0 for row in rows),
    }
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "scope": {
            "proved": (
                "weak-log physical-coloring emptiness for the pinned standard "
                "I_(1,n) wall modules of P17, P18, P19"
            ),
            "repair": (
                "independently confirms the frozen producer's repaired pre-MMP "
                "enumeration; an earlier draft invalidly reused canonical MMP "
                "orbits computed under the all-primitive physical scaling"
            ),
            "not_proved": [
                "absolute reflection maximality of P17, P18, or P19",
                "exhaustion of different commensurable integral wall hosts or similarities",
                "nonexistence of a nonrational or higher-dimensional carrier",
            ],
        },
        "producer_artifacts": producer_digests,
        "adjunction_dictionary": {
            "primitive_norm_one_physical_minus_one": "K.r=-1",
            "primitive_norm_one_doubled_boundary_B_2r": "K.r=+1",
            "primitive_norm_two_physical_minus_two": "K.r=0",
        },
        "rows": rows,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    payload["certificate_sha256"] = _canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not OUTPUT.is_file() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise AssertionError("archived P17--P19 coloring audit has drifted")
        print(f"archive check passed: {OUTPUT}")
        print(f"certificate_sha256={payload['certificate_sha256']}")
        return
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT)
        print(payload["certificate_sha256"])
        return
    print(rendered, end="")


if __name__ == "__main__":
    main()

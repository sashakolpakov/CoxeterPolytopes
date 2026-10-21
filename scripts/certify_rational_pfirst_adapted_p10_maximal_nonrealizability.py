#!/usr/bin/env python3
"""Exact certificate for the physically coloured P10 obstruction.

The finite calculation is deliberately small.  It checks the raw
I_(1,10) wall packet, its exact dual rays, reflection integrality, the
P10 support, and the two frozen AG/lineage comparisons.  The infinite
all-reflections descent is proved in the companion report; bounded reduced
root enumeration here is only a hostile regression test, not its proof.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt, lcm
from pathlib import Path
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/rational_pfirst_adapted_p10_maximal_nonrealizability.json"

AG_ARCHIVE = ROOT / "data/enumerations/coble_nodot_maxdegree3_closure.json"
AG_ARCHIVE_SHA256 = "56a0646809fff76514729a58887172eb916d28cbda699128ed2f96a23556e82d"
H11_ARCHIVE = ROOT / "data/enumerations/coble_type_i_two_boundary_full_ns_h11.json"
H11_ARCHIVE_SHA256 = "bcfa2de1e96bf4f41282345a64f8ee6f3a75e715287bab4962d7082584b4d985"

N = 10
DIM = N + 1


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dot(x: Sequence[int], y: Sequence[int]) -> int:
    return x[0] * y[0] - sum(a * b for a, b in zip(x[1:], y[1:]))


def primitive(v: Sequence[int]) -> tuple[int, ...]:
    g = 0
    for z in v:
        g = gcd(g, abs(z))
    if not g:
        raise ArithmeticError("zero vector has no primitive normalization")
    ans = tuple(z // g for z in v)
    first = next(z for z in ans if z)
    return tuple(-z for z in ans) if first < 0 else ans


def rank(rows: Sequence[Sequence[int | Fraction]]) -> int:
    if not rows:
        return 0
    a = [[Fraction(z) for z in row] for row in rows]
    r = 0
    for c in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        scale = a[r][c]
        a[r] = [z / scale for z in a[r]]
        for i in range(len(a)):
            if i != r and a[i][c]:
                scale = a[i][c]
                a[i] = [x - scale * y for x, y in zip(a[i], a[r])]
        r += 1
    return r


def null_line(rows: Sequence[Sequence[int]]) -> tuple[int, ...] | None:
    """Return a primitive generator when the exact nullspace is a line."""

    a = [[Fraction(z) for z in row] for row in rows]
    if not a:
        return None
    r = 0
    pivots: list[int] = []
    for c in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        scale = a[r][c]
        a[r] = [z / scale for z in a[r]]
        for i in range(len(a)):
            if i != r and a[i][c]:
                scale = a[i][c]
                a[i] = [x - scale * y for x, y in zip(a[i], a[r])]
        pivots.append(c)
        r += 1
    free = [c for c in range(len(a[0])) if c not in pivots]
    if len(free) != 1:
        return None
    ans = [Fraction(0)] * len(a[0])
    ans[free[0]] = 1
    for i, c in reversed(list(enumerate(pivots))):
        ans[c] = -sum(a[i][j] * ans[j] for j in free)
    den = 1
    for z in ans:
        den = lcm(den, z.denominator)
    return primitive([int(z * den) for z in ans])


def wall_packet() -> tuple[list[str], list[tuple[int, ...]]]:
    labels = ["E10"]
    roots: list[tuple[int, ...]] = []
    e10 = [0] * DIM
    e10[10] = 1
    roots.append(tuple(e10))
    for i in range(1, 10):
        v = [0] * DIM
        v[i] = 1
        v[i + 1] = -1
        labels.append(f"D{i}")
        roots.append(tuple(v))
    c = [0] * DIM
    c[0] = 1
    c[1] = c[2] = c[3] = -1
    labels.append("C")
    roots.append(tuple(c))
    f = (3,) + (-1,) * 10
    labels.append("F_boundary_half")
    roots.append(f)
    return labels, roots


def gram(vectors: Sequence[Sequence[int]]) -> list[list[int]]:
    return [[dot(x, y) for y in vectors] for x in vectors]


def determinant(rows: Sequence[Sequence[int]]) -> int:
    a = [[Fraction(z) for z in row] for row in rows]
    if len(a) != len(a[0]):
        raise ValueError("determinant requires a square matrix")
    ans = Fraction(1)
    for c in range(len(a)):
        pivot = next((i for i in range(c, len(a)) if a[i][c]), None)
        if pivot is None:
            return 0
        if pivot != c:
            a[c], a[pivot] = a[pivot], a[c]
            ans = -ans
        scale = a[c][c]
        ans *= scale
        for i in range(c + 1, len(a)):
            if a[i][c]:
                q = a[i][c] / scale
                a[i] = [x - q * y for x, y in zip(a[i], a[c])]
    if ans.denominator != 1:
        raise ArithmeticError("integral matrix acquired nonintegral determinant")
    return int(ans)


def reflection_matrix(root: Sequence[int]) -> list[list[int]]:
    square = dot(root, root)
    out: list[list[int]] = []
    for j in range(DIM):
        basis = [0] * DIM
        basis[j] = 1
        pairing = dot(basis, root)
        coeff = Fraction(-2 * pairing, square)
        col = [Fraction(basis[i]) + coeff * root[i] for i in range(DIM)]
        if any(z.denominator != 1 for z in col):
            raise ArithmeticError("nonintegral side reflection")
        out.append([int(z) for z in col])
    # Convert stored columns to ordinary row-major form.
    return [[out[j][i] for j in range(DIM)] for i in range(DIM)]


def matmul(a: Sequence[Sequence[int]], b: Sequence[Sequence[int]]) -> list[list[int]]:
    return [
        [sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def transpose(a: Sequence[Sequence[int]]) -> list[list[int]]:
    return [list(row) for row in zip(*a)]


def extreme_rays(roots: Sequence[Sequence[int]]) -> dict[tuple[int, ...], tuple[int, ...]]:
    # x.r >= 0 has ordinary coefficient row (r0,-r1,...,-r10).
    inequalities = [(r[0],) + tuple(-z for z in r[1:]) for r in roots]
    rays: dict[tuple[int, ...], tuple[int, ...]] = {}
    # A pointed full-dimensional cone in dimension 11 has every extreme ray
    # on at least ten linearly independent facets.
    for active in combinations(range(len(roots)), DIM - 1):
        rows = [inequalities[i] for i in active]
        if rank(rows) != DIM - 1:
            continue
        ray = null_line(rows)
        if ray is None:
            continue
        slacks = [dot(ray, root) for root in roots]
        if not all(z >= 0 for z in slacks):
            ray = tuple(-z for z in ray)
            slacks = [-z for z in slacks]
        if all(z >= 0 for z in slacks):
            ray = primitive(ray)
            rays[ray] = tuple(i for i, z in enumerate(slacks) if z == 0)
    return dict(sorted(rays.items()))


def partitions_with_square_sum(
    length: int, maximum: int, target: int
) -> Iterable[tuple[int, ...]]:
    if length == 0:
        if target == 0:
            yield ()
        return
    high = min(maximum, isqrt(target))
    for x in range(high, -1, -1):
        for tail in partitions_with_square_sum(length - 1, x, target - x * x):
            yield (x,) + tail


def bounded_reduced_root_audit(max_a: int = 12) -> dict[str, object]:
    rows: list[dict[str, object]] = []
    for a in range(1, max_a + 1):
        for q in (1, 2):
            for m in partitions_with_square_sum(10, a, a * a + q):
                if sum(m[:3]) > a:
                    continue
                d = sum(m) - 3 * a
                if d <= 0:
                    raise ArithmeticError("negative reduced root entered the wall cone")
                is_f = a == 3 and q == 1 and m == (1,) * 10
                reflected_a = a - 6 * d
                if not is_f and not abs(reflected_a) < a:
                    raise ArithmeticError("anticanonical descent failed")
                rows.append(
                    {
                        "a": a,
                        "q": q,
                        "m": list(m),
                        "d": d,
                        "reflected_a": reflected_a,
                        "simple_F": is_f,
                    }
                )
    return {
        "max_a": max_a,
        "reduced_root_count": len(rows),
        "rows": rows,
        "role": "bounded regression only; the report proves the infinite descent",
    }


def build() -> dict[str, object]:
    labels, roots = wall_packet()
    f = roots[-1]
    k = tuple(-z for z in f)
    b = tuple(2 * z for z in f)
    physical = [roots[0], *roots[1:11], b]

    expected_squares = [-1] + [-2] * 10 + [-1]
    if [dot(r, r) for r in roots] != expected_squares:
        raise ArithmeticError("primitive wall norms changed")
    if [dot(k, c) for c in physical] != [-1] + [0] * 10 + [2]:
        raise ArithmeticError("physical adjunction colours changed")
    if dot(k, k) != -1 or any(z % 2 != 1 for z in k):
        raise ArithmeticError("canonical characteristic data changed")
    if tuple(f[i] - b[i] // 2 for i in range(DIM)) != (0,) * DIM:
        raise ArithmeticError("wrong-log correction is no longer zero")

    g = gram(roots)
    expected_edges = {
        ("E10", "D9", 4),
        ("E10", "F_boundary_half", 0),
        ("D3", "C", 3),
        *{(f"D{i}", f"D{i + 1}", 3) for i in range(1, 9)},
    }
    edges: set[tuple[str, str, int]] = set()
    for i in range(len(roots)):
        for j in range(i + 1, len(roots)):
            p = g[i][j]
            if not p:
                continue
            norms = g[i][i] * g[j][j]
            if p * p == norms:
                order = 0
            elif 4 * p * p == norms:
                order = 3
            elif 2 * p * p == norms:
                order = 4
            else:
                raise ArithmeticError("unexpected Coxeter pairing")
            edges.add((labels[i], labels[j], order))
    if edges != expected_edges:
        raise ArithmeticError("P10 Coxeter support changed")

    form = [[0] * DIM for _ in range(DIM)]
    form[0][0] = 1
    for i in range(1, DIM):
        form[i][i] = -1
    reflection_checks = []
    for label, root in zip(labels, roots):
        matrix = reflection_matrix(root)
        if matmul(matmul(transpose(matrix), form), matrix) != form:
            raise ArithmeticError("reflection does not preserve I_(1,10)")
        if matmul(matrix, matrix) != [[int(i == j) for j in range(DIM)] for i in range(DIM)]:
            raise ArithmeticError("reflection is not involutive")
        det = determinant(matrix)
        if det != -1:
            raise ArithmeticError("reflection determinant changed")
        reflection_checks.append({"label": label, "determinant": det, "integral": True})

    # The first eleven roots (delete F) form a unimodular basis.  This also
    # proves pointedness of the twelve-inequality cone.
    basis_det = determinant(roots[:11])
    if abs(basis_det) != 1:
        raise ArithmeticError("standard root basis ceased to be unimodular")

    rays = extreme_rays(roots)
    if len(rays) != 19:
        raise ArithmeticError("P10 extreme-ray count changed")
    squares = [dot(v, v) for v in rays]
    if any(v[0] <= 0 for v in rays) or any(q < 0 for q in squares):
        raise ArithmeticError("hyperideal ray or wrong time orientation")
    if sum(q == 0 for q in squares) != 2:
        raise ArithmeticError("P10 cusp count changed")
    if min(dot(x, y) for x in rays for y in rays) < 0:
        raise ArithmeticError("dual rays do not share one causal sheet")
    for i in range(len(roots)):
        active = [v for v in rays if dot(v, roots[i]) == 0]
        if rank(active) != DIM - 1:
            raise ArithmeticError(f"facet {labels[i]} became redundant")

    interior = tuple(sum(v[j] for v in rays) for j in range(DIM))
    interior_slacks = [dot(interior, root) for root in roots]
    if dot(interior, interior) <= 0 or not all(z > 0 for z in interior_slacks):
        raise ArithmeticError("strict timelike interior certificate failed")

    if sha256(AG_ARCHIVE) != AG_ARCHIVE_SHA256:
        raise ArithmeticError("frozen complex Coble obstruction archive changed")
    ag = json.loads(AG_ARCHIVE.read_text())
    ag_no_go = ag["maximal_automaticity_counterexample"]
    if not (
        ag_no_go["raw_chamber"].startswith("classical P10")
        and ag_no_go["root_count"] == 10
        and ag_no_go["complete_cusp_lift_count"] == 1
        and ag_no_go["D_CobVin_arithmetic_gates_pass"]
        and not ag_no_go["exact_complex_one_boundary_strong_Coble_full_nef_carrier_exists"]
    ):
        raise ArithmeticError("frozen P10 AG obstruction changed")

    if sha256(H11_ARCHIVE) != H11_ARCHIVE_SHA256:
        raise ArithmeticError("frozen H11/26 comparison archive changed")
    h11 = json.loads(H11_ARCHIVE.read_text())
    if not (
        h11["full_NS_lattice"]["NS_identification"] == "odd unimodular I_(1,11)"
        and len(h11["facets"]["labels"]) == 26
        and h11["facets"]["labels"][:2] == ["B1", "B2"]
    ):
        raise ArithmeticError("H11/26 comparison invariants changed")

    reduced = bounded_reduced_root_audit()

    ray_rows = []
    for v, active in rays.items():
        ray_rows.append(
            {
                "vector": list(v),
                "square": dot(v, v),
                "active_facets": [labels[i] for i in active],
                "slacks": [dot(v, root) for root in roots],
            }
        )

    return {
        "schema": "rational-pfirst-adapted-p10-maximal-nonrealizability-v1",
        "status": "PASS_EXACT_RAW_MAXIMAL_PLUS_COMPLEX_PHYSICAL_CARRIER_NO_GO",
        "all_checks_pass": True,
        "scope": {
            "base_field_for_nonrealizability": "C",
            "raw_saturation_uses_classification": False,
            "forbidden_imports_used_in_raw_proof": [],
            "AG_terminal_input": "frozen Kondo/Dolgachev-Kondo one-boundary Coble obstruction",
            "not_claimed": [
                "positive-characteristic nonrealizability",
                "nonexistence of the abstract P10 polytope",
                "nonexistence of rational carriers with a different geometric colouring",
                "a universal adapted-host dimension bound",
            ],
        },
        "lattice": {
            "identification": "I_(1,10)",
            "rank": DIM,
            "hyperbolic_dimension": N,
            "K": list(k),
            "K_square": dot(k, k),
            "K_characteristic": True,
            "standard_root_basis_determinant": basis_det,
        },
        "walls": {
            "labels": labels,
            "primitive_normals": [list(r) for r in roots],
            "physical_normals": [list(r) for r in physical],
            "primitive_gram": g,
            "physical_square_split": {"minus_1": 1, "minus_2": 10, "minus_4": 1},
            "canonical_pairings_physical": [dot(k, c) for c in physical],
            "Coxeter_edges_order_0_means_parallel": [list(e) for e in sorted(edges)],
            "adapted_host": "M=<L,B/2>=L with doubled form 2q",
            "reflection_checks": reflection_checks,
        },
        "canonical_escape_and_log": {
            "F_minus_K": list(f),
            "B": list(b),
            "K_fails_E10_by": dot(k, roots[0]),
            "F_fails_boundary_half_by": dot(f, f),
            "D_equals_F_minus_B_over_2": [0] * DIM,
            "boundary_incidence": {
                "B_dot_each_minus_2_wall": [dot(b, r) for r in roots[1:11]],
                "B_dot_E10": dot(b, roots[0]),
                "weak_H1_H2_H3_gate": "PASS (one boundary; B.R=0; B.E=2)",
            },
        },
        "cone": {
            "strict_interior": list(interior),
            "strict_interior_square": dot(interior, interior),
            "strict_interior_slacks": interior_slacks,
            "extreme_ray_count": len(rays),
            "finite_ray_count": sum(q > 0 for q in squares),
            "ideal_ray_count": sum(q == 0 for q in squares),
            "hyperideal_ray_count": sum(q < 0 for q in squares),
            # JSON object keys are strings; normalize before archive comparison.
            "square_histogram": {
                str(square): count for square, count in sorted(Counter(squares).items())
            },
            "all_facets_irredundant": True,
            "rays": ray_rows,
        },
        "classification_free_reflection_saturation": {
            "primitive_integral_reflector_norms": [-1, -2],
            "descent_measure": "absolute H-coordinate a",
            "generators": "signed permutations, C=H-E1-E2-E3, F=3H-sum(Ei)",
            "infinite_proof_location": "companion report Section 3",
            "bounded_regression": reduced,
            "conclusion": "side group contains every integral hyperbolic reflection of I_(1,10)",
        },
        "lineage_and_AG_no_go": {
            "identified_raw_polytope": "classical P10 = R12-E8A1x2-C2-H10-12",
            "complex_physical_full_nef_carrier": False,
            "frozen_AG_reason": ag_no_go["AG_obstruction"],
            "AG_archive": str(AG_ARCHIVE.relative_to(ROOT)),
            "AG_archive_sha256": AG_ARCHIVE_SHA256,
            "H11_26_comparison": {
                "same_polytope": False,
                "reason": "P10 is H10/12 on I_(1,10) with one boundary; Type-I is H11/26 on I_(1,11) with two boundaries",
                "archive": str(H11_ARCHIVE.relative_to(ROOT)),
                "archive_sha256": H11_ARCHIVE_SHA256,
            },
        },
        "property_4_missing_gate": (
            "a P-native simultaneous reduced-irreducible effectivity/incidence certificate, "
            "especially distinguishing an effective smooth rational divisor B=2b from "
            "a merely numerical or nonreduced double of the primitive mirror class b"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    payload = build()
    if args.check_archive:
        frozen = json.loads(ARCHIVE.read_text())
        if payload != frozen:
            raise SystemExit("live certificate differs from frozen archive")
        print("archived adapted P10 certificate matches")
    elif args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(payload["status"])
        print(
            "rays",
            payload["cone"]["finite_ray_count"],
            payload["cone"]["ideal_ray_count"],
            payload["cone"]["hyperideal_ray_count"],
        )


if __name__ == "__main__":
    main()

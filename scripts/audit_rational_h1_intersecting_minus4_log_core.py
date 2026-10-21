#!/usr/bin/env python3
"""Independent exact hostile audit for the minimal H1 log-core note."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT
    / "data/enumerations/rational_h1_intersecting_minus4_log_core_hostile_audit.json"
)
REPORT = ROOT / "results/rational_h1_intersecting_minus4_log_core.md"

EXPECTED_LI_ARCHIVE_SHA256 = (
    "20991983ad5b4c4ab8fbd3340f3b438e91068f8e7e81cc851aaad9e2bc85d107"
)
EXPECTED_LI_TEX_SHA256 = (
    "da8ba69b8f96359f75fa95bfd5da44aa86f28dd98337f4f24030ed86cb0e4892"
)


def dot(
    left: list[Fraction],
    gram: list[list[Fraction]],
    right: list[Fraction],
) -> Fraction:
    return sum(
        (
            left[i] * gram[i][j] * right[j]
            for i in range(len(left))
            for j in range(len(right))
        ),
        Fraction(0),
    )


def determinant(matrix: list[list[Fraction]]) -> Fraction:
    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise ValueError("expected a nonempty square matrix")
    work = [row[:] for row in matrix]
    answer = Fraction(1)
    for column in range(n):
        pivot = next(
            (row for row in range(column, n) if work[row][column]), None
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            answer = -answer
        pivot_entry = work[column][column]
        answer *= pivot_entry
        for row in range(column + 1, n):
            scale = work[row][column] / pivot_entry
            for j in range(column + 1, n):
                work[row][j] -= scale * work[column][j]
    return answer


def inverse(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise ValueError("expected a nonempty square matrix")
    work = [
        row[:]
        + [Fraction(int(i == j)) for j in range(n)]
        for i, row in enumerate(matrix)
    ]
    for column in range(n):
        pivot = next(
            (row for row in range(column, n) if work[row][column]), None
        )
        if pivot is None:
            raise ValueError("singular matrix")
        work[column], work[pivot] = work[pivot], work[column]
        pivot_entry = work[column][column]
        work[column] = [entry / pivot_entry for entry in work[column]]
        for row in range(n):
            if row == column:
                continue
            scale = work[row][column]
            if scale:
                work[row] = [
                    work[row][j] - scale * work[column][j]
                    for j in range(2 * n)
                ]
    return [row[n:] for row in work]


def mat_vec(
    matrix: list[list[Fraction]], vector: list[Fraction]
) -> list[Fraction]:
    return [
        sum(
            (entry * value for entry, value in zip(row, vector, strict=True)),
            Fraction(0),
        )
        for row in matrix
    ]


def transpose(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    return [list(column) for column in zip(*matrix, strict=True)]


def mat_mul(
    left: list[list[Fraction]], right: list[list[Fraction]]
) -> list[list[Fraction]]:
    right_t = transpose(right)
    return [
        [
            sum(
                (a * b for a, b in zip(row, column, strict=True)),
                Fraction(0),
            )
            for column in right_t
        ]
        for row in left
    ]


def solve(matrix: list[list[int]], rhs: list[int]) -> list[Fraction]:
    rational = [[Fraction(entry) for entry in row] for row in matrix]
    return mat_vec(inverse(rational), [Fraction(entry) for entry in rhs])


def _basis_vector(rank: int, index: int) -> list[Fraction]:
    return [Fraction(int(i == index)) for i in range(rank)]


def build_li_h1_geometry() -> dict[str, Any]:
    """Rebuild Li's n=2, p1=p2=6 geometry in total-transform coordinates."""
    # NS basis: s,f,e11,...,e16,e21,...,e26 on F_2 and its 12 blowups.
    rank = 14
    ns_gram = [[Fraction(0) for _ in range(rank)] for _ in range(rank)]
    ns_gram[0][0] = -2
    ns_gram[0][1] = ns_gram[1][0] = 1
    for i in range(2, rank):
        ns_gram[i][i] = -1

    s = _basis_vector(rank, 0)
    f = _basis_vector(rank, 1)
    e = {
        (arm, j): _basis_vector(rank, 2 + (arm - 1) * 6 + (j - 1))
        for arm in (1, 2)
        for j in range(1, 7)
    }

    def add(*vectors: list[Fraction]) -> list[Fraction]:
        return [sum(entries, Fraction(0)) for entries in zip(*vectors, strict=True)]

    def scale(value: int | Fraction, vector: list[Fraction]) -> list[Fraction]:
        return [Fraction(value) * entry for entry in vector]

    base_boundary = add(s, scale(2, f))
    boundaries = {
        arm: add(base_boundary, *(scale(-1, e[arm, j]) for j in range(1, 7)))
        for arm in (1, 2)
    }

    physical: list[tuple[str, str, list[Fraction]]] = []
    for arm in (1, 2):
        physical.append((f"B{arm}", "B", boundaries[arm]))
    for arm in (1, 2):
        for j in range(1, 6):
            physical.append(
                (
                    f"E{arm}{j}",
                    "R",
                    add(e[arm, j], scale(-1, e[arm, j + 1])),
                )
            )
        physical.append((f"E{arm}6", "E", e[arm, 6]))
    physical.append(("F1", "R", add(f, scale(-1, e[1, 1]), scale(-1, e[2, 1]))))
    physical.append(("F2", "R", s))

    physical_names = [name for name, _, _ in physical]
    physical_squares = {
        name: int(dot(vector, ns_gram, vector))
        for name, _, vector in physical
    }
    physical_intersections = [
        [int(dot(left, ns_gram, right)) for _, _, right in physical]
        for _, _, left in physical
    ]

    # M basis: replace e16,e26 by B1/2,B2/2. This has index four over NS.
    retained_indices = [0, 1] + list(range(2, 7)) + list(range(8, 13))
    m_columns = [_basis_vector(rank, index) for index in retained_indices]
    m_columns += [scale(Fraction(1, 2), boundaries[1])]
    m_columns += [scale(Fraction(1, 2), boundaries[2])]
    change = transpose(m_columns)
    change_inverse = inverse(change)
    m_gram = mat_mul(
        mat_mul(transpose(change), ns_gram),
        change,
    )
    m_gram = [[2 * entry for entry in row] for row in m_gram]

    root_ns_vectors = [
        scale(Fraction(1, 2), vector) if species == "B" else vector
        for _, species, vector in physical
    ]
    root_m_vectors = [mat_vec(change_inverse, vector) for vector in root_ns_vectors]
    root_gram = [
        [dot(left, m_gram, right) for right in root_m_vectors]
        for left in root_m_vectors
    ]
    root_squares = [root_gram[i][i] for i in range(len(root_gram))]

    reflection_integral = []
    reflection_isometric = []
    for root, square in zip(root_m_vectors, root_squares, strict=True):
        covector = [
            sum(
                (m_gram[i][j] * root[j] for j in range(rank)),
                Fraction(0),
            )
            for i in range(rank)
        ]
        action = [
            [
                Fraction(int(i == j)) - Fraction(2) * root[i] * covector[j] / square
                for j in range(rank)
            ]
            for i in range(rank)
        ]
        reflection_integral.append(
            all(entry.denominator == 1 for row in action for entry in row)
        )
        reflection_isometric.append(
            mat_mul(mat_mul(transpose(action), m_gram), action) == m_gram
        )

    coxeter_orders: dict[str, int] = {"2": 0, "3": 0, "4": 0}
    unsupported_pairs: list[dict[str, Any]] = []
    for i in range(len(root_gram)):
        for j in range(i + 1, len(root_gram)):
            pairing = root_gram[i][j]
            if pairing == 0:
                coxeter_orders["2"] += 1
                continue
            ratio = pairing * pairing / (root_squares[i] * root_squares[j])
            order = {Fraction(1, 4): 3, Fraction(1, 2): 4}.get(ratio)
            if order is None:
                unsupported_pairs.append(
                    {
                        "left": physical_names[i],
                        "right": physical_names[j],
                        "pairing": str(pairing),
                        "normalized_square": str(ratio),
                    }
                )
            else:
                coxeter_orders[str(order)] += 1

    return {
        "picard_rank": rank,
        "physical_names": physical_names,
        "physical_species_counts": {
            "B(-4)": sum(value == -4 for value in physical_squares.values()),
            "R(-2)": sum(value == -2 for value in physical_squares.values()),
            "E(-1)": sum(value == -1 for value in physical_squares.values()),
        },
        "physical_squares": physical_squares,
        "physical_intersection_matrix": physical_intersections,
        "boundary_intersection": physical_intersections[0][1],
        "host_rank": rank,
        "host_gram_integral": all(
            entry.denominator == 1 for row in m_gram for entry in row
        ),
        "host_determinant_absolute": str(abs(determinant(m_gram))),
        "host_root_norm_counts": {
            "-2": sum(square == -2 for square in root_squares),
            "-4": sum(square == -4 for square in root_squares),
        },
        "all_wall_reflections_integral": all(reflection_integral),
        "all_wall_reflections_isometric": all(reflection_isometric),
        "coxeter_order_pair_counts": coxeter_orders,
        "unsupported_wall_pairs": unsupported_pairs,
        "all_physical_generators_negative": all(
            value < 0 for value in physical_squares.values()
        ),
    }


def free_tail(length: int) -> dict[str, Any]:
    if isinstance(length, bool) or not isinstance(length, int) or length < 1:
        raise ValueError("tail length must be a positive integer")
    roots = []
    for i in range(length - 1):
        root = [0] * length
        root[i], root[i + 1] = 1, -1
        roots.append(root)
    root = [0] * length
    root[-1] = 1
    roots.append(root)
    gram = [
        [
            -2 * sum(left[k] * right[k] for k in range(length))
            for right in roots
        ]
        for left in roots
    ]
    return {
        "length": length,
        "qM_gram": gram,
        "determinant": str(determinant([[Fraction(x) for x in row] for row in gram])),
        "expected_determinant": str((-2) ** length),
        "simple_actions": "adjacent coordinate swaps and terminal sign change",
        "generated_group": "W(B_m)=signed permutations",
        "generated_group_order": str(2**length * math.factorial(length)),
    }


def build_certificate() -> dict[str, Any]:
    pair = [[-4, 2], [2, -4]]
    discrepancies = solve(pair, [2, 2])
    geometry = build_li_h1_geometry()
    payload = {
        "schema_version": 1,
        "verdict": "PASS_with_scope_warning",
        "production_report_sha256": hashlib.sha256(REPORT.read_bytes()).hexdigest(),
        "core_arithmetic": {
            "pair_matrix": pair,
            "pair_determinant": str(
                determinant([[Fraction(x) for x in row] for row in pair])
            ),
            "Z_square": -4,
            "K_dot_Z": 4,
            "arithmetic_genus": 1,
            "adjoint_multidegrees": [0, 0],
            "adjoint_square_formula": "14-rho_core",
            "core_rank_bound": 14,
            "discrepancies": [str(entry) for entry in discrepancies],
        },
        "li_source_provenance": {
            "arxiv": "2207.12483v1",
            "archive_sha256": EXPECTED_LI_ARCHIVE_SHA256,
            "main_tex_sha256": EXPECTED_LI_TEX_SHA256,
            "source_scope": "Section n=2, dual-basis cone calculation",
            "parameters": {"p1": 6, "p2": 6},
        },
        "li_equality_geometry": geometry,
        "tail_controls": [free_tail(length) for length in (1, 2, 3, 8, 16)],
        "conditional_consequences": {
            "total_marked_facet_degree": "m<=T and rho(S)<=14+T",
            "trace_descent": "d(n+1)<=rho(core)<=14",
            "rational_field": "n<=13",
            "neither_hypothesis_is_automatic": True,
        },
        "scope_warning": (
            "Negative definiteness alone is used only for the local complex-"
            "analytic cusp contraction/discrepancy audit; the core theorem is "
            "independent of that contraction. No total-rank, finite-volume "
            "maximality, or unbounded maximal counterfamily theorem is proved."
        ),
    }
    certify(payload)
    payload["all_checks_pass"] = True
    return payload


def certify(payload: dict[str, Any]) -> None:
    core = payload["core_arithmetic"]
    if core != {
        "pair_matrix": [[-4, 2], [2, -4]],
        "pair_determinant": "12",
        "Z_square": -4,
        "K_dot_Z": 4,
        "arithmetic_genus": 1,
        "adjoint_multidegrees": [0, 0],
        "adjoint_square_formula": "14-rho_core",
        "core_rank_bound": 14,
        "discrepancies": ["-1", "-1"],
    }:
        raise AssertionError("the H1 core arithmetic changed")

    geometry = payload["li_equality_geometry"]
    if geometry["picard_rank"] != 14:
        raise AssertionError("Li p1=p2=6 has the wrong Picard rank")
    if geometry["physical_species_counts"] != {
        "B(-4)": 2,
        "R(-2)": 12,
        "E(-1)": 2,
    }:
        raise AssertionError("Li physical alphabet counts changed")
    if geometry["boundary_intersection"] != 2:
        raise AssertionError("Li boundary incidence changed")
    if geometry["host_root_norm_counts"] != {"-2": 4, "-4": 12}:
        raise AssertionError("Li host root split changed")
    if geometry["host_determinant_absolute"] != str(2**10):
        raise AssertionError("Li host determinant changed")
    if not geometry["host_gram_integral"]:
        raise AssertionError("Li host is not integral")
    if not geometry["all_wall_reflections_integral"]:
        raise AssertionError("a Li wall reflection is not integral")
    if not geometry["all_wall_reflections_isometric"]:
        raise AssertionError("a Li wall reflection is not isometric")
    if geometry["unsupported_wall_pairs"]:
        raise AssertionError("Li wall pair is not Coxeter crystallographic")
    if set(geometry["coxeter_order_pair_counts"]) != {"2", "3", "4"}:
        raise AssertionError("unexpected Coxeter order")

    for row in payload["tail_controls"]:
        length = row["length"]
        expected = [[0] * length for _ in range(length)]
        for i in range(length):
            expected[i][i] = -4 if i < length - 1 else -2
            if i + 1 < length:
                expected[i][i + 1] = expected[i + 1][i] = 2
        if row["qM_gram"] != expected:
            raise AssertionError(f"wrong free B_{length} tail")
        if row["determinant"] != str((-2) ** length):
            raise AssertionError(f"wrong free B_{length} determinant")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--li-archive", type=Path)
    parser.add_argument("--li-tex", type=Path)
    args = parser.parse_args()
    payload = build_certificate()
    if args.li_archive is not None:
        actual = hashlib.sha256(args.li_archive.read_bytes()).hexdigest()
        if actual != EXPECTED_LI_ARCHIVE_SHA256:
            raise SystemExit("Jennifer Li source-archive hash mismatch")
    if args.li_tex is not None:
        actual = hashlib.sha256(args.li_tex.read_bytes()).hexdigest()
        if actual != EXPECTED_LI_TEX_SHA256:
            raise SystemExit("Jennifer Li main-TeX hash mismatch")
    if args.write:
        ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        ARCHIVE.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if args.check_archive:
        frozen = json.loads(ARCHIVE.read_text(encoding="utf-8"))
        if frozen != payload:
            raise SystemExit("archive mismatch")
    print("PASS: hostile H1 core/construction/tail audit (scope warning retained)")


if __name__ == "__main__":
    main()

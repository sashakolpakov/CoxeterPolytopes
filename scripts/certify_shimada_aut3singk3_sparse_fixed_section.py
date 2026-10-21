#!/usr/bin/env python3
"""Sparse exact fixed-section gate for Shimada's X_0, X_1, and X_2 chambers.

The older X_1/X_2 certificates use SymPy matrices throughout and eventually
materialize the facet Gram matrix.  That is perfectly reasonable for 37 or
162 facets, but it is the wrong execution model for the 3,278-wall X_0 stock
while another memory-intensive search is live.  This module instead:

* parses the wall list twice and keeps only globally scaled integer vectors;
* extracts only the requested involution from the published group list;
* forms the involution's one- and two-wall orbits and oriented primitive sums;
* proves strictness with exact integral multiples of the projection witnesses;
* streams the complete pair-angle census and positive-Gram hash one row at a
  time, never retaining a quadratic pairing table.

The X_0 command is deliberately guarded.  ``--estimate-only`` is always safe,
but a full X_0 pass requires both ``--authorize-x0-full`` and at least 512 MiB
of currently free memory by default.  Passing the authorization flag records
an operational decision only; it does not enlarge the mathematical
denominator beyond the genuine source group ``Aut(X_0,a_0)``.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
import resource
import subprocess
import sys
import time
from typing import Any, Iterator, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = Path("/private/tmp/Aut3SingK3s_compdata/compdata")
X0_OUTPUT = (
    ROOT
    / "data/enumerations/shimada_aut3singk3_x0_symplectic_h11_sparse_screen.json"
)
SCHEMA = "coxeter.shimada.aut3singk3.sparse_fixed_section.v1"

FILE_SHA256 = {
    "GramS.txt": "b9c829e2253bdb62312d0484e622804d40b5bcf7ca6e2f1a96487d5729c84c00",
    "ample.txt": "33e29f7c6d2fc1cdb63fd10e07a5cc3192f1378c3cbceca362edbc1497b498ee",
    "walls.txt": "d9e63603667d69c71f760b3f73a887391d324c24c0c04103e14c3f0e2480eb55",
    "AutXa.txt": "d0844407e8c041177655b99b2d6af064a9077a93fa73dba33cd7a07f41111e80",
}

RANK = 20
_TOKEN = re.compile(r"\[|\]|-?\d+(?:\s*/\s*\d+)?")


@dataclass(frozen=True)
class SurfaceCase:
    key: str
    surface_index: int
    representative: int
    involution_type: str
    wall_count: int
    wall_orbit_sizes: tuple[int, ...]
    fixed_rank: int
    expected: dict[str, Any] | None


CASES = {
    "x0": SurfaceCase(
        key="x0",
        surface_index=0,
        representative=3,
        involution_type="unique source-group rank-12 involution class",
        wall_count=3278,
        wall_orbit_sizes=(
            60, 40, 180, 10, 144, 144, 240, 360, 180, 240, 240, 720, 720
        ),
        fixed_rank=12,
        expected=None,
    ),
    "x1-symplectic": SurfaceCase(
        key="x1-symplectic",
        surface_index=1,
        representative=3,
        involution_type="symplectic",
        wall_count=1094,
        wall_orbit_sizes=(45, 45, 10, 30, 72, 60, 12, 40, 120, 120, 120, 180, 120, 120),
        fixed_rank=12,
        expected={
            "wall_orbit_count": 594,
            "unique_restriction_count": 580,
            "restriction_square_sign_histogram": {
                "negative": 162, "zero": 28, "positive": 390,
            },
            "strict_facet_count": 162,
            "positive_wall_gram_sha256": (
                "cbd8b22f90a2dfdfa82cc5968a3b7e1dccbae681dc0e47b9b370bdc1b1eb2fbe"
            ),
            "coxeter_edge_counts": {
                "2": 2498, "3": 114, "4": 64, "6": 152,
                "infinity": 758, "ultraparallel": 8540, "forbidden": 915,
            },
        },
    ),
    "x2-rational": SurfaceCase(
        key="x2-rational",
        surface_index=2,
        representative=596,
        involution_type="rational",
        wall_count=478,
        wall_orbit_sizes=(36, 12, 40, 90, 30, 30, 120, 120),
        fixed_rank=11,
        expected={
            "wall_orbit_count": 242,
            "unique_restriction_count": 207,
            "restriction_square_sign_histogram": {
                "negative": 37, "zero": 45, "positive": 125,
            },
            "strict_facet_count": 37,
            "positive_wall_gram_sha256": (
                "48a86bd712574d5fdd057f09d5f56a5c8ec41e8128f8ad82040575064202f303"
            ),
            "coxeter_edge_counts": {
                "2": 405, "3": 15, "4": 70, "6": 0,
                "infinity": 121, "ultraparallel": 55, "forbidden": 0,
            },
        },
    ),
    "x2-symplectic": SurfaceCase(
        key="x2-symplectic",
        surface_index=2,
        representative=470,
        involution_type="symplectic",
        wall_count=478,
        wall_orbit_sizes=(36, 12, 40, 90, 30, 30, 120, 120),
        fixed_rank=12,
        expected={
            "wall_orbit_count": 266,
            "unique_restriction_count": 266,
            "restriction_square_sign_histogram": {
                "negative": 90, "zero": 0, "positive": 176,
            },
            "strict_facet_count": 90,
            "positive_wall_gram_sha256": (
                "023f4e46a6c582afdac7590e2db2eb7b537784e1d007107dd99317ba0d7a9bb7"
            ),
            "coxeter_edge_counts": {
                "2": 1152, "3": 48, "4": 8, "6": 32,
                "infinity": 155, "ultraparallel": 2258, "forbidden": 352,
            },
        },
    ),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _assignment(path: Path, name: str, index: int) -> Any:
    """Parse a small Magma-style rational assignment without ``eval``."""
    raw = path.read_text(encoding="utf-8")
    marker = f"{name}[{index}]:="
    start = raw.find(marker)
    if start < 0:
        raise ValueError(f"missing {marker} in {path}")
    start += len(marker)
    finish = raw.find(";", start)
    if finish < 0:
        raise ValueError(f"unterminated {marker} in {path}")
    tree = ast.parse(raw[start:finish].strip(), mode="eval")

    def visit(node: ast.AST) -> Any:
        if isinstance(node, ast.List):
            return [visit(child) for child in node.elts]
        if isinstance(node, ast.Constant) and isinstance(node.value, int):
            return Fraction(node.value)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -visit(node.operand)
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
            return visit(node.left) / visit(node.right)
        raise ValueError(f"unsupported syntax in {path}: {ast.dump(node)}")

    return visit(tree.body)


def _iter_rank3_vectors(
    path: Path, name: str, index: int
) -> Iterator[tuple[int, tuple[Fraction, ...]]]:
    """Yield ``(outer_index, vector)`` from a rank-three list assignment.

    Only the source text and one 20-entry vector are live in the parser.  This
    is used for both ``walls[surface][orbit][wall]`` and
    ``AutXa[surface][matrix][row]``.
    """
    raw = path.read_text(encoding="utf-8")
    marker = f"{name}[{index}]:="
    start = raw.find(marker)
    if start < 0:
        raise ValueError(f"missing {marker} in {path}")
    start += len(marker)
    finish = raw.find(";", start)
    if finish < 0:
        raise ValueError(f"unterminated {marker} in {path}")

    depth = 0
    outer = -1
    vector: list[Fraction] | None = None
    for match in _TOKEN.finditer(raw, start, finish):
        token = match.group(0)
        if token == "[":
            depth += 1
            if depth == 2:
                outer += 1
            elif depth == 3:
                vector = []
        elif token == "]":
            if depth == 3:
                if vector is None:
                    raise ValueError("rank-three parser lost its current vector")
                yield outer, tuple(vector)
                vector = None
            depth -= 1
            if depth < 0:
                raise ValueError("unbalanced closing bracket")
        elif vector is not None:
            vector.append(Fraction(token.replace(" ", "")))
    if depth != 0 or vector is not None:
        raise ValueError("unbalanced rank-three assignment")


def _integral_vector(values: Sequence[Fraction], label: str) -> tuple[int, ...]:
    if len(values) != RANK or any(value.denominator != 1 for value in values):
        raise ArithmeticError(f"{label} is not an integral rank-{RANK} vector")
    return tuple(value.numerator for value in values)


def _primitive_integer(values: Sequence[int]) -> tuple[int, ...]:
    content = math.gcd(*(abs(value) for value in values))
    if content == 0:
        raise ArithmeticError("zero restricted wall")
    # Do not normalize the sign: it is the oriented half-space sign.
    return tuple(value // content for value in values)


def _row_times_matrix(
    row: Sequence[int], matrix: Sequence[Sequence[int]]
) -> tuple[int, ...]:
    return tuple(
        sum(row[i] * matrix[i][j] for i in range(RANK))
        for j in range(RANK)
    )


def _dot(left: Sequence[int], right: Sequence[int]) -> int:
    return sum(a * b for a, b in zip(left, right, strict=True))


def _matrix_product(
    left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]
) -> tuple[tuple[int, ...], ...]:
    return tuple(_row_times_matrix(row, right) for row in left)


def _transpose(matrix: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(matrix[i][j] for i in range(RANK)) for j in range(RANK))


def _rank_fraction_free(rows: Sequence[Sequence[int]]) -> int:
    """Exact rank of a small integer matrix using rational row reduction."""
    work = [[Fraction(value) for value in row] for row in rows]
    if not work:
        return 0
    rank = 0
    for col in range(len(work[0])):
        pivot = next((i for i in range(rank, len(work)) if work[i][col]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][col]
        for i in range(rank + 1, len(work)):
            if not work[i][col]:
                continue
            multiple = work[i][col] / pivot_value
            for j in range(col, len(work[0])):
                work[i][j] -= multiple * work[rank][j]
        rank += 1
        if rank == len(work):
            break
    return rank


def _load_representative(
    path: Path, surface_index: int, representative: int
) -> tuple[tuple[int, ...], ...]:
    rows: list[tuple[int, ...]] = []
    iterator = _iter_rank3_vectors(path, "AutXa", surface_index)
    try:
        for outer, row in iterator:
            if outer < representative:
                continue
            if outer > representative:
                break
            rows.append(_integral_vector(row, f"AutXa[{surface_index}][{representative}] row"))
    finally:
        iterator.close()
    if len(rows) != RANK:
        raise ArithmeticError(
            f"representative {representative} has {len(rows)} rows, expected {RANK}"
        )
    return tuple(rows)


def _load_scaled_walls(
    path: Path,
    surface_index: int,
    expected_orbit_sizes: Sequence[int],
) -> tuple[list[tuple[int, ...]], int]:
    denominator = 1
    first_counts: Counter[int] = Counter()
    for orbit, wall in _iter_rank3_vectors(path, "walls", surface_index):
        if len(wall) != RANK:
            raise ArithmeticError(f"wall in orbit {orbit} does not have rank {RANK}")
        first_counts[orbit] += 1
        denominator = math.lcm(denominator, *(value.denominator for value in wall))
    if [first_counts[i] for i in range(len(first_counts))] != list(expected_orbit_sizes):
        raise ArithmeticError("published wall-orbit sizes changed")

    walls: list[tuple[int, ...]] = []
    second_counts: Counter[int] = Counter()
    for orbit, wall in _iter_rank3_vectors(path, "walls", surface_index):
        scaled = tuple(int(value * denominator) for value in wall)
        if any(Fraction(value, denominator) != source for value, source in zip(scaled, wall, strict=True)):
            raise ArithmeticError("global wall denominator did not clear a coordinate")
        walls.append(scaled)
        second_counts[orbit] += 1
    if second_counts != first_counts:
        raise ArithmeticError("the two wall-parser passes disagree")
    return walls, denominator


def _classify_pair(left_norm: int, right_norm: int, pairing: int) -> str:
    """Classify a crystallographic rational Coxeter product exactly."""
    if left_norm >= 0 or right_norm >= 0:
        raise ArithmeticError("pair classification requires negative facet norms")
    if pairing == 0:
        return "2"
    if pairing < 0:
        return "forbidden"
    denominator = left_norm * right_norm
    numerator = pairing * pairing
    if 4 * numerator == denominator:
        return "3"
    if 2 * numerator == denominator:
        return "4"
    if 4 * numerator == 3 * denominator:
        return "6"
    if numerator == denominator:
        return "infinity"
    if numerator > denominator:
        return "ultraparallel"
    return "forbidden"


def _positive_gram_and_pair_gate(
    facets: Sequence[tuple[int, ...]],
    covectors: Sequence[tuple[int, ...]],
    norms: Sequence[int],
    stop_on_first_forbidden: bool,
) -> tuple[dict[str, Any], str | None]:
    labels = ("2", "3", "4", "6", "infinity", "ultraparallel", "forbidden")
    counts: Counter[str] = Counter()
    forbidden: list[dict[str, Any]] = []

    if stop_on_first_forbidden:
        visited = 0
        for right in range(len(facets)):
            for left in range(right):
                pairing = _dot(covectors[left], facets[right])
                label = _classify_pair(norms[left], norms[right], pairing)
                counts[label] += 1
                visited += 1
                if label == "forbidden":
                    forbidden.append(
                        _forbidden_record(left, right, norms, pairing, facets)
                    )
                    return (
                        {
                            "exhaustive": False,
                            "visited_pair_count": visited,
                            "total_pair_count": len(facets) * (len(facets) - 1) // 2,
                            "coxeter_edge_counts_on_visited_prefix": {
                                key: counts[key] for key in labels
                            },
                            "all_pairs_coxeter": False,
                            "first_forbidden_pairs": forbidden,
                        },
                        None,
                    )
        return (
            {
                "exhaustive": True,
                "visited_pair_count": visited,
                "total_pair_count": visited,
                "coxeter_edge_counts": {key: counts[key] for key in labels},
                "all_pairs_coxeter": True,
                "first_forbidden_pairs": [],
            },
            None,
        )

    # Hash the JSON list-of-list of string-valued positive-Gram entries in the
    # exact format used by the frozen X_1/X_2 certificates.  Only one scalar is
    # live at a time.  The lower triangle simultaneously supplies the pair gate.
    digest = hashlib.sha256()
    digest.update(b"[")
    for right in range(len(facets)):
        if right:
            digest.update(b",")
        digest.update(b"[")
        for left in range(len(facets)):
            if left:
                digest.update(b",")
            pairing = _dot(covectors[right], facets[left])
            digest.update(json.dumps(str(-pairing)).encode())
            if left >= right:
                continue
            label = _classify_pair(norms[left], norms[right], pairing)
            counts[label] += 1
            if label == "forbidden" and len(forbidden) < 12:
                forbidden.append(
                    _forbidden_record(left, right, norms, pairing, facets)
                )
        digest.update(b"]")
    digest.update(b"]")
    pair_count = len(facets) * (len(facets) - 1) // 2
    return (
        {
            "exhaustive": True,
            "visited_pair_count": pair_count,
            "total_pair_count": pair_count,
            "coxeter_edge_counts": {key: counts[key] for key in labels},
            "all_pairs_coxeter": counts["forbidden"] == 0,
            "first_forbidden_pairs": forbidden,
        },
        digest.hexdigest(),
    )


def _forbidden_record(
    left: int,
    right: int,
    norms: Sequence[int],
    pairing: int,
    facets: Sequence[tuple[int, ...]],
) -> dict[str, Any]:
    ratio = Fraction(pairing * pairing, norms[left] * norms[right])
    return {
        "facet_pair_one_based": [left + 1, right + 1],
        "inner_product": pairing,
        "squares": [norms[left], norms[right]],
        "normalized_square_ratio": str(ratio),
        "primitive_facets_in_published_NS_basis": [
            list(facets[left]), list(facets[right]),
        ],
    }


def _strict_projection_witnesses(
    rays: Sequence[tuple[int, ...]],
    covectors: Sequence[tuple[int, ...]],
    norms: Sequence[int],
    ample: tuple[int, ...],
    ample_square: int,
    gram: Sequence[Sequence[int]],
) -> tuple[list[int], list[int], str]:
    negative = [i for i, norm in enumerate(norms) if norm < 0]
    strict: list[int] = []
    failed: list[int] = []
    records = []
    ample_pairings = [_dot(covector, ample) for covector in covectors]
    if ample_square <= 0 or min(ample_pairings) <= 0:
        raise ArithmeticError("the ample orientation is not strictly positive")

    for position in negative:
        record = _projection_witness_record(
            position, rays, norms, ample, ample_pairings, gram
        )
        ok = record["strict"]
        (strict if ok else failed).append(position)
        records.append(record)
    return strict, failed, _canonical_sha256(records)


def _projection_witness_record(
    position: int,
    rays: Sequence[tuple[int, ...]],
    norms: Sequence[int],
    ample: tuple[int, ...],
    ample_pairings: Sequence[int],
    gram: Sequence[Sequence[int]],
) -> dict[str, Any]:
    target = rays[position]
    target_norm = norms[position]
    target_ample = ample_pairings[position]
    # This is (-target_norm) times the usual orthogonal projection
    # a - (a,t)/(t,t) t, so its scale is positive.
    witness = tuple(
        (-target_norm) * ample[i] + target_ample * target[i]
        for i in range(RANK)
    )
    witness_covector = _row_times_matrix(witness, gram)
    witness_square = _dot(witness_covector, witness)
    slacks = [_dot(witness_covector, ray) for ray in rays]
    others = [value for i, value in enumerate(slacks) if i != position]
    ok = (
        target_norm < 0
        and witness_square > 0
        and slacks[position] == 0
        and bool(others)
        and min(others) > 0
    )
    return {
        "restriction_index_one_based": position + 1,
        "strict": ok,
        "scaled_integral_projection_witness": list(witness),
        "witness_square": witness_square,
        "target_slack": slacks[position],
        "minimum_other_slack": min(others) if others else None,
    }


def analyze_scaled(
    *,
    gram: Sequence[Sequence[int]],
    ample: Sequence[int],
    walls: Sequence[tuple[int, ...]],
    action: Sequence[Sequence[int]],
    global_wall_denominator: int = 1,
    stop_on_first_forbidden: bool = False,
) -> dict[str, Any]:
    """Analyze already scaled walls; public so small exact fixtures can replay it."""
    gram = tuple(tuple(int(value) for value in row) for row in gram)
    ample = tuple(int(value) for value in ample)
    action = tuple(tuple(int(value) for value in row) for row in action)
    if not (
        len(gram) == len(action) == len(ample) == RANK
        and all(len(row) == RANK for row in gram)
        and all(len(row) == RANK for row in action)
        and all(gram[i][j] == gram[j][i] for i in range(RANK) for j in range(RANK))
    ):
        raise ArithmeticError("ambient data are not rank-20 symmetric data")
    identity = tuple(
        tuple(int(i == j) for j in range(RANK)) for i in range(RANK)
    )
    if _matrix_product(action, action) != identity:
        raise ArithmeticError("the selected source element is not an involution")
    # Shimada's matrices act on row vectors on the right.  Thus the isometry
    # identity is A G A^T = G, not A^T G A = G.
    if _matrix_product(_matrix_product(action, gram), _transpose(action)) != gram:
        raise ArithmeticError("the selected source element does not preserve GramS")
    if _row_times_matrix(ample, action) != ample:
        raise ArithmeticError("the selected source element does not fix the ample class")
    equations = [
        [action[j][i] - identity[j][i] for j in range(RANK)]
        for i in range(RANK)
    ]
    fixed_rank = RANK - _rank_fraction_free(equations)

    wall_index = {wall: i for i, wall in enumerate(walls)}
    if len(wall_index) != len(walls):
        raise ArithmeticError("the complete source wall list contains duplicates")
    permutation = []
    for wall in walls:
        image = _row_times_matrix(wall, action)
        if image not in wall_index:
            raise ArithmeticError("the involution does not permute the source walls")
        permutation.append(wall_index[image])

    seen: set[int] = set()
    source_orbits: list[tuple[int, ...]] = []
    for start in range(len(walls)):
        if start in seen:
            continue
        image = permutation[start]
        if permutation[image] != start:
            raise ArithmeticError("the source-wall action is not involutive")
        orbit = (start,) if image == start else (start, image)
        seen.update(orbit)
        source_orbits.append(orbit)

    unique: dict[tuple[int, ...], dict[str, Any]] = {}
    for orbit in source_orbits:
        summed = tuple(sum(walls[index][j] for index in orbit) for j in range(RANK))
        ray = _primitive_integer(summed)
        if _row_times_matrix(ray, action) != ray:
            raise ArithmeticError("a restricted ray is not fixed by the involution")
        record = unique.setdefault(
            ray,
            {
                "primitive": ray,
                "source_orbit_multiplicity": 0,
                "first_source_orbit_one_based": [index + 1 for index in orbit],
            },
        )
        record["source_orbit_multiplicity"] += 1

    rays = list(unique)
    covectors = [_row_times_matrix(ray, gram) for ray in rays]
    norms = [_dot(covector, ray) for covector, ray in zip(covectors, rays, strict=True)]
    ample_covector = _row_times_matrix(ample, gram)
    ample_square = _dot(ample_covector, ample)
    strict_positions, failed_positions, witness_digest = _strict_projection_witnesses(
        rays, covectors, norms, ample, ample_square, gram
    )

    facets = [rays[i] for i in strict_positions]
    facet_covectors = [covectors[i] for i in strict_positions]
    facet_norms = [norms[i] for i in strict_positions]
    pair_gate, positive_gram_hash = _positive_gram_and_pair_gate(
        facets, facet_covectors, facet_norms, stop_on_first_forbidden
    )
    first_forbidden_witnesses = []
    if pair_gate["first_forbidden_pairs"]:
        ample_pairings = [_dot(covector, ample) for covector in covectors]
        facet_pair = pair_gate["first_forbidden_pairs"][0]["facet_pair_one_based"]
        for facet_one_based in facet_pair:
            ray_position = strict_positions[facet_one_based - 1]
            witness = _projection_witness_record(
                ray_position, rays, norms, ample, ample_pairings, gram
            )
            if witness["strict"] is not True:
                raise ArithmeticError("a forbidden-pair facet lost its strict witness")
            first_forbidden_witnesses.append(witness)

    sign_histogram = Counter(
        "negative" if norm < 0 else "zero" if norm == 0 else "positive"
        for norm in norms
    )
    for sign in ("negative", "zero", "positive"):
        sign_histogram[sign] += 0
    restriction_digest = _canonical_sha256(
        [
            {
                "primitive": list(ray),
                "norm": norm,
                "source_orbit_multiplicity": unique[ray]["source_orbit_multiplicity"],
                "first_source_orbit_one_based": unique[ray]["first_source_orbit_one_based"],
            }
            for ray, norm in zip(rays, norms, strict=True)
        ]
    )
    if failed_positions:
        classification = "irredundancy_unresolved"
    elif not pair_gate["all_pairs_coxeter"]:
        classification = "non_Coxeter_fixed_section"
    elif not pair_gate["exhaustive"]:
        classification = "pair_gate_incomplete"
    else:
        classification = "Coxeter_pair_survivor_volume_unchecked"
    return {
        "fixed_rank": fixed_rank,
        "hyperbolic_dimension": fixed_rank - 1,
        "source_action_checks": {
            "involutive": True,
            "preserves_ambient_gram": True,
            "fixes_ample_class": True,
            "permutes_complete_wall_set": True,
        },
        "complete_source_wall_count": len(walls),
        "global_source_wall_denominator": global_wall_denominator,
        "wall_orbit_count": len(source_orbits),
        "wall_orbit_size_histogram": {
            str(key): value
            for key, value in sorted(Counter(map(len, source_orbits)).items())
        },
        "unique_restriction_count": len(rays),
        "restriction_square_sign_histogram": dict(sign_histogram),
        "negative_restriction_count": len(strict_positions) + len(failed_positions),
        "strict_facet_count": len(strict_positions),
        "projection_witness_failure_count": len(failed_positions),
        "projection_witness_failure_indices_one_based": [i + 1 for i in failed_positions],
        "strict_projection_witnesses_sha256": witness_digest,
        "unique_restrictions_sha256": restriction_digest,
        "positive_wall_gram_sha256": positive_gram_hash,
        "pair_gate": pair_gate,
        "first_forbidden_pair_strict_witnesses": first_forbidden_witnesses,
        "classification": classification,
        "volume_status": (
            "not tested; Vinberg volume is authorized only after a complete "
            "Coxeter pair survivor"
        ),
    }


def _expected_checks(case: SurfaceCase, result: dict[str, Any]) -> dict[str, bool]:
    if case.expected is None:
        return {"no_frozen_downstream_result_exists": True}
    expected = case.expected
    return {
        "fixed_rank": result["fixed_rank"] == case.fixed_rank,
        "wall_orbit_count": result["wall_orbit_count"] == expected["wall_orbit_count"],
        "unique_restriction_count": (
            result["unique_restriction_count"] == expected["unique_restriction_count"]
        ),
        "restriction_square_sign_histogram": (
            result["restriction_square_sign_histogram"]
            == expected["restriction_square_sign_histogram"]
        ),
        "strict_facet_count": result["strict_facet_count"] == expected["strict_facet_count"],
        "all_negative_restrictions_have_projection_witnesses": (
            result["projection_witness_failure_count"] == 0
        ),
        "positive_wall_gram_sha256": (
            result["positive_wall_gram_sha256"]
            == expected["positive_wall_gram_sha256"]
        ),
        "complete_pair_relation_census": (
            result["pair_gate"].get("coxeter_edge_counts")
            == expected["coxeter_edge_counts"]
        ),
    }


def resource_estimate(case: SurfaceCase) -> dict[str, Any]:
    """Conservative transparent bound before knowing the involution orbit count."""
    count = case.wall_count
    rank = RANK
    # CPython upper model: 40-byte tuple header, 8-byte references, and a
    # conservative 36 bytes for each non-shared integer object.
    bytes_per_vector = 40 + rank * (8 + 36)
    retained_vector_collections = 3  # walls, primitive rays, Gram covectors
    vector_bytes = retained_vector_collections * count * bytes_per_vector
    # Exact byte size of the largest pinned parser input, AutXa.txt.
    source_text_bytes = 5_028_554
    structural_bytes = 8_000_000  # dict/set/list keys, records, parser overhead.
    conservative_incremental = vector_bytes + source_text_bytes + structural_bytes
    return {
        "source_wall_count": count,
        "ambient_rank": rank,
        "maximum_involution_wall_orbits": count,
        "maximum_unique_restrictions": count,
        "maximum_pair_count": count * (count - 1) // 2,
        "dense_source_gram_entry_count_avoided": count * count,
        "maximum_wall_action_integer_muladds": count * rank * rank,
        "maximum_strict_witness_integer_muladds": count * count * rank,
        "maximum_streamed_gram_integer_muladds": count * count * rank,
        "retained_vector_coordinate_slots_upper_bound": (
            retained_vector_collections * count * rank
        ),
        "cpython_vector_storage_model_bytes": vector_bytes,
        "pinned_source_text_buffer_bound_bytes": source_text_bytes,
        "structural_overhead_allowance_bytes": structural_bytes,
        "conservative_incremental_peak_bytes": conservative_incremental,
        "conservative_incremental_peak_mib": round(
            conservative_incremental / (1 << 20), 2
        ),
        "hard_operational_incremental_budget_mib": 64,
        "warning": (
            "This is an allocation model, not a measured RSS. The full X0 run "
            "remains gated until memory pressure clears."
        ),
    }


def _available_free_mib() -> float | None:
    """Return immediately free memory, deliberately excluding compressed pages."""
    try:
        completed = subprocess.run(
            ["vm_stat"], text=True, capture_output=True, check=True, timeout=5
        )
        page_match = re.search(r"page size of (\d+) bytes", completed.stdout)
        free_match = re.search(r"Pages free:\s+(\d+)\.", completed.stdout)
        if page_match and free_match:
            return int(page_match.group(1)) * int(free_match.group(1)) / (1 << 20)
    except (OSError, subprocess.SubprocessError):
        pass
    try:
        raw = Path("/proc/meminfo").read_text(encoding="utf-8")
        match = re.search(r"^MemFree:\s+(\d+) kB", raw, re.MULTILINE)
        if match:
            return int(match.group(1)) / 1024
    except OSError:
        pass
    return None


def build_case(
    case: SurfaceCase,
    source_dir: Path,
    *,
    stop_on_first_forbidden: bool = False,
) -> dict[str, Any]:
    for filename, expected in FILE_SHA256.items():
        path = source_dir / filename
        if _sha256(path) != expected:
            raise ArithmeticError(f"pinned source changed: {path}")

    gram_raw = _assignment(source_dir / "GramS.txt", "GramS", case.surface_index)
    ample_raw = _assignment(source_dir / "ample.txt", "ample", case.surface_index)
    gram = tuple(
        _integral_vector(row, f"GramS[{case.surface_index}] row")
        for row in gram_raw
    )
    ample = _integral_vector(ample_raw, f"ample[{case.surface_index}]")
    walls, denominator = _load_scaled_walls(
        source_dir / "walls.txt", case.surface_index, case.wall_orbit_sizes
    )
    if len(walls) != case.wall_count:
        raise ArithmeticError("complete source wall count changed")
    action = _load_representative(
        source_dir / "AutXa.txt", case.surface_index, case.representative
    )
    result = analyze_scaled(
        gram=gram,
        ample=ample,
        walls=walls,
        action=action,
        global_wall_denominator=denominator,
        stop_on_first_forbidden=stop_on_first_forbidden,
    )
    checks = _expected_checks(case, result)
    if case.expected is not None and not all(checks.values()):
        raise ArithmeticError(f"frozen {case.key} replay changed: {checks}")
    return {
        "schema": SCHEMA,
        "id": f"SHIMADA-AUT3SINGK3-{case.key.upper()}-SPARSE-FIXED-SECTION",
        "case": {
            "key": case.key,
            "surface": f"X_{case.surface_index}",
            "source_group_representative_zero_based": case.representative,
            "involution_type": case.involution_type,
        },
        "denominator_scope": {
            "included": (
                "the complete published induced-chamber wall set restricted by "
                "one representative of the unique stated conjugacy class in the "
                "genuine ample stabilizer"
            ),
            "excluded": (
                "the full abstract oriented weighted-wall graph group, if larger"
            ),
        },
        "algorithm": {
            "arithmetic": "exact Python integers after one global wall denominator",
            "wall_restriction": "oriented primitive sums over involution wall orbits",
            "irredundancy": "exact strict timelike orthogonal-projection witnesses",
            "pair_gate": "streamed exact rational crystallographic Coxeter test",
            "dense_pairing_matrix_materialized": False,
        },
        "resource_estimate": resource_estimate(case),
        "fixed_section": result,
        "frozen_replay_checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=sorted(CASES), required=True)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--estimate-only", action="store_true")
    parser.add_argument("--stop-on-first-forbidden", action="store_true")
    parser.add_argument("--authorize-x0-full", action="store_true")
    parser.add_argument("--minimum-free-mib", type=float, default=512.0)
    parser.add_argument("--measure-resources", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    case = CASES[args.case]
    if args.estimate_only:
        print(json.dumps(resource_estimate(case), indent=2, sort_keys=True))
        return
    if case.key == "x0":
        if not args.authorize_x0_full:
            raise SystemExit(
                "refusing full X0 pass without --authorize-x0-full; use "
                "--estimate-only while the compact search owns memory"
            )
        free_mib = _available_free_mib()
        if free_mib is None or free_mib < args.minimum_free_mib:
            raise SystemExit(
                f"refusing full X0 pass: immediately free memory is {free_mib!r} MiB, "
                f"below the {args.minimum_free_mib:.1f} MiB gate"
            )
    started = time.perf_counter()
    payload = build_case(
        case, args.source_dir, stop_on_first_forbidden=args.stop_on_first_forbidden
    )
    if args.measure_resources:
        raw_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        # Darwin reports bytes; Linux and the BSDs used in CI report KiB.
        max_rss_bytes = int(raw_rss if sys.platform == "darwin" else raw_rss * 1024)
        payload["measured_process_resources"] = {
            "elapsed_seconds": round(time.perf_counter() - started, 6),
            "maximum_resident_set_size_bytes": max_rss_bytes,
            "maximum_resident_set_size_mib": round(max_rss_bytes / (1 << 20), 3),
            "measurement_scope": "whole fresh Python process through certificate construction",
        }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = args.output or (X0_OUTPUT if case.key == "x0" else None)
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded, encoding="utf-8")
        print(json.dumps({
            "output": str(output),
            "classification": payload["fixed_section"]["classification"],
            "all_checks_pass": payload["all_checks_pass"],
        }, indent=2, sort_keys=True))
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()

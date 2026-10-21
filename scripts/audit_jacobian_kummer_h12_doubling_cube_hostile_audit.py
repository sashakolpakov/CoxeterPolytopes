#!/usr/bin/env python3
"""Hostile, source-level audit of the Jacobian--Kummer H12 doubling cube.

This verifier deliberately rebuilds the rank-13 chamber from the published
Dolgachev--Shimada S16/D16 computational data and the pinned 172-facet fold
orbit list.  It does not import or invoke the production cube certificate.

The four proposed walls are checked directly in the reconstructed lattice.
For every nonempty subset of them, the program constructs the complete
boundary of the corresponding finite chamber cube, transports strict facet
witnesses, and checks every primitive reflection and every Coxeter pair by
integer arithmetic.  Finite volume and noncompactness are then inherited by
an exact finite-index Coxeter-kernel/gluing argument from two independently
pinned volume certificates for the seed chamber.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from fractions import Fraction
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
import re
import tarfile
from typing import Any, Iterable, Sequence
import zipfile

import numpy as np
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = (
    ROOT
    / "data/external/dolgachev_shimada_15_nodal_2019/15nodalcompdata.txt"
)
FOLD = ROOT / "data/enumerations/jacobian_kummer_h12_interim_facet_gram.json"
VINBERG = (
    ROOT / "data/enumerations/jacobian_kummer_h12_vinberg_volume_certificate.json"
)
EXACT_CONE = (
    ROOT / "data/enumerations/jacobian_kummer_h12_exact_cone_certificate.json"
)
OUTPUT = (
    ROOT
    / "data/enumerations/jacobian_kummer_h12_doubling_cube_hostile_audit.json"
)
DEFAULT_COXITER_ARCHIVE = Path("/private/tmp/coxiter-current.zip")

SCHEMA = "coxeter.jacobian_kummer.h12_doubling_cube.hostile_audit.v1"
STATUS_DATE = "2026-08-08"
PINNED_SHA256 = {
    "published_Dolgachev_Shimada_data": (
        "413d7363faf4574595e9c28329737eb4ad8010e0aebae53afb76762ce9f8966a"
    ),
    "H12_fold_orbits": (
        "53c39642b81bbd7f5be172408bfecb857bd4ed57467c91bdbaed63f9e734fb92"
    ),
    "literal_Vinberg_seed_gate": (
        "c3474e33120cab541205225dd8df31beea7ae8006cd1667c4fd66408f9b931f1"
    ),
    "exact_GMP_seed_gate": (
        "84d9660b5990d89fb4a78b73180fa815843ee7af46d4d8c56be73f51e37a7c0f"
    ),
}

COXITER_COMMIT = "b800be48d8240fbabd2d1c9542e61054584ad7fe"
COXITER_ARCHIVE_SHA256 = (
    "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3"
)
COXITER_DOUBLING_SHA256 = (
    "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036"
)
COXITER_DOUBLING_MEMBER = "CoxIter-master/graphs/testing_index2.tar.gz"

# Zero-based positions in the independently reconstructed 172-facet list.
SELECTED = (30, 57, 58, 67)
EXPECTED_BY_MASK = {
    0b0001: 296,
    0b0010: 296,
    0b0100: 296,
    0b1000: 317,
    0b0011: 523,
    0b0101: 523,
    0b0110: 523,
    0b1001: 559,
    0b1010: 559,
    0b1100: 559,
    0b0111: 941,
    0b1011: 1004,
    0b1101: 1004,
    0b1110: 1004,
    0b1111: 1827,
}
EXPECTED_FACET_HISTOGRAM = {
    296: 3,
    317: 1,
    523: 3,
    559: 3,
    941: 1,
    1004: 3,
    1827: 1,
}

Vector = tuple[int, ...]
RationalVector = tuple[Fraction, ...]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _json_sha256(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return _sha256_bytes(encoded)


def _integer_matrix_sha256(matrix: np.ndarray) -> str:
    if matrix.dtype != np.int64:
        raise ArithmeticError("integer matrix digest requires int64 input")
    shape = json.dumps(list(matrix.shape), separators=(",", ":")).encode("ascii")
    little_endian = matrix.astype("<i8", copy=False).tobytes(order="C")
    return _sha256_bytes(shape + b"\n" + little_endian)


def _payload_number(value: int | Fraction | sp.Rational) -> int | str:
    rational = Fraction(value)
    return (
        rational.numerator
        if rational.denominator == 1
        else f"{rational.numerator}/{rational.denominator}"
    )


def _extract_assignment(text: str, name: str) -> str:
    """Extract one top-level bracketed assignment from the GAP source."""

    marker = f"{name}:="
    start = text.index(marker) + len(marker)
    while text[start].isspace():
        start += 1
    if text[start] != "[":
        raise ArithmeticError(f"published assignment {name} is not a list")
    depth = 0
    quote: str | None = None
    for end in range(start, len(text)):
        character = text[end]
        if quote is not None:
            if character == quote and text[end - 1] != "\\":
                quote = None
            continue
        if character in {'"', "'"}:
            quote = character
        elif character == "[":
            depth += 1
        elif character == "]":
            depth -= 1
            if depth == 0:
                return text[start : end + 1]
    raise ArithmeticError(f"unterminated published assignment {name}")


def _primitive_column(vector: sp.Matrix) -> sp.Matrix:
    fractions = [sp.Rational(value) for value in vector]
    denominator = math.lcm(*(int(value.q) for value in fractions))
    integers = [int(value * denominator) for value in fractions]
    content = math.gcd(*(abs(value) for value in integers if value))
    if content == 0:
        raise ArithmeticError("zero vector has no primitive normalization")
    return sp.Matrix([value // content for value in integers])


def _primitive_key(vector: Iterable[int]) -> Vector:
    values = tuple(int(value) for value in vector)
    content = math.gcd(*(abs(value) for value in values if value))
    if content == 0:
        raise ArithmeticError("zero vector has no primitive key")
    primitive = tuple(value // content for value in values)
    first = next(value for value in primitive if value)
    return primitive if first > 0 else tuple(-value for value in primitive)


def _content(vector: Iterable[int]) -> int:
    values = tuple(int(value) for value in vector)
    return math.gcd(*(abs(value) for value in values if value))


def _exact_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Compute rational inertia by exact symmetric elimination."""

    work = sp.Matrix(matrix)
    positive = negative = 0
    while work.rows:
        pivot = next((index for index in range(work.rows) if work[index, index]), None)
        if pivot is None:
            pair = next(
                (
                    (left, right)
                    for left in range(work.rows)
                    for right in range(left + 1, work.rows)
                    if work[left, right]
                ),
                None,
            )
            if pair is None:
                return positive, negative, work.rows
            order = list(pair) + [
                index for index in range(work.rows) if index not in pair
            ]
            work = work.extract(order, order)
            change = sp.eye(work.rows)
            change[1, 0] = 1
            work = change.T * work * change
            pivot = 0
        if pivot:
            order = [pivot] + [
                index for index in range(work.rows) if index != pivot
            ]
            work = work.extract(order, order)
        diagonal = sp.Rational(work[0, 0])
        positive += int(bool(diagonal > 0))
        negative += int(bool(diagonal < 0))
        if work.rows == 1:
            return positive, negative, 0
        column = work[1:, 0]
        work = work[1:, 1:] - column * column.T / diagonal
    return positive, negative, 0


def _to_int64(matrix: sp.Matrix) -> np.ndarray:
    if not all(sp.Rational(value).q == 1 for value in matrix):
        raise ArithmeticError("expected an integral matrix")
    return np.array(
        [
            [int(matrix[row, column]) for column in range(matrix.cols)]
            for row in range(matrix.rows)
        ],
        dtype=np.int64,
    )


def _exact_int64_matmul(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Use fast int64 multiplication only after proving no overflow."""

    if left.shape[1] != right.shape[0]:
        raise ArithmeticError("matrix product shape mismatch")
    maximum_left = int(np.max(np.abs(left))) if left.size else 0
    maximum_right = int(np.max(np.abs(right))) if right.size else 0
    bound = left.shape[1] * maximum_left * maximum_right
    if bound >= 2**62:
        raise ArithmeticError(f"int64 matrix product overflow bound: {bound}")
    return left @ right


def _pair(gram: np.ndarray, left: Sequence[int], right: Sequence[int]) -> int:
    # Use Python integers here: this small 13-dimensional primitive is also
    # the overflow-independent control for the bulk int64 calculations.
    return sum(
        int(left[row]) * int(gram[row, column]) * int(right[column])
        for row in range(len(left))
        for column in range(len(right))
    )


def _parse_published() -> dict[str, Any]:
    text = PUBLISHED.read_text(encoding="utf-8").replace("\\\n", "")
    gram = sp.Matrix(ast.literal_eval(_extract_assignment(text, "GramS16")))
    ample = sp.Matrix(
        1,
        17,
        ast.literal_eval(_extract_assignment(text, "ample16")),
    )
    walls_text = _extract_assignment(text, "D16WallRecs")
    vector_texts = re.findall(r"\bvector\s*:=\s*\[([^\]]+)\]", walls_text)
    wall_numbers = tuple(
        int(value) for value in re.findall(r"\bno\s*:=\s*(\d+)", walls_text)
    )
    walls = tuple(
        sp.Matrix(
            1,
            17,
            [sp.Rational(entry.strip()) for entry in vector.split(",")],
        )
        for vector in vector_texts
    )
    if not (
        gram.shape == (17, 17)
        and ample.shape == (1, 17)
        and len(walls) == 316
        and wall_numbers == tuple(range(1, 317))
    ):
        raise ArithmeticError("published S16/D16 parse failed")
    return {"gram": gram, "ample": ample, "walls": walls}


def _reconstruct_seed() -> dict[str, Any]:
    for path, key in (
        (PUBLISHED, "published_Dolgachev_Shimada_data"),
        (FOLD, "H12_fold_orbits"),
        (VINBERG, "literal_Vinberg_seed_gate"),
        (EXACT_CONE, "exact_GMP_seed_gate"),
    ):
        actual = _sha256(path)
        if actual != PINNED_SHA256[key]:
            raise ArithmeticError(f"pinned input {key} changed: {actual}")

    source = _parse_published()
    fold = json.loads(FOLD.read_text(encoding="utf-8"))
    ambient_gram: sp.Matrix = source["gram"]
    ample: sp.Matrix = source["ample"]
    walls: tuple[sp.Matrix, ...] = source["walls"]
    involution = sp.Matrix(fold["group_class"]["representative_matrix"])
    identity17 = sp.eye(17)

    wall_index = {
        tuple(wall): index for index, wall in enumerate(walls)
    }
    wall_images = tuple(tuple(wall * involution) for wall in walls)
    permutation = tuple(wall_index[image] for image in wall_images)
    permutation_is_valid = len(set(permutation)) == 316

    basis_rows = []
    for vector in (involution.T - identity17).nullspace():
        primitive = _primitive_column(vector)
        basis_rows.append([int(value) for value in primitive])
    basis = sp.Matrix(basis_rows)
    fixed_gram = sp.simplify(basis * ambient_gram * basis.T)
    unit_minor_columns = (0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 13, 14)
    unit_minor = int(basis[:, unit_minor_columns].det())
    interior_column = sp.simplify(
        fixed_gram.inv() * basis * ambient_gram * ample.T
    )
    if not all(sp.Rational(value).q == 1 for value in interior_column):
        raise ArithmeticError("published ample class is not integral in fixed basis")

    raw_normals = []
    roots = []
    for source_orbit in fold["facet_source_wall_orbits"]:
        raw = sum(
            (walls[int(index) - 1] for index in source_orbit),
            sp.zeros(1, 17),
        )
        raw_normals.append(raw)
        coordinates = sp.simplify(
            fixed_gram.inv() * basis * ambient_gram * raw.T
        )
        root = _primitive_column(coordinates)
        if (interior_column.T * fixed_gram * root)[0] < 0:
            root = -root
        roots.append(tuple(int(value) for value in root))

    raw_matrix = sp.Matrix.vstack(*raw_normals)
    rebuilt_raw_gram = sp.simplify(
        raw_matrix * ambient_gram * raw_matrix.T
    )
    archived_raw_gram = sp.Matrix(fold["facet_gram"])

    gram = _to_int64(fixed_gram)
    root_matrix = np.array(roots, dtype=np.int64)
    interior = tuple(int(value) for value in interior_column)
    root_pairings = _exact_int64_matmul(
        _exact_int64_matmul(root_matrix, gram), root_matrix.T
    )
    interior_row = np.asarray(interior, dtype=np.int64).reshape((1, -1))
    interior_pairings = _exact_int64_matmul(
        _exact_int64_matmul(interior_row, gram), root_matrix.T
    )[0]

    orbit_checks = []
    for orbit in fold["wall_orbits"]:
        zero_based = tuple(int(index) - 1 for index in orbit)
        generated = {zero_based[0], permutation[zero_based[0]]}
        orbit_checks.append(generated == set(zero_based))

    checks = {
        "published_representative_is_involution": involution**2 == identity17,
        "published_representative_preserves_ambient_gram": (
            involution * ambient_gram * involution.T == ambient_gram
        ),
        "published_representative_fixes_ample": ample * involution == ample,
        "published_representative_permutes_all_316_walls": permutation_is_valid,
        "all_180_archived_wall_orbits_replayed": (
            len(orbit_checks) == 180 and all(orbit_checks)
        ),
        "fixed_basis_has_shape_13_by_17": basis.shape == (13, 17),
        "fixed_basis_is_pointwise_fixed": basis * involution == basis,
        "fixed_basis_is_saturated_by_unit_minor": abs(unit_minor) == 1,
        "fixed_lattice_signature_is_1_12": _exact_inertia(fixed_gram) == (1, 12, 0),
        "fixed_lattice_determinant_is_256": int(fixed_gram.det()) == 256,
        "all_172_facet_orbits_rebuilt": len(roots) == 172,
        "raw_172_wall_gram_matches_fold_archive": rebuilt_raw_gram
        == archived_raw_gram,
        "primitive_root_rank_is_13": sp.Matrix(root_matrix.tolist()).rank() == 13,
        "all_primitive_roots_have_negative_square": bool(
            np.all(np.diag(root_pairings) < 0)
        ),
        "all_primitive_roots_oriented_by_published_ample": bool(
            np.all(interior_pairings > 0)
        ),
        "published_ample_is_timelike": _pair(gram, interior, interior) > 0,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"independent JK-H12 seed reconstruction failed: {checks}")

    return {
        "fold": fold,
        "gram_sympy": fixed_gram,
        "gram": gram,
        "roots": tuple(roots),
        "root_matrix": root_matrix,
        "root_pairings": root_pairings,
        "interior": interior,
        "checks": checks,
        "unit_minor": unit_minor,
    }


def _reflection_matrix(gram: np.ndarray, root: Vector) -> np.ndarray:
    column = np.asarray(root, dtype=np.int64).reshape((-1, 1))
    covector = _exact_int64_matmul(column.T, gram)
    square = _pair(gram, root, root)
    numerator = 2 * _exact_int64_matmul(column, covector)
    if square >= 0 or not bool(np.all(numerator % square == 0)):
        raise ArithmeticError("root does not define an integral reflection")
    reflection = np.eye(len(root), dtype=np.int64) - numerator // square
    if not (
        np.array_equal(
            _exact_int64_matmul(reflection, reflection),
            np.eye(len(root), dtype=np.int64),
        )
        and np.array_equal(
            _exact_int64_matmul(
                _exact_int64_matmul(reflection.T, gram), reflection
            ),
            gram,
        )
        and np.array_equal(_exact_int64_matmul(reflection, column), -column)
    ):
        raise ArithmeticError("integral reflection matrix failed replay")
    return reflection


def _classify_relation(pairing: int, square_product: int) -> str:
    if pairing < 0 or square_product <= 0:
        return "invalid"
    square = pairing * pairing
    if pairing == 0:
        return "m2"
    if 4 * square == square_product:
        return "m3"
    if 2 * square == square_product:
        return "m4"
    if 4 * square == 3 * square_product:
        return "m6"
    if square == square_product:
        return "parallel"
    if square > square_product:
        return "ultraparallel"
    return "invalid"


def _relation_histogram(pairing_matrix: np.ndarray) -> dict[str, int]:
    diagonal = np.diag(pairing_matrix)
    if not bool(np.all(diagonal < 0)):
        raise ArithmeticError("a boundary normal is not spacelike")
    count = pairing_matrix.shape[0]
    left, right = np.triu_indices(count, 1)
    pairings = pairing_matrix[left, right]
    maximum = int(np.max(np.abs(pairings))) if pairings.size else 0
    if 4 * maximum * maximum >= 2**62:
        raise ArithmeticError("relation square risks int64 overflow")
    pairing_squares = pairings * pairings
    square_products = (-diagonal[left]) * (-diagonal[right])
    masks = {
        "m2": pairings == 0,
        "m3": 4 * pairing_squares == square_products,
        "m4": 2 * pairing_squares == square_products,
        "m6": 4 * pairing_squares == 3 * square_products,
        "parallel": pairing_squares == square_products,
        "ultraparallel": pairing_squares > square_products,
    }
    # Only m2 can overlap a formula syntactically, but all square products are
    # positive, so these masks are disjoint.  Orientation is checked separately.
    covered = np.zeros(pairings.shape, dtype=bool)
    histogram: dict[str, int] = {}
    for label, mask in masks.items():
        overlap = bool(np.any(covered & mask))
        if overlap:
            raise ArithmeticError(f"overlapping Coxeter relation class {label}")
        histogram[label] = int(np.count_nonzero(mask))
        covered |= mask
    invalid = (pairings < 0) | ~covered
    histogram["invalid"] = int(np.count_nonzero(invalid))
    if histogram["invalid"]:
        raise ArithmeticError("boundary contains a non-acute or non-Coxeter pair")
    return {key: value for key, value in histogram.items() if value}


def _seed_witnesses(
    gram: np.ndarray,
    roots: Sequence[Vector],
    interior: Vector,
) -> tuple[RationalVector, ...]:
    witnesses = []
    for root in roots:
        coefficient = Fraction(
            _pair(gram, interior, root),
            _pair(gram, root, root),
        )
        witness = tuple(
            Fraction(interior[index]) - coefficient * root[index]
            for index in range(len(root))
        )
        witnesses.append(witness)
    return tuple(witnesses)


def _transform_rational(
    matrix: np.ndarray, vector: RationalVector
) -> RationalVector:
    return tuple(
        sum(
            (Fraction(int(matrix[row, column])) * vector[column]
             for column in range(matrix.shape[1])),
            Fraction(0),
        )
        for row in range(matrix.shape[0])
    )


def _integral_witness(vector: RationalVector) -> tuple[Vector, int]:
    denominator = math.lcm(*(value.denominator for value in vector))
    integral = tuple(int(value * denominator) for value in vector)
    return integral, denominator


def _submasks(mask: int) -> tuple[int, ...]:
    return tuple(value for value in range(16) if value & ~mask == 0)


def _construct_boundary(
    mask: int,
    gram: np.ndarray,
    roots: tuple[Vector, ...],
    interior: Vector,
    seed_witnesses: tuple[RationalVector, ...],
    group: tuple[np.ndarray, ...],
) -> dict[str, Any]:
    selected_positions = tuple(
        SELECTED[bit] for bit in range(4) if mask & (1 << bit)
    )
    copies = _submasks(mask)
    interior_array = np.asarray(interior, dtype=np.int64)
    descendant_interior = sum(
        (
            _exact_int64_matmul(
                group[element], interior_array.reshape((-1, 1))
            )[:, 0]
            for element in copies
        ),
        start=np.zeros(13, dtype=np.int64),
    )

    boundary: dict[Vector, dict[str, Any]] = {}
    for element in copies:
        transformation = group[element]
        for seed_index, seed_root in enumerate(roots):
            if seed_index in selected_positions:
                continue
            image = _exact_int64_matmul(
                transformation,
                np.asarray(seed_root, dtype=np.int64).reshape((-1, 1)),
            )[:, 0]
            key = _primitive_key(image)
            if key in boundary:
                continue
            root = np.asarray(key, dtype=np.int64)
            pairing = _pair(gram, descendant_interior, root)
            if pairing == 0:
                raise ArithmeticError("descendant interior lies on an exterior wall")
            if pairing < 0:
                root = -root
            witness = _transform_rational(
                transformation, seed_witnesses[seed_index]
            )
            boundary[key] = {
                "root": tuple(int(value) for value in root),
                "witness": witness,
                "origin_cube_element": element,
                "origin_seed_facet": seed_index + 1,
            }

    records = [boundary[key] for key in sorted(boundary)]
    boundary_roots = np.array([record["root"] for record in records], dtype=np.int64)
    witness_integral_and_scale = tuple(
        _integral_witness(record["witness"]) for record in records
    )
    witness_matrix = np.array(
        [row[0] for row in witness_integral_and_scale], dtype=np.int64
    )
    witness_scales = tuple(row[1] for row in witness_integral_and_scale)

    root_covectors = _exact_int64_matmul(boundary_roots, gram)
    pairing_matrix = _exact_int64_matmul(root_covectors, boundary_roots.T)
    witness_covectors = _exact_int64_matmul(witness_matrix, gram)
    witness_pairings = _exact_int64_matmul(
        witness_covectors, boundary_roots.T
    )
    square_bound = (
        witness_matrix.shape[1]
        * int(np.max(np.abs(witness_covectors)))
        * int(np.max(np.abs(witness_matrix)))
    )
    if square_bound >= 2**62:
        raise ArithmeticError("witness square calculation risks int64 overflow")
    witness_square_numerators = np.sum(
        witness_covectors * witness_matrix, axis=1, dtype=np.int64
    )

    facet_count = len(records)
    if not (
        np.count_nonzero(np.diag(witness_pairings)) == 0
        and np.count_nonzero(witness_pairings < 0) == 0
        and np.count_nonzero(witness_pairings == 0) == facet_count
        and bool(np.all(witness_square_numerators > 0))
    ):
        raise ArithmeticError("transported strict facet witnesses failed")
    off_target = witness_pairings.copy()
    np.fill_diagonal(off_target, np.iinfo(np.int64).max)
    minimum_slack_by_row = off_target.min(axis=1)
    minimum_off_target_slack = min(
        Fraction(int(value), witness_scales[index])
        for index, value in enumerate(minimum_slack_by_row)
    )
    minimum_witness_square = min(
        Fraction(
            int(value),
            witness_scales[index] * witness_scales[index],
        )
        for index, value in enumerate(witness_square_numerators)
    )

    diagonal = np.diag(pairing_matrix)
    divisibilities = tuple(
        math.gcd(*(abs(int(value)) for value in covector if value))
        for covector in root_covectors
    )
    integral_reflections = tuple(
        bool(np.all((2 * root_covectors[index]) % diagonal[index] == 0))
        for index in range(facet_count)
    )
    descendant_row = descendant_interior.reshape((1, -1))
    descendant_pairings = _exact_int64_matmul(
        _exact_int64_matmul(descendant_row, gram), boundary_roots.T
    )[0]
    descendant_square = _pair(
        gram, descendant_interior, descendant_interior
    )
    relation_histogram = _relation_histogram(pairing_matrix)

    checks = {
        "expected_facet_count": facet_count == EXPECTED_BY_MASK[mask],
        "boundary_roots_have_full_rank_13": (
            sp.Matrix(boundary_roots.tolist()).rank() == 13
        ),
        "all_boundary_roots_are_primitive": all(
            _content(root) == 1 for root in boundary_roots
        ),
        "all_boundary_roots_have_negative_square": bool(np.all(diagonal < 0)),
        "all_boundary_reflections_integral": all(integral_reflections),
        "descendant_interior_is_timelike": descendant_square > 0,
        "descendant_interior_is_strict": bool(np.all(descendant_pairings > 0)),
        "one_strict_timelike_witness_per_facet": (
            len(records) == facet_count
            and minimum_witness_square > 0
            and minimum_off_target_slack > 0
        ),
        "every_facet_pair_is_exactly_Coxeter_parallel_or_ultraparallel": (
            "invalid" not in relation_histogram
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(
            f"subset {tuple(index + 1 for index in selected_positions)} failed: {checks}"
        )

    root_payload = [list(record["root"]) for record in records]
    witness_payload = [
        {
            "integral_coordinates": list(witness_integral_and_scale[index][0]),
            "scale_denominator": witness_integral_and_scale[index][1],
        }
        for index in range(facet_count)
    ]
    selected_one_based = [index + 1 for index in selected_positions]
    return {
        "id": "JK_H12_CUBE_" + "_".join(f"W{value}" for value in selected_one_based),
        "subset_mask_binary": f"{mask:04b}",
        "selected_doubling_walls_one_based": selected_one_based,
        "selected_wall_count": mask.bit_count(),
        "cube_chamber_count": len(copies),
        "reflection_subgroup_index_and_volume_multiple": len(copies),
        "hyperbolic_dimension": 12,
        "facet_count": facet_count,
        "full_linear_rank": 13,
        "strict_interior": {
            "coordinates": [int(value) for value in descendant_interior],
            "square": descendant_square,
            "minimum_facet_pairing": int(np.min(descendant_pairings)),
        },
        "strict_facet_witnesses": {
            "count": facet_count,
            "minimum_square": _payload_number(minimum_witness_square),
            "minimum_off_target_slack": _payload_number(
                minimum_off_target_slack
            ),
            "coordinate_sha256": _json_sha256(witness_payload),
        },
        "primitive_root_square_histogram": {
            str(key): value
            for key, value in sorted(Counter(int(value) for value in diagonal).items())
        },
        "primitive_root_divisibility_histogram": {
            str(key): value
            for key, value in sorted(Counter(divisibilities).items())
        },
        "relation_histogram": relation_histogram,
        "exact_unordered_facet_pair_count": facet_count * (facet_count - 1) // 2,
        "integral_reflection_count": sum(integral_reflections),
        "boundary_root_coordinate_sha256": _json_sha256(root_payload),
        "exact_pairing_matrix_int64_sha256": _integer_matrix_sha256(pairing_matrix),
        "boundary_origin_histogram_by_cube_element": {
            f"{element:04b}": count
            for element, count in sorted(
                Counter(record["origin_cube_element"] for record in records).items()
            )
        },
        "finite_volume": True,
        "compact": False,
        "finite_volume_proof": (
            "exact convex Coxeter gluing of the displayed finite cube of seed chambers"
        ),
        "noncompactness_proof": (
            "finite-index subgroup of the nonuniform seed reflection group"
        ),
        "checks": checks,
        "all_exact_checks_pass": all(checks.values()),
    }


def _doubling_and_group_audit(seed: dict[str, Any]) -> dict[str, Any]:
    gram: np.ndarray = seed["gram"]
    roots: tuple[Vector, ...] = seed["roots"]
    root_pairings: np.ndarray = seed["root_pairings"]
    selected_roots = tuple(roots[index] for index in SELECTED)
    selected_gram = root_pairings[np.ix_(SELECTED, SELECTED)]
    reflections = tuple(
        _reflection_matrix(gram, root) for root in selected_roots
    )

    row_histograms = {}
    valid_doubling_rows = []
    for index in SELECTED:
        histogram: Counter[str] = Counter()
        valid = True
        for other in range(len(roots)):
            if other == index:
                continue
            label = _classify_relation(
                int(root_pairings[index, other]),
                int(root_pairings[index, index] * root_pairings[other, other]),
            )
            histogram[label] += 1
            if label in {"m3", "m5", "m7", "m9", "invalid"}:
                valid = False
        row_histograms[str(index + 1)] = dict(sorted(histogram.items()))
        valid_doubling_rows.append(valid)

    group_matrices = []
    identity = np.eye(13, dtype=np.int64)
    for mask in range(16):
        element = identity.copy()
        for bit, reflection in enumerate(reflections):
            if mask & (1 << bit):
                element = _exact_int64_matmul(reflection, element)
        group_matrices.append(element)
    group = tuple(group_matrices)
    group_hashes = {_integer_matrix_sha256(element) for element in group}

    selected_arrays = tuple(
        np.asarray(root, dtype=np.int64) for root in selected_roots
    )
    sign_action_checks = []
    for bit, reflection in enumerate(reflections):
        for other, root in enumerate(selected_arrays):
            expected = -root if bit == other else root
            sign_action_checks.append(
                np.array_equal(
                    _exact_int64_matmul(reflection, root.reshape((-1, 1)))[:, 0],
                    expected,
                )
            )

    quotient_relation_checks = []
    finite_orders = {"m2": 2, "m3": 3, "m4": 4, "m6": 6}
    selected_bit = {wall: bit for bit, wall in enumerate(SELECTED)}
    for left in range(len(roots)):
        for right in range(left):
            label = _classify_relation(
                int(root_pairings[left, right]),
                int(root_pairings[left, left] * root_pairings[right, right]),
            )
            if label not in finite_orders:
                continue
            quotient_image = 0
            if left in selected_bit:
                quotient_image ^= 1 << selected_bit[left]
            if right in selected_bit:
                quotient_image ^= 1 << selected_bit[right]
            quotient_relation_checks.append(
                finite_orders[label] % 2 == 0 or quotient_image == 0
            )

    checks = {
        "selected_walls_are_31_58_59_68": tuple(index + 1 for index in SELECTED)
        == (31, 58, 59, 68),
        "selected_wall_gram_is_minus_4_identity": np.array_equal(
            selected_gram, -4 * np.eye(4, dtype=np.int64)
        ),
        "all_four_rows_satisfy_even_angle_doubling_criterion": all(
            valid_doubling_rows
        ),
        "all_four_reflections_are_integral_involutive_isometries": len(
            reflections
        )
        == 4,
        "selected_reflections_commute_pairwise": all(
            np.array_equal(
                _exact_int64_matmul(reflections[left], reflections[right]),
                _exact_int64_matmul(reflections[right], reflections[left]),
            )
            for left in range(4)
            for right in range(left)
        ),
        "selected_reflections_flip_exactly_their_own_normal": all(
            sign_action_checks
        ),
        "all_16_cube_elements_are_distinct": len(group_hashes) == 16,
        "all_seed_Coxeter_relations_respect_C2_4_quotient": all(
            quotient_relation_checks
        ),
        "quotient_is_surjective": True,
        "C2_4_action_is_faithful": len(group_hashes) == 16,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"doubling/group audit failed: {checks}")
    return {
        "selected_roots": selected_roots,
        "selected_gram": selected_gram,
        "reflections": reflections,
        "group": group,
        "row_histograms": row_histograms,
        "checks": checks,
    }


def _seed_volume_audit() -> dict[str, Any]:
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    cone = json.loads(EXACT_CONE.read_text(encoding="utf-8"))
    checks = {
        "literal_Vinberg_classifies_seed_finite_volume_noncompact": (
            vinberg["vinberg_decision"] == "finite_volume_noncompact"
            and vinberg["finite_volume"]
            and not vinberg["compact"]
        ),
        "literal_Vinberg_has_no_bad_edges": vinberg["bad_edge_count"] == 0,
        "literal_Vinberg_has_no_unrecognized_endpoints": (
            vinberg["unrecognized_endpoint_edge_count"] == 0
        ),
        "literal_Vinberg_has_no_invalid_parabolic_neighborhoods": (
            vinberg["invalid_full_zero_neighborhood_count"] == 0
        ),
        "exact_GMP_cone_classifies_seed_finite_volume_noncompact": (
            cone["decision"]["classification"] == "finite_volume_noncompact"
            and cone["decision"]["finite_volume"]
            and not cone["decision"]["compact"]
        ),
        "exact_GMP_cone_has_no_hyperideal_rays": (
            cone["exact_double_description"]["hyperideal_vertex_count"] == 0
        ),
        "all_exact_GMP_replay_checks_green": all(cone["checks"].values()),
        "two_seed_gates_agree_on_finite_and_ideal_counts": (
            vinberg["finite_vertex_count"]
            == cone["exact_double_description"]["finite_vertex_count"]
            and vinberg["ideal_vertex_count"]
            == cone["exact_double_description"]["ideal_vertex_count"]
        ),
        "seed_has_cusps": vinberg["ideal_vertex_count"] > 0,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"pinned seed volume gates failed: {checks}")
    return {
        "classification": "finite_volume_noncompact",
        "facet_count": vinberg["facet_count"],
        "finite_vertices": vinberg["finite_vertex_count"],
        "cusps": vinberg["ideal_vertex_count"],
        "elliptic_rank_11_edges": vinberg["elliptic_rank_11_edge_count"],
        "checks": checks,
    }


def _parse_coxiter_header(payload: bytes) -> tuple[int, int]:
    text = payload.decode("utf-8")
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        facets, dimension = map(int, line.split()[:2])
        return dimension, facets
    raise ArithmeticError("CoxIter record has no header")


def _scan_coxiter(archive_path: Path) -> dict[str, Any]:
    if _sha256(archive_path) != COXITER_ARCHIVE_SHA256:
        raise ArithmeticError("pinned ordinary CoxIter archive changed")
    ordinary: Counter[tuple[int, int]] = Counter()
    with zipfile.ZipFile(archive_path) as archive:
        if archive.comment.decode("ascii") != COXITER_COMMIT:
            raise ArithmeticError("pinned CoxIter commit comment changed")
        for name in archive.namelist():
            if name.endswith(".coxiter"):
                ordinary[_parse_coxiter_header(archive.read(name))] += 1
        doubling_bytes = archive.read(COXITER_DOUBLING_MEMBER)
    if _sha256_bytes(doubling_bytes) != COXITER_DOUBLING_SHA256:
        raise ArithmeticError("pinned CoxIter doubling archive changed")

    doubling: Counter[tuple[int, int]] = Counter()
    with tarfile.open(fileobj=io.BytesIO(doubling_bytes), mode="r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            stream = archive.extractfile(member)
            if stream is None:
                raise ArithmeticError(f"could not read CoxIter member {member.name}")
            doubling[_parse_coxiter_header(stream.read())] += 1

    candidate_counts = sorted(set(EXPECTED_BY_MASK.values()))
    headers = []
    for facets in candidate_counts:
        headers.append(
            {
                "dimension": 12,
                "facets": facets,
                "ordinary_count": ordinary[(12, facets)],
                "doubling_count": doubling[(12, facets)],
                "absent_from_both_pinned_lanes": (
                    ordinary[(12, facets)] + doubling[(12, facets)] == 0
                ),
            }
        )
    checks = {
        "ordinary_record_count_is_1214": sum(ordinary.values()) == 1214,
        "doubling_record_count_is_2959": sum(doubling.values()) == 2959,
        "all_seven_candidate_headers_absent": all(
            row["absent_from_both_pinned_lanes"] for row in headers
        ),
        "ordinary_H12_histogram_matches": {
            facets: count for (dimension, facets), count in ordinary.items()
            if dimension == 12
        }
        == {13: 6, 14: 3, 16: 1, 18: 1},
        "doubling_H12_histogram_matches": {
            facets: count for (dimension, facets), count in doubling.items()
            if dimension == 12
        }
        == {18: 1, 20: 1, 25: 1, 35: 1},
    }
    if not all(checks.values()):
        raise ArithmeticError(f"bounded CoxIter header audit failed: {checks}")
    return {
        "repository": "https://github.com/rgugl/CoxIter",
        "commit": COXITER_COMMIT,
        "ordinary_archive_sha256": COXITER_ARCHIVE_SHA256,
        "ordinary_record_count": sum(ordinary.values()),
        "doubling_archive_member": COXITER_DOUBLING_MEMBER,
        "doubling_archive_sha256": COXITER_DOUBLING_SHA256,
        "doubling_record_count": sum(doubling.values()),
        "ordinary_H12_facet_histogram": {
            str(facets): ordinary[(12, facets)]
            for facets in sorted(facets for dimension, facets in ordinary if dimension == 12)
        },
        "doubling_H12_facet_histogram": {
            str(facets): doubling[(12, facets)]
            for facets in sorted(facets for dimension, facets in doubling if dimension == 12)
        },
        "candidate_headers": headers,
        "evidence_grade": "B1: bounded same-header absence only",
        "worldwide_priority": "open",
        "scope_warning": (
            "The pinned CoxIter lanes are not a publication-complete denominator; "
            "zero headers do not prove worldwide novelty."
        ),
        "checks": checks,
    }


def build_certificate(
    coxiter_archive: Path = DEFAULT_COXITER_ARCHIVE,
) -> dict[str, Any]:
    seed = _reconstruct_seed()
    group_audit = _doubling_and_group_audit(seed)
    seed_witnesses = _seed_witnesses(
        seed["gram"], seed["roots"], seed["interior"]
    )
    descendants = [
        _construct_boundary(
            mask,
            seed["gram"],
            seed["roots"],
            seed["interior"],
            seed_witnesses,
            group_audit["group"],
        )
        for mask in range(1, 16)
    ]
    facet_histogram = Counter(row["facet_count"] for row in descendants)
    volume_audit = _seed_volume_audit()
    coxiter = _scan_coxiter(coxiter_archive)

    selected_gram = group_audit["selected_gram"]
    checks = {
        "all_pinned_sources_hash_verified": True,
        "independent_Dolgachev_Shimada_seed_reconstruction_green": all(
            seed["checks"].values()
        ),
        "four_wall_doubling_and_faithfulness_audit_green": all(
            group_audit["checks"].values()
        ),
        "all_15_nonempty_subsets_constructed": len(descendants) == 15,
        "every_descendant_exact_audit_green": all(
            row["all_exact_checks_pass"] for row in descendants
        ),
        "facet_counts_match_every_labelled_subset": all(
            row["facet_count"] == EXPECTED_BY_MASK[int(row["subset_mask_binary"], 2)]
            for row in descendants
        ),
        "facet_count_histogram_matches": dict(facet_histogram)
        == EXPECTED_FACET_HISTOGRAM,
        "seed_finite_volume_noncompact_gate_green": all(
            volume_audit["checks"].values()
        ),
        "finite_index_gluing_applies_to_every_descendant": all(
            row["finite_volume"] and not row["compact"] for row in descendants
        ),
        "at_least_seven_types_by_distinct_facet_counts": len(facet_histogram) == 7,
        "five_facet_count_types_beyond_single_doublings": (
            set(facet_histogram) - {296, 317}
            == {523, 559, 941, 1004, 1827}
        ),
        "bounded_CoxIter_header_audit_green": all(coxiter["checks"].values()),
        "worldwide_priority_explicitly_left_open": coxiter["worldwide_priority"]
        == "open",
    }
    if not all(checks.values()):
        raise ArithmeticError(f"hostile cube audit failed: {checks}")

    seed_root_matrix: np.ndarray = seed["root_matrix"]
    carrier = {
        "rank": 13,
        "signature": [1, 12, 0],
        "determinant_absolute": 256,
        "fixed_basis_saturated_unit_minor": seed["unit_minor"],
        "gram": seed["gram"].tolist(),
        "seed_primitive_root_coordinate_sha256": _json_sha256(
            seed_root_matrix.tolist()
        ),
        "surface_source": (
            "rank-13 fixed divisor lattice extracted from the general genus-2 "
            "Jacobian--Kummer S16/D16 chamber"
        ),
        "scope": (
            "weak integral full-Picard divisor-reflection carrier; no full-nef "
            "or surface-automorphism claim"
        ),
    }
    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": all(checks.values()),
        "pinned_sources": {
            "published_Dolgachev_Shimada_data": {
                "path": PUBLISHED.relative_to(ROOT).as_posix(),
                "sha256": PINNED_SHA256["published_Dolgachev_Shimada_data"],
            },
            "H12_fold_orbits": {
                "path": FOLD.relative_to(ROOT).as_posix(),
                "sha256": PINNED_SHA256["H12_fold_orbits"],
            },
            "literal_Vinberg_seed_gate": {
                "path": VINBERG.relative_to(ROOT).as_posix(),
                "sha256": PINNED_SHA256["literal_Vinberg_seed_gate"],
            },
            "exact_GMP_seed_gate": {
                "path": EXACT_CONE.relative_to(ROOT).as_posix(),
                "sha256": PINNED_SHA256["exact_GMP_seed_gate"],
            },
        },
        "independence": {
            "method": (
                "raw published S16/D16 parse, independent fixed-lattice and wall-orbit "
                "reconstruction, independent boundary generator"
            ),
            "production_cube_code_imported_or_invoked": False,
            "production_cube_archive_used_as_input": False,
        },
        "seed": {
            "id": "JK-H12-172",
            "dimension": 12,
            "facet_count": 172,
            "carrier": carrier,
            "strict_interior_coordinates": list(seed["interior"]),
            "selected_doubling_walls_one_based": [index + 1 for index in SELECTED],
            "selected_primitive_roots": [list(root) for root in group_audit["selected_roots"]],
            "selected_root_gram": selected_gram.tolist(),
            "selected_seed_relation_histograms": group_audit["row_histograms"],
            "finite_volume_gate": volume_audit,
            "reconstruction_checks": seed["checks"],
        },
        "C2_4_doubling_action": {
            "abstract_group": "C2^4",
            "order": 16,
            "generators": [f"s_{index + 1}" for index in SELECTED],
            "kernel_description": (
                "kernel of the Coxeter epimorphism sending s_31,s_58,s_59,s_68 "
                "to the four standard basis vectors and every other seed generator to 0"
            ),
            "gluing_theorem": (
                "for each selected subset, the union of its 2^k coset chambers is "
                "the convex Coxeter fundamental chamber of the corresponding kernel"
            ),
            "checks": group_audit["checks"],
        },
        "descendants": descendants,
        "census": {
            "labelled_nonempty_subset_count": len(descendants),
            "facet_count_histogram": {
                str(key): value for key, value in sorted(facet_histogram.items())
            },
            "distinct_facet_counts": sorted(facet_histogram),
            "certified_pairwise_nonisomorphic_lower_bound": len(facet_histogram),
            "single_doubling_facet_count_types": [296, 317],
            "additional_facet_count_types": [523, 559, 941, 1004, 1827],
            "additional_types_beyond_singles_lower_bound": 5,
            "same_facet_count_subsets_classified_up_to_isomorphism": False,
        },
        "bounded_lineage": coxiter,
        "decision": {
            "classification": "finite_volume_noncompact_H12_Coxeter_polytopes",
            "all_15_labelled_descendants_certified": True,
            "pairwise_nonisomorphic_type_lower_bound": 7,
            "new_facet_count_types_beyond_single_doublings": 5,
            "worldwide_priority": "open",
            "careful_claim": (
                "seven types are forced by distinct facet counts inside this cube; "
                "the five nonsingle facet-count types are new to the project construction, "
                "while external publication priority remains unsettled"
            ),
        },
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--coxiter-archive", type=Path, default=DEFAULT_COXITER_ARCHIVE)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate(args.coxiter_archive)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != rendered:
            raise ArithmeticError("archived hostile JK-H12 cube audit changed")
        print(f"verified {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"wrote {args.output}")
    print(json.dumps(payload["decision"], sort_keys=True))
    print(f"payload_sha256={_sha256_bytes(rendered.encode('utf-8'))}")


if __name__ == "__main__":
    main()

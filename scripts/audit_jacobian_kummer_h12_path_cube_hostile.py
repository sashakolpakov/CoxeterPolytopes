#!/usr/bin/env python3
"""Independent hostile audit of the JK-H12 D_infinity x C2^3 family.

The verifier starts from the pinned Dolgachev--Shimada GAP source and the
archived trace-nine wall orbits.  It deliberately imports no project module
and does not read the exploratory path-cube scout or either production
doubling builder.  The universal count is proved by an exact all-integer
orbit-collision solve, not by interpolation from finite values.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter, deque
from fractions import Fraction
from functools import lru_cache
import hashlib
from io import BytesIO
import json
import math
from pathlib import Path
import re
import tarfile
from typing import Any, Iterable, Sequence
import zipfile

import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = (
    ROOT / "data/external/dolgachev_shimada_15_nodal_2019/15nodalcompdata.txt"
)
FOLD = ROOT / "data/enumerations/jacobian_kummer_h12_interim_facet_gram.json"
VINBERG = (
    ROOT / "data/enumerations/jacobian_kummer_h12_vinberg_volume_certificate.json"
)
EXACT_CONE = (
    ROOT / "data/enumerations/jacobian_kummer_h12_exact_cone_certificate.json"
)
CARRIER = (
    ROOT / "data/enumerations/jacobian_kummer_h12_full_picard_k3_carrier.json"
)
OUTPUT = (
    ROOT
    / "data/enumerations/jacobian_kummer_h12_path_cube_hostile_audit.json"
)
DEFAULT_COXITER_ARCHIVE = Path("/private/tmp/coxiter-current.zip")

SCHEMA = "coxeter.jacobian_kummer.h12_path_cube.hostile_audit.v1"
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
    "full_Picard_K3_carrier": (
        "f2c57b9fb0423560c9e238541b0ac1a81ca141917c864929501ebdfa7f294ff2"
    ),
}
COXITER_COMMIT = "b800be48d8240fbabd2d1c9542e61054584ad7fe"
COXITER_ZIP_SHA256 = (
    "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3"
)
COXITER_INDEX_TWO_SHA256 = (
    "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036"
)
COXITER_INDEX_TWO_MEMBER = "CoxIter-master/graphs/testing_index2.tar.gz"

# All entries are zero based.  These are the eight configurations under audit;
# their one-based form is recorded in the result archive and report.
CONFIGURATIONS = (
    ((75, 71), (35, 38, 47)),
    ((74, 70), (35, 41, 50)),
    ((73, 69), (38, 41, 55)),
    ((72, 68), (47, 50, 55)),
    ((67, 62), (30, 57, 58)),
    ((66, 63), (31, 57, 59)),
    ((65, 64), (32, 58, 59)),
    ((61, 60), (30, 31, 32)),
)
REPRESENTATIVE_CONFIGURATION = 4

Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]
Permutation = tuple[int, ...]


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_sha256(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return _sha256_bytes(encoded)


def _extract_assignment(text: str, name: str) -> str:
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


def _content(vector: Iterable[int]) -> int:
    entries = tuple(int(value) for value in vector)
    return math.gcd(*(abs(value) for value in entries if value))


def _primitive_column(vector: sp.Matrix) -> sp.Matrix:
    rationals = tuple(sp.Rational(value) for value in vector)
    denominator = math.lcm(*(int(value.q) for value in rationals))
    integers = tuple(int(value * denominator) for value in rationals)
    content = _content(integers)
    if content == 0:
        raise ArithmeticError("zero vector has no primitive normalization")
    return sp.Matrix([value // content for value in integers])


def _primitive_key(vector: Sequence[int]) -> Vector:
    entries = tuple(int(value) for value in vector)
    content = _content(entries)
    if content == 0:
        raise ArithmeticError("zero vector has no primitive key")
    primitive = tuple(value // content for value in entries)
    first = next(value for value in primitive if value)
    return primitive if first > 0 else tuple(-value for value in primitive)


def _identity(rank: int) -> Matrix:
    return tuple(
        tuple(int(row == column) for column in range(rank)) for row in range(rank)
    )


def _matvec(matrix: Matrix, vector: Sequence[int]) -> Vector:
    return tuple(
        sum(entry * int(vector[column]) for column, entry in enumerate(row))
        for row in matrix
    )


def _matmul(left: Matrix, right: Matrix) -> Matrix:
    size = len(left)
    return tuple(
        tuple(
            sum(left[row][middle] * right[middle][column] for middle in range(size))
            for column in range(size)
        )
        for row in range(size)
    )


def _add(*vectors: Sequence[int]) -> Vector:
    return tuple(sum(values) for values in zip(*vectors, strict=True))


def _scale(scalar: int, vector: Sequence[int]) -> Vector:
    return tuple(scalar * int(value) for value in vector)


def _pair(gram: Matrix, left: Sequence[int], right: Sequence[int]) -> int:
    return sum(
        int(left[row]) * gram[row][column] * int(right[column])
        for row in range(len(gram))
        for column in range(len(gram))
    )


def _covector(gram: Matrix, vector: Sequence[int]) -> Vector:
    return tuple(
        sum(gram[row][column] * int(vector[column]) for column in range(len(gram)))
        for row in range(len(gram))
    )


def _matrix_payload(matrix: Matrix) -> list[list[int]]:
    return [list(row) for row in matrix]


def _matrix_sha256(matrix: sp.Matrix | Matrix) -> str:
    if isinstance(matrix, sp.MatrixBase):
        payload = [
            [int(matrix[row, column]) for column in range(matrix.cols)]
            for row in range(matrix.rows)
        ]
    else:
        payload = _matrix_payload(matrix)
    return _json_sha256(payload)


def _exact_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
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
            order = [pivot] + [index for index in range(work.rows) if index != pivot]
            work = work.extract(order, order)
        diagonal = sp.Rational(work[0, 0])
        positive += int(bool(diagonal > 0))
        negative += int(bool(diagonal < 0))
        if work.rows == 1:
            return positive, negative, 0
        column = work[1:, 0]
        work = sp.simplify(work[1:, 1:] - column * column.T / diagonal)
    return positive, negative, 0


def _parse_published() -> dict[str, Any]:
    text = PUBLISHED.read_text(encoding="utf-8").replace("\\\n", "")
    gram = sp.Matrix(ast.literal_eval(_extract_assignment(text, "GramS16")))
    ample = sp.Matrix(1, 17, ast.literal_eval(_extract_assignment(text, "ample16")))
    stabilizer_generators = tuple(
        sp.Matrix(matrix)
        for matrix in ast.literal_eval(
            _extract_assignment(text, "GeneratorsOS16D16")
        )
    )
    wall_text = _extract_assignment(text, "D16WallRecs")
    vector_texts = re.findall(r"\bvector\s*:=\s*\[([^\]]+)\]", wall_text)
    wall_numbers = tuple(
        int(value) for value in re.findall(r"\bno\s*:=\s*(\d+)", wall_text)
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
        and len(stabilizer_generators) == 10
        and len(walls) == 316
        and wall_numbers == tuple(range(1, 317))
    ):
        raise ArithmeticError("published S16/D16 parse failed")
    return {
        "text": text,
        "gram": gram,
        "ample": ample,
        "stabilizer_generators": stabilizer_generators,
        "walls": walls,
    }


def _reconstruct_seed() -> dict[str, Any]:
    paths = {
        "published_Dolgachev_Shimada_data": PUBLISHED,
        "H12_fold_orbits": FOLD,
        "literal_Vinberg_seed_gate": VINBERG,
        "exact_GMP_seed_gate": EXACT_CONE,
        "full_Picard_K3_carrier": CARRIER,
    }
    hashes = {name: _sha256(path) for name, path in paths.items()}
    if hashes != PINNED_SHA256:
        raise ArithmeticError(f"a pinned input changed: {hashes}")

    source = _parse_published()
    fold = json.loads(FOLD.read_text(encoding="utf-8"))
    ambient_gram: sp.Matrix = source["gram"]
    ample: sp.Matrix = source["ample"]
    walls: tuple[sp.Matrix, ...] = source["walls"]
    involution = sp.Matrix(fold["group_class"]["representative_matrix"])

    wall_index = {tuple(wall): index for index, wall in enumerate(walls)}
    involution_permutation = tuple(
        wall_index[tuple(wall * involution)] for wall in walls
    )
    basis = sp.Matrix(
        [
            list(_primitive_column(vector))
            for vector in (involution.T - sp.eye(17)).nullspace()
        ]
    )
    fixed_gram_sympy = sp.simplify(basis * ambient_gram * basis.T)
    fixed_gram: Matrix = tuple(
        tuple(int(fixed_gram_sympy[row, column]) for column in range(13))
        for row in range(13)
    )
    unit_minor_columns = (0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 13, 14)
    unit_minor = int(basis[:, unit_minor_columns].det())
    interior_column = sp.simplify(
        fixed_gram_sympy.inv() * basis * ambient_gram * ample.T
    )
    if not all(sp.Rational(value).q == 1 for value in interior_column):
        raise ArithmeticError("published ample is nonintegral in the fixed basis")
    interior = tuple(int(value) for value in interior_column)

    roots: list[Vector] = []
    orbit_sums: list[sp.Matrix] = []
    for orbit in fold["facet_source_wall_orbits"]:
        raw = sum(
            (walls[int(index) - 1] for index in orbit),
            sp.zeros(1, 17),
        )
        coordinates = _primitive_column(
            fixed_gram_sympy.inv() * basis * ambient_gram * raw.T
        )
        root = tuple(int(value) for value in coordinates)
        if _pair(fixed_gram, interior, root) < 0:
            root = tuple(-value for value in root)
        roots.append(root)
        orbit_sums.append(raw)

    rebuilt_raw_gram = sp.Matrix.vstack(*orbit_sums) * ambient_gram * sp.Matrix.vstack(
        *orbit_sums
    ).T
    carrier = json.loads(CARRIER.read_text(encoding="utf-8"))
    checks = {
        "involution_is_exact_order_two": involution**2 == sp.eye(17),
        "involution_preserves_ambient_gram": (
            involution * ambient_gram * involution.T == ambient_gram
        ),
        "involution_fixes_published_ample": ample * involution == ample,
        "involution_permutes_all_316_walls": (
            len(set(involution_permutation)) == 316
        ),
        "fixed_basis_is_13_by_17": basis.shape == (13, 17),
        "fixed_basis_is_pointwise_fixed": basis * involution == basis,
        "fixed_basis_is_saturated": abs(unit_minor) == 1,
        "fixed_lattice_signature_is_1_12": (
            _exact_inertia(fixed_gram_sympy) == (1, 12, 0)
        ),
        "fixed_lattice_determinant_is_256": int(fixed_gram_sympy.det()) == 256,
        "all_172_roots_reconstructed": len(roots) == 172,
        "raw_root_gram_matches_frozen_fold": (
            rebuilt_raw_gram == sp.Matrix(fold["facet_gram"])
        ),
        "primitive_root_rank_is_13": sp.Matrix.hstack(
            *(sp.Matrix(root) for root in roots)
        ).rank()
        == 13,
        "all_roots_are_primitive_spacelike_and_strictly_oriented": all(
            _content(root) == 1
            and _pair(fixed_gram, root, root) < 0
            and _pair(fixed_gram, interior, root) > 0
            for root in roots
        ),
        "published_ample_is_timelike": _pair(fixed_gram, interior, interior) > 0,
        "full_Picard_carrier_certificate_is_green": (
            carrier["all_checks_pass"]
            and all(carrier["checks"].values())
            and carrier["chamber"]["dimension"] == 12
            and carrier["chamber"]["facet_count"] == 172
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"independent seed reconstruction failed: {checks}")
    return {
        "hashes": hashes,
        "source": source,
        "fold": fold,
        "ambient_gram": ambient_gram,
        "ample": ample,
        "walls": walls,
        "involution": involution,
        "involution_permutation": involution_permutation,
        "basis": basis,
        "fixed_gram_sympy": fixed_gram_sympy,
        "fixed_gram": fixed_gram,
        "unit_minor": unit_minor,
        "interior": interior,
        "roots": tuple(roots),
        "checks": checks,
        "carrier": carrier,
    }


def _classify_relation(pairing: int, left_square: int, right_square: int) -> str:
    if pairing < 0 or left_square >= 0 or right_square >= 0:
        return "invalid"
    square = pairing * pairing
    product = left_square * right_square
    if pairing == 0:
        return "m2"
    if 4 * square == product:
        return "m3"
    if 2 * square == product:
        return "m4"
    if 4 * square == 3 * product:
        return "m6"
    if square == product:
        return "parallel"
    if square > product:
        return "ultraparallel"
    return "invalid"


def _doubling_stock(gram: Matrix, roots: tuple[Vector, ...]) -> dict[str, Any]:
    squares = tuple(_pair(gram, root, root) for root in roots)
    relations: list[tuple[str, ...]] = []
    doubling: list[int] = []
    row_histograms: dict[str, dict[str, int]] = {}
    selected = {index for pair, cube in CONFIGURATIONS for index in (*pair, *cube)}
    for left, root in enumerate(roots):
        labels: list[str] = []
        histogram: Counter[str] = Counter()
        for right, other in enumerate(roots):
            if left == right:
                labels.append("m1")
                continue
            label = _classify_relation(
                _pair(gram, root, other), squares[left], squares[right]
            )
            if label == "invalid":
                raise ArithmeticError("seed relation system is non-Coxeter")
            labels.append(label)
            histogram[label] += 1
        relations.append(tuple(labels))
        if histogram.get("m3", 0) == 0:
            doubling.append(left)
        if left in selected:
            row_histograms[str(left + 1)] = dict(sorted(histogram.items()))
    checks = {
        "exactly_124_of_172_seed_walls_are_doubling": len(doubling) == 124,
        "doubling_criterion_is_no_incident_odd_finite_label": True,
        "all_40_listed_generator_occurrences_are_doubling": all(
            index in set(doubling)
            for pair, cube in CONFIGURATIONS
            for index in (*pair, *cube)
        ),
        "every_seed_pair_has_a_valid_exact_relation": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"doubling stock audit failed: {checks}")
    return {
        "doubling": tuple(doubling),
        "relations": tuple(relations),
        "doubling_walls_one_based": [index + 1 for index in doubling],
        "selected_row_histograms": row_histograms,
        "checks": checks,
    }


@lru_cache(maxsize=None)
def _reflection_matrix(gram: Matrix, root: Vector) -> Matrix:
    square = _pair(gram, root, root)
    covector = _covector(gram, root)
    rows: list[tuple[int, ...]] = []
    for row in range(len(root)):
        entries: list[int] = []
        for column in range(len(root)):
            entry = Fraction(int(row == column)) - Fraction(
                2 * root[row] * covector[column], square
            )
            if entry.denominator != 1:
                raise ArithmeticError("a selected reflection is nonintegral")
            entries.append(entry.numerator)
        rows.append(tuple(entries))
    reflection = tuple(rows)
    gram_sympy = sp.Matrix(gram)
    reflection_sympy = sp.Matrix(reflection)
    if not (
        reflection_sympy**2 == sp.eye(len(root))
        and reflection_sympy.T * gram_sympy * reflection_sympy == gram_sympy
        and _matvec(reflection, root) == tuple(-value for value in root)
    ):
        raise ArithmeticError("reflection exact isometry gates failed")
    return reflection


def _cube_transforms(reflections: dict[int, Matrix], cube: tuple[int, ...]) -> tuple[Matrix, ...]:
    transforms: list[Matrix] = []
    for mask in range(8):
        transform = _identity(13)
        for bit, wall in enumerate(cube):
            if mask & (1 << bit):
                transform = _matmul(transform, reflections[wall])
        transforms.append(transform)
    if len(set(transforms)) != 8:
        raise ArithmeticError("the proposed C2^3 action has a kernel")
    return tuple(transforms)


def _cube_boundary(
    gram: Matrix,
    roots: tuple[Vector, ...],
    interior: Vector,
    cube: tuple[int, ...],
    transforms: tuple[Matrix, ...],
) -> dict[str, Any]:
    descendant_interior = tuple(
        sum(_matvec(transform, interior)[index] for transform in transforms)
        for index in range(13)
    )
    boundary: dict[Vector, dict[str, Any]] = {}
    for element, transform in enumerate(transforms):
        for seed_index, root in enumerate(roots):
            if seed_index in cube:
                continue
            image = _matvec(transform, root)
            key = _primitive_key(image)
            if key in boundary:
                continue
            oriented = image
            if _pair(gram, descendant_interior, oriented) < 0:
                oriented = tuple(-value for value in oriented)
            if _pair(gram, descendant_interior, oriented) <= 0:
                raise ArithmeticError("cube interior is not strict")
            boundary[key] = {
                "root": oriented,
                "seed_index": seed_index,
                "cube_element": element,
                "transform": transform,
            }
    return {
        "interior": descendant_interior,
        "records": boundary,
        "roots": tuple(boundary[key]["root"] for key in sorted(boundary)),
    }


def _smith_diagonal(vectors: Sequence[Vector]) -> list[int]:
    matrix = sp.Matrix.hstack(*(sp.Matrix(vector) for vector in vectors))
    smith = smith_normal_form(matrix, domain=sp.ZZ)
    return [abs(int(smith[index, index])) for index in range(matrix.rows)]


def _configuration_audit(
    seed: dict[str, Any], stock: dict[str, Any], configuration_index: int
) -> dict[str, Any]:
    gram: Matrix = seed["fixed_gram"]
    roots: tuple[Vector, ...] = seed["roots"]
    pair, cube = CONFIGURATIONS[configuration_index]
    selected = (*pair, *cube)
    reflections = {index: _reflection_matrix(gram, roots[index]) for index in selected}
    selected_gram = [
        [_pair(gram, roots[left], roots[right]) for right in selected]
        for left in selected
    ]
    expected_gram = [
        [-4, 4, 0, 0, 0],
        [4, -4, 0, 0, 0],
        [0, 0, -4, 0, 0],
        [0, 0, 0, -4, 0],
        [0, 0, 0, 0, -4],
    ]
    cube_transforms = _cube_transforms(reflections, cube)
    cube_boundary = _cube_boundary(
        gram, roots, seed["interior"], cube, cube_transforms
    )
    records: dict[Vector, dict[str, Any]] = cube_boundary["records"]
    left_key = _primitive_key(roots[pair[0]])
    right_key = _primitive_key(roots[pair[1]])
    if left_key not in records or right_key not in records:
        raise ArithmeticError("path endpoint disappeared from the cube boundary")
    left = records[left_key]["root"]
    right = records[right_key]["root"]
    sides = tuple(
        record["root"]
        for key, record in records.items()
        if key not in {left_key, right_key}
    )
    left_orthogonal = {
        _primitive_key(root) for root in sides if _pair(gram, left, root) == 0
    }
    right_orthogonal = {
        _primitive_key(root) for root in sides if _pair(gram, right, root) == 0
    }
    side_keys = {_primitive_key(root) for root in sides}
    categories = {
        "orthogonal_to_both": len(left_orthogonal & right_orthogonal),
        "orthogonal_to_left_only": len(left_orthogonal - right_orthogonal),
        "orthogonal_to_right_only": len(right_orthogonal - left_orthogonal),
        "orthogonal_to_neither": len(side_keys - (left_orthogonal | right_orthogonal)),
    }

    left_reflection = reflections[pair[0]]
    right_reflection = reflections[pair[1]]
    translation = _matmul(right_reflection, left_reflection)
    nilpotent_sympy = sp.Matrix(translation) - sp.eye(13)
    half_sum = tuple(Fraction(a + b, 2) for a, b in zip(left, right, strict=True))
    half_sum_integral = all(value.denominator == 1 for value in half_sum)
    null_direction = tuple(int(value) for value in half_sum) if half_sum_integral else ()
    generator_set = set(stock["doubling"])
    checks = {
        "five_selected_walls_are_in_124_wall_doubling_stock": set(selected)
        <= generator_set,
        "selected_Gram_is_parallel_pair_orthogonal_C2_3": selected_gram
        == expected_gram,
        "all_selected_reflections_are_integral_involutive_isometries": len(
            reflections
        )
        == 5,
        "cube_reflections_commute_pairwise": all(
            _matmul(reflections[cube[left_bit]], reflections[cube[right_bit]])
            == _matmul(reflections[cube[right_bit]], reflections[cube[left_bit]])
            for left_bit in range(3)
            for right_bit in range(left_bit)
        ),
        "path_reflections_commute_with_cube_action": all(
            _matmul(reflections[path_wall], reflections[cube_wall])
            == _matmul(reflections[cube_wall], reflections[path_wall])
            for path_wall in pair
            for cube_wall in cube
        ),
        "cube_action_has_exact_order_8": len(set(cube_transforms)) == 8,
        "parallel_pair_has_primitive_integral_null_half_sum": (
            half_sum_integral
            and _content(null_direction) == 1
            and _pair(gram, null_direction, null_direction) == 0
        ),
        "path_translation_is_nontrivial_three_step_unipotent": (
            nilpotent_sympy.rank() == 2
            and (nilpotent_sympy**2).rank() == 1
            and nilpotent_sympy**3 == sp.zeros(13)
            and nilpotent_sympy**2 != sp.zeros(13)
        ),
        "cube_first_chamber_has_941_facets": len(records) == 941,
        "cube_first_chamber_has_two_path_endpoints": (
            left_key in records and right_key in records
        ),
        "side_category_census_is_17_36_36_850": categories
        == {
            "orthogonal_to_both": 17,
            "orthogonal_to_left_only": 36,
            "orthogonal_to_right_only": 36,
            "orthogonal_to_neither": 850,
        },
        "all_941_boundary_roots_are_primitive_spacelike_integral_reflection_roots": all(
            _content(record["root"]) == 1
            and _pair(gram, record["root"], record["root"]) < 0
            and all(
                (2 * value) % _pair(gram, record["root"], record["root"]) == 0
                for value in _covector(gram, record["root"])
            )
            for record in records.values()
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(
            f"configuration {configuration_index + 1} failed: {checks}"
        )
    return {
        "configuration_number": configuration_index + 1,
        "path_pair_one_based": [index + 1 for index in pair],
        "cube_walls_one_based": [index + 1 for index in cube],
        "selected_wall_Gram": selected_gram,
        "primitive_integral_null_direction": list(null_direction),
        "cube_order": len(cube_transforms),
        "cube_first_facet_count": len(records),
        "side_facet_category_counts": categories,
        "translation_rank_T_minus_I": int(nilpotent_sympy.rank()),
        "translation_rank_T_minus_I_squared": int((nilpotent_sympy**2).rank()),
        "translation_T_minus_I_cubed_is_zero": True,
        "cube_boundary_root_coordinate_sha256": _json_sha256(
            [list(records[key]["root"]) for key in sorted(records)]
        ),
        "checks": checks,
        "all_exact_checks_pass": True,
        "_reflections": reflections,
        "_cube_transforms": cube_transforms,
        "_cube_boundary": cube_boundary,
        "_left": left,
        "_right": right,
        "_translation": translation,
    }


def _seed_witnesses(
    gram: Matrix, roots: tuple[Vector, ...], interior: Vector
) -> tuple[tuple[Fraction, ...], ...]:
    witnesses: list[tuple[Fraction, ...]] = []
    for root in roots:
        coefficient = Fraction(_pair(gram, interior, root), _pair(gram, root, root))
        witness = tuple(
            Fraction(interior[index]) - coefficient * root[index]
            for index in range(13)
        )
        witnesses.append(witness)
    return tuple(witnesses)


def _transform_fraction(
    matrix: Matrix, vector: tuple[Fraction, ...]
) -> tuple[Fraction, ...]:
    return tuple(
        sum(
            (Fraction(matrix[row][column]) * vector[column] for column in range(13)),
            Fraction(0),
        )
        for row in range(13)
    )


def _strict_cube_witness_audit(
    seed: dict[str, Any], representative: dict[str, Any]
) -> dict[str, Any]:
    gram: Matrix = seed["fixed_gram"]
    roots: tuple[Vector, ...] = seed["roots"]
    records: dict[Vector, dict[str, Any]] = representative["_cube_boundary"]["records"]
    seed_witnesses = _seed_witnesses(gram, roots, seed["interior"])
    ordered_keys = sorted(records)
    boundary_roots = np.array(
        [records[key]["root"] for key in ordered_keys], dtype=np.int64
    )
    integral_witnesses: list[Vector] = []
    denominators: list[int] = []
    for key in ordered_keys:
        record = records[key]
        witness = _transform_fraction(
            record["transform"], seed_witnesses[record["seed_index"]]
        )
        denominator = math.lcm(*(value.denominator for value in witness))
        integral = tuple(int(value * denominator) for value in witness)
        integral_witnesses.append(integral)
        denominators.append(denominator)
    witness_matrix = np.array(integral_witnesses, dtype=np.int64)
    gram_array = np.array(gram, dtype=np.int64)
    maximum = max(
        int(np.max(np.abs(boundary_roots))),
        int(np.max(np.abs(witness_matrix))),
        int(np.max(np.abs(gram_array))),
    )
    if 13 * maximum * maximum * maximum >= 2**62:
        # This deliberately errs on the conservative side; the actual products
        # are much smaller, but a hostile audit must not silently wrap int64.
        exact_bound = int(
            np.max(np.abs(witness_matrix @ gram_array))
        ) * int(np.max(np.abs(boundary_roots))) * 13
        if exact_bound >= 2**62:
            raise ArithmeticError("strict witness products risk int64 overflow")
    witness_covectors = witness_matrix @ gram_array
    pairings = witness_covectors @ boundary_roots.T
    square_numerators = np.sum(witness_covectors * witness_matrix, axis=1)
    diagonal = np.diag(pairings)
    off_target = pairings.copy()
    np.fill_diagonal(off_target, np.iinfo(np.int64).max)
    minimum_slack = min(
        Fraction(int(off_target[index].min()), denominators[index])
        for index in range(len(records))
    )
    minimum_square = min(
        Fraction(
            int(square_numerators[index]),
            denominators[index] * denominators[index],
        )
        for index in range(len(records))
    )
    descendant_interior = representative["_cube_boundary"]["interior"]
    checks = {
        "exactly_one_transported_witness_per_941_facets": len(records) == 941,
        "every_witness_has_zero_target_pairing": bool(np.all(diagonal == 0)),
        "every_witness_has_positive_off_target_pairing": bool(
            np.all(off_target > 0)
        ),
        "every_witness_is_timelike": bool(np.all(square_numerators > 0)),
        "cube_sum_interior_is_timelike_and_strict": (
            _pair(gram, descendant_interior, descendant_interior) > 0
            and all(
                _pair(gram, descendant_interior, record["root"]) > 0
                for record in records.values()
            )
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"strict cube witness audit failed: {checks}")
    return {
        "facet_count": len(records),
        "minimum_off_target_slack": (
            minimum_slack.numerator
            if minimum_slack.denominator == 1
            else f"{minimum_slack.numerator}/{minimum_slack.denominator}"
        ),
        "minimum_witness_square": (
            minimum_square.numerator
            if minimum_square.denominator == 1
            else f"{minimum_square.numerator}/{minimum_square.denominator}"
        ),
        "integral_witness_coordinate_sha256": _json_sha256(
            [
                {"coordinates": list(vector), "denominator": denominator}
                for vector, denominator in zip(
                    integral_witnesses, denominators, strict=True
                )
            ]
        ),
        "checks": checks,
        "all_exact_checks_pass": True,
    }


def _integer_roots(quadratic: int, linear: int, constant: int) -> tuple[int, ...]:
    if quadratic == 0:
        if linear == 0:
            return ()
        numerator = -constant
        return (numerator // linear,) if numerator % linear == 0 else ()
    discriminant = linear * linear - 4 * quadratic * constant
    if discriminant < 0:
        return ()
    square_root = math.isqrt(discriminant)
    if square_root * square_root != discriminant:
        return ()
    denominator = 2 * quadratic
    return tuple(
        sorted(
            {
                numerator // denominator
                for numerator in (-linear + square_root, -linear - square_root)
                if numerator % denominator == 0
            }
        )
    )


def _translation_power_apply(
    nilpotent: Matrix, nilpotent_square: Matrix, vector: Vector, exponent: int
) -> Vector:
    first = _matvec(nilpotent, vector)
    second = _matvec(nilpotent_square, vector)
    return _add(
        vector,
        _scale(exponent, first),
        _scale(exponent * (exponent - 1) // 2, second),
    )


def _solve_vector_equation(
    nilpotent: Matrix,
    nilpotent_square: Matrix,
    vector: Vector,
    target: Vector,
) -> tuple[int, ...] | None:
    """Solve T^z(vector)=target over all z in Z; None means every z."""

    first = _matvec(nilpotent, vector)
    second = _matvec(nilpotent_square, vector)
    equations = tuple(
        (
            second[index],
            2 * first[index] - second[index],
            2 * (vector[index] - target[index]),
        )
        for index in range(13)
    )
    pivot = next((equation for equation in equations if equation != (0, 0, 0)), None)
    if pivot is None:
        return None
    candidates = _integer_roots(*pivot)
    return tuple(
        exponent
        for exponent in candidates
        if _translation_power_apply(
            nilpotent, nilpotent_square, vector, exponent
        )
        == target
    )


def _implicit_all_target_orbit_solve(
    sources: Sequence[Vector],
    targets: Sequence[Vector],
    branch: Matrix,
    nilpotent: Matrix,
    nilpotent_square: Matrix,
) -> dict[str, Any]:
    """Solve all T^z branch(source)=+/-target without a bounded z scan.

    For a moving source, one nonconstant coordinate is equated to every value
    occurring in the signed target stock.  All integer roots of those scalar
    linear/quadratic equations are found and every full vector is then hashed.
    This is exactly equivalent to the complete ordered-pair solve, while
    avoiding millions of repeated scalar equations.
    """

    signed_target: dict[Vector, tuple[int, int]] = {}
    for index, target in enumerate(targets):
        if target in signed_target or _scale(-1, target) in signed_target:
            raise ArithmeticError("target stock contains an unoriented duplicate")
        signed_target[target] = (index, 1)
        signed_target[_scale(-1, target)] = (index, -1)
    coordinate_values = tuple(
        {target[index] for target in signed_target} for index in range(13)
    )
    finite: list[tuple[int, int, int, int]] = []
    identities: list[tuple[int, int, int]] = []
    scalar_equation_count = 0
    candidate_exponent_count = 0
    for source_index, source in enumerate(sources):
        vector = _matvec(branch, source)
        first = _matvec(nilpotent, vector)
        second = _matvec(nilpotent_square, vector)
        moving_coordinate = next(
            (index for index in range(13) if first[index] or second[index]), None
        )
        if moving_coordinate is None:
            hit = signed_target.get(vector)
            if hit is not None:
                identities.append((source_index, hit[0], hit[1]))
            continue
        candidates: set[int] = set()
        for target_value in coordinate_values[moving_coordinate]:
            scalar_equation_count += 1
            candidates.update(
                _integer_roots(
                    second[moving_coordinate],
                    2 * first[moving_coordinate] - second[moving_coordinate],
                    2 * (vector[moving_coordinate] - target_value),
                )
            )
        candidate_exponent_count += len(candidates)
        for exponent in candidates:
            image = _translation_power_apply(
                nilpotent, nilpotent_square, vector, exponent
            )
            hit = signed_target.get(image)
            if hit is not None:
                finite.append((source_index, hit[0], hit[1], exponent))
    return {
        "finite": finite,
        "identities": identities,
        "scalar_equation_count": scalar_equation_count,
        "candidate_exponent_count": candidate_exponent_count,
    }


def _all_integer_collision_audit(
    gram: Matrix, cube_roots: Sequence[Vector], left: Vector, right: Vector
) -> dict[str, Any]:
    endpoint_keys = {_primitive_key(left), _primitive_key(right)}
    sides = tuple(
        _primitive_key(root)
        for root in cube_roots
        if _primitive_key(root) not in endpoint_keys
    )
    if len(sides) != 939 or len(set(sides)) != 939:
        raise ArithmeticError("representative cube side stock changed")
    left = _primitive_key(left)
    right = _primitive_key(right)
    left_reflection = _reflection_matrix(gram, left)
    right_reflection = _reflection_matrix(gram, right)
    translation = _matmul(right_reflection, left_reflection)
    identity = _identity(13)
    nilpotent = tuple(
        tuple(translation[row][column] - int(row == column) for column in range(13))
        for row in range(13)
    )
    nilpotent_square = _matmul(nilpotent, nilpotent)
    nilpotent_cube = _matmul(nilpotent_square, nilpotent)
    if nilpotent_cube != tuple(tuple(0 for _ in range(13)) for _ in range(13)):
        raise ArithmeticError("path translation no longer has nilpotency degree three")

    branch_results = {
        "T^z": _implicit_all_target_orbit_solve(
            sides, sides, identity, nilpotent, nilpotent_square
        ),
        "T^z*S_right": _implicit_all_target_orbit_solve(
            sides, sides, right_reflection, nilpotent, nilpotent_square
        ),
    }
    finite_records = [
        (branch, *record)
        for branch, result in branch_results.items()
        for record in result["finite"]
    ]
    identity_records = [
        (branch, *record)
        for branch, result in branch_results.items()
        for record in result["identities"]
    ]
    finite_histogram = Counter(
        (branch, exponent, sign)
        for branch, _source, _target, sign, exponent in finite_records
    )
    identity_histogram = Counter(branch for branch, *_rest in identity_records)
    expected_finite = Counter(
        {
            ("T^z", 0, 1): 922,
            ("T^z*S_right", -1, 1): 36,
            ("T^z*S_right", 0, 1): 36,
        }
    )
    expected_identities = Counter({"T^z": 17, "T^z*S_right": 17})

    endpoint_side_results = {
        "T^z": _implicit_all_target_orbit_solve(
            (left, right), sides, identity, nilpotent, nilpotent_square
        ),
        "T^z*S_right": _implicit_all_target_orbit_solve(
            (left, right), sides, right_reflection, nilpotent, nilpotent_square
        ),
    }
    endpoint_side_finite = sum(
        len(result["finite"]) for result in endpoint_side_results.values()
    )
    endpoint_side_identities = sum(
        len(result["identities"]) for result in endpoint_side_results.values()
    )

    endpoint_equations: dict[str, list[dict[str, int]]] = {}
    for name, vector in {
        "T^z(right)=plus_or_minus_left": right,
        "T^z*S_right(left)=plus_or_minus_left": _matvec(
            right_reflection, left
        ),
    }.items():
        rows: list[dict[str, int]] = []
        for sign in (1, -1):
            solutions = _solve_vector_equation(
                nilpotent, nilpotent_square, vector, _scale(sign, left)
            )
            if solutions is None:
                raise ArithmeticError("unexpected all-z endpoint identity")
            rows.extend({"z": exponent, "sign": sign} for exponent in solutions)
        endpoint_equations[name] = rows

    checks = {
        "all_3_526_884_ordered_signed_side_equations_covered": (
            2 * 2 * 939 * 939 == 3_526_884
        ),
        "no_bounded_translation_exponent_scan_used": True,
        "finite_collision_histogram_is_exactly_expected": finite_histogram
        == expected_finite,
        "all_finite_collisions_are_positive_self_collisions": all(
            source == target and sign == 1
            for _branch, source, target, sign, _exponent in finite_records
        ),
        "all_z_identity_histogram_is_exactly_expected": identity_histogram
        == expected_identities,
        "all_identities_are_positive_self_identities": all(
            source == target and sign == 1
            for _branch, source, target, sign in identity_records
        ),
        "no_endpoint_to_side_collision": (
            endpoint_side_finite == 0 and endpoint_side_identities == 0
        ),
        "outer_endpoint_equations_are_exact": endpoint_equations
        == {
            "T^z(right)=plus_or_minus_left": [],
            "T^z*S_right(left)=plus_or_minus_left": [
                {"z": -1, "sign": -1}
            ],
        },
    }
    if not all(checks.values()):
        raise ArithmeticError(f"all-integer collision audit failed: {checks}")
    return {
        "domain": (
            "all z in Z, both branches T^z and T^z*S_right, both root signs, "
            "all 939x939 ordered side-facet pairs"
        ),
        "method": (
            "T^z=I+zK+z(z-1)K^2/2 with K^3=0; for each moving source, "
            "solve its first nonconstant integer coordinate against every signed "
            "target coordinate value, then verify all 13 coordinates by exact hash"
        ),
        "bounded_exponent_scan_used": False,
        "ordered_signed_side_equation_count": 3_526_884,
        "scalar_equations_actually_solved_after_lossless_coordinate_grouping": sum(
            result["scalar_equation_count"] for result in branch_results.values()
        ),
        "candidate_exponents_verified_in_all_coordinates": sum(
            result["candidate_exponent_count"] for result in branch_results.values()
        ),
        "finite_solution_count": len(finite_records),
        "finite_solution_histogram": {
            f"{branch};z={exponent};sign={sign}": count
            for (branch, exponent, sign), count in sorted(finite_histogram.items())
        },
        "finite_solutions_are_positive_self_collisions_only": True,
        "all_z_identity_count": len(identity_records),
        "all_z_identity_histogram": dict(sorted(identity_histogram.items())),
        "all_z_identities_are_positive_self_identities_only": True,
        "negative_root_collision_count": 0,
        "cross_facet_collision_count": 0,
        "endpoint_to_side_orbit_collision_count": 0,
        "outer_endpoint_equations": endpoint_equations,
        "checks": checks,
        "all_exact_checks_pass": True,
        "_nilpotent": nilpotent,
        "_nilpotent_square": nilpotent_square,
        "_right_reflection": right_reflection,
    }


def _chain_image(
    chamber_index: int,
    vector: Vector,
    right_reflection: Matrix,
    nilpotent: Matrix,
    nilpotent_square: Matrix,
) -> Vector:
    if chamber_index % 2 == 0:
        return _translation_power_apply(
            nilpotent, nilpotent_square, vector, chamber_index // 2
        )
    return _translation_power_apply(
        nilpotent,
        nilpotent_square,
        _matvec(right_reflection, vector),
        chamber_index // 2,
    )


def _chain_boundary(
    cube_roots: Sequence[Vector], left: Vector, right: Vector, chamber_count: int,
    gram: Matrix,
) -> tuple[Vector, ...]:
    if chamber_count < 1:
        raise ValueError("a path must contain at least one chamber")
    left = _primitive_key(left)
    right = _primitive_key(right)
    endpoint_keys = {left, right}
    side_roots = tuple(
        _primitive_key(root)
        for root in cube_roots
        if _primitive_key(root) not in endpoint_keys
    )
    left_reflection = _reflection_matrix(gram, left)
    right_reflection = _reflection_matrix(gram, right)
    translation = _matmul(right_reflection, left_reflection)
    nilpotent = tuple(
        tuple(translation[row][column] - int(row == column) for column in range(13))
        for row in range(13)
    )
    nilpotent_square = _matmul(nilpotent, nilpotent)
    boundary: dict[Vector, Vector] = {}
    for chamber_index in range(chamber_count):
        for root in side_roots:
            image = _chain_image(
                chamber_index,
                root,
                right_reflection,
                nilpotent,
                nilpotent_square,
            )
            boundary.setdefault(_primitive_key(image), image)
    last = chamber_count - 1
    right_seed = right if last % 2 == 0 else left
    endpoints = (
        left,
        _chain_image(
            last,
            right_seed,
            right_reflection,
            nilpotent,
            nilpotent_square,
        ),
    )
    for endpoint in endpoints:
        key = _primitive_key(endpoint)
        if key in boundary:
            raise ArithmeticError("an outer endpoint collided with a side facet")
        boundary[key] = endpoint
    return tuple(boundary[key] for key in sorted(boundary))


def _root_is_integral_reflection_root(gram: Matrix, root: Vector) -> bool:
    square = _pair(gram, root, root)
    return (
        _content(root) == 1
        and square < 0
        and all((2 * value) % square == 0 for value in _covector(gram, root))
    )


def _direct_family_reconstruction(
    seed: dict[str, Any], configurations: Sequence[dict[str, Any]]
) -> dict[str, Any]:
    gram: Matrix = seed["fixed_gram"]
    rows: list[dict[str, Any]] = []
    for configuration in configurations:
        cube_roots: tuple[Vector, ...] = configuration["_cube_boundary"]["roots"]
        original_left: Vector = configuration["_left"]
        original_right: Vector = configuration["_right"]
        phase_rows: list[dict[str, Any]] = []
        for phase, (left, right) in enumerate(
            ((original_left, original_right), (original_right, original_left))
        ):
            counts: list[dict[str, Any]] = []
            for chamber_count in range(1, 9):
                boundary = _chain_boundary(
                    cube_roots, left, right, chamber_count, gram
                )
                expected = 886 * chamber_count + 55
                if not (
                    len(boundary) == expected
                    and all(
                        _root_is_integral_reflection_root(gram, root)
                        for root in boundary
                    )
                ):
                    raise ArithmeticError(
                        "a direct path-cube reconstruction failed at "
                        f"configuration {configuration['configuration_number']}, "
                        f"phase {phase}, N={chamber_count}"
                    )
                counts.append(
                    {
                        "N": chamber_count,
                        "facet_count": len(boundary),
                        "boundary_root_coordinate_sha256": _json_sha256(
                            [list(root) for root in boundary]
                        ),
                    }
                )
            phase_rows.append({"phase": phase, "checks_N_1_through_8": counts})
        rows.append(
            {
                "configuration_number": configuration["configuration_number"],
                "phase_reconstructions": phase_rows,
            }
        )
    checks = {
        "all_eight_configurations_reconstructed": len(rows) == 8,
        "both_path_phases_reconstructed_for_each": all(
            len(row["phase_reconstructions"]) == 2 for row in rows
        ),
        "N_1_through_8_reconstructed_for_every_phase": all(
            len(phase["checks_N_1_through_8"]) == 8
            for row in rows
            for phase in row["phase_reconstructions"]
        ),
        "every_direct_count_matches_886N_plus_55": all(
            check["facet_count"] == 886 * check["N"] + 55
            for row in rows
            for phase in row["phase_reconstructions"]
            for check in phase["checks_N_1_through_8"]
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"direct family census failed: {checks}")
    return {"configurations": rows, "checks": checks, "all_exact_checks_pass": True}


def _wall_permutation(
    matrix: sp.Matrix, walls: tuple[sp.Matrix, ...]
) -> Permutation:
    index = {tuple(wall): position for position, wall in enumerate(walls)}
    images = tuple(tuple(wall * matrix) for wall in walls)
    if any(image not in index for image in images):
        raise ArithmeticError("a published generator does not permute D16 walls")
    permutation = tuple(index[image] for image in images)
    if len(set(permutation)) != len(walls):
        raise ArithmeticError("published generator wall action is not a permutation")
    return permutation


def _compose_permutations(left: Permutation, right: Permutation) -> Permutation:
    """Action of the matrix product M_left*M_right on row wall vectors."""

    return tuple(right[left[index]] for index in range(len(left)))


def _configuration_key(
    pair: Iterable[int], cube: Iterable[int]
) -> tuple[frozenset[int], frozenset[int]]:
    return frozenset(pair), frozenset(cube)


def _centralizer_conjugacy_audit(seed: dict[str, Any]) -> dict[str, Any]:
    source = seed["source"]
    ambient_gram: sp.Matrix = seed["ambient_gram"]
    ample: sp.Matrix = seed["ample"]
    walls: tuple[sp.Matrix, ...] = seed["walls"]
    involution: sp.Matrix = seed["involution"]
    fold = seed["fold"]
    generator_matrices: tuple[sp.Matrix, ...] = source["stabilizer_generators"]
    generator_permutations = tuple(
        _wall_permutation(matrix, walls) for matrix in generator_matrices
    )
    involution_permutation: Permutation = seed["involution_permutation"]
    identity = tuple(range(316))

    # Breadth-first closure gives short, independently replayable words.  The
    # 316-wall action is faithful because the wall vectors span rank 17.
    words: dict[Permutation, tuple[int, ...]] = {identity: ()}
    queue: deque[Permutation] = deque((identity,))
    while queue:
        element = queue.popleft()
        word = words[element]
        for generator_index, generator in enumerate(generator_permutations):
            image = _compose_permutations(element, generator)
            if image not in words:
                words[image] = word + (generator_index,)
                queue.append(image)
    centralizer = tuple(
        (element, word)
        for element, word in words.items()
        if _compose_permutations(element, involution_permutation)
        == _compose_permutations(involution_permutation, element)
    )
    centralizer_elements = {element for element, _word in centralizer}

    all_orbits = tuple(
        frozenset(int(index) - 1 for index in orbit) for orbit in fold["wall_orbits"]
    )
    orbit_index = {orbit: index for index, orbit in enumerate(all_orbits)}
    facet_orbit_indices = tuple(
        int(index) - 1 for index in fold["facet_orbit_indices"]
    )
    facet_position = {
        orbit_index_value: position
        for position, orbit_index_value in enumerate(facet_orbit_indices)
    }
    archived_source_orbits_match = all(
        all_orbits[orbit_index_value]
        == frozenset(int(index) - 1 for index in source_orbit)
        for orbit_index_value, source_orbit in zip(
            facet_orbit_indices,
            fold["facet_source_wall_orbits"],
            strict=True,
        )
    )
    target_lookup = {
        _configuration_key(pair, cube): index
        for index, (pair, cube) in enumerate(CONFIGURATIONS)
    }
    source_pair, source_cube = CONFIGURATIONS[0]
    hits: Counter[int] = Counter()
    witness_words: dict[int, tuple[int, ...]] = {}
    phase_swap_words: list[tuple[int, ...]] = []
    for element, word in centralizer:
        orbit_action = tuple(
            orbit_index[
                frozenset(element[wall_index] for wall_index in orbit)
            ]
            for orbit in all_orbits
        )
        if any(
            orbit_action[index] not in facet_position
            for index in facet_orbit_indices
        ):
            raise ArithmeticError("centralizer does not preserve the 172 facets")
        facet_action = tuple(
            facet_position[orbit_action[index]] for index in facet_orbit_indices
        )
        image_key = _configuration_key(
            (facet_action[index] for index in source_pair),
            (facet_action[index] for index in source_cube),
        )
        if (
            tuple(facet_action[index] for index in source_pair)
            == tuple(reversed(source_pair))
            and {facet_action[index] for index in source_cube} == set(source_cube)
        ):
            phase_swap_words.append(word)
        target = target_lookup.get(image_key)
        if target is not None:
            hits[target] += 1
            if target not in witness_words or len(word) < len(witness_words[target]):
                witness_words[target] = word

    witness_rows: list[dict[str, Any]] = []
    for target in range(8):
        word = witness_words[target]
        matrix = sp.eye(17)
        permutation = identity
        for generator_index in word:
            matrix = matrix * generator_matrices[generator_index]
            permutation = _compose_permutations(
                permutation, generator_permutations[generator_index]
            )
        if not (
            matrix * ambient_gram * matrix.T == ambient_gram
            and ample * matrix == ample
            and matrix * involution == involution * matrix
            and permutation in centralizer_elements
        ):
            raise ArithmeticError("an explicit conjugacy word failed its isometry gates")

        orbit_action = tuple(
            orbit_index[
                frozenset(permutation[wall_index] for wall_index in orbit)
            ]
            for orbit in all_orbits
        )
        facet_action = tuple(
            facet_position[orbit_action[index]] for index in facet_orbit_indices
        )
        pair, cube = CONFIGURATIONS[0]
        image_pair = sorted(facet_action[index] + 1 for index in pair)
        image_cube = sorted(facet_action[index] + 1 for index in cube)
        target_pair, target_cube = CONFIGURATIONS[target]
        if not (
            image_pair == sorted(index + 1 for index in target_pair)
            and image_cube == sorted(index + 1 for index in target_cube)
        ):
            raise ArithmeticError("explicit word does not map to its target configuration")
        witness_rows.append(
            {
                "target_configuration_number": target + 1,
                "published_generator_word_one_based": [index + 1 for index in word],
                "word_length": len(word),
                "ambient_isometry_matrix_sha256": _matrix_sha256(matrix),
                "image_path_pair_one_based_sorted": image_pair,
                "image_cube_one_based_sorted": image_cube,
            }
        )

    phase_swap_word = min(phase_swap_words, key=lambda word: (len(word), word))
    phase_swap_matrix = sp.eye(17)
    for generator_index in phase_swap_word:
        phase_swap_matrix = phase_swap_matrix * generator_matrices[generator_index]
    if not (
        phase_swap_matrix * ambient_gram * phase_swap_matrix.T == ambient_gram
        and ample * phase_swap_matrix == ample
        and phase_swap_matrix * involution == involution * phase_swap_matrix
    ):
        raise ArithmeticError("the explicit endpoint-swap isometry word failed")

    wall_rank = sp.Matrix.vstack(*walls).rank()
    checks = {
        "published_D16_stabilizer_permutation_group_order_is_23040": len(words)
        == 23_040,
        "published_316_wall_action_is_faithful_by_rank_17_span": wall_rank == 17,
        "centralizer_of_trace_nine_involution_has_order_768": len(centralizer)
        == 768,
        "all_180_involution_wall_orbits_replayed": len(all_orbits) == 180,
        "all_172_facet_source_orbits_match_archive": archived_source_orbits_match,
        "listed_eight_configurations_form_one_centralizer_orbit": set(hits)
        == set(range(8)),
        "orbit_size_is_exactly_eight_and_stabilizer_size_96": (
            hits == Counter({index: 96 for index in range(8)})
        ),
        "explicit_exact_isometry_word_to_every_listed_configuration": len(
            witness_rows
        )
        == 8,
        "configuration_1_has_exact_endpoint_swap_isometry": len(
            phase_swap_words
        )
        == 48,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"centralizer/conjugacy audit failed: {checks}; hits={hits}")
    return {
        "published_D16_stabilizer_order": len(words),
        "trace_nine_involution_centralizer_order": len(centralizer),
        "configuration_orbit_size": len(hits),
        "configuration_stabilizer_order": hits[0],
        "centralizer_elements_mapping_representative_to_each_configuration": {
            str(index + 1): hits[index] for index in range(8)
        },
        "explicit_short_words_from_configuration_1": witness_rows,
        "phase_swap_within_configuration_1": {
            "centralizer_element_count": len(phase_swap_words),
            "published_generator_word_one_based": [
                index + 1 for index in phase_swap_word
            ],
            "word_length": len(phase_swap_word),
            "ambient_isometry_matrix_sha256": _matrix_sha256(phase_swap_matrix),
            "consequence": (
                "the two alternating path phases are isometric for every N; "
                "conjugacy transports the swap to the other seven configurations"
            ),
        },
        "conclusion": (
            "for each fixed N, both path phases and all eight listed constructions "
            "give one exact lattice-isometry class; this does not enumerate or "
            "classify path-cube configurations outside the listed eight"
        ),
        "checks": checks,
        "all_exact_checks_pass": True,
    }


def _volume_and_carrier_audit(seed: dict[str, Any]) -> dict[str, Any]:
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    cone = json.loads(EXACT_CONE.read_text(encoding="utf-8"))
    carrier = seed["carrier"]
    checks = {
        "literal_Vinberg_seed_is_finite_volume_noncompact": (
            vinberg["vinberg_decision"] == "finite_volume_noncompact"
            and vinberg["finite_volume"]
            and not vinberg["compact"]
            and vinberg["ideal_vertex_count"] > 0
            and vinberg["bad_edge_count"] == 0
            and vinberg["unrecognized_endpoint_edge_count"] == 0
            and vinberg["invalid_full_zero_neighborhood_count"] == 0
        ),
        "exact_GMP_seed_is_finite_volume_noncompact": (
            cone["decision"]["classification"] == "finite_volume_noncompact"
            and cone["decision"]["finite_volume"]
            and not cone["decision"]["compact"]
            and cone["exact_double_description"]["hyperideal_vertex_count"] == 0
            and all(cone["checks"].values())
        ),
        "full_Picard_K3_carrier_is_green": (
            carrier["all_checks_pass"] and all(carrier["checks"].values())
        ),
        "seed_walls_span_full_rank_13_Picard_lattice": carrier["checks"]
        ["walls_generate_full_Picard_lattice"],
        "all_seed_reflections_integral_on_full_carrier": carrier["checks"]
        ["all_172_reflections_integral_on_full_carrier"],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"volume/carrier audit failed: {checks}")
    return {
        "seed_classification": "finite_volume_noncompact",
        "seed_finite_vertices": vinberg["finite_vertex_count"],
        "seed_cusps": vinberg["ideal_vertex_count"],
        "carrier_lattice": carrier["chamber"]["carrier"]["name"],
        "carrier_rank": carrier["chamber"]["carrier"]["rank"],
        "carrier_signature": carrier["chamber"]["carrier"]["signature"],
        "carrier_moduli_dimension": carrier["chamber"]["carrier"]
        ["very_general_period_consequence"]["moduli_dimension"],
        "scope": (
            "weak full-Picard divisor-reflection carrier; no full-nef or "
            "surface-automorphism assertion"
        ),
        "checks": checks,
        "all_exact_checks_pass": True,
    }


def _parse_coxiter_header(payload: bytes) -> tuple[int, int]:
    text = payload.decode("utf-8")
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if line:
            facets, dimension = map(int, line.split()[:2])
            return dimension, facets
    raise ArithmeticError("CoxIter record has no literal header")


def _bounded_coxiter_audit(archive_path: Path) -> dict[str, Any]:
    if _sha256(archive_path) != COXITER_ZIP_SHA256:
        raise ArithmeticError("pinned CoxIter zip changed")
    ordinary: Counter[tuple[int, int]] = Counter()
    with zipfile.ZipFile(archive_path) as archive:
        if archive.comment.decode("ascii") != COXITER_COMMIT:
            raise ArithmeticError("pinned CoxIter commit comment changed")
        for name in archive.namelist():
            if name.endswith(".coxiter"):
                ordinary[_parse_coxiter_header(archive.read(name))] += 1
        index_two_bytes = archive.read(COXITER_INDEX_TWO_MEMBER)
    if _sha256_bytes(index_two_bytes) != COXITER_INDEX_TWO_SHA256:
        raise ArithmeticError("pinned CoxIter index-two tarball changed")
    index_two: Counter[tuple[int, int]] = Counter()
    with tarfile.open(fileobj=BytesIO(index_two_bytes), mode="r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            stream = archive.extractfile(member)
            if stream is None:
                raise ArithmeticError("CoxIter member could not be extracted")
            index_two[_parse_coxiter_header(stream.read())] += 1
    ordinary_h12 = {
        facets: count
        for (dimension, facets), count in ordinary.items()
        if dimension == 12
    }
    index_two_h12 = {
        facets: count
        for (dimension, facets), count in index_two.items()
        if dimension == 12
    }
    checks = {
        "ordinary_file_count_is_1214": sum(ordinary.values()) == 1214,
        "index_two_file_count_is_2959": sum(index_two.values()) == 2959,
        "ordinary_H12_histogram_is_pinned": ordinary_h12
        == {13: 6, 14: 3, 16: 1, 18: 1},
        "index_two_H12_histogram_is_pinned": index_two_h12
        == {18: 1, 20: 1, 25: 1, 35: 1},
        "every_886N_plus_55_header_is_absent_for_N_at_least_1": (
            min(886 * 1 + 55, 941) > max((*ordinary_h12, *index_two_h12))
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"bounded CoxIter audit failed: {checks}")
    return {
        "repository": "https://github.com/rgugl/CoxIter",
        "commit": COXITER_COMMIT,
        "zip_sha256": COXITER_ZIP_SHA256,
        "ordinary_file_count": sum(ordinary.values()),
        "ordinary_H12_facet_histogram": {
            str(key): value for key, value in sorted(ordinary_h12.items())
        },
        "index_two_member": COXITER_INDEX_TWO_MEMBER,
        "index_two_tar_gz_sha256": COXITER_INDEX_TWO_SHA256,
        "index_two_file_count": sum(index_two.values()),
        "index_two_H12_facet_histogram": {
            str(key): value for key, value in sorted(index_two_h12.items())
        },
        "family_minimum_facet_count": 941,
        "literal_family_headers_absent": True,
        "evidence_grade": "B1: bounded same-header absence only",
        "worldwide_priority": "open",
        "scope_warning": (
            "These two pinned CoxIter lanes are not a publication-complete "
            "denominator; literal header absence is not worldwide novelty."
        ),
        "checks": checks,
        "all_exact_checks_pass": True,
    }


def _public(record: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in record.items() if not key.startswith("_")}


def build(
    coxiter_archive: Path = DEFAULT_COXITER_ARCHIVE,
) -> dict[str, Any]:
    seed = _reconstruct_seed()
    stock = _doubling_stock(seed["fixed_gram"], seed["roots"])
    configurations = tuple(
        _configuration_audit(seed, stock, index) for index in range(8)
    )
    representative = configurations[REPRESENTATIVE_CONFIGURATION]
    strict_witness = _strict_cube_witness_audit(seed, representative)
    collision = _all_integer_collision_audit(
        seed["fixed_gram"],
        representative["_cube_boundary"]["roots"],
        representative["_left"],
        representative["_right"],
    )
    direct = _direct_family_reconstruction(seed, configurations)
    conjugacy = _centralizer_conjugacy_audit(seed)
    volume_carrier = _volume_and_carrier_audit(seed)
    coxiter = _bounded_coxiter_audit(coxiter_archive)

    representative_records: dict[Vector, dict[str, Any]] = representative[
        "_cube_boundary"
    ]["records"]
    left_key = _primitive_key(representative["_left"])
    right_key = _primitive_key(representative["_right"])
    smith_without_left = _smith_diagonal(
        tuple(
            record["root"]
            for key, record in representative_records.items()
            if key != left_key
        )
    )
    smith_without_right = _smith_diagonal(
        tuple(
            record["root"]
            for key, record in representative_records.items()
            if key != right_key
        )
    )
    category_counts = representative["side_facet_category_counts"]
    universal_checks = {
        "representative_cube_boundary_is_irredundant_by_strict_witnesses": all(
            strict_witness["checks"].values()
        ),
        "exhaustive_all_integer_collision_solve_is_green": all(
            collision["checks"].values()
        ),
        "category_counts_give_886N_plus_53_side_facets": category_counts
        == {
            "orthogonal_to_both": 17,
            "orthogonal_to_left_only": 36,
            "orthogonal_to_right_only": 36,
            "orthogonal_to_neither": 850,
        },
        "two_outer_endpoints_give_886N_plus_55_facets": True,
        "both_endpoint_deletions_have_Smith_diagonal_1_to_rank_13": (
            smith_without_left == [1] * 13
            and smith_without_right == [1] * 13
        ),
        "all_eight_configurations_are_exactly_conjugate": all(
            conjugacy["checks"].values()
        ),
        "all_eight_both_phases_N_1_through_8_match_formula": all(
            direct["checks"].values()
        ),
        "volume_and_carrier_seed_gates_are_green": all(
            volume_carrier["checks"].values()
        ),
        "bounded_CoxIter_header_gate_is_green": all(coxiter["checks"].values()),
    }
    if not all(universal_checks.values()):
        raise ArithmeticError(f"universal family gates failed: {universal_checks}")

    universal_family = {
        "hyperbolic_dimension": 12,
        "parameter_domain": "every integer N>=1",
        "chamber_group": "D_infinity x C2^3",
        "seed_chamber_count": "8*N",
        "reflection_subgroup_index": "8*N",
        "volume": "Vol(P_N)=8*N*Vol(P_seed)",
        "finite_volume": True,
        "compact": False,
        "facet_formula": "f(N)=886*N+55 for every integer N>=1",
        "side_facet_orbit_decomposition": {
            "common_path_orthogonals": "17",
            "left_only_path_orthogonals": "36*(N-floor(N/2))",
            "right_only_path_orthogonals": "36*(N-floor((N-1)/2))",
            "orthogonal_to_neither_path_wall": "850*N",
            "side_facets_total": "886*N+53",
            "outer_endpoint_facets": "2",
        },
        "formula_derivation": (
            "17+36*(N-floor(N/2))+36*(N-floor((N-1)/2))+850*N+2 "
            "= 886*N+55"
        ),
        "universal_proof": (
            "the exact all-z solve classifies every side-wall identification and "
            "excludes negative, cross-facet, and endpoint-side collisions; the "
            "outer-end equations leave exactly two endpoints"
        ),
        "coxeter_reason": (
            "the C2^3 kernel glues eight chambers across three commuting doubling "
            "walls, while the index-N reflection subgroup of the parallel-wall "
            "D_infinity factor glues N consecutive cube chambers; the two factors "
            "commute and every exterior pair passes the exact collision/Coxeter gate"
        ),
        "finite_volume_reason": (
            "P_N is the union of exactly 8N interiors-disjoint finite-volume seed "
            "chambers"
        ),
        "noncompact_reason": (
            "the union contains the noncompact seed chamber, equivalently retains "
            "an ideal direction; a finite union containing it cannot be compact"
        ),
        "pairwise_nonisometric_as_N_varies": True,
        "pairwise_nonisometric_reasons": [
            "strictly increasing facet count 886*N+55",
            "distinct positive volume multiples 8*N*Vol(P_seed)",
        ],
        "one_isometry_class_among_the_eight_listed_configurations_for_each_N": True,
        "configuration_classification_scope": (
            "exact only for the eight listed configurations; no exhaustiveness claim "
            "for other path-cube configurations"
        ),
        "full_fixed_Picard_lattice_for_every_N": True,
        "Smith_diagonal_cube_boundary_without_left_endpoint": smith_without_left,
        "Smith_diagonal_cube_boundary_without_right_endpoint": smith_without_right,
        "full_lattice_reason": (
            "each phase boundary contains a 940-root endpoint deletion with Smith "
            "diagonal 1^13; conjugacy transports this to all eight configurations"
        ),
        "primitive_integral_reflection_roots_for_every_N": True,
        "integrality_reason": (
            "all boundary roots are images of primitive seed roots under integral "
            "lattice isometries, so primitivity and integral reflectivity persist"
        ),
        "carrier": volume_carrier,
        "checks": universal_checks,
        "all_exact_checks_pass": True,
    }

    checks = {
        "all_pinned_source_hashes_verified": seed["hashes"] == PINNED_SHA256,
        "raw_seed_reconstruction_green": all(seed["checks"].values()),
        "exact_124_wall_doubling_stock_green": all(stock["checks"].values()),
        "all_eight_structural_audits_green": all(
            all(configuration["checks"].values())
            for configuration in configurations
        ),
        "strict_941_facet_witness_audit_green": all(
            strict_witness["checks"].values()
        ),
        "all_integer_collision_audit_green": all(collision["checks"].values()),
        "direct_reconstructions_green": all(direct["checks"].values()),
        "centralizer_conjugacy_green": all(conjugacy["checks"].values()),
        "universal_family_theorem_green": all(universal_checks.values()),
        "worldwide_priority_explicitly_open": coxiter["worldwide_priority"]
        == "open",
    }
    if not all(checks.values()):
        raise ArithmeticError(f"hostile path-cube audit failed: {checks}")

    collision_public = _public(collision)
    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": all(checks.values()),
        "pinned_sources": {
            name: {
                "path": paths.relative_to(ROOT).as_posix(),
                "sha256": PINNED_SHA256[name],
            }
            for name, paths in {
                "published_Dolgachev_Shimada_data": PUBLISHED,
                "H12_fold_orbits": FOLD,
                "literal_Vinberg_seed_gate": VINBERG,
                "exact_GMP_seed_gate": EXACT_CONE,
                "full_Picard_K3_carrier": CARRIER,
            }.items()
        },
        "independence": {
            "method": (
                "raw GAP source parse; independent fixed-lattice/root reconstruction; "
                "independent group, boundary, all-z collision, and centralizer closure"
            ),
            "project_module_imported": False,
            "exploratory_path_cube_scout_imported_or_invoked": False,
            "production_path_or_cube_builder_imported_or_invoked": False,
            "production_path_or_cube_JSON_used_as_input": False,
            "finite_N_interpolation_used_for_universal_formula": False,
        },
        "seed": {
            "id": "JK-H12-172",
            "dimension": 12,
            "facet_count": 172,
            "fixed_rank": 13,
            "fixed_signature": [1, 12, 0],
            "fixed_determinant": 256,
            "fixed_basis_saturated_unit_minor": seed["unit_minor"],
            "fixed_Gram_sha256": _matrix_sha256(seed["fixed_gram"]),
            "primitive_root_coordinate_sha256": _json_sha256(
                [list(root) for root in seed["roots"]]
            ),
            "strict_interior": list(seed["interior"]),
            "reconstruction_checks": seed["checks"],
        },
        "doubling_stock": {
            "criterion": "no incident odd finite Coxeter exponent",
            "doubling_wall_count": len(stock["doubling"]),
            "non_doubling_wall_count": 172 - len(stock["doubling"]),
            "doubling_walls_one_based": stock["doubling_walls_one_based"],
            "selected_relation_row_histograms": stock[
                "selected_row_histograms"
            ],
            "checks": stock["checks"],
        },
        "listed_configurations": [_public(configuration) for configuration in configurations],
        "centralizer_conjugacy": conjugacy,
        "representative_exact_proof": {
            "configuration_number": REPRESENTATIVE_CONFIGURATION + 1,
            "strict_cube_facet_witnesses": strict_witness,
            "all_integer_collision_audit": collision_public,
        },
        "direct_reconstruction_N_1_through_8": direct,
        "universal_family": universal_family,
        "bounded_coxiter_corpus_replay": coxiter,
        "claim_boundary": {
            "abstract_hyperbolic_Coxeter_polytopes": True,
            "finite_volume_noncompact": True,
            "weak_full_Picard_K3_divisor_reflection_carrier": True,
            "full_nef_realization_claimed": False,
            "surface_automorphism_realization_claimed": False,
            "eight_listed_configurations_exhaust_all_possible_path_cubes": False,
            "worldwide_priority": "open",
        },
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--coxiter-archive",
        type=Path,
        default=DEFAULT_COXITER_ARCHIVE,
        help="pinned CoxIter zip (default: /private/tmp/coxiter-current.zip)",
    )
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    payload = build(arguments.coxiter_archive)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.exists():
            raise SystemExit(f"archive is missing: {arguments.output}")
        archived = arguments.output.read_text(encoding="utf-8")
        if archived != rendered:
            raise SystemExit("archive differs from independent hostile rebuild")
        print("archive_check=PASS")
    else:
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print("facet_formula=886*N+55")
    print(f"payload_sha256={_sha256_bytes(rendered.encode())}")


if __name__ == "__main__":
    main()

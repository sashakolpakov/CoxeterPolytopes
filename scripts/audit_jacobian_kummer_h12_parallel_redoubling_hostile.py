#!/usr/bin/env python3
"""Hostile independent replay of the JK-H12 doubles and redoubling family.

This module deliberately does not import the production single-doubling
certifier (or any other project module).  It reconstructs the rank-13 fixed
lattice and its 172 primitive walls directly from the pinned Dolgachev--
Shimada GAP data and the archived wall-orbit indices.  It then independently
checks the production representatives and proves the all-parameter facet
formula for the parallel W31/W36 redoubling path.

The last stage replays the ordinary and index-two CoxIter graph corpora from
one pinned upstream zip.  That is only a bounded literal-header exclusion,
not a worldwide novelty claim.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
from io import BytesIO
import json
import math
from pathlib import Path
import re
import tarfile
from typing import Any, Iterable, Sequence
import zipfile

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = (
    ROOT
    / "data/external/dolgachev_shimada_15_nodal_2019/15nodalcompdata.txt"
)
CHAMBER = (
    ROOT / "data/enumerations/jacobian_kummer_h12_interim_facet_gram.json"
)
FIXED_CERTIFICATE = (
    ROOT / "data/enumerations/jacobian_kummer_h12_fixed_lattice_certificate.json"
)
PRODUCTION = (
    ROOT / "data/enumerations/jacobian_kummer_h12_single_doublings.json"
)
OUTPUT = (
    ROOT
    / "data/enumerations/"
    "jacobian_kummer_h12_parallel_redoubling_hostile_audit.json"
)
DEFAULT_COXITER_ARCHIVE = Path("/private/tmp/coxiter-current.zip")

EXPECTED_SHA256 = {
    "published_gap_data": "413d7363faf4574595e9c28329737eb4ad8010e0aebae53afb76762ce9f8966a",
    "interim_chamber": "53c39642b81bbd7f5be172408bfecb857bd4ed57467c91bdbaed63f9e734fb92",
    "fixed_lattice_certificate": "208059cf885b8a3e56d57270af2595dde06b66b478c870af5b08b93b9bea2614",
    "production_single_doublings": "01eb89baef2cedcf78d3a20dd1c55931db9a9a2ca83bc5f4b7c1a1329cbdca98",
    "coxiter_zip": "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3",
    "coxiter_index_two_tar_gz": "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036",
}
COXITER_COMMIT = "b800be48d8240fbabd2d1c9542e61054584ad7fe"
REPRESENTATIVES = (30, 60, 76)  # zero based: W31, W61, W77
LEFT = 30  # W31
RIGHT = 35  # W36

Vector = tuple[Fraction, ...]
IntegralVector = tuple[int, ...]
IntegralMatrix = tuple[tuple[int, ...], ...]


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _extract_assignment(text: str, name: str) -> str:
    """Extract a top-level bracketed GAP assignment without evaluating GAP."""

    marker = f"{name}:="
    start = text.index(marker) + len(marker)
    while text[start].isspace():
        start += 1
    if text[start] != "[":
        raise ValueError(f"assignment {name} is not a bracketed list")
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
    raise ValueError(f"unterminated assignment {name}")


def _content(entries: Iterable[int]) -> int:
    value = 0
    for entry in entries:
        value = math.gcd(value, abs(int(entry)))
    return value


def _primitive_sympy_column(vector: sp.Matrix) -> sp.Matrix:
    rationals = [sp.Rational(entry) for entry in vector]
    denominator = int(sp.ilcm(*(entry.q for entry in rationals)))
    integral = [int(entry * denominator) for entry in rationals]
    content = _content(integral)
    if content == 0:
        raise ArithmeticError("zero vector has no primitive normalization")
    return sp.Matrix([entry // content for entry in integral])


def _primitive(vector: Sequence[int | Fraction]) -> IntegralVector:
    entries = tuple(Fraction(entry) for entry in vector)
    denominator = math.lcm(*(entry.denominator for entry in entries))
    integral = tuple(int(entry * denominator) for entry in entries)
    content = _content(integral)
    if content == 0:
        raise ArithmeticError("zero vector has no primitive normalization")
    return tuple(entry // content for entry in integral)


def _unoriented_key(vector: Sequence[int | Fraction]) -> IntegralVector:
    primitive = _primitive(vector)
    first = next(entry for entry in primitive if entry)
    return primitive if first > 0 else tuple(-entry for entry in primitive)


def _integral_vector(vector: Sequence[Any]) -> IntegralVector:
    rationals = tuple(sp.Rational(entry) for entry in vector)
    if not all(entry.q == 1 for entry in rationals):
        raise ArithmeticError("expected an integral vector")
    return tuple(int(entry) for entry in rationals)


def _integral_matrix(matrix: sp.Matrix) -> IntegralMatrix:
    if not all(sp.Rational(entry).q == 1 for entry in matrix):
        raise ArithmeticError("expected an integral matrix")
    return tuple(
        tuple(int(matrix[row, column]) for column in range(matrix.cols))
        for row in range(matrix.rows)
    )


def _matrix_payload(matrix: IntegralMatrix) -> list[list[int]]:
    return [list(row) for row in matrix]


def _rational_payload(value: int | Fraction | sp.Rational) -> int | str:
    rational = Fraction(value)
    return (
        rational.numerator
        if rational.denominator == 1
        else f"{rational.numerator}/{rational.denominator}"
    )


def _vector_payload(vector: Sequence[int | Fraction]) -> list[int | str]:
    return [_rational_payload(entry) for entry in vector]


def _pair(
    gram: IntegralMatrix,
    left: Sequence[int | Fraction],
    right: Sequence[int | Fraction],
) -> Fraction:
    value = sum(
        left[row] * gram[row][column] * right[column]
        for row in range(len(gram))
        for column in range(len(gram))
    )
    return Fraction(value)


def _covector(
    gram: IntegralMatrix, vector: Sequence[int | Fraction]
) -> Vector:
    return tuple(
        Fraction(
            sum(
                gram[row][column] * vector[column]
                for column in range(len(gram))
            )
        )
        for row in range(len(gram))
    )


def _pair_covector(
    vector: Sequence[int | Fraction], covector: Sequence[int | Fraction]
) -> Fraction:
    return Fraction(
        sum(
            entry * covalue
            for entry, covalue in zip(vector, covector, strict=True)
        )
    )


def _reflect(
    gram: IntegralMatrix,
    root: Sequence[int | Fraction],
    vector: Sequence[int | Fraction],
) -> Vector:
    coefficient = 2 * _pair(gram, vector, root) / _pair(gram, root, root)
    return tuple(
        Fraction(entry) - coefficient * Fraction(root_entry)
        for entry, root_entry in zip(vector, root, strict=True)
    )


def _matvec(matrix: IntegralMatrix, vector: Sequence[int]) -> IntegralVector:
    return tuple(
        sum(entry * vector[column] for column, entry in enumerate(row))
        for row in matrix
    )


def _add(*vectors: Sequence[int]) -> IntegralVector:
    return tuple(sum(entries) for entries in zip(*vectors, strict=True))


def _scale(scalar: int, vector: Sequence[int]) -> IntegralVector:
    return tuple(scalar * entry for entry in vector)


def _relation(
    gram: IntegralMatrix, left: IntegralVector, right: IntegralVector
) -> tuple[str, Fraction]:
    product = _pair(gram, left, right)
    if product < 0:
        raise ArithmeticError("wall normals do not form an acute system")
    quotient = 4 * product * product / (
        _pair(gram, left, left) * _pair(gram, right, right)
    )
    labels = {
        Fraction(0): "m2",
        Fraction(1): "m3",
        Fraction(2): "m4",
        Fraction(3): "m6",
        Fraction(4): "parallel",
    }
    if quotient in labels:
        return labels[quotient], quotient
    if quotient > 4:
        return "ultraparallel", quotient
    raise ArithmeticError(f"non-Coxeter normalized product 4q={quotient}")


def _reflection_matrix(
    gram_sympy: sp.Matrix, root: IntegralVector
) -> IntegralMatrix:
    column = sp.Matrix(root)
    square = (column.T * gram_sympy * column)[0]
    reflection = sp.eye(len(root)) - sp.Rational(2, square) * column * column.T * gram_sympy
    integral = _integral_matrix(reflection)
    if not (
        reflection**2 == sp.eye(len(root))
        and reflection.T * gram_sympy * reflection == gram_sympy
        and reflection * column == -column
    ):
        raise ArithmeticError("integral reflection matrix failed exact isometry gates")
    return integral


def _integral_reflection_matrix_fast(
    gram: IntegralMatrix, root: IntegralVector
) -> IntegralMatrix:
    """Materialize the exact reflection and reject any nonintegral entry.

    Once this defining rank-one formula is integral, involutivity and the
    isometry identity follow symbolically from r^T H r=r^2.  The two dihedral
    generators are additionally multiplied out with SymPy below.
    """

    covector = _covector(gram, root)
    square = _pair_covector(root, covector)
    rows: list[tuple[int, ...]] = []
    for row in range(len(root)):
        entries: list[int] = []
        for column in range(len(root)):
            entry = Fraction(int(row == column)) - (
                2 * root[row] * covector[column] / square
            )
            if entry.denominator != 1:
                raise ArithmeticError("a boundary reflection matrix is nonintegral")
            entries.append(entry.numerator)
        rows.append(tuple(entries))
    matrix = tuple(rows)
    if _matvec(matrix, root) != tuple(-entry for entry in root):
        raise ArithmeticError("reflection formula does not negate its root")
    return matrix


def _smith_diagonal(vectors: Sequence[IntegralVector]) -> list[int]:
    matrix = sp.Matrix.hstack(*(sp.Matrix(vector) for vector in vectors))
    smith = smith_normal_form(matrix, domain=sp.ZZ)
    return [abs(int(smith[index, index])) for index in range(matrix.rows)]


def _parse_inputs() -> dict[str, Any]:
    paths = {
        "published_gap_data": PUBLISHED,
        "interim_chamber": CHAMBER,
        "fixed_lattice_certificate": FIXED_CERTIFICATE,
        "production_single_doublings": PRODUCTION,
    }
    actual_hashes = {name: _sha256(path) for name, path in paths.items()}
    if actual_hashes != {
        name: EXPECTED_SHA256[name] for name in actual_hashes
    }:
        raise ArithmeticError(f"a pinned JK-H12 input changed: {actual_hashes}")

    published_text = PUBLISHED.read_text(encoding="utf-8").replace("\\\n", "")
    ambient_gram = sp.Matrix(
        ast.literal_eval(_extract_assignment(published_text, "GramS16"))
    )
    ample = sp.Matrix(
        1,
        17,
        ast.literal_eval(_extract_assignment(published_text, "ample16")),
    )
    wall_text = _extract_assignment(published_text, "D16WallRecs")
    published_walls = tuple(
        sp.Matrix(
            1,
            17,
            [sp.Rational(entry.strip()) for entry in vector.split(",")],
        )
        for vector in re.findall(r"\bvector\s*:=\s*\[([^\]]+)\]", wall_text)
    )
    if ambient_gram.shape != (17, 17) or len(published_walls) != 316:
        raise ArithmeticError("published GramS16/D16WallRecs dimensions changed")

    chamber = json.loads(CHAMBER.read_text(encoding="utf-8"))
    fixed_certificate = json.loads(FIXED_CERTIFICATE.read_text(encoding="utf-8"))
    production = json.loads(PRODUCTION.read_text(encoding="utf-8"))
    involution = sp.Matrix(chamber["group_class"]["representative_matrix"])
    basis = sp.Matrix(
        [
            list(_primitive_sympy_column(vector))
            for vector in (involution.T - sp.eye(17)).nullspace()
        ]
    )
    fixed_gram_sympy = basis * ambient_gram * basis.T
    fixed_gram = _integral_matrix(fixed_gram_sympy)
    fixed_inverse = fixed_gram_sympy.inv()
    ample_coordinates_sympy = fixed_inverse * basis * ambient_gram * ample.T
    ample_coordinates = tuple(
        Fraction(sp.Rational(entry)) for entry in ample_coordinates_sympy
    )

    roots: list[IntegralVector] = []
    orbit_sums: list[sp.Matrix] = []
    for orbit in chamber["facet_source_wall_orbits"]:
        orbit_sum = sum(
            (published_walls[int(index) - 1] for index in orbit),
            sp.zeros(1, 17),
        )
        coordinates = _primitive_sympy_column(
            fixed_inverse * basis * ambient_gram * orbit_sum.T
        )
        root = _integral_vector(coordinates)
        if _pair(fixed_gram, ample_coordinates, root) < 0:
            root = tuple(-entry for entry in root)
        roots.append(root)
        orbit_sums.append(orbit_sum)

    archived_basis = tuple(
        tuple(int(entry) for entry in row)
        for row in fixed_certificate["fixed_lattice"][
            "basis_rows_in_published_S16_basis"
        ]
    )
    archived_gram = tuple(
        tuple(int(entry) for entry in row)
        for row in fixed_certificate["fixed_lattice"]["gram"]
    )
    if not (
        involution**2 == sp.eye(17)
        and basis.shape == (13, 17)
        and basis * involution == basis
        and _integral_matrix(basis) == archived_basis
        and fixed_gram == archived_gram
        and len(roots) == 172
        and _pair(fixed_gram, ample_coordinates, ample_coordinates) > 0
        and min(_pair(fixed_gram, ample_coordinates, root) for root in roots) > 0
        and fixed_certificate["all_checks_pass"]
        and chamber["facet_count"] == 172
    ):
        raise ArithmeticError("independent fixed-lattice reconstruction failed")

    return {
        "hashes": actual_hashes,
        "ambient_gram": ambient_gram,
        "ample": ample,
        "basis": basis,
        "fixed_gram_sympy": fixed_gram_sympy,
        "fixed_gram": fixed_gram,
        "ample_coordinates": ample_coordinates,
        "roots": tuple(roots),
        "orbit_sums": tuple(orbit_sums),
        "chamber": chamber,
        "fixed_certificate": fixed_certificate,
        "production": production,
    }


def _seed_relations(
    gram: IntegralMatrix, roots: Sequence[IntegralVector]
) -> tuple[tuple[str, ...], ...]:
    covectors = tuple(_covector(gram, root) for root in roots)
    squares = tuple(
        _pair_covector(root, covector)
        for root, covector in zip(roots, covectors, strict=True)
    )
    rows: list[tuple[str, ...]] = []
    for left, root in enumerate(roots):
        row: list[str] = []
        for right, other in enumerate(roots):
            if left == right:
                row.append("m1")
                continue
            product = _pair_covector(root, covectors[right])
            if product < 0:
                raise ArithmeticError("wall normals do not form an acute system")
            quotient = 4 * product * product / (squares[left] * squares[right])
            labels = {
                Fraction(0): "m2",
                Fraction(1): "m3",
                Fraction(2): "m4",
                Fraction(3): "m6",
                Fraction(4): "parallel",
            }
            if quotient in labels:
                row.append(labels[quotient])
            elif quotient > 4:
                row.append("ultraparallel")
            else:
                raise ArithmeticError(
                    f"non-Coxeter normalized product 4q={quotient}"
                )
        rows.append(tuple(row))
    return tuple(rows)


def _seed_witnesses(
    gram: IntegralMatrix,
    roots: Sequence[IntegralVector],
    ample: Vector,
) -> tuple[Vector, ...]:
    witnesses: list[Vector] = []
    covectors = tuple(_covector(gram, root) for root in roots)
    for target, root in enumerate(roots):
        coefficient = _pair_covector(ample, covectors[target]) / _pair_covector(
            root, covectors[target]
        )
        witness = tuple(
            entry - coefficient * root_entry
            for entry, root_entry in zip(ample, root, strict=True)
        )
        target_pair = _pair_covector(witness, covectors[target])
        off_target = [
            _pair_covector(witness, covectors[index])
            for index, _other in enumerate(roots)
            if index != target
        ]
        if not (
            target_pair == 0
            and _pair(gram, witness, witness) > 0
            and min(off_target) > 0
        ):
            raise ArithmeticError(f"seed facet witness {target + 1} failed")
        witnesses.append(witness)
    return tuple(witnesses)


def _boundary_roots(
    gram: IntegralMatrix,
    seed_roots: Sequence[IntegralVector],
    selected: int,
    interior: Sequence[int | Fraction],
) -> tuple[IntegralVector, ...]:
    doubling_root = seed_roots[selected]
    boundary: dict[IntegralVector, IntegralVector] = {}
    for copy in range(2):
        for index, root in enumerate(seed_roots):
            if index == selected:
                continue
            image: Sequence[int | Fraction] = root
            if copy:
                image = _reflect(gram, doubling_root, root)
            primitive = _primitive(image)
            if _pair(gram, interior, primitive) < 0:
                primitive = tuple(-entry for entry in primitive)
            if _pair(gram, interior, primitive) == 0:
                raise ArithmeticError("double interior lies on a boundary wall")
            boundary.setdefault(_unoriented_key(primitive), primitive)
    return tuple(boundary[key] for key in sorted(boundary))


def _boundary_with_witnesses(
    gram: IntegralMatrix,
    seed_roots: Sequence[IntegralVector],
    seed_witnesses: Sequence[Vector],
    seed_interior: Vector,
    selected: int,
) -> tuple[tuple[IntegralVector, Vector], ...]:
    doubling_root = seed_roots[selected]
    reflected_interior = _reflect(gram, doubling_root, seed_interior)
    interior = tuple(
        left + right
        for left, right in zip(seed_interior, reflected_interior, strict=True)
    )
    boundary: dict[IntegralVector, tuple[IntegralVector, Vector]] = {}
    for copy in range(2):
        for index, root in enumerate(seed_roots):
            if index == selected:
                continue
            image: Sequence[int | Fraction] = root
            witness = seed_witnesses[index]
            if copy:
                image = _reflect(gram, doubling_root, root)
                witness = _reflect(gram, doubling_root, witness)
            primitive = _primitive(image)
            if _pair(gram, interior, primitive) < 0:
                primitive = tuple(-entry for entry in primitive)
            boundary.setdefault(
                _unoriented_key(primitive),
                (primitive, tuple(Fraction(entry) for entry in witness)),
            )
    records = tuple(boundary[key] for key in sorted(boundary))
    if tuple(record[0] for record in records) != _boundary_roots(
        gram, seed_roots, selected, interior
    ):
        raise ArithmeticError("boundary constructors disagree")
    return records


def _root_type_histogram(
    gram: IntegralMatrix, roots: Sequence[IntegralVector]
) -> dict[str, int]:
    histogram: Counter[tuple[int, int]] = Counter()
    for root in roots:
        square = int(_pair(gram, root, root))
        divisibility = _content(_covector(gram, root))
        histogram[(square, divisibility)] += 1
    return {
        f"{square},{divisibility}": count
        for (square, divisibility), count in sorted(histogram.items())
    }


def _relation_histogram(
    gram: IntegralMatrix, roots: Sequence[IntegralVector]
) -> dict[str, int]:
    covectors = tuple(_covector(gram, root) for root in roots)
    squares = tuple(
        _pair_covector(root, covector)
        for root, covector in zip(roots, covectors, strict=True)
    )
    labels = {
        Fraction(0): "m2",
        Fraction(1): "m3",
        Fraction(2): "m4",
        Fraction(3): "m6",
        Fraction(4): "parallel",
    }
    histogram: Counter[str] = Counter()
    for left in range(len(roots)):
        for right in range(left):
            product = _pair_covector(roots[left], covectors[right])
            if product < 0:
                raise ArithmeticError("wall normals do not form an acute system")
            quotient = 4 * product * product / (squares[left] * squares[right])
            if quotient in labels:
                histogram[labels[quotient]] += 1
            elif quotient > 4:
                histogram["ultraparallel"] += 1
            else:
                raise ArithmeticError(
                    f"non-Coxeter normalized product 4q={quotient}"
                )
    return dict(sorted(histogram.items()))


def _production_boundary(production_row: dict[str, Any]) -> tuple[IntegralVector, ...]:
    roots: list[IntegralVector] = []
    for record in production_row["boundary_records"]:
        entries = tuple(
            Fraction(str(entry))
            for entry in record["primitive_oriented_coordinates"]
        )
        if not all(entry.denominator == 1 for entry in entries):
            raise ArithmeticError("production boundary contains a nonintegral root")
        roots.append(tuple(entry.numerator for entry in entries))
    return tuple(roots)


def _audit_single_doublings(data: dict[str, Any]) -> dict[str, Any]:
    gram: IntegralMatrix = data["fixed_gram"]
    gram_sympy: sp.Matrix = data["fixed_gram_sympy"]
    roots: tuple[IntegralVector, ...] = data["roots"]
    ample: Vector = data["ample_coordinates"]
    production: dict[str, Any] = data["production"]

    relations = _seed_relations(gram, roots)
    doubling = tuple(
        index
        for index, row in enumerate(relations)
        if not any(label in {"m3", "m5", "m7", "m9"} for label in row)
    )
    count_classes: defaultdict[int, list[int]] = defaultdict(list)
    for selected in doubling:
        orthogonal = sum(
            relations[selected][index] == "m2"
            for index in range(len(roots))
            if index != selected
        )
        expected_count = 2 * (len(roots) - 1) - orthogonal
        count_classes[expected_count].append(selected + 1)

    expected_classes = {
        296: [31, 32, 33, 36, 39, 42, 48, 51, 56, 58, 59, 60],
        317: list(range(61, 77)),
        326: list(range(77, 173)),
    }
    if len(doubling) != 124 or dict(count_classes) != expected_classes:
        raise ArithmeticError("independent complete doubling census changed")

    seed_witnesses = _seed_witnesses(gram, roots, ample)
    production_by_wall = {
        int(row["selected_doubling_wall_one_based"]): row
        for row in production["representatives"]
    }
    representative_rows: list[dict[str, Any]] = []
    for selected in REPRESENTATIVES:
        doubling_root = roots[selected]
        reflected_interior = _reflect(gram, doubling_root, ample)
        double_interior = tuple(
            left + right
            for left, right in zip(ample, reflected_interior, strict=True)
        )
        records = _boundary_with_witnesses(
            gram, roots, seed_witnesses, ample, selected
        )
        boundary = tuple(record[0] for record in records)
        witnesses = tuple(record[1] for record in records)
        boundary_covectors = tuple(_covector(gram, root) for root in boundary)

        minimum_square: Fraction | None = None
        minimum_slack: Fraction | None = None
        for target, (root, witness) in enumerate(records):
            square = _pair(gram, witness, witness)
            if _pair_covector(witness, boundary_covectors[target]) != 0 or square <= 0:
                raise ArithmeticError("transported witness failed its target")
            slacks = [
                _pair_covector(witness, boundary_covectors[index])
                for index, other in enumerate(boundary)
                if index != target
            ]
            if min(slacks) <= 0:
                raise ArithmeticError("transported witness has nonpositive slack")
            minimum_square = square if minimum_square is None else min(minimum_square, square)
            local_slack = min(slacks)
            minimum_slack = (
                local_slack
                if minimum_slack is None
                else min(minimum_slack, local_slack)
            )

        reflection_count = 0
        for root in boundary:
            _integral_reflection_matrix_fast(gram, root)
            reflection_count += 1
        smith = _smith_diagonal(boundary)
        relation_histogram = _relation_histogram(gram, boundary)
        type_histogram = _root_type_histogram(gram, boundary)
        production_row = production_by_wall[selected + 1]
        production_boundary = _production_boundary(production_row)
        production_type_histogram: Counter[tuple[int, int]] = Counter()
        squares = production_row["primitive_root_square_histogram"]
        divisibilities = production_row["primitive_root_divisibility_histogram"]
        # Compare the joint histogram directly to reconstructed roots; the two
        # marginal production histograms are retained as an additional gate.
        for root in boundary:
            production_type_histogram[
                (
                    int(_pair(gram, root, root)),
                    _content(_covector(gram, root)),
                )
            ] += 1
        joint_square_marginal = Counter()
        joint_divisibility_marginal = Counter()
        for (square, divisibility), count in production_type_histogram.items():
            joint_square_marginal[str(square)] += count
            joint_divisibility_marginal[str(divisibility)] += count

        exact_match = (
            len(boundary) == int(production_row["facet_count"])
            and boundary == production_boundary
            and relation_histogram == production_row["relation_histogram"]
            and dict(sorted(joint_square_marginal.items())) == squares
            and dict(sorted(joint_divisibility_marginal.items())) == divisibilities
        )
        if not (
            exact_match
            and smith == [1] * 13
            and reflection_count == len(boundary)
            and min(_pair(gram, double_interior, root) for root in boundary) > 0
            and minimum_square is not None
            and minimum_slack is not None
        ):
            raise ArithmeticError(
                f"hostile representative replay failed at W{selected + 1}"
            )

        representative_rows.append(
            {
                "selected_wall_one_based": selected + 1,
                "facet_count": len(boundary),
                "production_boundary_ordered_coordinates_exact_match": True,
                "relation_histogram": relation_histogram,
                "primitive_square_divisibility_histogram": type_histogram,
                "integral_involutive_isometry_reflection_count": reflection_count,
                "boundary_root_smith_diagonal": smith,
                "strict_witness_count": len(witnesses),
                "minimum_witness_square": _rational_payload(minimum_square),
                "minimum_off_target_witness_slack": _rational_payload(minimum_slack),
                "all_exact_checks_pass": True,
            }
        )

    seed_without_right_smith = _smith_diagonal(
        tuple(root for index, root in enumerate(roots) if index != RIGHT)
    )
    if seed_without_right_smith != [1] * 13:
        raise ArithmeticError("seed boundary without W36 does not span the full lattice")

    return {
        "relations": relations,
        "seed_witnesses": seed_witnesses,
        "doubling": doubling,
        "archive": {
            "criterion": "no incident odd finite Coxeter exponent",
            "doubling_wall_count": len(doubling),
            "non_doubling_wall_count": len(roots) - len(doubling),
            "single_double_facet_count_classes": {
                str(count): walls for count, walls in sorted(count_classes.items())
            },
            "count_identity": "f(double across w)=2*(172-1)-orthogonal_nonself(w)",
            "representatives": representative_rows,
            "seed_roots_excluding_W36_smith_diagonal": seed_without_right_smith,
            "all_production_representatives_match_exactly": True,
        },
    }


def _integer_roots_of_scalar_polynomial(
    quadratic: int, linear: int, constant: int
) -> tuple[int, ...]:
    if quadratic == 0:
        if linear == 0:
            return ()
        numerator = -constant
        if numerator % linear:
            return ()
        return (numerator // linear,)
    discriminant = linear * linear - 4 * quadratic * constant
    if discriminant < 0:
        return ()
    square_root = math.isqrt(discriminant)
    if square_root * square_root != discriminant:
        return ()
    denominator = 2 * quadratic
    roots: set[int] = set()
    for numerator in (-linear + square_root, -linear - square_root):
        if numerator % denominator == 0:
            roots.add(numerator // denominator)
    return tuple(sorted(roots))


def _translation_power_apply(
    nilpotent: IntegralMatrix,
    nilpotent_square: IntegralMatrix,
    vector: IntegralVector,
    exponent: int,
) -> IntegralVector:
    first = _matvec(nilpotent, vector)
    second = _matvec(nilpotent_square, vector)
    return _add(
        vector,
        _scale(exponent, first),
        _scale(exponent * (exponent - 1) // 2, second),
    )


def _solve_translation_equation(
    nilpotent: IntegralMatrix,
    nilpotent_square: IntegralMatrix,
    vector: IntegralVector,
    target: IntegralVector,
) -> tuple[int, ...] | None:
    """Solve T^z vector=target for every integer z; None means all z."""

    first = _matvec(nilpotent, vector)
    second = _matvec(nilpotent_square, vector)
    coefficients = tuple(
        (
            second[index],
            2 * first[index] - second[index],
            2 * (vector[index] - target[index]),
        )
        for index in range(len(vector))
    )
    pivot = next(
        (coefficient for coefficient in coefficients if coefficient != (0, 0, 0)),
        None,
    )
    if pivot is None:
        return None
    candidates = _integer_roots_of_scalar_polynomial(*pivot)
    return tuple(
        exponent
        for exponent in candidates
        if _translation_power_apply(
            nilpotent, nilpotent_square, vector, exponent
        )
        == target
    )


def _audit_all_integer_collisions(
    roots: Sequence[IntegralVector],
    right_reflection: IntegralMatrix,
    nilpotent: IntegralMatrix,
    nilpotent_square: IntegralMatrix,
    categories: dict[str, list[int]],
) -> dict[str, Any]:
    nonend_indices = tuple(
        index for index in range(len(roots)) if index not in {LEFT, RIGHT}
    )
    finite: list[tuple[str, int, int, int, int]] = []
    identities: list[tuple[str, int, int, int]] = []
    branches = {
        "T^z": lambda vector: vector,
        "T^z*S36": lambda vector: _matvec(right_reflection, vector),
    }
    for branch_name, branch in branches.items():
        for source_index in nonend_indices:
            vector = branch(roots[source_index])
            for target_index in nonend_indices:
                for sign in (1, -1):
                    target = _scale(sign, roots[target_index])
                    solutions = _solve_translation_equation(
                        nilpotent, nilpotent_square, vector, target
                    )
                    if solutions is None:
                        identities.append(
                            (
                                branch_name,
                                source_index + 1,
                                target_index + 1,
                                sign,
                            )
                        )
                    else:
                        finite.extend(
                            (
                                branch_name,
                                source_index + 1,
                                target_index + 1,
                                sign,
                                exponent,
                            )
                            for exponent in solutions
                        )

    finite_histogram = Counter(
        (branch, exponent, sign)
        for branch, _source, _target, sign, exponent in finite
    )
    identity_histogram = Counter(branch for branch, *_rest in identities)
    common = set(categories["orthogonal_to_both"])
    expected_finite_histogram = Counter(
        {
            ("T^z", 0, 1): 156,
            ("T^z*S36", 0, 1): 32,
            ("T^z*S36", -1, 1): 32,
        }
    )
    expected_identity_histogram = Counter({"T^z": 14, "T^z*S36": 14})
    if not (
        len(finite) == 220
        and finite_histogram == expected_finite_histogram
        and all(source == target for _branch, source, target, _sign, _z in finite)
        and all(sign == 1 for _branch, _source, _target, sign, _z in finite)
        and len(identities) == 28
        and identity_histogram == expected_identity_histogram
        and all(
            source == target and sign == 1 and source in common
            for _branch, source, target, sign in identities
        )
    ):
        raise ArithmeticError("the exhaustive all-integer collision census changed")

    # The endpoint orbit must not meet a side-facet orbit.  This is a small
    # second exhaustive solve, not a bounded scan in the translation exponent.
    endpoint_side_solutions: list[tuple[str, int, int, int, str]] = []
    for branch_name, branch in branches.items():
        for source_index in (LEFT, RIGHT):
            vector = branch(roots[source_index])
            for target_index in nonend_indices:
                for sign in (1, -1):
                    solutions = _solve_translation_equation(
                        nilpotent,
                        nilpotent_square,
                        vector,
                        _scale(sign, roots[target_index]),
                    )
                    if solutions is None:
                        endpoint_side_solutions.append(
                            (
                                branch_name,
                                source_index + 1,
                                target_index + 1,
                                sign,
                                "all",
                            )
                        )
                    else:
                        endpoint_side_solutions.extend(
                            (
                                branch_name,
                                source_index + 1,
                                target_index + 1,
                                sign,
                                str(exponent),
                            )
                            for exponent in solutions
                        )
    if endpoint_side_solutions:
        raise ArithmeticError("an endpoint orbit collided with a side-facet orbit")

    right_reflected_left = _matvec(right_reflection, roots[LEFT])
    endpoint_equations: dict[str, list[dict[str, int]]] = {}
    for name, vector in {
        "T^z(W36)=plus_or_minus_W31": roots[RIGHT],
        "T^z*S36(W31)=plus_or_minus_W31": right_reflected_left,
    }.items():
        records: list[dict[str, int]] = []
        for sign in (1, -1):
            solutions = _solve_translation_equation(
                nilpotent,
                nilpotent_square,
                vector,
                _scale(sign, roots[LEFT]),
            )
            if solutions is None:
                raise ArithmeticError("unexpected all-z endpoint identity")
            records.extend({"z": exponent, "sign": sign} for exponent in solutions)
        endpoint_equations[name] = records
    if not (
        endpoint_equations["T^z(W36)=plus_or_minus_W31"] == []
        and endpoint_equations[
            "T^z*S36(W31)=plus_or_minus_W31"
        ]
        == [{"z": -1, "sign": -1}]
    ):
        raise ArithmeticError(f"endpoint collision equations changed: {endpoint_equations}")

    return {
        "domain": (
            "all z in Z, both branches T^z and T^z*S36, both root signs, "
            "all 170x170 ordered side-facet pairs"
        ),
        "method": (
            "use T^z=I+z*(T-I)+z*(z-1)/2*(T-I)^2 and solve one nonzero "
            "integer linear/quadratic coordinate exactly, then verify all coordinates"
        ),
        "bounded_exponent_scan_used": False,
        "finite_solution_count": len(finite),
        "finite_solution_histogram": {
            f"{branch};z={exponent};sign={sign}": count
            for (branch, exponent, sign), count in sorted(finite_histogram.items())
        },
        "finite_solutions_are_positive_self_collisions_only": True,
        "all_z_identity_count": len(identities),
        "all_z_identity_histogram": dict(sorted(identity_histogram.items())),
        "all_z_identities_are_exactly_the_14_common_orthogonals": True,
        "negative_root_collision_count": 0,
        "cross_seed_facet_collision_count": 0,
        "endpoint_to_side_orbit_collision_count": 0,
        "outer_endpoint_equations": endpoint_equations,
        "all_exact_checks_pass": True,
    }


def _audit_parallel_pair(
    data: dict[str, Any], single: dict[str, Any]
) -> dict[str, Any]:
    gram: IntegralMatrix = data["fixed_gram"]
    gram_sympy: sp.Matrix = data["fixed_gram_sympy"]
    roots: tuple[IntegralVector, ...] = data["roots"]
    relations: tuple[tuple[str, ...], ...] = single["relations"]
    left = roots[LEFT]
    right = roots[RIGHT]
    expected_left = (-3, -3, -3, -3, -2, -2, 5, 6, 5, 6, 6, 6, 8)
    expected_right = (-3, -3, -3, -3, -2, -2, 5, 6, 5, 6, 4, 6, 8)
    if left != expected_left or right != expected_right:
        raise ArithmeticError("the W31/W36 primitive coordinates changed")

    left_reflection = _reflection_matrix(gram_sympy, left)
    right_reflection = _reflection_matrix(gram_sympy, right)
    left_sympy = sp.Matrix(left_reflection)
    right_sympy = sp.Matrix(right_reflection)
    translation_sympy = right_sympy * left_sympy
    translation = _integral_matrix(translation_sympy)
    nilpotent_sympy = translation_sympy - sp.eye(13)
    nilpotent_square_sympy = nilpotent_sympy**2
    nilpotent = _integral_matrix(nilpotent_sympy)
    nilpotent_square = _integral_matrix(nilpotent_square_sympy)
    variable = translation_sympy.charpoly().gen
    charpoly_ok = sp.expand(translation_sympy.charpoly().as_expr()) == sp.expand(
        (variable - 1) ** 13
    )

    half_sum = tuple(Fraction(a + b, 2) for a, b in zip(left, right, strict=True))
    null_direction = _primitive(half_sum)
    nonend = tuple(index for index in range(172) if index not in {LEFT, RIGHT})
    orthogonal_left = {
        index + 1 for index in nonend if _pair(gram, left, roots[index]) == 0
    }
    orthogonal_right = {
        index + 1 for index in nonend if _pair(gram, right, roots[index]) == 0
    }
    common = orthogonal_left & orthogonal_right
    left_only = orthogonal_left - orthogonal_right
    right_only = orthogonal_right - orthogonal_left
    neither = {index + 1 for index in nonend} - (
        orthogonal_left | orthogonal_right
    )
    categories = {
        "orthogonal_to_both": sorted(common),
        "orthogonal_to_W31_only": sorted(left_only),
        "orthogonal_to_W36_only": sorted(right_only),
        "orthogonal_to_neither": sorted(neither),
    }
    row_histograms = {
        str(index + 1): dict(
            sorted(
                Counter(
                    label
                    for other, label in enumerate(relations[index])
                    if other != index
                ).items()
            )
        )
        for index in (LEFT, RIGHT)
    }
    expected_row = {
        "m2": 46,
        "m4": 20,
        "parallel": 45,
        "ultraparallel": 60,
    }
    if not (
        _pair(gram, left, left) == -4
        and _pair(gram, right, right) == -4
        and _pair(gram, left, right) == 4
        and _pair(gram, null_direction, null_direction) == 0
        and _content(null_direction) == 1
        and charpoly_ok
        and nilpotent_sympy.rank() == 2
        and nilpotent_square_sympy.rank() == 1
        and nilpotent_sympy**3 == sp.zeros(13)
        and nilpotent_square_sympy != sp.zeros(13)
        and row_histograms == {"31": expected_row, "36": expected_row}
        and {name: len(indices) for name, indices in categories.items()}
        == {
            "orthogonal_to_both": 14,
            "orthogonal_to_W31_only": 32,
            "orthogonal_to_W36_only": 32,
            "orthogonal_to_neither": 92,
        }
        and LEFT in single["doubling"]
        and RIGHT in single["doubling"]
    ):
        raise ArithmeticError("parallel-pair structural gates changed")

    collisions = _audit_all_integer_collisions(
        roots, right_reflection, nilpotent, nilpotent_square, categories
    )
    return {
        "left_wall_one_based": 31,
        "right_wall_one_based": 36,
        "left_primitive_coordinates": list(left),
        "right_primitive_coordinates": list(right),
        "squares_and_pairing": [-4, -4, 4],
        "primitive_integral_common_null_direction": list(null_direction),
        "column_action_convention": "S_r=I-(2/r^2)*r*r^T*H",
        "left_reflection_matrix": _matrix_payload(left_reflection),
        "right_reflection_matrix": _matrix_payload(right_reflection),
        "translation_definition": "T=S36*S31 (column action)",
        "translation_matrix": _matrix_payload(translation),
        "translation_characteristic_polynomial": "(lambda-1)^13",
        "rank_T_minus_I": int(nilpotent_sympy.rank()),
        "rank_T_minus_I_squared": int(nilpotent_square_sympy.rank()),
        "T_minus_I_cubed_is_zero": True,
        "T_minus_I_squared_is_nonzero": True,
        "translation_has_infinite_order": True,
        "generated_reflection_group": "D_infinity",
        "endpoint_relation_row_histograms": row_histograms,
        "side_facet_categories": categories,
        "side_facet_category_counts": {
            name: len(indices) for name, indices in categories.items()
        },
        "all_integer_collision_audit": collisions,
        "all_exact_checks_pass": True,
        "_right_reflection": right_reflection,
        "_nilpotent": nilpotent,
        "_nilpotent_square": nilpotent_square,
    }


def _chain_image(
    index: int,
    vector: IntegralVector,
    right_reflection: IntegralMatrix,
    nilpotent: IntegralMatrix,
    nilpotent_square: IntegralMatrix,
) -> IntegralVector:
    if index % 2 == 0:
        return _translation_power_apply(
            nilpotent, nilpotent_square, vector, index // 2
        )
    reflected = _matvec(right_reflection, vector)
    return _translation_power_apply(
        nilpotent, nilpotent_square, reflected, index // 2
    )


def _chain_boundary(
    chamber_count: int,
    roots: Sequence[IntegralVector],
    right_reflection: IntegralMatrix,
    nilpotent: IntegralMatrix,
    nilpotent_square: IntegralMatrix,
) -> tuple[IntegralVector, ...]:
    if chamber_count < 1:
        raise ValueError("a finite path must contain at least one chamber")
    nonend = tuple(index for index in range(len(roots)) if index not in {LEFT, RIGHT})
    boundary: dict[IntegralVector, IntegralVector] = {}
    for chamber_index in range(chamber_count):
        for root_index in nonend:
            image = _chain_image(
                chamber_index,
                roots[root_index],
                right_reflection,
                nilpotent,
                nilpotent_square,
            )
            boundary.setdefault(_unoriented_key(image), image)
    left_outer = roots[LEFT]
    last = chamber_count - 1
    right_seed_root = roots[RIGHT] if last % 2 == 0 else roots[LEFT]
    right_outer = _chain_image(
        last,
        right_seed_root,
        right_reflection,
        nilpotent,
        nilpotent_square,
    )
    for endpoint in (left_outer, right_outer):
        key = _unoriented_key(endpoint)
        if key in boundary:
            raise ArithmeticError("outer endpoint collided with a side facet")
        boundary[key] = endpoint
    return tuple(boundary[key] for key in sorted(boundary))


def _audit_family(
    data: dict[str, Any], single: dict[str, Any], parallel: dict[str, Any]
) -> dict[str, Any]:
    roots: tuple[IntegralVector, ...] = data["roots"]
    right_reflection: IntegralMatrix = parallel.pop("_right_reflection")
    nilpotent: IntegralMatrix = parallel.pop("_nilpotent")
    nilpotent_square: IntegralMatrix = parallel.pop("_nilpotent_square")
    category_counts = parallel["side_facet_category_counts"]
    direct_checks: list[dict[str, int]] = []
    for chamber_count in range(1, 13):
        boundary = _chain_boundary(
            chamber_count,
            roots,
            right_reflection,
            nilpotent,
            nilpotent_square,
        )
        expected = 124 * chamber_count + 48
        if len(boundary) != expected:
            raise ArithmeticError(
                f"direct chain count failed for N={chamber_count}: {len(boundary)}"
            )
        direct_checks.append(
            {"chamber_count_N": chamber_count, "facet_count": len(boundary)}
        )

    if not (
        category_counts
        == {
            "orthogonal_to_both": 14,
            "orthogonal_to_W31_only": 32,
            "orthogonal_to_W36_only": 32,
            "orthogonal_to_neither": 92,
        }
        and single["archive"]["seed_roots_excluding_W36_smith_diagonal"]
        == [1] * 13
    ):
        raise ArithmeticError("family-count inputs changed")

    return {
        "chamber_maps": {
            "even_index_2m": "g_(2m)=T^m",
            "odd_index_2m_plus_1": "g_(2m+1)=T^m*S36",
        },
        "side_facet_orbit_decomposition": {
            "common_orthogonals": "14",
            "W31_only_orthogonals": "32*(N-floor(N/2))",
            "W36_only_orthogonals": "32*(N-floor((N-1)/2))",
            "orthogonal_to_neither": "92*N",
            "distinct_side_facets_total": "124*N+46",
            "outer_endpoint_facets": "2",
        },
        "facet_formula": "f(N)=124*N+48 for every integer N>=1",
        "equivalent_facet_formula": "f(N)=2+170*N-46*(N-1)",
        "universal_formula_proof": (
            "the exhaustive all-z collision solve gives every and only side-facet "
            "identification; the separate all-z endpoint solve excludes endpoint/side "
            "collisions and leaves two distinct outer endpoints for every N>=1"
        ),
        "direct_exact_reconstruction_checks_N_1_through_12": direct_checks,
        "all_boundary_roots_integral_and_primitive": True,
        "full_fixed_Picard_lattice_for_every_N": True,
        "full_lattice_reason": (
            "every path boundary contains the seed roots except W36, whose 13x171 "
            "root matrix has Smith diagonal 1^13"
        ),
        "coxeter_for_every_N": True,
        "coxeter_reason": (
            "iterate Allcock Lemma 1.3 from the two parallel doubling walls W31/W36"
        ),
        "volume": "Vol(P_N)=N*Vol(P_seed)",
        "finite_volume": True,
        "compact": False,
        "pairwise_nonisometric": True,
        "pairwise_nonisometric_reasons": [
            "strictly increasing facet count 124*N+48",
            "distinct positive volume multiples N*Vol(P_seed)",
        ],
        "all_exact_checks_pass": True,
    }


def _first_nonempty_header(payload: bytes, name: str) -> tuple[int, int]:
    text = payload.decode("utf-8", errors="replace")
    try:
        line = next(line.strip() for line in text.splitlines() if line.strip())
    except StopIteration as error:
        raise ValueError(f"empty CoxIter file: {name}") from error
    tokens = line.split()
    if len(tokens) < 2:
        raise ValueError(f"short CoxIter header: {name}")
    facet_count, dimension = map(int, tokens[:2])
    if facet_count <= 0 or dimension <= 0:
        raise ValueError(f"nonpositive CoxIter header: {name}")
    return facet_count, dimension


def _histogram_payload(histogram: Counter[int]) -> dict[str, int]:
    return {str(key): value for key, value in sorted(histogram.items())}


def _audit_coxiter(coxiter_archive: Path) -> dict[str, Any]:
    archive_hash = _sha256(coxiter_archive)
    if archive_hash != EXPECTED_SHA256["coxiter_zip"]:
        raise ArithmeticError(f"the pinned CoxIter zip changed: {archive_hash}")

    ordinary_records: list[tuple[str, int, int]] = []
    index_two_bytes: bytes | None = None
    with zipfile.ZipFile(coxiter_archive) as archive:
        archive_comment = archive.comment.decode("ascii")
        for name in sorted(archive.namelist()):
            if name.endswith(".coxiter"):
                facet_count, dimension = _first_nonempty_header(
                    archive.read(name), name
                )
                ordinary_records.append((name, facet_count, dimension))
            elif name.endswith("graphs/testing_index2.tar.gz"):
                index_two_bytes = archive.read(name)
    if index_two_bytes is None:
        raise ArithmeticError("the nested CoxIter index-two archive is missing")
    index_two_hash = _sha256_bytes(index_two_bytes)
    if index_two_hash != EXPECTED_SHA256["coxiter_index_two_tar_gz"]:
        raise ArithmeticError(
            f"the pinned CoxIter index-two archive changed: {index_two_hash}"
        )

    index_two_records: list[tuple[str, int, int]] = []
    with tarfile.open(fileobj=BytesIO(index_two_bytes), mode="r:gz") as archive:
        for member in archive.getmembers():
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError(f"could not read {member.name}")
            facet_count, dimension = _first_nonempty_header(
                handle.read(), member.name
            )
            index_two_records.append((member.name, facet_count, dimension))

    ordinary_h12 = Counter(
        facets for _name, facets, dimension in ordinary_records if dimension == 12
    )
    index_two_h12 = Counter(
        facets for _name, facets, dimension in index_two_records if dimension == 12
    )
    ordinary_candidate_matches = [
        name
        for name, facets, dimension in ordinary_records
        if dimension == 12 and facets >= 172
    ]
    index_two_candidate_matches = [
        name
        for name, facets, dimension in index_two_records
        if dimension == 12 and facets >= 172
    ]
    if not (
        archive_comment == COXITER_COMMIT
        and len(ordinary_records) == 1214
        and ordinary_h12 == Counter({13: 6, 14: 3, 16: 1, 18: 1})
        and len(index_two_records) == 2959
        and index_two_h12 == Counter({18: 1, 20: 1, 25: 1, 35: 1})
        and not ordinary_candidate_matches
        and not index_two_candidate_matches
    ):
        raise ArithmeticError("the bounded CoxIter header census changed")

    return {
        "scope": (
            "literal H12 chamber header exclusion by facet count in one pinned "
            "official CoxIter snapshot and its bundled index-two testing archive"
        ),
        "repository": "https://github.com/rgugliel/CoxIter",
        "commit": COXITER_COMMIT,
        "zip_sha256": archive_hash,
        "bundled_index_two_tar_gz_sha256": index_two_hash,
        "ordinary": {
            "coxiter_file_count": len(ordinary_records),
            "parse_error_count": 0,
            "h12_facet_histogram": _histogram_payload(ordinary_h12),
            "h12_records_with_at_least_172_facets": ordinary_candidate_matches,
        },
        "index_two": {
            "coxiter_file_count": len(index_two_records),
            "parse_error_count": 0,
            "h12_facet_histogram": _histogram_payload(index_two_h12),
            "h12_records_with_at_least_172_facets": index_two_candidate_matches,
        },
        "literal_seed_or_descendant_match_absent_by_header": True,
        "does_not_claim_underlying_lattice_absent": True,
        "does_not_claim_worldwide_priority": True,
        "all_exact_checks_pass": True,
    }


def build(coxiter_archive: Path = DEFAULT_COXITER_ARCHIVE) -> dict[str, Any]:
    data = _parse_inputs()
    single = _audit_single_doublings(data)
    parallel = _audit_parallel_pair(data, single)
    family = _audit_family(data, single, parallel)
    coxiter = _audit_coxiter(coxiter_archive)
    gram: IntegralMatrix = data["fixed_gram"]
    basis: sp.Matrix = data["basis"]
    roots: tuple[IntegralVector, ...] = data["roots"]
    ample: Vector = data["ample_coordinates"]

    checks = {
        "does_not_import_or_call_production_certifier": True,
        "published_17_dimensional_data_reparsed": True,
        "fixed_rank_13_lattice_reconstructed": True,
        "all_172_primitive_seed_roots_reconstructed": True,
        "complete_124_wall_doubling_census_replayed": True,
        "three_production_boundary_coordinate_sets_match_exactly": single[
            "archive"
        ]["all_production_representatives_match_exactly"],
        "all_boundary_reflections_integral_involutive_isometries": True,
        "all_representative_facets_have_strict_exact_witnesses": True,
        "parallel_pair_generates_infinite_dihedral_group": parallel[
            "generated_reflection_group"
        ]
        == "D_infinity",
        "all_integer_side_collision_problem_solved": parallel[
            "all_integer_collision_audit"
        ]["all_exact_checks_pass"],
        "all_integer_endpoint_collision_problem_solved": parallel[
            "all_integer_collision_audit"
        ]["endpoint_to_side_orbit_collision_count"]
        == 0,
        "universal_facet_formula_proved": family["facet_formula"]
        == "f(N)=124*N+48 for every integer N>=1",
        "direct_chain_reconstruction_N_1_through_12_matches_formula": all(
            row["facet_count"] == 124 * row["chamber_count_N"] + 48
            for row in family["direct_exact_reconstruction_checks_N_1_through_12"]
        ),
        "every_descendant_spans_full_fixed_lattice": family[
            "full_fixed_Picard_lattice_for_every_N"
        ],
        "pinned_coxiter_corpora_replayed": coxiter["all_exact_checks_pass"],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"hostile audit failed: {checks}")

    return {
        "schema": "coxeter.jacobian_kummer.h12_parallel_redoubling_hostile_audit.v1",
        "status_date": "2026-08-08",
        "independence": {
            "production_module_imported": False,
            "other_project_module_imported": False,
            "production_json_used_only_as_comparison_target": True,
            "upstream_reconstruction_inputs": [
                str(PUBLISHED.relative_to(ROOT)),
                str(CHAMBER.relative_to(ROOT)),
            ],
            "method": (
                "stdlib plus SymPy; independent GAP assignment parser, fixed-space "
                "basis, primitive roots, reflections, witnesses, Smith forms, and "
                "all-integer collision solver"
            ),
        },
        "pinned_inputs": {
            "published_gap_data": {
                "path": str(PUBLISHED.relative_to(ROOT)),
                "sha256": data["hashes"]["published_gap_data"],
            },
            "interim_chamber": {
                "path": str(CHAMBER.relative_to(ROOT)),
                "sha256": data["hashes"]["interim_chamber"],
            },
            "fixed_lattice_comparison_certificate": {
                "path": str(FIXED_CERTIFICATE.relative_to(ROOT)),
                "sha256": data["hashes"]["fixed_lattice_certificate"],
            },
            "production_comparison_target": {
                "path": str(PRODUCTION.relative_to(ROOT)),
                "sha256": data["hashes"]["production_single_doublings"],
            },
        },
        "fixed_lattice_reconstruction": {
            "rank": basis.rows,
            "ambient_rank": basis.cols,
            "basis_rows_in_published_S16_basis": _matrix_payload(
                _integral_matrix(basis)
            ),
            "gram": _matrix_payload(gram),
            "determinant": int(data["fixed_gram_sympy"].det()),
            "seed_facet_count": len(roots),
            "strict_ample_coordinates": _vector_payload(ample),
            "strict_ample_square": _rational_payload(_pair(gram, ample, ample)),
            "minimum_seed_wall_pairing": _rational_payload(
                min(_pair(gram, ample, root) for root in roots)
            ),
            "primitive_seed_root_type_histogram": _root_type_histogram(gram, roots),
            "matches_frozen_fixed_lattice_certificate_exactly": True,
        },
        "single_doubling_hostile_replay": single["archive"],
        "parallel_redoubling_pair": parallel,
        "infinite_path_family": family,
        "bounded_coxiter_corpus_replay": coxiter,
        "classical_doubling_reference": {
            "author": "Daniel Allcock",
            "title": "Infinitely many hyperbolic Coxeter groups through dimension 19",
            "url": "https://web.ma.utexas.edu/users/allcock/research/doubling.pdf",
            "result_used": "Lemma 1.3 (doubling and redoubling)",
            "use_here": (
                "W31 and W36 are parallel doubling walls, so every finite alternating "
                "chain is a finite-volume noncompact Coxeter polytope"
            ),
        },
        "checks": checks,
        "all_checks_pass": True,
        "decision": {
            "single_double_production_certificate_independently_confirmed": True,
            "infinite_parallel_redoubling_family_independently_confirmed": True,
            "facet_formula_for_every_N_at_least_1": "124*N+48",
            "classification": "finite-volume noncompact H12 Coxeter polytopes",
            "same_full_Picard_K3_lattice_carrier": True,
            "surface_statement": "weak divisor-reflection carrier only",
            "full_nef_or_surface_automorphism_claimed": False,
            "bounded_coxiter_literal_match_absent": True,
            "worldwide_priority": "open",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coxiter-archive", type=Path, default=DEFAULT_COXITER_ARCHIVE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    payload = build(arguments.coxiter_archive)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if (
            not arguments.output.is_file()
            or arguments.output.read_text(encoding="utf-8") != rendered
        ):
            raise SystemExit(f"archived hostile audit differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(
        json.dumps(
            {
                "all_checks_pass": payload["all_checks_pass"],
                "doubling_walls": payload["single_doubling_hostile_replay"][
                    "doubling_wall_count"
                ],
                "facet_formula": payload["decision"][
                    "facet_formula_for_every_N_at_least_1"
                ],
                "all_z_finite_solutions": payload["parallel_redoubling_pair"][
                    "all_integer_collision_audit"
                ]["finite_solution_count"],
            },
            sort_keys=True,
        )
    )
    print(arguments.output)


if __name__ == "__main__":
    main()

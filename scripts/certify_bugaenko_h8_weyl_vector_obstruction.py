#!/usr/bin/env python3
"""Certify that Bugaenko's eleven walls have no generalized Weyl vector.

The only input is the frozen marked-polytope certificate
``data/known_examples/bugaenko_compact_h8.json``.  Its eleven normalized
walls have square two over ``F=Q(phi)``.  A generalized lattice Weyl vector
would therefore have to satisfy

    (rho, r_i) = -1,  i=0,...,10.

The first nine walls are a basis over ``F``.  We solve their nine equations
exactly, then evaluate the two remaining equations.  Both have the same
nonzero residual.  Thus no solution exists even in the ambient vector space
``L tensor_O F``; integrality is not the obstruction.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

from coxeter_polytopes.exact.golden import Golden, ZERO
from coxeter_polytopes.exact.matrix import determinant, rank


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "known_examples" / "bugaenko_compact_h8.json"
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "bugaenko_h8_weyl_vector_obstruction.json"
)

Vector = list[Golden]
Matrix = list[list[Golden]]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_sha256(payload: dict[str, Any]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _decode_matrix(payload: Sequence[Sequence[dict[str, str]]]) -> Matrix:
    return [[Golden.from_json(entry) for entry in row] for row in payload]


def _encode(value: Golden) -> dict[str, str]:
    return value.to_json()


def _encode_vector(vector: Sequence[Golden]) -> list[dict[str, str]]:
    return [_encode(entry) for entry in vector]


def _pair(left: Sequence[Golden], gram: Matrix, right: Sequence[Golden]) -> Golden:
    return sum(
        (
            left[i] * gram[i][j] * right[j]
            for i in range(len(left))
            for j in range(len(right))
        ),
        start=ZERO,
    )


def _equation_row(root: Sequence[Golden], gram: Matrix) -> Vector:
    """Return the row representing ``rho -> (rho, root)``."""

    return [
        sum((gram[i][j] * root[j] for j in range(len(root))), start=ZERO)
        for i in range(len(root))
    ]


def _solve_square(matrix: Matrix, right: Sequence[Golden]) -> Vector:
    """Solve a nonsingular square system by exact Gauss--Jordan reduction."""

    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise ValueError("the coefficient matrix must be nonempty and square")
    if len(right) != size:
        raise ValueError("the right-hand side has the wrong length")

    work = [
        [Golden.coerce(entry) for entry in row] + [Golden.coerce(right[index])]
        for index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]), None
        )
        if pivot is None:
            raise ArithmeticError("the coefficient matrix is singular")
        work[column], work[pivot] = work[pivot], work[column]
        pivot_value = work[column][column]
        work[column] = [entry / pivot_value for entry in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                work[row][index] - factor * work[column][index]
                for index in range(size + 1)
            ]
    return [work[index][-1] for index in range(size)]


def analyze() -> dict[str, Any]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    gram = _decode_matrix(source["lattice_gram"])
    roots = _decode_matrix(source["simple_roots"])
    dimension = int(source["dimension"])
    ambient_rank = dimension + 1
    pivot_indices = list(range(ambient_rank))
    pivot_roots = [roots[index] for index in pivot_indices]
    pivot_determinant = determinant(pivot_roots)

    root_norms = [_pair(root, gram, root) for root in roots]
    equations = [_equation_row(root, gram) for root in pivot_roots]
    candidate = _solve_square(equations, [Golden(-1)] * ambient_rank)
    pairings = [_pair(candidate, gram, root) for root in roots]
    residuals = [pairing + 1 for pairing in pairings]
    failing_indices = [index for index, residual in enumerate(residuals) if residual]
    augmented = [
        _equation_row(root, gram) + [Golden(-1)]
        for root in roots
    ]
    full_equations = [_equation_row(root, gram) for root in roots]

    checks = {
        "ambient_signature_rank_matches_dimension": len(gram) == ambient_rank,
        "eleven_normalized_square_two_walls": len(roots) == 11
        and all(norm == 2 for norm in root_norms),
        "first_nine_walls_are_a_field_basis": bool(pivot_determinant),
        "first_nine_weyl_equations_hold": all(
            residual == 0 for residual in residuals[:ambient_rank]
        ),
        "last_two_weyl_equations_fail": failing_indices == [9, 10],
        "last_two_residuals_equal_and_nonzero": residuals[9] == residuals[10]
        and bool(residuals[9]),
        "coefficient_rank_is_nine": rank(full_equations) == ambient_rank,
        "augmented_rank_is_ten": rank(augmented) == ambient_rank + 1,
    }
    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise ArithmeticError(f"Bugaenko Weyl-vector audit failed: {failed}")

    payload: dict[str, Any] = {
        "schema_version": "1.0.0",
        "certificate_kind": "generalized_weyl_vector_inconsistency",
        "source": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": _sha256(SOURCE),
            "status": source["completeness"]["evidence"]["status"],
            "dimension": dimension,
            "facet_count": len(roots),
        },
        "field": {
            "name": "Q(phi)",
            "relation": "phi^2=phi+1",
            "integer_ring": "Z[phi]",
        },
        "system": {
            "ambient_rank": ambient_rank,
            "wall_count": len(roots),
            "root_norms": _encode_vector(root_norms),
            "required_pairing": _encode(Golden(-1)),
            "equation": "(rho,r_i)=-r_i^2/2=-1 for every normalized wall",
            "coefficient_rank": rank(full_equations),
            "augmented_rank": rank(augmented),
        },
        "unique_nine_wall_candidate": {
            "pivot_wall_indices_zero_based": pivot_indices,
            "pivot_root_coordinate_determinant": _encode(pivot_determinant),
            "coordinates_in_archived_ambient_basis": _encode_vector(candidate),
            "square": _encode(_pair(candidate, gram, candidate)),
            "pairings_with_all_eleven_walls": _encode_vector(pairings),
            "residuals_pairing_plus_one": _encode_vector(residuals),
            "failing_wall_indices_zero_based": failing_indices,
            "common_nonzero_residual": _encode(residuals[9]),
            "common_residual_field_norm": str(residuals[9].norm()),
        },
        "checks": checks,
        "conclusion": {
            "generalized_weyl_vector_exists_over_Q_phi": False,
            "integral_weyl_vector_exists_over_Z_phi": False,
            "reason": (
                "the first nine wall equations have a unique solution, while "
                "each remaining wall has residual 4-8*phi of field norm -80"
            ),
            "scope": (
                "direct linear-algebra obstruction for the marked Bugaenko "
                "H8 chamber; not an application or extension of the integral "
                "Z-lattice classification of Gritsenko--Nikulin"
            ),
        },
        "all_checks_pass": True,
    }
    payload["certificate_sha256"] = _canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check-archive",
        action="store_true",
        help="compare the exact replay with the frozen JSON certificate",
    )
    args = parser.parse_args()
    payload = analyze()
    if args.check_archive:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if payload != archived:
            raise SystemExit("archive mismatch")
        print(f"archive replay passed: {payload['certificate_sha256']}")
        return
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

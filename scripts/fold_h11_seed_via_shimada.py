"""Fold the H11 Pell isometry into the K3 eigen-nef chamber exactly.

An infinite-order field-lattice isometry that stabilizes one complete
norm-two chamber proves that the reflection subgroup has infinite index.  We
start with the explicit H11 Pell isometry, use Shimada's complete induced-nef
atlas to locate its image point, and whenever a genuine ``(-2)`` boundary is
encountered reflect the image point back across the corresponding restricted
field mirror.  Auxiliary Borcherds walls merely change the physical induced
chamber and are crossed without modifying the arithmetic isometry.

The search is conclusive only if it terminates in a strict induced-chamber
interior and an exact trace-of-a-power test proves that the folded matrix has
infinite projective order.  State or reflection limits are reported as
incomplete.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.geometric.golden_halfspaces import (
    golden_dot,
    golden_quadratic_value,
    inverse_golden_matrix,
)
from scripts.coupled_discriminant11_candidate import field_gram
from scripts.explore_shimada_disc11_global_eigen_nef import (
    _Q_GRAM,
    _physical_state_key,
    _restricted_rows,
)
from scripts.h11_h19_symmetry_seeds import h11_seed
from scripts.shimada_disc11_atlas import ShimadaChamberState, ShimadaDisc11Atlas
from scripts.shimada_disc11_conjugate import (
    AMPLE_FIELD_COEFFICIENTS,
    DISC11_SOURCE_SHA256,
)


GoldenVector = tuple[Golden, ...]
GoldenMatrix = tuple[tuple[Golden, ...], ...]
WallLabel = tuple[int, int]


def _identity(size: int) -> GoldenMatrix:
    return tuple(
        tuple(Golden(int(i == j)) for j in range(size)) for i in range(size)
    )


def _transpose(matrix: GoldenMatrix) -> GoldenMatrix:
    return tuple(tuple(matrix[j][i] for j in range(len(matrix))) for i in range(len(matrix)))


def _matrix_product(left: GoldenMatrix, right: GoldenMatrix) -> GoldenMatrix:
    return tuple(
        tuple(
            sum(
                (left[i][k] * right[k][j] for k in range(len(right))),
                Golden(0),
            )
            for j in range(len(right[0]))
        )
        for i in range(len(left))
    )


def _row_product(row: GoldenVector, matrix: GoldenMatrix) -> GoldenVector:
    return tuple(
        sum((row[i] * matrix[i][j] for i in range(len(row))), Golden(0))
        for j in range(len(matrix[0]))
    )


def _matrix_power(matrix: GoldenMatrix, exponent: int) -> GoldenMatrix:
    if exponent < 0:
        raise ValueError("matrix exponent must be nonnegative")
    result = _identity(len(matrix))
    base = matrix
    power = exponent
    while power:
        if power & 1:
            result = _matrix_product(result, base)
        base = _matrix_product(base, base)
        power >>= 1
    return result


def _seed_row_action() -> GoldenMatrix:
    column_seed = h11_seed()
    block = tuple(tuple(entry for entry in row) for row in column_seed)
    column_action = tuple(
        tuple(
            block[i][j]
            if i < 2 and j < 2
            else Golden(int(i == j))
            for j in range(10)
        )
        for i in range(10)
    )
    # h11_seed satisfies A^T H A = H for column coordinates.  Its transpose
    # is therefore the corresponding right action on row coordinates.
    return _transpose(column_action)


def _isometry_verified(matrix: GoldenMatrix) -> bool:
    gram = tuple(tuple(entry for entry in row) for row in field_gram())
    return _matrix_product(_matrix_product(matrix, gram), _transpose(matrix)) == gram


def _field_pair(left: GoldenVector, right: GoldenVector) -> Golden:
    gram = tuple(tuple(entry for entry in row) for row in field_gram())
    return sum(
        (
            left[i] * gram[i][j] * right[j]
            for i in range(10)
            for j in range(10)
        ),
        Golden(0),
    )


def _reflection_in_functional(functional: GoldenVector) -> GoldenMatrix:
    gram = tuple(tuple(entry for entry in row) for row in field_gram())
    inverse = inverse_golden_matrix(gram)
    normal = tuple(
        sum(
            (inverse[i][j] * functional[j] for j in range(10)),
            Golden(0),
        )
        for i in range(10)
    )
    norm = golden_dot(functional, normal)
    if norm.sign() <= 0:
        raise ValueError("a K3 boundary functional has non-spacelike dual normal")
    reflection = tuple(
        tuple(
            Golden(int(i == j)) - 2 * functional[i] * normal[j] / norm
            for j in range(10)
        )
        for i in range(10)
    )
    if not _isometry_verified(reflection):
        raise ArithmeticError("the reconstructed field reflection is not an isometry")
    return reflection


@dataclass(frozen=True, slots=True)
class LocationResult:
    kind: str
    state: ShimadaChamberState
    gallery: tuple[tuple[int, int, int], ...]
    boundary_label: WallLabel | None = None
    boundary_functional: GoldenVector | None = None


@dataclass(frozen=True, slots=True)
class PointFoldResult:
    point: GoldenVector
    row_action: GoldenMatrix
    terminal: LocationResult
    records: tuple[dict[str, object], ...]


def _walk_until_boundary_or_interior(
    atlas: ShimadaDisc11Atlas,
    point: GoldenVector,
    *,
    max_gallery: int,
) -> LocationResult:
    current = atlas.representative_state(68)
    gallery = []
    seen = set()
    for _ in range(max_gallery):
        rows = _restricted_rows(atlas, current)
        signs = {label: golden_dot(row, point).sign() for label, row in rows.items()}
        if all(sign > 0 for sign in signs.values()):
            if golden_quadratic_value(_Q_GRAM, point).sign() <= 0:
                raise ValueError("the located point is not timelike")
            return LocationResult("interior", current, tuple(gallery))
        if any(sign == 0 for sign in signs.values()):
            return LocationResult("singular", current, tuple(gallery))
        key = _physical_state_key(atlas, current)
        if key in seen:
            return LocationResult("cycle", current, tuple(gallery))
        seen.add(key)
        negative = next(label for label, sign in signs.items() if sign < 0)
        orbit, member = negative
        neighbor = atlas.neighboring_state(current, orbit, member)
        if neighbor is None:
            return LocationResult(
                "root_boundary",
                current,
                tuple(gallery),
                negative,
                rows[negative],
            )
        gallery.append((current.index, orbit, member))
        current = neighbor
    return LocationResult("gallery_limit", current, tuple(gallery))


def _infinite_order_trace_power(
    matrix: GoldenMatrix,
    *,
    max_power: int = 40,
) -> tuple[int | None, Golden | None]:
    """Prove infinite order once ``|trace(g^k)| > rank`` exactly."""

    power = _identity(len(matrix))
    for exponent in range(1, max_power + 1):
        power = _matrix_product(power, matrix)
        trace = sum((power[i][i] for i in range(len(matrix))), Golden(0))
        if (trace - len(matrix)).sign() > 0 or (-trace - len(matrix)).sign() > 0:
            return exponent, trace
    return None, None


def fold_point_to_eigen_nef(
    atlas: ShimadaDisc11Atlas,
    initial_point: GoldenVector,
    *,
    initial_row_action: GoldenMatrix | None = None,
    max_root_reflections: int = 100,
    max_gallery: int = 2000,
) -> PointFoldResult:
    """Fold a timelike field point into the selected K3 eigen-nef chamber."""

    point = tuple(initial_point)
    folded = initial_row_action or _identity(len(point))
    records = []
    terminal: LocationResult | None = None
    last_location: LocationResult | None = None
    for reflection_index in range(max_root_reflections + 1):
        location = _walk_until_boundary_or_interior(
            atlas, point, max_gallery=max_gallery
        )
        last_location = location
        records.append(
            {
                "fold_index": reflection_index,
                "result": location.kind,
                "representative": location.state.index,
                "auxiliary_gallery_length": len(location.gallery),
                "boundary_label": (
                    list(location.boundary_label)
                    if location.boundary_label is not None
                    else None
                ),
                "negative_pairing_with_fixed_ample": (
                    -_field_pair(tuple(AMPLE_FIELD_COEFFICIENTS), point)
                ).to_json(),
            }
        )
        if location.kind == "interior":
            terminal = location
            break
        if location.kind != "root_boundary" or location.boundary_functional is None:
            terminal = location
            break
        reflection = _reflection_in_functional(location.boundary_functional)
        point = _row_product(point, reflection)
        folded = _matrix_product(folded, reflection)
    if terminal is None:
        assert last_location is not None
        terminal = LocationResult(
            "reflection_limit",
            last_location.state,
            last_location.gallery,
        )
    return PointFoldResult(point, folded, terminal, tuple(records))


def fold_h11_seed(
    path: Path,
    *,
    seed_power: int = 1,
    max_root_reflections: int = 100,
    max_gallery: int = 2000,
) -> dict[str, object]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")
    if seed_power < 1:
        raise ValueError("seed_power must be positive")
    atlas = ShimadaDisc11Atlas(path)
    seed_action = _matrix_power(_seed_row_action(), seed_power)
    initial_point = _row_product(tuple(AMPLE_FIELD_COEFFICIENTS), seed_action)
    gram = tuple(tuple(entry for entry in row) for row in field_gram())
    start_square = golden_quadratic_value(gram, AMPLE_FIELD_COEFFICIENTS)
    fold = fold_point_to_eigen_nef(
        atlas,
        initial_point,
        initial_row_action=seed_action,
        max_root_reflections=max_root_reflections,
        max_gallery=max_gallery,
    )
    point = fold.point
    folded = fold.row_action
    terminal = fold.terminal
    records = fold.records

    trace_power, trace_value = _infinite_order_trace_power(folded)
    same_chamber = terminal is not None and terminal.kind == "interior"
    integral_matrix = all(entry.is_integral for row in folded for entry in row)
    square_preserved = golden_quadratic_value(gram, point) == start_square
    verified = (
        same_chamber
        and _isometry_verified(folded)
        and integral_matrix
        and square_preserved
        and trace_power is not None
    )
    return {
        "verified_infinite_order_same_chamber_obstruction": verified,
        "source_sha256": digest,
        "seed_power": seed_power,
        "seed_isometry_verified": _isometry_verified(_seed_row_action()),
        "root_reflections_used": sum(
            record["result"] == "root_boundary" for record in records
        ),
        "fold_records": records,
        "terminal_same_chamber": same_chamber,
        "terminal_representative": terminal.state.index if terminal else None,
        "terminal_auxiliary_gallery": (
            [
                {
                    "representative": representative,
                    "orbit": orbit,
                    "member": member,
                }
                for representative, orbit, member in terminal.gallery
            ]
            if terminal and terminal.kind == "interior"
            else []
        ),
        "folded_matrix_integral_over_Zphi": integral_matrix,
        "folded_matrix_isometry_verified": _isometry_verified(folded),
        "point_square_preserved": square_preserved,
        "infinite_order_trace_power": trace_power,
        "infinite_order_trace_value": (
            trace_value.to_json() if trace_value is not None else None
        ),
        "folded_row_action": [
            [entry.to_json() for entry in row] for row in folded
        ],
        "conclusion": (
            "H11_norm_two_reflection_group_has_infinite_index"
            if verified
            else "bounded_folding_search_inconclusive"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    parser.add_argument("--seed-power", type=int, default=1)
    parser.add_argument("--max-root-reflections", type=int, default=100)
    parser.add_argument("--max-gallery", type=int, default=2000)
    args = parser.parse_args()
    print(
        json.dumps(
            fold_h11_seed(
                args.disc11,
                seed_power=args.seed_power,
                max_root_reflections=args.max_root_reflections,
                max_gallery=args.max_gallery,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()

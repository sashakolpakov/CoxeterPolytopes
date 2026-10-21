#!/usr/bin/env python3
"""Certify rigidity of the archived OG10 golden Vinberg prefixes.

For a finite oriented norm-two wall set, the multiset of exact inner
products in each Gram row is invariant under every setwise isometry.  The
archived 171-, 321-, and 457-wall sets have pairwise distinct row
multisets.  Since their first eleven roots form a ``Z[phi]``-basis, their
setwise stabilizers in ``O(H_OG, Z[phi])`` are trivial.

The certificate also exhausts extensions of the height-zero diagram
automorphisms and records the natural infinite-order binary-block isometry
obtained from two wall reflections.  The latter does not stabilize any of
the archived prefixes.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence, TypeAlias

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.search.og10_golden_vinberg import og10_field_gram


ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = (
    ROOT / "data/candidates/og10_golden_vinberg_bound8.json",
    ROOT / "data/candidates/og10_golden_vinberg_through_2_plus_4phi.json",
    ROOT / "data/candidates/og10_golden_vinberg_through_3_plus_4phi.json",
)

ZPhi: TypeAlias = tuple[int, int]
QPhi: TypeAlias = tuple[Fraction, Fraction]
Vector: TypeAlias = tuple[ZPhi, ...]


def _add(left: tuple, right: tuple) -> tuple:
    return left[0] + right[0], left[1] + right[1]


def _neg(value: tuple) -> tuple:
    return -value[0], -value[1]


def _scale(value: tuple, scalar: int | Fraction) -> tuple:
    return scalar * value[0], scalar * value[1]


def _multiply(left: tuple, right: tuple) -> tuple:
    """Multiply coefficients using phi^2=phi+1."""

    a, b = left
    c, d = right
    return a * c + b * d, a * d + b * c + b * d


def _digest(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode(
        "ascii"
    )
    return hashlib.sha256(encoded).hexdigest()


def _integer(value: Golden) -> ZPhi:
    if not value.is_integral:
        raise ValueError(f"expected a Z[phi] value, got {value}")
    return int(value.a), int(value.b)


def _roots(payload: dict[str, object]) -> tuple[Vector, ...]:
    serialized = payload["simple_roots"]
    if not isinstance(serialized, list):
        raise TypeError("simple_roots must be a list")
    return tuple(
        tuple(_integer(Golden.from_json(entry)) for entry in root)
        for root in serialized
    )


def _gram() -> tuple[tuple[ZPhi, ...], ...]:
    return tuple(tuple(_integer(entry) for entry in row) for row in og10_field_gram())


def _covector(
    root: Sequence[ZPhi], gram: Sequence[Sequence[ZPhi]]
) -> tuple[ZPhi, ...]:
    return tuple(
        _sum(_multiply(gram[i][j], root[j]) for j in range(len(root)))
        for i in range(len(root))
    )


def _sum(values: Iterable[tuple]) -> tuple:
    result = (0, 0)
    for value in values:
        result = _add(result, value)
    return result


def _dot(left: Sequence[ZPhi], covector: Sequence[ZPhi]) -> ZPhi:
    return _sum(_multiply(left[i], covector[i]) for i in range(len(left)))


def _inner_matrix(
    roots: Sequence[Vector], gram: Sequence[Sequence[ZPhi]]
) -> tuple[tuple[ZPhi, ...], ...]:
    covectors = tuple(_covector(root, gram) for root in roots)
    return tuple(tuple(_dot(root, covector) for covector in covectors) for root in roots)


def _fingerprints(
    inner: Sequence[Sequence[ZPhi]], size: int
) -> tuple[tuple[tuple[ZPhi, int], ...], ...]:
    return tuple(
        tuple(sorted(Counter(row[:size]).items())) for row in inner[:size]
    )


def _fingerprint_payload(
    fingerprints: Sequence[Sequence[tuple[ZPhi, int]]],
) -> list[list[list[int]]]:
    return [
        [[value[0], value[1], multiplicity] for value, multiplicity in fingerprint]
        for fingerprint in fingerprints
    ]


def _height_zero_automorphisms(
    height_zero_gram: Sequence[Sequence[ZPhi]],
) -> tuple[tuple[int, ...], ...]:
    """Enumerate every permutation preserving the ten-node exact Gram matrix."""

    size = len(height_zero_gram)
    row_types = [Counter(row) for row in height_zero_gram]
    candidates = [
        tuple(j for j in range(size) if row_types[i] == row_types[j])
        for i in range(size)
    ]
    solutions: list[tuple[int, ...]] = []
    assignment: dict[int, int] = {}
    used: set[int] = set()

    def extend() -> None:
        if len(assignment) == size:
            solutions.append(tuple(assignment[i] for i in range(size)))
            return
        unassigned = [i for i in range(size) if i not in assignment]
        source = min(
            unassigned,
            key=lambda i: sum(target not in used for target in candidates[i]),
        )
        for target in candidates[source]:
            if target in used:
                continue
            if any(
                height_zero_gram[source][other]
                != height_zero_gram[target][assignment[other]]
                for other in assignment
            ):
                continue
            assignment[source] = target
            used.add(target)
            extend()
            used.remove(target)
            del assignment[source]

    extend()
    return tuple(sorted(solutions))


def _height_zero_extension(
    roots: Sequence[Vector], permutation: Sequence[int], p_sign: int
) -> tuple[tuple[QPhi, ...], ...]:
    """Return columns of the unique field isometry extending diagram data.

    The orthogonal basis is ``p=(2,-1,0,...,0)`` followed by the ten
    height-zero simple roots.  Since ``e0=(p-r0)/2`` and
    ``ej=-r_(j-1)``, the extension is explicit.
    """

    p: tuple[QPhi, ...] = (
        (Fraction(2), Fraction(0)),
        (Fraction(-1), Fraction(0)),
        *((Fraction(0), Fraction(0)) for _ in range(9)),
    )
    height_zero: tuple[tuple[QPhi, ...], ...] = tuple(
        tuple((Fraction(a), Fraction(b)) for a, b in root) for root in roots[:10]
    )
    image_e0 = tuple(
        _scale(
            _add(_scale(p[i], p_sign), _neg(height_zero[permutation[0]][i])),
            Fraction(1, 2),
        )
        for i in range(11)
    )
    other_columns = tuple(
        tuple(_neg(entry) for entry in height_zero[permutation[j]])
        for j in range(10)
    )
    return (image_e0, *other_columns)


def _columns_integral(columns: Sequence[Sequence[QPhi]]) -> bool:
    return all(
        a.denominator == 1 and b.denominator == 1
        for column in columns
        for a, b in column
    )


def _apply_columns(columns: Sequence[Sequence[ZPhi]], vector: Vector) -> Vector:
    return tuple(
        _sum(_multiply(vector[j], columns[j][i]) for j in range(len(vector)))
        for i in range(len(vector))
    )


def _reflection(
    vector: Vector, root: Vector, root_covector: Sequence[ZPhi]
) -> Vector:
    # Every root used here has norm two.
    coefficient = _dot(vector, root_covector)
    return tuple(
        _add(vector[i], _neg(_multiply(coefficient, root[i])))
        for i in range(len(vector))
    )


def _binary_translation(
    roots: Sequence[Vector], gram: Sequence[Sequence[ZPhi]]
) -> tuple[tuple[ZPhi, ...], ...]:
    """Columns of s_(root 0) o s_(root 10)."""

    first = roots[0]
    second = roots[10]
    first_covector = _covector(first, gram)
    second_covector = _covector(second, gram)
    columns: list[Vector] = []
    for j in range(11):
        basis = tuple((1, 0) if i == j else (0, 0) for i in range(11))
        image = _reflection(basis, second, second_covector)
        columns.append(_reflection(image, first, first_covector))
    return tuple(columns)


def _matrix_isometry(
    columns: Sequence[Sequence[tuple]], gram: Sequence[Sequence[ZPhi]]
) -> bool:
    transformed_covectors = tuple(_covector(column, gram) for column in columns)
    return all(
        _dot(columns[i], transformed_covectors[j]) == gram[i][j]
        for i in range(11)
        for j in range(11)
    )


def build() -> dict[str, object]:
    payloads = [json.loads(path.read_text()) for path in ARTIFACTS]
    root_sets = [_roots(payload) for payload in payloads]
    sizes = tuple(len(roots) for roots in root_sets)
    if sizes != (171, 321, 457):
        raise ValueError(f"unexpected archived prefix sizes {sizes}")
    roots = root_sets[-1]
    if any(root_sets[i] != roots[: sizes[i]] for i in range(2)):
        raise ValueError("the archived root sets are not nested ordered prefixes")
    if any(
        _digest(payload["simple_roots"]) != payload["roots_sha256"]
        for payload in payloads
    ):
        raise ValueError("an archived root digest does not match")

    gram = _gram()
    inner = _inner_matrix(roots, gram)
    if any(inner[i][i] != (2, 0) for i in range(len(roots))):
        raise ValueError("an archived wall is not a norm-two root")

    prefix_records: list[dict[str, object]] = []
    all_fingerprints_unique = True
    for path, size in zip(ARTIFACTS, sizes, strict=True):
        fingerprints = _fingerprints(inner, size)
        multiplicities = Counter(fingerprints)
        unique = len(multiplicities) == size
        all_fingerprints_unique &= unique
        prefix_records.append(
            {
                "artifact": str(path.relative_to(ROOT)),
                "wall_count": size,
                "distinct_exact_gram_row_multisets": len(multiplicities),
                "largest_fingerprint_multiplicity": max(multiplicities.values()),
                "ordered_fingerprint_sha256": _digest(
                    _fingerprint_payload(fingerprints)
                ),
                "oriented_wall_set_stabilizer": "identity" if unique else "unresolved",
            }
        )

    expected_first_ten = tuple(
        tuple((-1, 0) if j == i + 1 else (0, 0) for j in range(11))
        for i in range(10)
    )
    leading_coefficient = roots[10][0]
    leading_coefficient_norm = (
        leading_coefficient[0] ** 2
        + leading_coefficient[0] * leading_coefficient[1]
        - leading_coefficient[1] ** 2
    )
    first_eleven_form_basis = (
        roots[:10] == expected_first_ten
        and roots[10] == ((-1, 1), (1, 0), *((0, 0) for _ in range(9)))
        and leading_coefficient_norm == -1
    )

    p = ((2, 0), (-1, 0), *((0, 0) for _ in range(9)))
    p_covector = _covector(p, gram)
    p_is_orthogonal_control_line = (
        _dot(p, p_covector) == (-2, -8)
        and all(_dot(root, p_covector) == (0, 0) for root in roots[:10])
    )

    height_zero_gram = tuple(tuple(row[:10]) for row in inner[:10])
    diagram_automorphisms = _height_zero_automorphisms(height_zero_gram)
    identity = tuple(range(10))
    isolated_swap = (9, *range(1, 9), 0)
    extension_records = []
    for permutation in diagram_automorphisms:
        for p_sign in (1, -1):
            columns = _height_zero_extension(roots, permutation, p_sign)
            extension_records.append(
                {
                    "diagram_permutation": list(permutation),
                    "p_sign": p_sign,
                    "isometry_over_Q_phi": _matrix_isometry(columns, gram),
                    "integral_over_Z_phi": _columns_integral(columns),
                    "preserves_selected_p_ray": p_sign == 1,
                }
            )

    translation = _binary_translation(roots, gram)
    translation_images = tuple(_apply_columns(translation, root) for root in roots)
    root_lookup = set(roots)
    signed_lookup = root_lookup | {
        tuple(_neg(entry) for entry in root) for root in roots
    }
    retained_oriented = tuple(
        i for i, image in enumerate(translation_images) if image in root_lookup
    )
    retained_unoriented = tuple(
        i for i, image in enumerate(translation_images) if image in signed_lookup
    )
    binary_block = tuple(tuple(translation[j][i] for j in range(2)) for i in range(2))
    expected_binary_block = (((0, 1), (0, -1)), ((-1, -1), (0, 2)))
    # The displayed matrix is row-major; columns above are transposed here.
    binary_trace = _add(binary_block[0][0], binary_block[1][1])
    binary_determinant = _add(
        _multiply(binary_block[0][0], binary_block[1][1]),
        _neg(_multiply(binary_block[0][1], binary_block[1][0])),
    )
    phi = (1 + math.sqrt(5)) / 2
    spectral_radius = (3 * phi + math.sqrt(9 * phi * phi - 4)) / 2

    checks = {
        "archived_prefixes_are_nested_171_321_457": sizes == (171, 321, 457),
        "all_roots_have_field_norm_two": all(
            inner[i][i] == (2, 0) for i in range(len(roots))
        ),
        "every_prefix_has_unique_exact_gram_row_multisets": all_fingerprints_unique,
        "first_eleven_roots_form_a_Z_phi_basis": first_eleven_form_basis,
        "p_has_norm_minus_2_minus_8phi_and_height_zero_orthogonal": (
            p_is_orthogonal_control_line
        ),
        "height_zero_diagram_group_has_order_two": diagram_automorphisms
        == (identity, isolated_swap),
        "all_height_zero_diagram_extensions_are_field_isometries": all(
            record["isometry_over_Q_phi"] for record in extension_records
        ),
        "isolated_A1_swap_has_no_integral_extension": all(
            not record["integral_over_Z_phi"]
            for record in extension_records
            if tuple(record["diagram_permutation"]) == isolated_swap
        ),
        "selected_p_ray_height_zero_stabilizer_is_identity": [
            record
            for record in extension_records
            if record["integral_over_Z_phi"]
            and record["preserves_selected_p_ray"]
        ]
        == [
            {
                "diagram_permutation": list(identity),
                "p_sign": 1,
                "isometry_over_Q_phi": True,
                "integral_over_Z_phi": True,
                "preserves_selected_p_ray": True,
            }
        ],
        "binary_translation_is_integral_isometry": _matrix_isometry(
            translation, gram
        ),
        "binary_translation_has_expected_block": binary_block
        == expected_binary_block,
        "binary_translation_characteristic_polynomial_t2_minus_3phi_t_plus_1": (
            binary_trace == (0, 3) and binary_determinant == (1, 0)
        ),
        "binary_translation_is_infinite_order": Golden(-2, 3).sign() > 0,
        "binary_translation_retains_only_nine_archived_walls": (
            retained_oriented == tuple(range(1, 10))
            and retained_unoriented == retained_oriented
        ),
    }

    return {
        "schema_version": "0.1.0",
        "certificate_kind": "exact_bounded_OG10_golden_prefix_stabilizer",
        "candidate_id": "og10_golden_h10_carrier",
        "prefixes": prefix_records,
        "basis": {
            "indices": list(range(11)),
            "determinant_up_to_sign": "phi-1",
            "determinant_field_norm": leading_coefficient_norm,
        },
        "height_zero": {
            "diagram_automorphisms": [list(value) for value in diagram_automorphisms],
            "extension_checks": extension_records,
            "integral_time_reversing_involution_characteristic_polynomial": (
                "(t+1)*(t-1)^10"
            ),
        },
        "infinite_order_candidate": {
            "construction": "reflection(root_0) o reflection(root_10)",
            "binary_matrix": [
                [[value[0], value[1]] for value in row] for row in binary_block
            ],
            "characteristic_polynomial_over_Z_phi": "(t-1)^9*(t^2-3*phi*t+1)",
            "characteristic_polynomial_over_Z": (
                "(t-1)^18*(t^4-3*t^3-7*t^2-3*t+1)"
            ),
            "identity_embedding_spectral_radius_decimal": f"{spectral_radius:.15f}",
            "retained_oriented_wall_indices": list(retained_oriented),
            "retained_unoriented_wall_indices": list(retained_unoriented),
            "stabilizes_457_wall_prefix": False,
        },
        "checks": checks,
        "verified": all(checks.values()),
        "conclusion": (
            "the_setwise_stabilizer_in_O(H_OG,Z[phi])_of_each_archived_"
            "oriented_wall_prefix_is_trivial;_this_is_a_bounded_prefix_"
            "theorem_and_does_not_compute_the_stabilizer_of_the_unfinished_"
            "infinite_wall_chamber"
        ),
    }


def main() -> None:
    print(json.dumps(build(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

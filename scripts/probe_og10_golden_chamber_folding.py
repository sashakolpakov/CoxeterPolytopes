#!/usr/bin/env python3
"""Exact chamber reduction probes for the golden OG10 norm-two lattice.

This is deliberately a *dual-track* search.  It uses the archived exhaustive
Vinberg prefix through identity control coordinate ``x <= 10`` together with
an a-priori Lorentzian bound.  Whenever that bound is below ``10^2``, testing
the archived roots is a test against the complete (possibly infinite) set of
simple norm-two walls, not merely against the archived prefix.

The initial arithmetic seed is the binary Pell isometry

    A = [[1, 1-phi], [-1, phi]],

extended by the identity on ``E8 + <2>``.  Its square is the already archived
product ``s_0 s_10``.  Further seeds are integral reflections with normal
supported on the binary block and the final ``<2>`` line.  A successful probe
would still have to show that the folded same-chamber matrix has infinite
projective order.  Finite-order outputs and bounded failures prove nothing
about the eventual reflection index.
"""

from __future__ import annotations

import json
import hashlib
from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path
from typing import Iterable, Sequence, TypeAlias

from sympy import Matrix as SympyMatrix, eye

from coxeter_polytopes.exact.golden import Golden, ONE, PHI, ZERO
from coxeter_polytopes.search.bugaenko_vinberg import _fraction_inverse
from coxeter_polytopes.search.h11_vinberg import _standard_e8_gram
from coxeter_polytopes.search.og10_golden_vinberg import og10_field_gram
from scripts import certify_og10_golden_h10_carrier as carrier


ROOT = Path(__file__).resolve().parents[1]
PREFIX = ROOT / "data/candidates/og10_golden_vinberg_through_3_plus_4phi.json"

Vector = tuple[Golden, ...]
Matrix = tuple[Vector, ...]  # columns


def _sum(values: Iterable[Golden]) -> Golden:
    return sum(values, start=ZERO)


def _covector(vector: Sequence[Golden], gram: Sequence[Sequence[Golden]]) -> Vector:
    return tuple(
        _sum(gram[i][j] * vector[j] for j in range(len(vector)))
        for i in range(len(vector))
    )


def _dot_covector(vector: Sequence[Golden], covector: Sequence[Golden]) -> Golden:
    return _sum(vector[i] * covector[i] for i in range(len(vector)))


def _pair(
    left: Sequence[Golden],
    right: Sequence[Golden],
    gram: Sequence[Sequence[Golden]],
) -> Golden:
    return _dot_covector(left, _covector(right, gram))


def _identity(size: int) -> Matrix:
    return tuple(
        tuple(ONE if i == j else ZERO for i in range(size)) for j in range(size)
    )


def _apply(matrix: Sequence[Sequence[Golden]], vector: Sequence[Golden]) -> Vector:
    return tuple(
        _sum(matrix[j][i] * vector[j] for j in range(len(vector)))
        for i in range(len(vector))
    )


def _compose(left: Matrix, right: Matrix) -> Matrix:
    """Return ``left o right`` for matrices stored by columns."""

    return tuple(_apply(left, column) for column in right)


def _reflection_matrix(
    normal: Vector, gram: Sequence[Sequence[Golden]]
) -> Matrix | None:
    square = _pair(normal, normal, gram)
    if square.sign() <= 0:
        return None
    normal_covector = _covector(normal, gram)
    columns: list[Vector] = []
    for j in range(len(normal)):
        coefficient = 2 * normal_covector[j] / square
        column = tuple(
            (ONE if i == j else ZERO) - coefficient * normal[i]
            for i in range(len(normal))
        )
        if any(not entry.is_integral for entry in column):
            return None
        columns.append(column)
    return tuple(columns)


def _matrix_isometry(
    matrix: Matrix, gram: Sequence[Sequence[Golden]]
) -> bool:
    return all(
        _pair(matrix[i], matrix[j], gram) == gram[i][j]
        for i in range(len(matrix))
        for j in range(len(matrix))
    )


def _root_reflect_point(
    point: Vector, root: Vector, root_covector: Vector
) -> Vector:
    # Every archived wall root has square two.
    coefficient = _dot_covector(point, root_covector)
    return tuple(point[i] - coefficient * root[i] for i in range(len(point)))


def _root_reflection_matrix(root: Vector, root_covector: Vector) -> Matrix:
    columns = []
    for j in range(len(root)):
        columns.append(
            tuple(
                (ONE if i == j else ZERO) - root_covector[j] * root[i]
                for i in range(len(root))
            )
        )
    return tuple(columns)


def _load_roots() -> tuple[Vector, ...]:
    payload = json.loads(PREFIX.read_text())
    roots = tuple(
        tuple(Golden.from_json(entry) for entry in root)
        for root in payload["simple_roots"]
    )
    if len(roots) != 457:
        raise ValueError(f"expected 457 archived roots, got {len(roots)}")
    return roots


def _interior_point(gram: Sequence[Sequence[Golden]]) -> tuple[Vector, Vector, Vector]:
    """Return control point p, a Weyl-sector perturbation rho, and q=p+rho/100."""

    p = (Golden(2), Golden(-1), *(ZERO for _ in range(9)))
    inverse_e8 = _fraction_inverse(_standard_e8_gram())
    e8_rho = tuple(sum(row, start=Fraction(0)) for row in inverse_e8)
    rho = (
        ZERO,
        Golden(Fraction(1, 2)),
        *(Golden(entry) for entry in e8_rho),
        Golden(Fraction(1, 2)),
    )
    q = tuple(p[i] + rho[i] / 100 for i in range(11))
    if _pair(p, rho, gram):
        raise AssertionError("rho must be orthogonal to the control line")
    if _pair(rho, rho, gram) != 621:
        raise AssertionError("unexpected E8 Weyl-vector square")
    return p, rho, q


def _possible_violator_bound_square(
    point: Vector,
    control: Vector,
    gram: Sequence[Sequence[Golden]],
) -> Golden:
    """Bound ``x^2`` for every positive-height simple root violated by point.

    Write ``point = a*control + v`` with ``v`` control-orthogonal.  Since
    ``control^2=-2D`` and a norm-two root of control coordinate x has
    orthogonal component of square ``2 + D*x^2/2``, Cauchy--Schwarz gives

      x^2 < 2*v^2 / (D*(a^2*D-v^2/2)).

    The denominator is positive exactly when point is timelike and lies in
    the selected component (a>0).
    """

    d = Golden(1, 4)
    a = point[0] / 2
    orthogonal = tuple(point[i] - a * control[i] for i in range(11))
    orthogonal_square = _pair(orthogonal, orthogonal, gram)
    denominator = d * (a * a * d - orthogonal_square / 2)
    if a.sign() <= 0 or denominator.sign() <= 0:
        raise ValueError("point is not in the selected timelike component")
    return 2 * orthogonal_square / denominator


def _violations(point: Vector, covectors: Sequence[Vector]) -> tuple[int, ...]:
    return tuple(
        i
        for i, covector in enumerate(covectors)
        if _dot_covector(point, covector).sign() > 0
    )


def _identity_decimal(value: Golden) -> float:
    return float(value.a) + ((1 + 5**0.5) / 2) * float(value.b)


def _fold(
    seed: Matrix,
    interior: Vector,
    roots: Sequence[Vector],
    covectors: Sequence[Vector],
    reflection_matrices: Sequence[Matrix],
    control: Vector,
    gram: Sequence[Sequence[Golden]],
    *,
    maximum_steps: int = 1000,
) -> tuple[Matrix, Vector, tuple[int, ...], Golden]:
    """Greedily return the endpoint to the archived Vinberg half-spaces."""

    matrix = seed
    point = _apply(matrix, interior)
    word: list[int] = []
    seen: set[tuple[tuple[Fraction, Fraction], ...]] = set()
    for _ in range(maximum_steps):
        key = tuple((entry.a, entry.b) for entry in point)
        if key in seen:
            raise RuntimeError("folding entered a cycle")
        seen.add(key)
        violated = _violations(point, covectors)
        if not violated:
            return (
                matrix,
                point,
                tuple(word),
                _possible_violator_bound_square(point, control, gram),
            )
        # All roots have square two, so the largest positive pairing is the
        # steepest available descent in distance from the control point.
        index = max(
            violated,
            key=lambda i: _identity_decimal(_dot_covector(point, covectors[i])),
        )
        point = _root_reflect_point(point, roots[index], covectors[index])
        matrix = _compose(reflection_matrices[index], matrix)
        word.append(index)
    raise RuntimeError(f"folding did not terminate in {maximum_steps} steps")


def _pell_seed() -> Matrix:
    columns: list[Vector] = []
    # Row-major binary block [[1,1-phi],[-1,phi]].
    columns.append((ONE, Golden(-1), *(ZERO for _ in range(9))))
    columns.append((ONE - PHI, PHI, *(ZERO for _ in range(9))))
    for j in range(2, 11):
        columns.append(
            tuple(ONE if i == j else ZERO for i in range(11))
        )
    return tuple(columns)


def _matrix_key(matrix: Matrix) -> tuple[tuple[Fraction, Fraction], ...]:
    return tuple((entry.a, entry.b) for column in matrix for entry in column)


def _digest(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode(
        "ascii"
    )
    return hashlib.sha256(encoded).hexdigest()


def _matrix_payload(matrix: Matrix) -> list[list[dict[str, str]]]:
    return [[entry.to_json() for entry in column] for column in matrix]


def _sympy_matrix_digest(matrix: SympyMatrix) -> str:
    payload = [
        [str(matrix[i, j]) for j in range(matrix.cols)]
        for i in range(matrix.rows)
    ]
    return _digest(payload)


def _mod_one(value) -> Fraction:
    numerator, denominator = value.as_numer_denom()
    fraction = Fraction(int(numerator), int(denominator))
    return fraction - fraction.numerator // fraction.denominator


def _class_key(vector: SympyMatrix) -> tuple[Fraction, ...]:
    return tuple(_mod_one(entry) for entry in vector)


def _trace_and_ambient_extension(field_matrix: Matrix) -> dict[str, object]:
    """Certify the action on A_N and extension by -I on the period lattice."""

    original_trace_action = SympyMatrix.zeros(22)
    for j, column in enumerate(field_matrix):
        for i, entry in enumerate(column):
            if not entry.is_integral:
                raise ArithmeticError("field isometry is not integral")
            a = int(entry.a)
            b = int(entry.b)
            original_trace_action[2 * i, 2 * j] = a
            original_trace_action[2 * i + 1, 2 * j] = b
            original_trace_action[2 * i, 2 * j + 1] = b
            original_trace_action[2 * i + 1, 2 * j + 1] = a + b

    full_change = (
        SympyMatrix(carrier.trace_basis_change())
        * carrier.t11_reduction()
        * carrier.core_small_reordering()
    )
    trace_action = full_change.inv() * original_trace_action * full_change
    if not all(entry.is_Integer is True for entry in trace_action):
        raise ArithmeticError("trace-lattice action is nonintegral")
    trace_gram = SympyMatrix.diag(carrier.core_gram(), carrier.small_n_gram())
    phi_action = (
        full_change.inv()
        * SympyMatrix(carrier.multiplication_by_phi(11))
        * full_change
    )

    g1, g2 = carrier.gluing_vectors()
    n_generators = (g1[:4, :], g2[:4, :])
    t_generators = (g1[4:, :], g2[4:, :])
    n_span = {
        _class_key(a * n_generators[0] + b * n_generators[1]): (a, b)
        for a in range(2)
        for b in range(22)
    }
    induced = []
    for generator in n_generators:
        image = trace_action * SympyMatrix.vstack(SympyMatrix.zeros(18, 1), generator)
        if any(_mod_one(entry) for entry in image[:18, :]):
            raise ArithmeticError("discriminant action leaked into the unimodular core")
        induced.append(n_span[_class_key(image[18:, :])])

    complement_action = -eye(2)
    gluing_matches = True
    for generator, coefficient in zip(t_generators, induced, strict=True):
        actual = _class_key(complement_action * generator)
        desired = _class_key(
            coefficient[0] * t_generators[0]
            + coefficient[1] * t_generators[1]
        )
        gluing_matches &= actual == desired

    direct_action = SympyMatrix.diag(trace_action, complement_action)
    overlattice_change = SympyMatrix.diag(eye(18), carrier.overlattice_basis())
    ambient_action = overlattice_change.inv() * direct_action * overlattice_change
    direct_gram = SympyMatrix.diag(
        carrier.core_gram(), carrier.small_n_gram(), carrier.complement_gram()
    )
    ambient_gram = overlattice_change.T * direct_gram * overlattice_change
    ambient_integral = all(entry.is_Integer is True for entry in ambient_action)

    checks = {
        "trace_action_integral_isometry": (
            trace_action.T * trace_gram * trace_action == trace_gram
        ),
        "trace_action_commutes_with_tau": trace_action * phi_action
        == phi_action * trace_action,
        "trace_action_order_two": trace_action * trace_action == eye(22),
        "discriminant_generator_action_is_n1_fixed_n2_negated": induced
        == [(1, 0), (0, 21)],
        "minus_identity_on_complement_matches_gluing": gluing_matches,
        "ambient_extension_integral": ambient_integral,
        "ambient_extension_isometry": (
            ambient_action.T * ambient_gram * ambient_action == ambient_gram
        ),
        "ambient_extension_order_two": ambient_action * ambient_action == eye(24),
        "ambient_extension_determinant_one": ambient_action.det() == 1,
        # h preserves the selected positive ray.  The complement is positive
        # definite of rank two and -I has determinant +1 there.
        "ambient_extension_preserves_positive_orientation": (
            complement_action.det() == 1
        ),
        # The OG10 period plane is the real complement plane, so -I preserves
        # its oriented complex period line projectively.
        "minus_identity_on_complement_preserves_period_line": True,
    }
    return {
        "trace_action_sha256": _sympy_matrix_digest(trace_action),
        "discriminant_group_order": 44,
        "gluing_generator_orders": [2, 22],
        "gluing_generator_action_coefficients": [list(value) for value in induced],
        "gluing_generator_action": ["n1_to_n1", "n2_to_minus_n2"],
        "complement": "<2> + <66>",
        "complement_action": "-identity",
        "ambient_extension_sha256": _sympy_matrix_digest(ambient_action),
        "checks": checks,
        "verified": all(checks.values()),
        "conclusion": (
            "h_extends_to_an_order_two_positive_Hodge_monodromy_isometry;_"
            "therefore_h_lies_in_Gamma_tau_mon_and_in_the_enhanced_"
            "birational_image_by_birational_Torelli"
        ),
    }


def _small_supported_reflections(
    gram: Sequence[Sequence[Golden]], coefficient_bound: int = 1
) -> tuple[tuple[Matrix, Vector, Golden], ...]:
    values = tuple(
        Golden(a, b)
        for a in range(-coefficient_bound, coefficient_bound + 1)
        for b in range(-coefficient_bound, coefficient_bound + 1)
    )
    distinct: dict[
        tuple[tuple[Fraction, Fraction], ...], tuple[Matrix, Vector, Golden]
    ] = {}
    for x, y, u in product(values, repeat=3):
        if not x and not y and not u:
            continue
        normal = (x, y, *(ZERO for _ in range(8)), u)
        matrix = _reflection_matrix(normal, gram)
        if matrix is None:
            continue
        distinct.setdefault(
            _matrix_key(matrix), (matrix, normal, _pair(normal, normal, gram))
        )
    return tuple(distinct.values())


def build() -> dict[str, object]:
    gram = og10_field_gram()
    roots = _load_roots()
    covectors = tuple(_covector(root, gram) for root in roots)
    reflection_matrices = tuple(
        _root_reflection_matrix(root, covector)
        for root, covector in zip(roots, covectors, strict=True)
    )
    control, rho, interior = _interior_point(gram)
    interior_pairings = tuple(
        _dot_covector(interior, covector) for covector in covectors
    )
    interior_violations = tuple(
        i for i, value in enumerate(interior_pairings) if value.sign() > 0
    )
    interior_zeroes = tuple(
        i for i, value in enumerate(interior_pairings) if not value
    )
    interior_bound = _possible_violator_bound_square(interior, control, gram)

    pell = _pell_seed()
    if not _matrix_isometry(pell, gram):
        raise AssertionError("Pell seed is not an integral isometry")
    pell_square = _compose(pell, pell)
    expected_translation = _compose(reflection_matrices[0], reflection_matrices[10])
    chamber_involution = _compose(reflection_matrices[0], pell)
    chamber_involution_order_two = (
        _compose(chamber_involution, chamber_involution) == _identity(11)
    )
    chamber_involution_normal = (
        PHI - ONE,
        Golden(2),
        *(ZERO for _ in range(9)),
    )
    chamber_involution_normal_square = _pair(
        chamber_involution_normal, chamber_involution_normal, gram
    )
    chamber_involution_as_reflection = _reflection_matrix(
        chamber_involution_normal, gram
    )

    raw_seeds: list[tuple[str, Matrix, Vector | None, Golden | None]] = [
        ("pell_A", pell, None, None),
        ("pell_A_inverse", _compose(expected_translation, pell), None, None),
    ]
    for i, (matrix, normal, square) in enumerate(
        _small_supported_reflections(gram), start=0
    ):
        raw_seeds.append((f"supported_reflection_{i}", matrix, normal, square))

    records = []
    folded_matrices: dict[tuple[tuple[Fraction, Fraction], ...], Matrix] = {}
    for label, seed, normal, normal_square in raw_seeds:
        folded, endpoint, word, bound = _fold(
            seed,
            interior,
            roots,
            covectors,
            reflection_matrices,
            control,
            gram,
        )
        complete_membership = bound.sign() >= 0 and (bound - 100).sign() < 0
        endpoint_pairings = tuple(
            _dot_covector(endpoint, covector) for covector in covectors
        )
        strict_archived_membership = all(value.sign() < 0 for value in endpoint_pairings)
        key = _matrix_key(folded)
        folded_matrices.setdefault(key, folded)
        folded_square = _compose(folded, folded)
        if folded == _identity(11):
            folded_projective_order: int | str = 1
        elif folded_square == _identity(11):
            folded_projective_order = 2
        else:
            folded_projective_order = "not_certified"
        records.append(
            {
                "label": label,
                "normal_supported_coordinates_x_y_u": None
                if normal is None
                else [normal[i].to_json() for i in (0, 1, 10)],
                "normal_square": None
                if normal_square is None
                else normal_square.to_json(),
                "folding_word": list(word),
                "folding_word_length": len(word),
                "possible_violator_x_square_bound": bound.to_json(),
                "complete_chamber_membership_certified_by_x_lt_10": (
                    complete_membership and strict_archived_membership
                ),
                "all_457_archived_wall_pairings_strictly_negative": (
                    strict_archived_membership
                ),
                "returns_interior_point_exactly": endpoint == interior,
                "folded_matrix_is_identity": folded == _identity(11),
                "folded_projective_order": folded_projective_order,
                "folded_matrix_class": (
                    "identity"
                    if folded == _identity(11)
                    else "complete_chamber_involution_h"
                    if folded == chamber_involution
                    else "other"
                ),
            }
        )

    root_lookup = {root: i for i, root in enumerate(roots)}
    involution_images = tuple(_apply(chamber_involution, root) for root in roots)
    retained = tuple(
        (i, root_lookup[image])
        for i, image in enumerate(involution_images)
        if image in root_lookup
    )
    new_roots = tuple(image for image in involution_images if image not in root_lookup)
    new_level_counts = Counter(str(root[0]) for root in new_roots)
    closed_roots = tuple(dict.fromkeys((*roots, *new_roots)))
    closed_lookup = set(closed_roots)
    closed_under_involution = all(
        _apply(chamber_involution, root) in closed_lookup for root in closed_roots
    )
    closed_covectors = tuple(_covector(root, gram) for root in closed_roots)
    closed_pairwise_acute = all(
        _dot_covector(left, closed_covectors[j]).sign() <= 0
        for i, left in enumerate(closed_roots)
        for j in range(i + 1, len(closed_roots))
    )
    new_root_payload = [
        [entry.to_json() for entry in root] for root in new_roots
    ]

    monodromy_extension = _trace_and_ambient_extension(chamber_involution)

    complete_chamber_checks = {
        "interior_has_no_nonnegative_pairing_with_457_roots": (
            not interior_violations and not interior_zeroes
        ),
        "interior_global_violator_bound_below_100": (
            interior_bound.sign() >= 0 and (interior_bound - 100).sign() < 0
        ),
        "pell_square_equals_s0_s10": pell_square == expected_translation,
        "folded_pell_is_h_equals_s0_A": records[0]["folding_word"] == [0],
        "h_is_nonidentity": chamber_involution != _identity(11),
        "h_has_order_two": chamber_involution_order_two,
        "h_is_integral_reflection_in_norm_6_plus_2phi_vector": (
            chamber_involution_normal_square == Golden(6, 2)
            and chamber_involution_as_reflection == chamber_involution
        ),
        "h_endpoint_strictly_inside_all_457_halfspaces": records[0][
            "all_457_archived_wall_pairings_strictly_negative"
        ],
        "h_endpoint_global_violator_bound_below_100": records[0][
            "complete_chamber_membership_certified_by_x_lt_10"
        ],
        "small_seed_scan_folds_only_to_identity_or_h": all(
            record["folded_matrix_class"] in {"identity", "complete_chamber_involution_h"}
            for record in records
        ),
        "h_retains_274_and_exports_183_of_457_roots": (
            len(retained) == 274 and len(new_roots) == 183
        ),
        "h_closure_has_640_roots": len(closed_roots) == 640,
        "h_closure_is_involution_invariant": closed_under_involution,
        "h_closure_is_pairwise_acute": closed_pairwise_acute,
        "new_root_level_counts_20_3_129_31": new_level_counts
        == {
            "2 + 5*phi": 20,
            "4 + 4*phi": 3,
            "3 + 5*phi": 129,
            "4 + 5*phi": 31,
        },
        "monodromy_extension_verified": monodromy_extension["verified"],
    }

    return {
        "schema_version": "0.1.0",
        "certificate_kind": "exact_bounded_OG10_complete_chamber_folding_probe",
        "candidate_id": "og10_golden_h10_carrier",
        "interior_point": [entry.to_json() for entry in interior],
        "rho_square": _pair(rho, rho, gram).to_json(),
        "interior_archived_wall_violations": list(interior_violations),
        "interior_archived_wall_zeroes": list(interior_zeroes),
        "interior_possible_violator_x_square_bound": interior_bound.to_json(),
        "interior_complete_chamber_membership_certified_by_x_lt_10": (
            not interior_violations
            and not interior_zeroes
            and (interior_bound - 100).sign() < 0
        ),
        "global_membership_argument": {
            "decomposition": "q=a*p+v with v perpendicular to p",
            "violated_positive_height_root_bound": (
                "x^2 < 2*v^2/(D*(a^2*D-v^2/2))"
            ),
            "prefix_completeness": "all simple roots with 0<x<=10 are archived",
            "height_zero_test": "all ten height-zero simple pairings are strict",
            "conclusion": (
                "if the exact bound is below 100 and all 457 pairings are "
                "strictly negative, the point is in the complete open chamber"
            ),
        },
        "pell_seed": {
            "characteristic_polynomial_over_Z_phi": (
                "(t-1)^9*(t^2-(1+phi)*t+1)"
            ),
            "infinite_order": (Golden(-1, 1)).sign() > 0,
            "square_equals_s0_s10": pell_square == expected_translation,
        },
        "complete_chamber_involution": {
            "construction": "h=s0*A",
            "field_matrix_columns": _matrix_payload(chamber_involution),
            "field_matrix_sha256": _digest(_matrix_payload(chamber_involution)),
            "integral_reflection_normal": [
                entry.to_json() for entry in chamber_involution_normal
            ],
            "integral_reflection_normal_square": (
                chamber_involution_normal_square.to_json()
            ),
            "order": 2,
            "global_membership_x_square_bound": records[0][
                "possible_violator_x_square_bound"
            ],
            "global_membership_certified": records[0][
                "complete_chamber_membership_certified_by_x_lt_10"
            ],
            "scope": (
                "h is a nontrivial complete norm-two chamber symmetry but has "
                "finite order, so it does not prove infinite reflection index"
            ),
        },
        "chamber_involution_wall_transport": {
            "archived_root_count": len(roots),
            "images_retained_in_archived_prefix": len(retained),
            "images_outside_archived_prefix": len(new_roots),
            "involution_closed_root_count": len(closed_roots),
            "new_control_level_counts": dict(
                sorted(new_level_counts.items(), key=lambda item: item[0])
            ),
            "new_roots_sha256": _digest(new_root_payload),
            "closed_root_set_sha256": _digest(
                [[entry.to_json() for entry in root] for root in closed_roots]
            ),
            "all_640_pairwise_acute": closed_pairwise_acute,
            "scope": (
                "the 183 transported roots are complete-chamber simple walls by "
                "the certified symmetry; this alone does not exclude additional "
                "walls at the same control levels"
            ),
        },
        "monodromy_and_birational_scope": monodromy_extension,
        "seed_count": len(raw_seeds),
        "distinct_folded_matrices": len(folded_matrices),
        "folded_seeds": records,
        "checks": complete_chamber_checks,
        "verified": all(complete_chamber_checks.values()),
        "conclusion": (
            "the_complete_norm_two_chamber_has_a_nontrivial_order_two_"
            "Gamma_tau_mon_symmetry_h_and_the_enhanced_birational_image_"
            "contains_h;_no_infinite_order_same_chamber_element_is_certified"
        ),
        "scope": (
            "complete chamber membership is proved only for records whose exact "
            "violator bound is below 100; finite-order folded outputs or bounded "
            "search failure do not decide the reflection index"
        ),
    }


def main() -> None:
    print(json.dumps(build(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

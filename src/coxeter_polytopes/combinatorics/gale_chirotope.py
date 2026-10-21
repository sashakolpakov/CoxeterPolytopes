r"""Uniform rank-four Gale oriented-matroid constraints for H9/14 nerves.

A simplicial 9-polytope on 14 vertices can be perturbed inside its realization
space so that its affine oriented matroid is uniform of rank 10.  Its Gale
dual is then a uniform rank-four oriented matroid.  This gives a finite,
structurally independent necessary polytopality test:

* there are ``C(14,4)=1001`` chirotope signs;
* the uniform three-term Grassmann--Pluecker axiom gives 45,045 constraints;
* a five-set is a facet complement precisely when its circuit has one sign.
* total cyclicity and primal vertex extremality are imposed separately.

The module can either audit a concrete sign certificate without a solver or
ask Z3 for one.  A passing oriented matroid is still not a realizability
certificate: nonrealizable oriented matroids also satisfy these axioms.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import comb
from time import monotonic
from typing import Iterable, Mapping, Sequence


GROUND_SET_SIZE = 14
GALE_RANK = 4
EXPECTED_CHIROTOPE_SIGN_COUNT = comb(GROUND_SET_SIZE, GALE_RANK)
EXPECTED_GRASSMANN_PLUECKER_COUNT = (
    comb(GROUND_SET_SIZE, GALE_RANK - 2)
    * comb(GROUND_SET_SIZE - (GALE_RANK - 2), 4)
)


def rank_four_bases() -> tuple[tuple[int, int, int, int], ...]:
    """Return the canonical lexicographic basis order."""

    return tuple(combinations(range(GROUND_SET_SIZE), GALE_RANK))


RANK_FOUR_BASES = rank_four_bases()
RANK_FOUR_BASIS_INDEX = {
    basis: index for index, basis in enumerate(RANK_FOUR_BASES)
}


def _mask(vertices: Iterable[int]) -> int:
    result = 0
    for vertex in vertices:
        if not 0 <= vertex < GROUND_SET_SIZE:
            raise ValueError("Gale element lies outside the 14-element ground set")
        result |= 1 << vertex
    return result


def _permutation_is_odd(values: Sequence[int]) -> bool:
    return bool(
        sum(
            values[left] > values[right]
            for left in range(len(values))
            for right in range(left + 1, len(values))
        )
        % 2
    )


def _normalize_signs(signs: Sequence[int]) -> tuple[int, ...]:
    normalized = tuple(int(sign) for sign in signs)
    if len(normalized) != EXPECTED_CHIROTOPE_SIGN_COUNT:
        raise ValueError("a rank-four chirotope needs exactly 1001 basis signs")
    if any(sign not in (-1, 1) for sign in normalized):
        raise ValueError("a uniform chirotope sign must be +1 or -1")
    return normalized


def ordered_chirotope_sign(
    signs: Sequence[int], ordered_basis: Sequence[int]
) -> int:
    """Evaluate the alternating chirotope on one ordered four-tuple."""

    normalized = _normalize_signs(signs)
    return _ordered_chirotope_sign_normalized(normalized, ordered_basis)


def _ordered_chirotope_sign_normalized(
    normalized: tuple[int, ...], ordered_basis: Sequence[int]
) -> int:
    """Fast internal ordered evaluation after one global validation."""

    ordered = tuple(ordered_basis)
    if len(ordered) != GALE_RANK or len(set(ordered)) != GALE_RANK:
        raise ValueError("a uniform rank-four basis has four distinct elements")
    try:
        value = normalized[RANK_FOUR_BASIS_INDEX[tuple(sorted(ordered))]]
    except KeyError as exc:
        raise ValueError("basis element lies outside the ground set") from exc
    return -value if _permutation_is_odd(ordered) else value


def first_grassmann_pluecker_failure(
    signs: Sequence[int],
) -> tuple[tuple[int, int], tuple[int, int, int, int], tuple[int, int, int]] | None:
    """Return the first failed uniform three-term chirotope relation."""

    normalized = _normalize_signs(signs)
    ground = tuple(range(GROUND_SET_SIZE))
    for fixed in combinations(ground, GALE_RANK - 2):
        remaining = tuple(element for element in ground if element not in fixed)
        for first, second, third, fourth in combinations(remaining, 4):
            products = (
                _ordered_chirotope_sign_normalized(
                    normalized, (*fixed, first, second)
                )
                * _ordered_chirotope_sign_normalized(
                    normalized, (*fixed, third, fourth)
                ),
                -_ordered_chirotope_sign_normalized(
                    normalized, (*fixed, first, third)
                )
                * _ordered_chirotope_sign_normalized(
                    normalized, (*fixed, second, fourth)
                ),
                _ordered_chirotope_sign_normalized(
                    normalized, (*fixed, first, fourth)
                )
                * _ordered_chirotope_sign_normalized(
                    normalized, (*fixed, second, third)
                ),
            )
            if products[0] == products[1] == products[2]:
                return fixed, (first, second, third, fourth), products
    return None


def positive_circuit_masks(signs: Sequence[int]) -> tuple[int, ...]:
    """Return all positive five-element circuit supports.

    For ``B=(b_0<...<b_4)``, circuit coefficients have signs
    ``(-1)^i chi(B without b_i)``.  The circuit is positive, up to its global
    sign, exactly when all five coefficient signs agree.
    """

    normalized = _normalize_signs(signs)
    result = []
    for support in combinations(range(GROUND_SET_SIZE), GALE_RANK + 1):
        coefficient_signs = tuple(
            (-1 if index % 2 else 1)
            * normalized[RANK_FOUR_BASIS_INDEX[support[:index] + support[index + 1 :]]]
            for index in range(GALE_RANK + 1)
        )
        if len(set(coefficient_signs)) == 1:
            result.append(_mask(support))
    return tuple(sorted(result))


@dataclass(frozen=True)
class GaleChirotopeAudit:
    """Independent audit of a fixed Gale sign certificate."""

    sign_count: int
    first_grassmann_pluecker_failure: (
        tuple[tuple[int, int], tuple[int, int, int, int], tuple[int, int, int]]
        | None
    )
    positive_circuit_masks: tuple[int, ...]
    expected_facet_complement_masks: tuple[int, ...]
    positive_circuits_match_facets: bool
    every_element_in_positive_circuit: bool
    every_element_avoided_by_positive_circuit: bool

    @property
    def passes(self) -> bool:
        return (
            self.first_grassmann_pluecker_failure is None
            and self.positive_circuits_match_facets
            and self.every_element_in_positive_circuit
            and self.every_element_avoided_by_positive_circuit
        )


def audit_uniform_rank_four_gale(
    facet_complement_masks: Iterable[int], signs: Sequence[int]
) -> GaleChirotopeAudit:
    """Check all chirotope axioms and the facet/positive-circuit biconditional."""

    normalized = _normalize_signs(signs)
    expected = tuple(sorted(set(facet_complement_masks)))
    full_mask = (1 << GROUND_SET_SIZE) - 1
    if not expected:
        raise ValueError("a Gale facet certificate cannot be empty")
    if any(
        mask <= 0
        or mask & ~full_mask
        or mask.bit_count() != GALE_RANK + 1
        for mask in expected
    ):
        raise ValueError("facet complements must be five-set masks")
    positive = positive_circuit_masks(normalized)
    covered = 0
    avoided = 0
    for support in positive:
        covered |= support
        avoided |= full_mask ^ support
    return GaleChirotopeAudit(
        sign_count=len(normalized),
        first_grassmann_pluecker_failure=first_grassmann_pluecker_failure(
            normalized
        ),
        positive_circuit_masks=positive,
        expected_facet_complement_masks=expected,
        positive_circuits_match_facets=positive == expected,
        every_element_in_positive_circuit=covered == full_mask,
        every_element_avoided_by_positive_circuit=avoided == full_mask,
    )


@dataclass(frozen=True)
class GaleChirotopeSolveResult:
    """One exact oriented-matroid extension disposition."""

    status: str
    signs: tuple[int, ...]
    backend: str
    backend_version: str
    elapsed_seconds: float
    assertion_count: int
    reason_unknown: str | None


def _z3_ordered_chirotope(variables, ordered_basis, z3):
    ordered = tuple(ordered_basis)
    base = variables[RANK_FOUR_BASIS_INDEX[tuple(sorted(ordered))]]
    return z3.Not(base) if _permutation_is_odd(ordered) else base


def add_uniform_rank_four_gale_constraints(
    solver,
    facet_literals: Mapping[int, object],
):
    """Couple rank-four chirotope signs to five-set facet literals.

    ``facet_literals`` may omit five-sets; omitted sets are constrained not to
    be positive circuits.  The returned tuple contains the 1001 sign variables.
    This helper imports Z3 lazily so the independent audit has no solver
    dependency.
    """

    try:
        import z3  # type: ignore[import-not-found]
    except ImportError as exc:  # pragma: no cover - optional search backend
        raise RuntimeError("Gale chirotope search requires z3-solver") from exc

    variables = tuple(
        z3.Bool(f"gale_chi_{'_'.join(map(str, basis))}")
        for basis in RANK_FOUR_BASES
    )
    # Global sign reversal leaves all circuits and GP products unchanged.
    solver.add(variables[0])

    ground = tuple(range(GROUND_SET_SIZE))
    for fixed in combinations(ground, GALE_RANK - 2):
        remaining = tuple(element for element in ground if element not in fixed)
        for first, second, third, fourth in combinations(remaining, 4):
            first_product = (
                _z3_ordered_chirotope(
                    variables, (*fixed, first, second), z3
                )
                == _z3_ordered_chirotope(
                    variables, (*fixed, third, fourth), z3
                )
            )
            second_product = z3.Not(
                _z3_ordered_chirotope(
                    variables, (*fixed, first, third), z3
                )
                == _z3_ordered_chirotope(
                    variables, (*fixed, second, fourth), z3
                )
            )
            third_product = (
                _z3_ordered_chirotope(
                    variables, (*fixed, first, fourth), z3
                )
                == _z3_ordered_chirotope(
                    variables, (*fixed, second, third), z3
                )
            )
            solver.add(
                z3.Or(first_product, second_product, third_product),
                z3.Or(
                    z3.Not(first_product),
                    z3.Not(second_product),
                    z3.Not(third_product),
                ),
            )

    positive_by_element: list[list[object]] = [
        [] for _ in range(GROUND_SET_SIZE)
    ]
    positive_avoiding_element: list[list[object]] = [
        [] for _ in range(GROUND_SET_SIZE)
    ]
    for support in combinations(ground, GALE_RANK + 1):
        coefficient_signs = []
        for index in range(GALE_RANK + 1):
            sign = variables[
                RANK_FOUR_BASIS_INDEX[support[:index] + support[index + 1 :]]
            ]
            coefficient_signs.append(sign if index % 2 == 0 else z3.Not(sign))
        positive = z3.And(
            *[
                coefficient_signs[index] == coefficient_signs[0]
                for index in range(1, GALE_RANK + 1)
            ]
        )
        literal = facet_literals.get(_mask(support), z3.BoolVal(False))
        solver.add(positive == literal)
        for element in support:
            positive_by_element[element].append(positive)
        for element in set(ground) - set(support):
            positive_avoiding_element[element].append(positive)

    # The dual of an acyclic affine oriented matroid is totally cyclic:
    # every element lies in a positive circuit.  Geometrically, every vertex
    # of a bounded polytope is omitted by at least one facet.  This condition
    # is automatic for a genuine closed pseudomanifold facet selection.
    for element_circuits in positive_by_element:
        solver.add(z3.Or(*element_circuits))

    # Every labelled primal point must actually be an extreme vertex.  By the
    # Gale face criterion, {e} is a face exactly when a positive circuit
    # avoids e, equivalently when e occurs in at least one selected facet.
    # This too is automatic for the intended closed pseudomanifold input, but
    # is explicit here because the standalone helper accepts an arbitrary
    # facet-complement mapping.
    for avoiding_circuits in positive_avoiding_element:
        solver.add(z3.Or(*avoiding_circuits))
    return variables


def solve_uniform_rank_four_gale(
    facet_complement_masks: Iterable[int], *, timeout_ms: int | None = None
) -> GaleChirotopeSolveResult:
    """Find a uniform rank-four oriented-matroid extension of one nerve."""

    try:
        import z3  # type: ignore[import-not-found]
    except ImportError as exc:  # pragma: no cover - optional search backend
        raise RuntimeError("Gale chirotope search requires z3-solver") from exc

    selected = tuple(sorted(set(facet_complement_masks)))
    full_mask = (1 << GROUND_SET_SIZE) - 1
    if not selected or any(
        mask <= 0
        or mask & ~full_mask
        or mask.bit_count() != GALE_RANK + 1
        for mask in selected
    ):
        raise ValueError("selected facet complements must be nonempty five-sets")

    solver = z3.Solver()
    if timeout_ms is not None:
        if timeout_ms <= 0:
            raise ValueError("timeout must be positive")
        solver.set(timeout=timeout_ms)
    start = monotonic()
    signs = add_uniform_rank_four_gale_constraints(
        solver, {mask: z3.BoolVal(True) for mask in selected}
    )
    assertion_count = len(solver.assertions())
    disposition = solver.check()
    elapsed = monotonic() - start
    if disposition == z3.sat:
        model = solver.model()
        certificate = tuple(
            1 if z3.is_true(model.eval(sign)) else -1 for sign in signs
        )
        status = "sat"
        reason_unknown = None
    elif disposition == z3.unsat:
        certificate = ()
        status = "unsat"
        reason_unknown = None
    else:
        certificate = ()
        status = "unknown"
        reason_unknown = solver.reason_unknown()
    return GaleChirotopeSolveResult(
        status=status,
        signs=certificate,
        backend="z3",
        backend_version=z3.get_version_string(),
        elapsed_seconds=elapsed,
        assertion_count=assertion_count,
        reason_unknown=reason_unknown,
    )

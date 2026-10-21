r"""Exact certificates for polyhedral cones over :math:`\mathbf Q(\phi)`.

The routines in this module use finite exact phase-I and Lorentzian active-set
algorithms; they do not decide feasibility with floating-point tolerances.  A
finite half-space reduction is certified by

* exact positive-proportionality classes;
* one strict face point for every asserted facet; and
* one positive linear dependence for every conically redundant class; and
* an exact nonpositive Lorentzian maximum for any causally inactive class.

This is the certificate pattern used by the restricted Shimada chambers in
the discriminant-11 K3 calculation.  It is kept independent of that atlas so
that the final verifier need only know vectors in ``Q(phi)``.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Hashable, Mapping, Sequence, TypeAlias

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.exact.matrix import inertia, rank


Label: TypeAlias = Hashable
GoldenVector: TypeAlias = tuple[Golden, ...]
GoldenMatrix: TypeAlias = tuple[tuple[Golden, ...], ...]
DependenceTerms: TypeAlias = tuple[tuple[Golden, Label], ...]


@dataclass(frozen=True, slots=True)
class PhaseOneCertificate:
    """Exact phase-I certificate for a nonnegative linear system.

    ``objective`` is the minimum sum of the artificial variables in the
    standard-form system built by :func:`_phase_one_nonnegative`.  It is zero
    exactly when the requested system is feasible.  The dual vector verifies
    the same optimum by strong linear-programming duality over the ordered
    field ``Q(phi)``.
    """

    feasible: bool
    objective: Golden
    primal: tuple[Golden, ...] | None
    dual: tuple[Golden, ...]
    primal_objective_matches: bool
    dual_objective_matches: bool
    dual_feasible: bool

    @property
    def passes(self) -> bool:
        return (
            self.primal_objective_matches
            and self.dual_objective_matches
            and self.dual_feasible
            and (self.feasible == (self.objective == 0))
            and (self.primal is not None) == self.feasible
        )


@dataclass(frozen=True, slots=True)
class LorentzianFaceDecision:
    """Exact decision for one strict face of a Lorentzian half-space cone.

    A facet witness satisfies the target equality, every other inequality
    strictly, and has positive quadratic value in the time component selected
    by ``interior``.  If the first exact linear-programming witness is not
    timelike, ``maximum_point`` and its KKT data certify the global maximum of
    the quadratic form on the normalized closed face.
    """

    is_facet: bool
    reason: str
    witness: GoldenVector | None
    strict_linear_certificate: PhaseOneCertificate
    maximum_point: GoldenVector | None = None
    maximum_quadratic_value: Golden | None = None
    active_labels: tuple[Label, ...] = ()
    active_multipliers: tuple[Golden, ...] = ()
    equality_multipliers: tuple[Golden, Golden] | None = None


@dataclass(frozen=True, slots=True)
class ClosedLorentzianConeDecision:
    """Exact future-timelike decision for a closed homogeneous polyhedron.

    The input inequalities are ``row*x >= 0``.  The homogeneous cone is
    normalized by

    ``B(future_reference, x) = B(future_reference, future_reference)``.

    Every future-timelike ray meets this affine slice exactly once.  The
    phase-I certificate proves whether the *closed* slice is nonempty, even
    when it has empty relative interior in the ambient space.  A positive
    ``witness`` proves retention.  Otherwise ``maximum_point`` and the KKT
    fields certify the global Lorentzian maximum on the slice.
    """

    has_future_timelike_point: bool
    reason: str
    witness: GoldenVector | None
    closed_linear_certificate: PhaseOneCertificate
    maximum_point: GoldenVector | None = None
    maximum_quadratic_value: Golden | None = None
    active_labels: tuple[Label, ...] = ()
    active_multipliers: tuple[Golden, ...] = ()
    equality_multiplier: Golden | None = None
    retained_linear_labels: tuple[Label, ...] = ()
    linear_redundancies: tuple[tuple[Label, DependenceTerms], ...] = ()



def _vector(values: Sequence[Golden | int]) -> GoldenVector:
    return tuple(Golden.coerce(value) for value in values)


def _matrix(values: Sequence[Sequence[Golden | int]]) -> GoldenMatrix:
    rows = tuple(_vector(row) for row in values)
    if not rows or any(len(row) != len(rows) for row in rows):
        raise ValueError("the quadratic Gram matrix must be nonempty and square")
    return rows


def golden_dot(left: Sequence[Golden | int], right: Sequence[Golden | int]) -> Golden:
    """Return the exact coordinate pairing of two equally sized vectors."""

    if len(left) != len(right):
        raise ValueError("vectors have different dimensions")
    return sum(
        (Golden.coerce(a) * Golden.coerce(b) for a, b in zip(left, right)),
        Golden(0),
    )


def golden_quadratic_value(
    gram: Sequence[Sequence[Golden | int]],
    point: Sequence[Golden | int],
) -> Golden:
    """Evaluate ``point * gram * point^T`` exactly."""

    matrix = _matrix(gram)
    vector = _vector(point)
    if len(matrix) != len(vector):
        raise ValueError("the point and quadratic Gram matrix have different dimensions")
    return sum(
        (
            vector[i] * matrix[i][j] * vector[j]
            for i in range(len(vector))
            for j in range(len(vector))
        ),
        Golden(0),
    )


def inverse_golden_matrix(
    matrix: Sequence[Sequence[Golden | int]],
) -> GoldenMatrix:
    """Invert a square matrix over ``Q(phi)`` by exact row reduction."""

    source = _matrix(matrix)
    size = len(source)
    work = [
        list(source[row])
        + [Golden(int(row == column)) for column in range(size)]
        for row in range(size)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if pivot is None:
            raise ZeroDivisionError("the matrix is singular")
        work[column], work[pivot] = work[pivot], work[column]
        pivot_value = work[column][column]
        work[column] = [entry / pivot_value for entry in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                work[row][index] - factor * work[column][index]
                for index in range(2 * size)
            ]
    return tuple(tuple(row[size:]) for row in work)


def _mat_vec(
    matrix: Sequence[Sequence[Golden]], vector: Sequence[Golden]
) -> tuple[Golden, ...]:
    return tuple(
        sum((entry * value for entry, value in zip(row, vector)), Golden(0))
        for row in matrix
    )


def _phase_one_standard(
    matrix: Sequence[Sequence[Golden]],
    right_hand_side: Sequence[Golden],
    objective: Sequence[Golden],
    initial_basis: Sequence[int],
    *,
    original_variable_count: int,
) -> PhaseOneCertificate:
    """Minimize a phase-I objective by exact Bland revised simplex.

    The supplied basis must be feasible.  Bland's entering and leaving rules
    make the finite pivot algorithm non-cycling, including in degenerate
    cases.  This is deliberately implemented over :class:`Golden`, rather
    than delegated to a floating solver or SymPy's rational-only simplex.
    """

    rows = tuple(tuple(entry for entry in row) for row in matrix)
    rhs = tuple(right_hand_side)
    costs = tuple(objective)
    row_count = len(rows)
    if row_count == 0:
        primal = tuple(Golden(0) for _ in range(original_variable_count))
        return PhaseOneCertificate(
            True,
            Golden(0),
            primal,
            (),
            True,
            True,
            True,
        )
    width = len(rows[0])
    if (
        any(len(row) != width for row in rows)
        or len(rhs) != row_count
        or len(costs) != width
        or len(initial_basis) != row_count
        or original_variable_count > width
    ):
        raise ValueError("inconsistent phase-I standard-form dimensions")
    if any(value.sign() < 0 for value in rhs):
        raise ValueError("the initial phase-I basis must have nonnegative RHS")

    basis = list(initial_basis)
    if len(set(basis)) != row_count:
        raise ValueError("phase-I basis columns must be distinct")

    # Canonical primal-simplex tableau.  The initial auxiliary columns form an
    # identity basis, so only the objective row needs elimination.
    tableau = [list(rows[row]) + [rhs[row]] for row in range(row_count)]
    objective_row = list(costs) + [Golden(0)]
    for row, column in enumerate(basis):
        coefficient = costs[column]
        if not coefficient:
            continue
        objective_row = [
            objective_row[index] - coefficient * tableau[row][index]
            for index in range(width + 1)
        ]
    tableau.append(objective_row)

    # A Bland simplex has finitely many bases.  The guard is an assertion of
    # the theorem, not a search cutoff: reaching a repeated basis indicates an
    # implementation error rather than returning an undecided result.
    seen_bases: set[tuple[int, ...]] = set()
    while True:
        basis_key = tuple(basis)
        if basis_key in seen_bases:
            raise ArithmeticError("Bland phase-I simplex repeated a basis")
        seen_bases.add(basis_key)
        basis_set = set(basis)
        entering = next(
            (
                column
                for column in range(width)
                if column not in basis_set and tableau[-1][column].sign() < 0
            ),
            None,
        )
        if entering is None:
            break
        candidates = [
            row
            for row in range(row_count)
            if tableau[row][entering].sign() > 0
        ]
        if not candidates:
            raise ArithmeticError("the nonnegative phase-I objective became unbounded")
        minimum_ratio = min(
            tableau[row][-1] / tableau[row][entering] for row in candidates
        )
        tied = [
            row
            for row in candidates
            if tableau[row][-1] / tableau[row][entering] == minimum_ratio
        ]
        leaving_row = min(tied, key=lambda row: basis[row])
        pivot = tableau[leaving_row][entering]
        tableau[leaving_row] = [entry / pivot for entry in tableau[leaving_row]]
        pivot_row = tableau[leaving_row]
        for row in range(row_count + 1):
            if row == leaving_row or not tableau[row][entering]:
                continue
            coefficient = tableau[row][entering]
            tableau[row] = [
                tableau[row][column] - coefficient * pivot_row[column]
                for column in range(width + 1)
            ]
        basis[leaving_row] = entering

    basic_values = tuple(tableau[row][-1] for row in range(row_count))
    if any(value.sign() < 0 for value in basic_values):
        raise ArithmeticError("phase-I pivot lost primal feasibility")
    full_primal = [Golden(0) for _ in range(width)]
    for row, column in enumerate(basis):
        full_primal[column] = basic_values[row]

    # Reconstruct the exact dual once at the terminal basis.  This gives an
    # independently auditable Farkas/strong-duality certificate without paying
    # for a matrix inverse at every simplex pivot.
    basis_matrix = tuple(
        tuple(rows[row][column] for column in basis)
        for row in range(row_count)
    )
    inverse_basis = inverse_golden_matrix(basis_matrix)
    basis_costs = tuple(costs[column] for column in basis)
    dual = tuple(
        sum(
            (basis_costs[row] * inverse_basis[row][column] for row in range(row_count)),
            Golden(0),
        )
        for column in range(row_count)
    )
    reduced_costs = tuple(
        costs[column]
        - sum(
            (dual[row] * rows[row][column] for row in range(row_count)),
            Golden(0),
        )
        for column in range(width)
    )
    optimum = sum(
        (costs[column] * full_primal[column] for column in range(width)),
        Golden(0),
    )
    primal_standard_feasible = (
        all(value.sign() >= 0 for value in full_primal)
        and all(
            sum(
                (rows[row][column] * full_primal[column] for column in range(width)),
                Golden(0),
            )
            == rhs[row]
            for row in range(row_count)
        )
    )
    dual_objective = sum(
        (dual[row] * rhs[row] for row in range(row_count)), Golden(0)
    )
    dual_feasible = all(value.sign() >= 0 for value in reduced_costs)
    feasible = optimum == 0
    primal = tuple(full_primal[:original_variable_count]) if feasible else None
    return PhaseOneCertificate(
        feasible,
        optimum,
        primal,
        dual,
        primal_standard_feasible
        and sum(
            (costs[column] * full_primal[column] for column in range(width)),
            Golden(0),
        ) == optimum,
        dual_objective == optimum,
        dual_feasible,
    )


def _phase_one_nonnegative(
    variable_count: int,
    *,
    equalities: Sequence[tuple[Sequence[Golden], Golden]] = (),
    inequalities_le: Sequence[tuple[Sequence[Golden], Golden]] = (),
) -> PhaseOneCertificate:
    """Decide an exact linear system in nonnegative variables.

    Equalities have the form ``row*x = rhs`` and inequalities the form
    ``row*x <= rhs``.  Slack, surplus, and artificial variables are generated
    canonically; the returned primal contains only the requested variables.
    """

    if variable_count < 0:
        raise ValueError("variable_count must be nonnegative")
    constraints: list[tuple[str, tuple[Golden, ...], Golden]] = []
    for kind, source in (("eq", equalities), ("le", inequalities_le)):
        for raw_row, raw_rhs in source:
            row = tuple(Golden.coerce(entry) for entry in raw_row)
            rhs = Golden.coerce(raw_rhs)
            if len(row) != variable_count:
                raise ValueError("a phase-I constraint has the wrong width")
            constraints.append((kind, row, rhs))

    # First assign the canonical auxiliary columns.
    next_column = variable_count
    records: list[tuple[str, tuple[Golden, ...], Golden, int | None, int | None]] = []
    artificial_columns: set[int] = set()
    for kind, row, rhs in constraints:
        if rhs.sign() < 0:
            row = tuple(-entry for entry in row)
            rhs = -rhs
            if kind == "le":
                # The flipped inequality is >=: surplus -1 plus artificial +1.
                surplus = next_column
                artificial = next_column + 1
                next_column += 2
                artificial_columns.add(artificial)
                records.append(("ge", row, rhs, surplus, artificial))
            else:
                artificial = next_column
                next_column += 1
                artificial_columns.add(artificial)
                records.append(("eq", row, rhs, None, artificial))
        elif kind == "le":
            slack = next_column
            next_column += 1
            records.append(("le", row, rhs, slack, None))
        else:
            artificial = next_column
            next_column += 1
            artificial_columns.add(artificial)
            records.append(("eq", row, rhs, None, artificial))

    width = next_column
    standard_rows: list[list[Golden]] = []
    standard_rhs: list[Golden] = []
    basis: list[int] = []
    for kind, row, rhs, auxiliary, artificial in records:
        standard = list(row) + [Golden(0) for _ in range(width - variable_count)]
        if kind == "le":
            assert auxiliary is not None and artificial is None
            standard[auxiliary] = Golden(1)
            basis.append(auxiliary)
        elif kind == "ge":
            assert auxiliary is not None and artificial is not None
            standard[auxiliary] = Golden(-1)
            standard[artificial] = Golden(1)
            basis.append(artificial)
        else:
            assert artificial is not None
            standard[artificial] = Golden(1)
            basis.append(artificial)
        standard_rows.append(standard)
        standard_rhs.append(rhs)
    objective = [
        Golden(1) if column in artificial_columns else Golden(0)
        for column in range(width)
    ]
    return _phase_one_standard(
        standard_rows,
        standard_rhs,
        objective,
        basis,
        original_variable_count=variable_count,
    )


def project_to_halfspace_wall(
    interior: Sequence[Golden | int],
    functional: Sequence[Golden | int],
    gram: Sequence[Sequence[Golden | int]],
) -> GoldenVector:
    """Orthogonally project an interior point to a functional's hyperplane.

    The functional is a covector.  Its metric-dual normal is therefore
    ``gram^{-1} functional``.  An isotropic dual normal is rejected because
    the projection formula would be undefined.
    """

    point = _vector(interior)
    wall = _vector(functional)
    if len(point) != len(wall):
        raise ValueError("the point and functional have different dimensions")
    inverse = inverse_golden_matrix(gram)
    if len(inverse) != len(point):
        raise ValueError("the vectors and quadratic Gram matrix have different dimensions")
    normal = tuple(
        sum(
            (inverse[row][column] * wall[column] for column in range(len(wall))),
            Golden(0),
        )
        for row in range(len(wall))
    )
    normal_square = golden_dot(wall, normal)
    if not normal_square:
        raise ValueError("the functional has isotropic metric-dual normal")
    scale = golden_dot(wall, point) / normal_square
    result = tuple(point[index] - scale * normal[index] for index in range(len(point)))
    if golden_dot(wall, result):
        raise ArithmeticError("exact wall projection failed")
    return result


@dataclass(frozen=True, slots=True)
class PositiveProportionalClass:
    """One oriented ray of nonzero half-space functionals."""

    members: tuple[Label, ...]
    ratios_to_first: tuple[Golden, ...]

    @property
    def passes(self) -> bool:
        return bool(self.members) and all(ratio.sign() > 0 for ratio in self.ratios_to_first)


def positive_proportional_ratio(
    base: Sequence[Golden | int],
    candidate: Sequence[Golden | int],
) -> Golden | None:
    """Return the positive ratio ``candidate / base``, or ``None``."""

    left = _vector(base)
    right = _vector(candidate)
    if len(left) != len(right):
        raise ValueError("functionals have different dimensions")
    pivot = next((index for index, entry in enumerate(left) if entry), None)
    if pivot is None:
        raise ValueError("a half-space functional cannot be zero")
    if not any(right):
        raise ValueError("a half-space functional cannot be zero")
    ratio = right[pivot] / left[pivot]
    if ratio.sign() <= 0:
        return None
    if right != tuple(ratio * entry for entry in left):
        return None
    return ratio


def positive_proportional_classes(
    rows: Mapping[Label, Sequence[Golden | int]],
) -> tuple[PositiveProportionalClass, ...]:
    """Partition oriented functionals into exact positive rays.

    Mapping insertion order is retained, which makes certificate generation
    deterministic for an atlas with a fixed wall order.
    """

    if not rows:
        raise ValueError("at least one half-space functional is required")
    vectors = {label: _vector(row) for label, row in rows.items()}
    dimensions = {len(row) for row in vectors.values()}
    if len(dimensions) != 1:
        raise ValueError("half-space functionals have different dimensions")
    if 0 in dimensions or any(not any(row) for row in vectors.values()):
        raise ValueError("a half-space functional cannot be zero")

    mutable: list[list[Label]] = []
    for label, row in vectors.items():
        for wall_class in mutable:
            if positive_proportional_ratio(vectors[wall_class[0]], row) is not None:
                wall_class.append(label)
                break
        else:
            mutable.append([label])

    result = []
    for wall_class in mutable:
        first = vectors[wall_class[0]]
        ratios = tuple(
            positive_proportional_ratio(first, vectors[label]) for label in wall_class
        )
        if any(ratio is None for ratio in ratios):
            raise ArithmeticError("positive-proportionality partition is inconsistent")
        result.append(
            PositiveProportionalClass(
                tuple(wall_class),
                tuple(ratio for ratio in ratios if ratio is not None),
            )
        )
    return tuple(result)


def _unique_solution_for_support(
    target: GoldenVector,
    generators: Sequence[GoldenVector],
) -> tuple[Golden, ...] | None:
    """Solve ``sum coefficient_i * generator_i = target`` for a fixed support.

    A dependent support is skipped.  This loses no positive-span witness:
    a conic representation always has a linearly independent sub-support.
    """

    variable_count = len(generators)
    if variable_count == 0:
        return None
    dimension = len(target)
    if any(len(row) != dimension for row in generators):
        raise ValueError("target and generators have different dimensions")
    work = [
        [generators[column][coordinate] for column in range(variable_count)]
        + [target[coordinate]]
        for coordinate in range(dimension)
    ]

    pivot_row = 0
    pivot_rows: dict[int, int] = {}
    for column in range(variable_count):
        pivot = next(
            (row for row in range(pivot_row, dimension) if work[row][column]),
            None,
        )
        if pivot is None:
            return None
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        work[pivot_row] = [entry / pivot_value for entry in work[pivot_row]]
        for row in range(dimension):
            if row == pivot_row or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                work[row][index] - factor * work[pivot_row][index]
                for index in range(variable_count + 1)
            ]
        pivot_rows[column] = pivot_row
        pivot_row += 1

    if any(
        not any(row[:variable_count]) and row[variable_count]
        for row in work
    ):
        return None
    solution = tuple(work[pivot_rows[column]][variable_count] for column in range(variable_count))
    reconstructed = tuple(
        sum(
            (solution[index] * generators[index][coordinate] for index in range(variable_count)),
            Golden(0),
        )
        for coordinate in range(dimension)
    )
    return solution if reconstructed == target else None


def find_positive_span(
    target: Sequence[Golden | int],
    generators: Mapping[Label, Sequence[Golden | int]],
    *,
    max_support: int | None = None,
) -> DependenceTerms | None:
    """Find an exact strictly positive conic expression for ``target``.

    This is a small deterministic certificate generator, not a numerical LP.
    It enumerates linearly independent supports, as justified by the conic
    Caratheodory theorem, and verifies the returned identity exactly.
    """

    goal = _vector(target)
    labels = tuple(generators)
    vectors = {label: _vector(generators[label]) for label in labels}
    if any(len(vector) != len(goal) for vector in vectors.values()):
        raise ValueError("target and generators have different dimensions")
    support_limit = min(len(goal), len(labels))
    if max_support is not None:
        if max_support < 1:
            raise ValueError("max_support must be positive")
        support_limit = min(support_limit, max_support)

    for size in range(1, support_limit + 1):
        for support in combinations(labels, size):
            solution = _unique_solution_for_support(
                goal, tuple(vectors[label] for label in support)
            )
            if solution is not None and all(value.sign() > 0 for value in solution):
                return tuple(zip(solution, support))
    return None


def positive_span_on_support(
    target: Sequence[Golden | int],
    generators: Mapping[Label, Sequence[Golden | int]],
) -> DependenceTerms | None:
    """Reconstruct a positive dependence on one proposed full support.

    Numerical cone solvers usually return a linearly independent basic
    support.  This helper solves that support once over ``Q(phi)`` and checks
    positivity, avoiding the exponential subset enumeration performed by
    :func:`find_positive_span` when the support is already known.
    """

    goal = _vector(target)
    labels = tuple(generators)
    vectors = tuple(_vector(generators[label]) for label in labels)
    solution = _unique_solution_for_support(goal, vectors)
    if solution is None or not all(value.sign() > 0 for value in solution):
        return None
    return tuple(zip(solution, labels))


def find_nonnegative_span_phase_one(
    target: Sequence[Golden | int],
    generators: Mapping[Label, Sequence[Golden | int]],
) -> tuple[DependenceTerms | None, PhaseOneCertificate]:
    """Decide cone membership by an exact finite phase-I simplex.

    Unlike :func:`find_positive_span`, this routine does not enumerate
    supports.  Zero basic coefficients are discarded, so every term in a
    returned dependence is strictly positive and the identity can be checked
    directly by :func:`audit_golden_halfspace_certificate`.
    """

    goal = _vector(target)
    labels = tuple(generators)
    vectors = tuple(_vector(generators[label]) for label in labels)
    if any(len(vector) != len(goal) for vector in vectors):
        raise ValueError("target and generators have different dimensions")
    equalities = tuple(
        (
            tuple(vector[coordinate] for vector in vectors),
            goal[coordinate],
        )
        for coordinate in range(len(goal))
    )
    certificate = _phase_one_nonnegative(
        len(labels),
        equalities=equalities,
    )
    if not certificate.passes:
        raise ArithmeticError("exact cone-membership phase-I audit failed")
    if not certificate.feasible:
        return None, certificate
    assert certificate.primal is not None
    terms = tuple(
        (coefficient, label)
        for coefficient, label in zip(certificate.primal, labels)
        if coefficient.sign() > 0
    )
    if not terms:
        return None, certificate
    reconstructed = tuple(
        sum(
            (coefficient * _vector(generators[label])[coordinate] for coefficient, label in terms),
            Golden(0),
        )
        for coordinate in range(len(goal))
    )
    if reconstructed != goal:
        raise ArithmeticError("phase-I cone-membership reconstruction failed")
    return terms, certificate


def _closed_free_slice_feasibility(
    rows: Mapping[Label, GoldenVector],
    time_functional: GoldenVector,
    time_value: Golden,
) -> tuple[GoldenVector | None, PhaseOneCertificate]:
    """Solve ``rows*x >= 0`` and ``time*x = time_value`` exactly.

    Coordinates are free, so each one is split as ``x^+ - x^-``.  Unlike the
    strict-face helper below, zero row pairings are deliberately retained:
    carrier cells can meet the eigenspace only in a lower-dimensional face.
    """

    dimension = len(time_functional)
    if any(len(row) != dimension for row in rows.values()):
        raise ValueError("closed slice rows have inconsistent dimensions")

    def split(row: GoldenVector) -> tuple[Golden, ...]:
        return (*row, *(-entry for entry in row))

    certificate = _phase_one_nonnegative(
        2 * dimension,
        equalities=((split(time_functional), time_value),),
        inequalities_le=tuple(
            (
                tuple(-entry for entry in split(row)),
                Golden(0),
            )
            for row in rows.values()
        ),
    )
    if not certificate.passes:
        raise ArithmeticError("closed-slice phase-I audit failed")
    if not certificate.feasible:
        return None, certificate
    assert certificate.primal is not None
    point = tuple(
        certificate.primal[index]
        - certificate.primal[index + dimension]
        for index in range(dimension)
    )
    if golden_dot(time_functional, point) != time_value:
        raise ArithmeticError("closed-slice phase I missed its normalization")
    if any(golden_dot(row, point).sign() < 0 for row in rows.values()):
        raise ArithmeticError("closed-slice phase I violated an inequality")
    return point, certificate


def _sequential_dual_cone_reduction(
    rows: Mapping[Label, GoldenVector],
) -> tuple[
    dict[Label, GoldenVector],
    tuple[tuple[Label, DependenceTerms], ...],
]:
    """Preserve a row cone while deleting exact conic redundancies.

    A closed slice need not admit a strict point, so its dual row cone need
    not be pointed.  Deletions are consequently performed sequentially rather
    than deleting every non-extreme input ray simultaneously.  Each deletion
    preserves the generated cone and carries an exact phase-I dependence.
    """

    nonzero_rows = {label: row for label, row in rows.items() if any(row)}
    zero_redundancies: list[tuple[Label, DependenceTerms]] = [
        (label, ()) for label, row in rows.items() if not any(row)
    ]
    if not nonzero_rows:
        return {}, tuple(zero_redundancies)

    classes = positive_proportional_classes(nonzero_rows)
    representatives = {
        item.members[0]: nonzero_rows[item.members[0]] for item in classes
    }
    redundancies = list(zero_redundancies)
    for item in classes:
        representative = item.members[0]
        for member, ratio in zip(item.members[1:], item.ratios_to_first[1:]):
            redundancies.append(
                (member, ((ratio, representative),))
            )

    retained = dict(representatives)
    for target in tuple(representatives):
        generators = {
            label: row for label, row in retained.items() if label != target
        }
        if not generators:
            continue
        dependence, certificate = find_nonnegative_span_phase_one(
            retained[target], generators
        )
        if not certificate.passes:
            raise ArithmeticError("closed-cone redundancy phase I failed")
        if dependence is None:
            continue
        redundancies.append((target, dependence))
        del retained[target]
    return retained, tuple(redundancies)


def decide_lorentzian_closed_cone(
    rows: Mapping[Label, Sequence[Golden | int]],
    *,
    quadratic_gram: Sequence[Sequence[Golden | int]],
    future_reference: Sequence[Golden | int],
) -> ClosedLorentzianConeDecision:
    """Decide whether a closed half-space cone contains a future timelike ray.

    This is the carrier-cell feasibility test.  It is finite and exact over
    ``Q(phi)`` and remains complete when the feasible section is contained in
    several walls.  After exact phase I, a non-timelike basic point triggers
    maximization of the Lorentzian form on the normalized closed polyhedron.

    On the affine time hyperplane the form is strictly concave and coercive.
    A maximum therefore exists whenever the slice is nonempty.  The normal
    cone of a polyhedron is generated by its active inequalities without a
    Slater hypothesis.  Conic Caratheodory then gives an independent active
    support of size at most ``dimension - 1``, which is enumerated below.
    """

    vectors = {label: _vector(row) for label, row in rows.items()}
    gram = _matrix(quadratic_gram)
    reference = _vector(future_reference)
    dimension = len(reference)
    if len(gram) != dimension or any(
        len(row) != dimension for row in vectors.values()
    ):
        raise ValueError("closed Lorentzian cone data have inconsistent dimensions")
    gram_inertia = inertia(gram)
    if (gram_inertia.positive, gram_inertia.negative, gram_inertia.zero) != (
        1,
        dimension - 1,
        0,
    ):
        raise ValueError("quadratic_gram must have Lorentzian inertia (1,n-1)")
    reference_square = golden_quadratic_value(gram, reference)
    if reference_square.sign() <= 0:
        raise ValueError("future_reference must be timelike")
    time_functional = tuple(
        sum(
            (
                reference[row] * gram[row][column]
                for row in range(dimension)
            ),
            Golden(0),
        )
        for column in range(dimension)
    )

    linear_point, phase_one = _closed_free_slice_feasibility(
        vectors, time_functional, reference_square
    )
    if linear_point is None:
        return ClosedLorentzianConeDecision(
            False,
            "normalized_closed_slice_infeasible",
            None,
            phase_one,
        )
    if golden_quadratic_value(gram, linear_point).sign() > 0:
        return ClosedLorentzianConeDecision(
            True,
            "phase_one_witness_timelike",
            linear_point,
            phase_one,
        )

    reduced_rows, redundancies = _sequential_dual_cone_reduction(vectors)
    labels = tuple(reduced_rows)
    maximum_point: GoldenVector | None = None
    maximum_labels: tuple[Label, ...] = ()
    maximum_multipliers: tuple[Golden, ...] = ()
    maximum_equality_multiplier: Golden | None = None

    for support_size in range(dimension):
        for support in combinations(labels, support_size):
            active_rows = tuple(reduced_rows[label] for label in support)
            constraint_rows = (*active_rows, time_functional)
            if rank(constraint_rows) != len(constraint_rows):
                continue
            multiplier_count = support_size + 1
            system_size = dimension + multiplier_count
            system = [
                [Golden(0) for _ in range(system_size)]
                for _ in range(system_size)
            ]
            rhs = [Golden(0) for _ in range(system_size)]
            for row in range(dimension):
                for column in range(dimension):
                    system[row][column] = Golden(2) * gram[row][column]
                for offset, functional in enumerate(constraint_rows):
                    system[row][dimension + offset] = functional[row]
            for offset, functional in enumerate(constraint_rows):
                equation_row = dimension + offset
                for column in range(dimension):
                    system[equation_row][column] = functional[column]
                if offset == support_size:
                    rhs[equation_row] = reference_square

            solution = _solve_square_linear_system(system, rhs)
            if solution is None:
                continue
            candidate = tuple(solution[:dimension])
            multipliers = tuple(
                solution[dimension : dimension + support_size]
            )
            if any(value.sign() < 0 for value in multipliers):
                continue
            if golden_dot(time_functional, candidate) != reference_square:
                continue
            if any(
                golden_dot(row, candidate).sign() < 0
                for row in reduced_rows.values()
            ):
                continue
            equality_multiplier = solution[-1]
            stationarity = tuple(
                Golden(2)
                * sum(
                    (
                        gram[coordinate][column] * candidate[column]
                        for column in range(dimension)
                    ),
                    Golden(0),
                )
                + sum(
                    (
                        multipliers[offset]
                        * active_rows[offset][coordinate]
                        for offset in range(support_size)
                    ),
                    Golden(0),
                )
                + equality_multiplier * time_functional[coordinate]
                for coordinate in range(dimension)
            )
            if any(stationarity):
                continue
            maximum_point = candidate
            maximum_labels = support
            maximum_multipliers = multipliers
            maximum_equality_multiplier = equality_multiplier
            break
        if maximum_point is not None:
            break
    if maximum_point is None or maximum_equality_multiplier is None:
        raise ArithmeticError("closed Lorentzian KKT enumeration found no maximum")

    maximum_square = golden_quadratic_value(gram, maximum_point)
    common = {
        "maximum_point": maximum_point,
        "maximum_quadratic_value": maximum_square,
        "active_labels": maximum_labels,
        "active_multipliers": maximum_multipliers,
        "equality_multiplier": maximum_equality_multiplier,
        "retained_linear_labels": labels,
        "linear_redundancies": redundancies,
    }
    if maximum_square.sign() <= 0:
        return ClosedLorentzianConeDecision(
            False,
            "normalized_closed_slice_has_nonpositive_global_quadratic_maximum",
            None,
            phase_one,
            **common,
        )
    if any(golden_dot(row, maximum_point).sign() < 0 for row in vectors.values()):
        raise ArithmeticError("a reduced closed-cone maximum failed the full rows")
    return ClosedLorentzianConeDecision(
        True,
        "positive_normalized_closed_quadratic_maximum",
        maximum_point,
        phase_one,
        **common,
    )


def _strict_free_face_feasibility(
    target: GoldenVector,
    others: Mapping[Label, GoldenVector],
    time_functional: GoldenVector,
) -> tuple[GoldenVector | None, PhaseOneCertificate]:
    """Solve ``target*x=0``, ``others*x>0``, and ``time*x>0`` exactly."""

    dimension = len(target)
    if len(time_functional) != dimension or any(
        len(row) != dimension for row in others.values()
    ):
        raise ValueError("strict face data have inconsistent dimensions")

    # A free coordinate x is x^+ - x^- with both parts nonnegative.  Because
    # the system is homogeneous and finite, strict positivity is equivalent
    # to simultaneous lower bounds >= 1 after one positive rescaling.
    def split(row: GoldenVector) -> tuple[Golden, ...]:
        return (*row, *(-entry for entry in row))

    inequalities = []
    for row in (*others.values(), time_functional):
        free_row = split(row)
        inequalities.append((tuple(-entry for entry in free_row), Golden(-1)))
    equality = (split(target), Golden(0))
    certificate = _phase_one_nonnegative(
        2 * dimension,
        equalities=(equality,),
        inequalities_le=tuple(inequalities),
    )
    if not certificate.passes:
        raise ArithmeticError("strict face phase-I audit failed")
    if not certificate.feasible:
        return None, certificate
    assert certificate.primal is not None
    point = tuple(
        certificate.primal[index] - certificate.primal[index + dimension]
        for index in range(dimension)
    )
    if golden_dot(target, point) != 0:
        raise ArithmeticError("strict face phase-I missed its target equality")
    if any(golden_dot(row, point).sign() <= 0 for row in others.values()):
        raise ArithmeticError("strict face phase-I returned a non-strict point")
    if golden_dot(time_functional, point).sign() <= 0:
        raise ArithmeticError("strict face phase-I returned the wrong time component")
    return point, certificate


def _solve_square_linear_system(
    matrix: Sequence[Sequence[Golden]], right_hand_side: Sequence[Golden]
) -> tuple[Golden, ...] | None:
    try:
        inverse = inverse_golden_matrix(matrix)
    except ZeroDivisionError:
        return None
    return _mat_vec(inverse, tuple(right_hand_side))


def decide_lorentzian_strict_face(
    target: Sequence[Golden | int],
    others: Mapping[Label, Sequence[Golden | int]],
    *,
    quadratic_gram: Sequence[Sequence[Golden | int]],
    interior: Sequence[Golden | int],
) -> LorentzianFaceDecision:
    """Decide strict timelike feasibility of one homogeneous face.

    The algorithm is finite and exact over ``Q(phi)``.  Linear phase I first
    decides whether the target hyperplane contains a point strictly satisfying
    every other inequality in the future component.  If its basic witness is
    not timelike, the routine maximizes the Lorentzian quadratic form on

    ``B(interior, x) = B(interior, interior)``.

    The form is strictly concave and coercive on that affine hyperplane.  A
    strict phase-I point supplies Slater's condition, so an exact enumeration
    of linearly independent KKT active supports is both necessary and
    sufficient.  Thus a nonpositive maximum is an exact causal-inactivity
    certificate, not a failed witness search.
    """

    wall = _vector(target)
    rows = {label: _vector(row) for label, row in others.items()}
    gram = _matrix(quadratic_gram)
    point_inside = _vector(interior)
    dimension = len(wall)
    if len(gram) != dimension or len(point_inside) != dimension:
        raise ValueError("Lorentzian face data have inconsistent dimensions")
    gram_inertia = inertia(gram)
    if (gram_inertia.positive, gram_inertia.negative, gram_inertia.zero) != (
        1,
        dimension - 1,
        0,
    ):
        raise ValueError("quadratic_gram must have Lorentzian inertia (1,n-1)")
    interior_square = golden_quadratic_value(gram, point_inside)
    if interior_square.sign() <= 0:
        raise ValueError("interior must be timelike")
    if golden_dot(wall, point_inside).sign() <= 0 or any(
        golden_dot(row, point_inside).sign() <= 0 for row in rows.values()
    ):
        raise ValueError("interior must satisfy every half-space strictly")

    time_functional = tuple(
        sum(
            (point_inside[row] * gram[row][column] for row in range(dimension)),
            Golden(0),
        )
        for column in range(dimension)
    )
    linear_point, phase_one = _strict_free_face_feasibility(
        wall, rows, time_functional
    )
    if linear_point is None:
        return LorentzianFaceDecision(
            False,
            "strict_future_linear_face_infeasible",
            None,
            phase_one,
        )
    if golden_quadratic_value(gram, linear_point).sign() > 0:
        return LorentzianFaceDecision(
            True,
            "phase_one_witness_timelike",
            linear_point,
            phase_one,
        )

    time_value = golden_dot(time_functional, linear_point)
    slater = tuple(
        interior_square * entry / time_value for entry in linear_point
    )
    equality_rows = (wall, time_functional)
    equality_rhs = (Golden(0), interior_square)
    labels = tuple(rows)
    maximum_point: GoldenVector | None = None
    maximum_labels: tuple[Label, ...] = ()
    maximum_multipliers: tuple[Golden, ...] = ()
    equality_multipliers: tuple[Golden, Golden] | None = None

    # In the quotient by the two equality rows, a minimal normal-cone support
    # has at most n-2 independent active inequalities.
    for support_size in range(dimension - 1):
        for support in combinations(labels, support_size):
            active_rows = tuple(rows[label] for label in support)
            constraint_rows = (*active_rows, *equality_rows)
            if rank(constraint_rows) != len(constraint_rows):
                continue
            multiplier_count = support_size + 2
            system_size = dimension + multiplier_count
            system = [
                [Golden(0) for _ in range(system_size)]
                for _ in range(system_size)
            ]
            rhs = [Golden(0) for _ in range(system_size)]
            for row in range(dimension):
                for column in range(dimension):
                    system[row][column] = Golden(2) * gram[row][column]
                for offset, functional in enumerate(constraint_rows):
                    system[row][dimension + offset] = functional[row]
            for offset, functional in enumerate(constraint_rows):
                equation_row = dimension + offset
                for column in range(dimension):
                    system[equation_row][column] = functional[column]
                if offset >= support_size:
                    rhs[equation_row] = equality_rhs[offset - support_size]
            solution = _solve_square_linear_system(system, rhs)
            if solution is None:
                continue
            candidate = tuple(solution[:dimension])
            active_multipliers = tuple(
                solution[dimension : dimension + support_size]
            )
            if any(value.sign() < 0 for value in active_multipliers):
                continue
            if golden_dot(wall, candidate) != 0:
                continue
            if golden_dot(time_functional, candidate) != interior_square:
                continue
            if any(golden_dot(row, candidate).sign() < 0 for row in rows.values()):
                continue
            candidate_equality_multipliers = tuple(
                solution[dimension + support_size :]
            )
            stationarity = tuple(
                Golden(2)
                * sum(
                    (gram[coordinate][column] * candidate[column] for column in range(dimension)),
                    Golden(0),
                )
                + sum(
                    (
                        active_multipliers[offset]
                        * active_rows[offset][coordinate]
                        for offset in range(support_size)
                    ),
                    Golden(0),
                )
                + sum(
                    (
                        candidate_equality_multipliers[offset]
                        * equality_rows[offset][coordinate]
                        for offset in range(2)
                    ),
                    Golden(0),
                )
                for coordinate in range(dimension)
            )
            if any(stationarity):
                continue
            maximum_point = candidate
            maximum_labels = support
            maximum_multipliers = active_multipliers
            equality_multipliers = candidate_equality_multipliers  # type: ignore[assignment]
            break
        if maximum_point is not None:
            break
    if maximum_point is None or equality_multipliers is None:
        raise ArithmeticError("exact Lorentzian KKT enumeration found no maximum")

    maximum_square = golden_quadratic_value(gram, maximum_point)
    if maximum_square.sign() <= 0:
        return LorentzianFaceDecision(
            False,
            "normalized_face_has_nonpositive_global_quadratic_maximum",
            None,
            phase_one,
            maximum_point,
            maximum_square,
            maximum_labels,
            maximum_multipliers,
            equality_multipliers,
        )

    if all(golden_dot(row, maximum_point).sign() > 0 for row in rows.values()):
        witness = maximum_point
    else:
        # Convex mixing with Slater makes every inequality strict.  Powers of
        # two are unbounded here; positivity at epsilon=0 proves termination.
        denominator = 2
        while True:
            epsilon = Golden(1) / denominator
            witness = tuple(
                (Golden(1) - epsilon) * maximum_point[index]
                + epsilon * slater[index]
                for index in range(dimension)
            )
            if golden_quadratic_value(gram, witness).sign() > 0:
                break
            denominator *= 2
    if (
        golden_dot(wall, witness) != 0
        or any(golden_dot(row, witness).sign() <= 0 for row in rows.values())
        or golden_dot(time_functional, witness).sign() <= 0
        or golden_quadratic_value(gram, witness).sign() <= 0
    ):
        raise ArithmeticError("Lorentzian active-set witness audit failed")
    return LorentzianFaceDecision(
        True,
        "positive_normalized_quadratic_maximum",
        witness,
        phase_one,
        maximum_point,
        maximum_square,
        maximum_labels,
        maximum_multipliers,
        equality_multipliers,
    )


@dataclass(frozen=True, slots=True)
class FacetWitnessAudit:
    wall_class: tuple[Label, ...]
    target_pairings_zero: bool
    other_pairings_strictly_positive: bool
    quadratic_value: Golden | None
    timelike: bool
    same_time_component: bool

    @property
    def passes(self) -> bool:
        return (
            self.target_pairings_zero
            and self.other_pairings_strictly_positive
            and self.timelike
            and self.same_time_component
        )


@dataclass(frozen=True, slots=True)
class PositiveDependenceAudit:
    target: Label
    terms: DependenceTerms
    coefficients_positive: bool
    sources_are_facets: bool
    exact_identity: bool

    @property
    def passes(self) -> bool:
        return (
            self.coefficients_positive
            and self.sources_are_facets
            and self.exact_identity
        )


@dataclass(frozen=True, slots=True)
class CausalInactivityAudit:
    """Audit that a linearly exposed class never reaches the time cone."""

    target: Label
    reason: str
    exact_nonfacet_decision: bool

    @property
    def passes(self) -> bool:
        return self.exact_nonfacet_decision


@dataclass(frozen=True, slots=True)
class GoldenHalfspaceCertificateAudit:
    """Complete exact audit of an asserted irredundant half-space system."""

    primitive_wall_count: int
    proportional_classes: tuple[PositiveProportionalClass, ...]
    facet_class_count: int
    redundant_class_count: int
    facet_partition_exact: bool
    redundancy_partition_exact: bool
    strict_interior: bool
    interior_quadratic_value: Golden | None
    interior_timelike: bool
    facet_witnesses: tuple[FacetWitnessAudit, ...]
    positive_dependences: tuple[PositiveDependenceAudit, ...]
    causal_inactivities: tuple[CausalInactivityAudit, ...]

    @property
    def passes(self) -> bool:
        return (
            self.facet_partition_exact
            and self.redundancy_partition_exact
            and self.strict_interior
            and self.interior_timelike
            and all(witness.passes for witness in self.facet_witnesses)
            and all(dependence.passes for dependence in self.positive_dependences)
            and all(inactivity.passes for inactivity in self.causal_inactivities)
        )


def audit_golden_halfspace_certificate(
    rows: Mapping[Label, Sequence[Golden | int]],
    *,
    facet_classes: Sequence[Sequence[Label]],
    facet_witnesses: Mapping[Label, Sequence[Golden | int]],
    redundancies: Mapping[Label, Sequence[tuple[Golden | int, Label]]],
    causal_inactivities: Mapping[Label, LorentzianFaceDecision] | None = None,
    interior: Sequence[Golden | int],
    quadratic_gram: Sequence[Sequence[Golden | int]] | None = None,
) -> GoldenHalfspaceCertificateAudit:
    """Verify a complete irredundant reduction without numerical decisions.

    ``facet_witnesses`` is keyed by any member of its asserted facet class.
    Every nonfacet positive-proportionality class must occur exactly once,
    either in ``redundancies`` or in ``causal_inactivities``.  Dependence
    sources are required to lie in asserted facet classes, preventing circular
    redundancy proofs.  A causal-inactivity decision is recomputed exactly and
    covers the distinct possibility that an inequality is linearly exposed
    but its face never reaches the selected Lorentzian time cone.
    If no quadratic form is supplied, timelikeness is deliberately left
    unproved and the overall audit does not pass.
    """

    vectors = {label: _vector(row) for label, row in rows.items()}
    interior_point = _vector(interior)
    classes = positive_proportional_classes(vectors)
    computed_by_set = {frozenset(item.members): item for item in classes}
    asserted_sets = tuple(frozenset(wall_class) for wall_class in facet_classes)
    asserted_set_collection = set(asserted_sets)
    facet_partition_exact = (
        len(asserted_sets) == len(asserted_set_collection)
        and all(wall_class in computed_by_set for wall_class in asserted_sets)
    )
    facet_labels = set().union(*asserted_sets) if asserted_sets else set()

    witness_audits = []
    witness_classes: set[frozenset[Label]] = set()
    for key, raw_point in facet_witnesses.items():
        containing = next((item for item in classes if key in item.members), None)
        if containing is None:
            wall_class: tuple[Label, ...] = (key,)
            target_zero = False
            other_positive = False
        else:
            wall_class = containing.members
            witness_classes.add(frozenset(wall_class))
            point = _vector(raw_point)
            target = set(wall_class)
            pairings = {
                label: golden_dot(functional, point)
                for label, functional in vectors.items()
            }
            target_zero = all(pairings[label] == 0 for label in target)
            other_positive = all(
                value.sign() > 0
                for label, value in pairings.items()
                if label not in target
            )
        point = _vector(raw_point)
        square = (
            golden_quadratic_value(quadratic_gram, point)
            if quadratic_gram is not None
            else None
        )
        time_pairing = (
            sum(
                (
                    interior_point[i]
                    * Golden.coerce(quadratic_gram[i][j])
                    * point[j]
                    for i in range(len(point))
                    for j in range(len(point))
                ),
                Golden(0),
            )
            if quadratic_gram is not None
            else None
        )
        witness_audits.append(
            FacetWitnessAudit(
                wall_class,
                target_zero,
                other_positive,
                square,
                square is not None and square.sign() > 0,
                time_pairing is not None and time_pairing.sign() > 0,
            )
        )
    facet_partition_exact = (
        facet_partition_exact
        and witness_classes == asserted_set_collection
        and len(facet_witnesses) == len(asserted_sets)
    )

    dependence_audits = []
    redundancy_classes: set[frozenset[Label]] = set()
    for target, raw_terms in redundancies.items():
        containing = next((item for item in classes if target in item.members), None)
        if containing is not None:
            redundancy_classes.add(frozenset(containing.members))
        terms: DependenceTerms = tuple(
            (Golden.coerce(coefficient), source) for coefficient, source in raw_terms
        )
        coefficients_positive = bool(terms) and all(
            coefficient.sign() > 0 for coefficient, _ in terms
        )
        sources_are_facets = all(source in facet_labels for _, source in terms)
        exact_identity = target in vectors and all(source in vectors for _, source in terms)
        if exact_identity:
            value = tuple(
                sum(
                    (
                        coefficient * vectors[source][coordinate]
                        for coefficient, source in terms
                    ),
                    Golden(0),
                )
                for coordinate in range(len(vectors[target]))
            )
            exact_identity = value == vectors[target]
        dependence_audits.append(
            PositiveDependenceAudit(
                target,
                terms,
                coefficients_positive,
                sources_are_facets,
                exact_identity,
            )
        )

    inactivity_audits = []
    inactivity_classes: set[frozenset[Label]] = set()
    for target, asserted_decision in (causal_inactivities or {}).items():
        containing = next((item for item in classes if target in item.members), None)
        if containing is None or quadratic_gram is None:
            inactivity_audits.append(
                CausalInactivityAudit(target, asserted_decision.reason, False)
            )
            continue
        wall_class = frozenset(containing.members)
        inactivity_classes.add(wall_class)
        other_rows = {
            label: row
            for label, row in vectors.items()
            if label not in wall_class
        }
        recomputed = decide_lorentzian_strict_face(
            vectors[target],
            other_rows,
            quadratic_gram=quadratic_gram,
            interior=interior,
        )
        agrees = (
            not recomputed.is_facet
            and not asserted_decision.is_facet
            and recomputed.reason == asserted_decision.reason
            and recomputed.maximum_quadratic_value
            == asserted_decision.maximum_quadratic_value
        )
        inactivity_audits.append(
            CausalInactivityAudit(target, recomputed.reason, agrees)
        )

    expected_redundancy_classes = set(computed_by_set) - asserted_set_collection
    redundancy_partition_exact = (
        not (redundancy_classes & inactivity_classes)
        and redundancy_classes | inactivity_classes == expected_redundancy_classes
        and len(redundancies) == len(redundancy_classes)
        and len(causal_inactivities or {}) == len(inactivity_classes)
    )

    strict_interior = all(
        golden_dot(functional, interior_point).sign() > 0
        for functional in vectors.values()
    )
    interior_square = (
        golden_quadratic_value(quadratic_gram, interior_point)
        if quadratic_gram is not None
        else None
    )
    return GoldenHalfspaceCertificateAudit(
        primitive_wall_count=len(vectors),
        proportional_classes=classes,
        facet_class_count=len(asserted_sets),
        redundant_class_count=len(expected_redundancy_classes),
        facet_partition_exact=facet_partition_exact,
        redundancy_partition_exact=redundancy_partition_exact,
        strict_interior=strict_interior,
        interior_quadratic_value=interior_square,
        interior_timelike=interior_square is not None and interior_square.sign() > 0,
        facet_witnesses=tuple(witness_audits),
        positive_dependences=tuple(dependence_audits),
        causal_inactivities=tuple(inactivity_audits),
    )

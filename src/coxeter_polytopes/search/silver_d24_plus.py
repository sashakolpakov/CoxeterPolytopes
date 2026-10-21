"""First-packet exact silver search on ``<-2sqrt(2)> + D24+``.

The positive block is the even unimodular Niemeier lattice with root system
``D24``.  We use a glue-aware integral basis: the simple roots
``alpha_2,...,alpha_24`` together with the half-spinor
``g=(1/2,...,1/2)``.  Its Euclidean Gram matrix is integral, even, and has
determinant one.  A dominant Dynkin-weight vector is retained exactly when
its reconstructed coordinates in this basis lie in ``Z[sqrt(2)]``; this
enumerates both D24 root/spinor cosets without pretending the Cartan root
lattice is unimodular.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import groupby
from typing import Sequence, TypeAlias

from coxeter_polytopes.exact.silver import SILVER_ONE, SILVER_ZERO, SQRT2, Silver
from coxeter_polytopes.search.bugaenko_vinberg import (
    _fraction_inverse,
    _prefix_schur_complements,
)
from coxeter_polytopes.search.silver_e8_cubes import (
    SilverE8CubeSeed,
    SilverH24ControlLevel,
    TWO_SQRT2_E8_CUBE_SEED,
    _bounded_for_caps,
    _mat_vec,
    _quadratic,
    coefficient_is_arithmetically_admissible,
)


Vector: TypeAlias = tuple[Silver, ...]
RationalMatrix: TypeAlias = tuple[tuple[Fraction, ...], ...]
IntegerMatrix: TypeAlias = tuple[tuple[int, ...], ...]


UNIT_D24_PLUS_NORM_PROOF = (
    "Put pi=sqrt(2) and epsilon=1+pi. Since D24+ is even unimodular and "
    "a=epsilon is a unit, q(r)|2I(r) has no support away from (pi). If "
    "the control x is a pi-adic unit, q=-epsilon*x^2 mod pi^2, hence "
    "v_pi(q)=0. Otherwise primitivity and unimodularity make the spatial "
    "pairing ideal a unit; evenness gives v_pi(q)>=2 and q|2I gives "
    "v_pi(q)<=2. Thus the norm ideals are O and (2). Every totally "
    "positive unit of Z[sqrt(2)] is a square, leaving canonical norms 1,2."
)


UNIT_D24_PLUS_SEED = SilverE8CubeSeed(
    "silver-H24-D24+-a-1+sqrt2",
    SILVER_ONE + SQRT2,
    (SILVER_ONE, Silver(2)),
    UNIT_D24_PLUS_NORM_PROOF,
)


def _dot(left: Sequence[Fraction], right: Sequence[Fraction]) -> Fraction:
    return sum((x * y for x, y in zip(left, right)), start=Fraction(0))


@lru_cache(maxsize=1)
def d24_simple_roots_euclidean() -> tuple[tuple[Fraction, ...], ...]:
    roots = []
    for index in range(23):
        root = [Fraction(0)] * 24
        root[index] = Fraction(1)
        root[index + 1] = Fraction(-1)
        roots.append(tuple(root))
    final = [Fraction(0)] * 24
    final[-2] = final[-1] = Fraction(1)
    roots.append(tuple(final))
    return tuple(roots)


@lru_cache(maxsize=1)
def d24_cartan_gram() -> RationalMatrix:
    roots = d24_simple_roots_euclidean()
    return tuple(tuple(_dot(left, right) for right in roots) for left in roots)


@lru_cache(maxsize=1)
def d24_plus_basis_euclidean() -> tuple[tuple[Fraction, ...], ...]:
    roots = d24_simple_roots_euclidean()
    spinor = tuple(Fraction(1, 2) for _ in range(24))
    return (*roots[1:], spinor)


@lru_cache(maxsize=1)
def d24_plus_gram() -> RationalMatrix:
    basis = d24_plus_basis_euclidean()
    return tuple(tuple(_dot(left, right) for right in basis) for left in basis)


@lru_cache(maxsize=1)
def _d24_plus_euclidean_basis_inverse() -> RationalMatrix:
    basis = d24_plus_basis_euclidean()
    columns = tuple(
        tuple(basis[column][row] for column in range(24)) for row in range(24)
    )
    return _fraction_inverse(columns)


def d24_plus_integral_basis_coordinates(
    coordinates: Sequence[Silver],
) -> Vector:
    if len(coordinates) != 24:
        raise ValueError("D24+ has 24 Euclidean coordinates")
    inverse = _d24_plus_euclidean_basis_inverse()
    return tuple(
        sum(
            (
                Silver(inverse[row][column]) * coordinates[column]
                for column in range(24)
                if inverse[row][column]
            ),
            start=SILVER_ZERO,
        )
        for row in range(24)
    )


@lru_cache(maxsize=1)
def _d24_dominant_data() -> tuple[RationalMatrix, tuple[RationalMatrix, ...]]:
    roots = d24_simple_roots_euclidean()
    basis = d24_plus_basis_euclidean()
    # M_ij=<alpha_i,basis_j>; for coordinates c in the D24+ basis,
    # the dominant weights are z=M c and c=M^{-1}z.
    pairing = tuple(tuple(_dot(root, vector) for vector in basis) for root in roots)
    inverse_pairing = _fraction_inverse(pairing)
    gram = d24_plus_gram()
    energy = tuple(
        tuple(
            sum(
                inverse_pairing[k][i] * gram[k][ell] * inverse_pairing[ell][j]
                for k in range(24)
                for ell in range(24)
            )
            for j in range(24)
        )
        for i in range(24)
    )
    # This must be C_D24^{-1}; check rather than relying on a diagram convention.
    cartan_inverse = _fraction_inverse(d24_cartan_gram())
    if energy != cartan_inverse:
        raise AssertionError("D24+ dominant energy differs from inverse Cartan")
    return inverse_pairing, _prefix_schur_complements(energy)


def _integer_determinant(matrix: RationalMatrix) -> int:
    size = len(matrix)
    work = [list(row) for row in matrix]
    determinant = Fraction(1)
    sign = 1
    for column in range(size):
        pivot = next(row for row in range(column, size) if work[row][column])
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
        value = work[column][column]
        determinant *= value
        for row in range(column + 1, size):
            factor = work[row][column] / value
            for index in range(column + 1, size):
                work[row][index] -= factor * work[column][index]
    result = determinant * sign
    if result.denominator != 1:
        raise AssertionError("integral Gram determinant was fractional")
    return result.numerator


@dataclass(frozen=True, slots=True)
class SilverDominantD24PlusWeight:
    square: Silver
    weights: Vector
    coordinates: Vector


@dataclass(frozen=True, slots=True)
class SilverD24PlusEnumeration:
    target: Silver
    root_norm: Silver
    alphabet: tuple[Silver, ...]
    visited_by_depth: tuple[int, ...]
    pruned_by_depth: tuple[int, ...]
    nonintegral_leaf_count: int
    vectors: tuple[SilverDominantD24PlusWeight, ...]


@lru_cache(maxsize=None)
def enumerate_dominant_d24_plus_weights(
    target: Silver, root_norm: Silver
) -> SilverD24PlusEnumeration:
    inverse_pairing, schur = _d24_dominant_data()
    inverse_cartan = _fraction_inverse(d24_cartan_gram())
    twice_target = 2 * target
    alphabet = tuple(
        value
        for value in _bounded_for_caps((twice_target, twice_target))
        if value.sign(0) >= 0
        and all((twice_target - value * value).sign(place) >= 0 for place in (0, 1))
        and (2 * value / root_norm).is_integral
    )
    visited = [0] * 25
    pruned = [0] * 25
    nonintegral = 0
    output: list[SilverDominantD24PlusWeight] = []

    def extend(prefix: Vector) -> None:
        nonlocal nonintegral
        depth = len(prefix)
        visited[depth] += 1
        if depth:
            minimum = _quadratic(prefix, schur[depth - 1])
            if any((target - minimum).sign(place) < 0 for place in (0, 1)):
                pruned[depth] += 1
                return
        if depth == 24:
            square = _quadratic(prefix, inverse_cartan)
            coordinates = _mat_vec(inverse_pairing, prefix)
            if not all(entry.is_integral for entry in coordinates):
                nonintegral += 1
                return
            output.append(SilverDominantD24PlusWeight(square, prefix, coordinates))
            return
        for entry in alphabet:
            extend((*prefix, entry))

    extend(())
    vectors = tuple(sorted(output, key=lambda item: (item.square, item.weights)))
    return SilverD24PlusEnumeration(
        target,
        root_norm,
        alphabet,
        tuple(visited),
        tuple(pruned),
        nonintegral,
        vectors,
    )


def _silver_dominant_euclidean(vector: Sequence[Silver]) -> bool:
    return all(
        (vector[index] - vector[index + 1]).sign(0) >= 0
        for index in range(len(vector) - 2)
    ) and (vector[-2] - vector[-1]).sign(0) >= 0 and (
        vector[-2] + vector[-1]
    ).sign(0) >= 0


def _in_silver_d24_root_coset(vector: Sequence[Silver]) -> bool:
    if len(vector) != 24 or not all(value.is_integral for value in vector):
        return False
    return all(
        sum(getattr(value, coefficient) for value in vector).numerator % 2 == 0
        for coefficient in ("a", "b")
    )


def _in_silver_d24_spinor_coset(vector: Sequence[Silver]) -> bool:
    if len(vector) != 24:
        return False
    parities = []
    translated_sums = []
    for coefficient in ("a", "b"):
        doubled = [2 * getattr(value, coefficient) for value in vector]
        if any(value.denominator != 1 for value in doubled):
            return False
        numerators = [value.numerator for value in doubled]
        parity = numerators[0] % 2
        if any(value % 2 != parity for value in numerators):
            return False
        parities.append(parity)
        translated_sums.append(sum((value - parity) // 2 for value in numerators))
    return any(parities) and all(value % 2 == 0 for value in translated_sums)


def _silver_coefficient_box(target: Silver) -> tuple[Silver, ...]:
    return tuple(
        value
        for value in _bounded_for_caps((target, target))
        if all((target - value * value).sign(place) >= 0 for place in (0, 1))
    )


def _silver_half_spinor_alphabet(
    target: Silver, residue: Silver
) -> tuple[Silver, ...]:
    return tuple(
        sorted(
            value / 2
            for value in _silver_coefficient_box(4 * target)
            if (value - residue).is_integral
            and int((value - residue).a) % 2 == 0
            and int((value - residue).b) % 2 == 0
            and all(
                (target - (value / 2) ** 2).sign(place) >= 0
                for place in (0, 1)
            )
        )
    )


def _silver_sum_squares_from_alphabet(
    target: Silver,
    *,
    alphabet: Sequence[Silver],
    membership,
    coset: str,
) -> tuple[tuple[Vector, str], ...]:
    positive = tuple(
        sorted(
            {
                value if value.sign(0) >= 0 else -value
                for value in alphabet
                if value
            },
            reverse=True,
        )
    )
    squares = tuple(value * value for value in positive)

    @lru_cache(maxsize=None)
    def solve(
        index: int, slots: int, remainder: Silver
    ) -> tuple[tuple[int, ...], ...]:
        if any(remainder.sign(place) < 0 for place in (0, 1)):
            return ()
        if index == len(positive):
            return ((),) if not remainder else ()
        output = []
        running = remainder
        for multiplicity in range(slots + 1):
            if all(running.sign(place) >= 0 for place in (0, 1)):
                for suffix in solve(index + 1, slots - multiplicity, running):
                    output.append((multiplicity, *suffix))
            running -= squares[index]
        return tuple(output)

    output = []
    for multiplicities in solve(0, 24, target):
        nonnegative = tuple(
            value
            for value, multiplicity in zip(positive, multiplicities)
            for _ in range(multiplicity)
        ) + (SILVER_ZERO,) * (24 - sum(multiplicities))
        representatives = (nonnegative,)
        if nonnegative[-1]:
            representatives = (
                nonnegative,
                (*nonnegative[:-1], -nonnegative[-1]),
            )
        for representative in representatives:
            if _silver_dominant_euclidean(representative) and membership(
                representative
            ):
                output.append((representative, coset))
    return tuple(output)


@lru_cache(maxsize=None)
def enumerate_dominant_d24_plus_euclidean(
    target: Silver,
) -> tuple[tuple[Vector, str], ...]:
    output = list(
        _silver_sum_squares_from_alphabet(
            target,
            alphabet=_silver_coefficient_box(target),
            membership=_in_silver_d24_root_coset,
            coset="D24",
        )
    )
    for residue in (SILVER_ONE, SQRT2, SILVER_ONE + SQRT2):
        output.extend(
            _silver_sum_squares_from_alphabet(
                target,
                alphabet=_silver_half_spinor_alphabet(target, residue),
                membership=_in_silver_d24_spinor_coset,
                coset=f"spinor-{residue}",
            )
        )
    return tuple(sorted(set(output)))


def d24_plus_pair(
    seed: SilverE8CubeSeed, left: Sequence[Silver], right: Sequence[Silver]
) -> Silver:
    if len(left) != 25 or len(right) != 25:
        raise ValueError("a silver H24 D24+ vector has rank 25")
    gram = d24_plus_gram()
    return -seed.timelike_coefficient * left[0] * right[0] + sum(
        (
            left[1 + i] * Silver(gram[i][j]) * right[1 + j]
            for i in range(24)
            for j in range(24)
            if gram[i][j]
        ),
        start=SILVER_ZERO,
    )


@dataclass(frozen=True, slots=True)
class SilverD24PlusCandidate:
    norm: Silver
    height: Silver
    control: Silver
    target: Silver
    vector: Vector


@lru_cache(maxsize=None)
def enumerate_d24_plus_control_level(
    seed: SilverE8CubeSeed, level: SilverH24ControlLevel
) -> tuple[SilverD24PlusCandidate, ...]:
    enumeration = enumerate_dominant_d24_plus_weights(level.target, level.norm)
    roots = []
    for spatial in enumeration.vectors:
        if spatial.square != level.target:
            continue
        vector = (level.control, *spatial.coordinates)
        if d24_plus_pair(seed, vector, vector) != level.norm:
            raise AssertionError("D24+ reconstruction returned the wrong norm")
        roots.append(
            SilverD24PlusCandidate(
                level.norm, level.height, level.control, level.target, vector
            )
        )
    return tuple(sorted(roots, key=lambda root: root.vector))


@lru_cache(maxsize=None)
def enumerate_d24_plus_unit_control_level(
    level: SilverH24ControlLevel,
) -> tuple[SilverD24PlusCandidate, ...]:
    seed = UNIT_D24_PLUS_SEED
    roots = []
    for spatial, _coset in enumerate_dominant_d24_plus_euclidean(level.target):
        coordinates = d24_plus_integral_basis_coordinates(spatial)
        if not all(value.is_integral for value in coordinates):
            raise AssertionError("a D24+ Euclidean vector left its integral basis")
        vector = (level.control, *coordinates)
        if d24_plus_pair(seed, vector, vector) != level.norm:
            raise AssertionError("the unit D24+ reconstruction has wrong norm")
        roots.append(
            SilverD24PlusCandidate(
                level.norm,
                level.height,
                level.control,
                level.target,
                vector,
            )
        )
    return tuple(sorted(roots, key=lambda root: root.vector))


@lru_cache(maxsize=None)
def _control_levels(
    seed: SilverE8CubeSeed, height_bound: int
) -> tuple[SilverH24ControlLevel, ...]:
    if seed not in (TWO_SQRT2_E8_CUBE_SEED, UNIT_D24_PLUS_SEED):
        raise ValueError("unsupported D24+ timelike coefficient")
    levels = set()
    coefficient = seed.timelike_coefficient
    for norm in seed.canonical_mirror_norms:
        caps = (norm * height_bound / coefficient, norm / (-coefficient))
        for control in _bounded_for_caps(caps):
            if control.sign(0) <= 0 or not (2 * coefficient * control / norm).is_integral:
                continue
            target = norm + coefficient * control * control
            if target.sign(1) < 0:
                continue
            height = coefficient * control * control / norm
            if (Silver(height_bound) - height).sign(0) < 0:
                continue
            levels.add(SilverH24ControlLevel(norm, height, control, target))
    return tuple(sorted(levels, key=lambda item: (item.height, item.norm, item.control)))


@dataclass(frozen=True, slots=True)
class SilverD24PlusPacket:
    height: Silver
    candidate_count: int
    roots: tuple[SilverD24PlusCandidate, ...]


@dataclass(frozen=True, slots=True)
class SilverD24PlusSearch:
    seed: SilverE8CubeSeed
    height_bound: int
    control_levels: tuple[SilverH24ControlLevel, ...]
    packets: tuple[SilverD24PlusPacket, ...]
    accepted_roots: tuple[SilverD24PlusCandidate, ...]
    walls: tuple[Vector, ...]
    wall_count: int


@lru_cache(maxsize=1)
def d24_height_zero_simple_roots() -> tuple[Vector, ...]:
    """Return the genuine negative D24 simple roots in the glue basis.

    If ``z_i=<alpha_i,y>`` and ``c`` are the coordinates of ``y`` in the
    displayed D24+ basis, then ``c=M^{-1}z``.  For ``y=alpha_j`` the vector
    ``z`` is the j-th Cartan column.  Using unit vectors here would instead
    construct fundamental weights; this Cartan multiplication is essential.
    """

    inverse_pairing, _schur = _d24_dominant_data()
    cartan = d24_cartan_gram()
    roots = tuple(
        (
            SILVER_ZERO,
            *tuple(
                -value
                for value in _mat_vec(
                    inverse_pairing,
                    tuple(Silver(cartan[row][column]) for row in range(24)),
                )
            ),
        )
        for column in range(24)
    )
    seed = TWO_SQRT2_E8_CUBE_SEED
    gram = tuple(
        tuple(d24_plus_pair(seed, left, right) for right in roots)
        for left in roots
    )
    expected = tuple(
        tuple(Silver(value) for value in row) for row in d24_cartan_gram()
    )
    if gram != expected:
        raise AssertionError("height-zero glue-basis roots do not have D24 Cartan Gram")
    return roots


@lru_cache(maxsize=None)
def run_silver_d24_plus_first_shell(
    height_bound: int = 10,
) -> SilverD24PlusSearch:
    seed = TWO_SQRT2_E8_CUBE_SEED
    if not coefficient_is_arithmetically_admissible(seed.timelike_coefficient):
        raise AssertionError("the silver D24+ seed is not admissible")
    gram = d24_plus_gram()
    if _integer_determinant(gram) != 1 or any(
        entry.denominator != 1 for row in gram for entry in row
    ) or any(gram[index][index].numerator % 2 for index in range(24)):
        raise AssertionError("the displayed D24+ basis is not even unimodular")
    levels = _control_levels(seed, height_bound)
    simple_roots = d24_height_zero_simple_roots()
    walls = list(simple_roots)
    accepted = []
    packets = []
    for height, iterator in groupby(levels, key=lambda item: item.height):
        packet_levels = tuple(iterator)
        candidates = tuple(
            root
            for level in packet_levels
            for root in enumerate_d24_plus_control_level(seed, level)
        )
        survivors = tuple(
            root
            for root in candidates
            if all(d24_plus_pair(seed, root.vector, wall).sign(0) <= 0 for wall in walls)
        )
        if any(
            d24_plus_pair(seed, left.vector, right.vector).sign(0) > 0
            for index, left in enumerate(survivors)
            for right in survivors[index + 1 :]
        ):
            raise AssertionError("same-height D24+ survivors are not mutually acute")
        accepted.extend(survivors)
        walls.extend(root.vector for root in survivors)
        packets.append(SilverD24PlusPacket(height, len(candidates), survivors))
        if survivors:
            break
    return SilverD24PlusSearch(
        seed,
        height_bound,
        levels,
        tuple(packets),
        tuple(accepted),
        tuple(walls),
        len(walls),
    )


@lru_cache(maxsize=None)
def run_silver_d24_plus_unit_prefix(
    height_bound: int = 30,
) -> SilverD24PlusSearch:
    """Run the glue-complete unit-coefficient D24+ prefix exactly."""

    seed = UNIT_D24_PLUS_SEED
    if not coefficient_is_arithmetically_admissible(seed.timelike_coefficient):
        raise AssertionError("the unit silver D24+ seed is not admissible")
    gram = d24_plus_gram()
    if _integer_determinant(gram) != 1 or any(
        entry.denominator != 1 for row in gram for entry in row
    ) or any(gram[index][index].numerator % 2 for index in range(24)):
        raise AssertionError("the displayed D24+ basis is not even unimodular")
    levels = _control_levels(seed, height_bound)
    walls = list(d24_height_zero_simple_roots())
    accepted = []
    packets = []
    for height, iterator in groupby(levels, key=lambda item: item.height):
        packet_levels = tuple(iterator)
        candidates = tuple(
            root
            for level in packet_levels
            for root in enumerate_d24_plus_unit_control_level(level)
        )
        survivors = tuple(
            root
            for root in candidates
            if all(
                d24_plus_pair(seed, root.vector, wall).sign(0) <= 0
                for wall in walls
            )
        )
        if any(
            d24_plus_pair(seed, left.vector, right.vector).sign(0) > 0
            for index, left in enumerate(survivors)
            for right in survivors[index + 1 :]
        ):
            raise AssertionError("same-height unit D24+ survivors are not acute")
        accepted.extend(survivors)
        walls.extend(root.vector for root in survivors)
        packets.append(SilverD24PlusPacket(height, len(candidates), survivors))
    return SilverD24PlusSearch(
        seed,
        height_bound,
        levels,
        tuple(packets),
        tuple(accepted),
        tuple(walls),
        len(walls),
    )


@dataclass(frozen=True, slots=True)
class SilverD24PlusExtremeRay:
    vector: Vector
    active_facets: tuple[int, ...]
    norm: Silver
    kind: str


@dataclass(frozen=True, slots=True)
class SilverD24PlusInsertionProfile:
    facet_index: int
    negative_ray_count: int
    zero_ray_count: int
    positive_ray_count: int
    output_ray_count: int


@dataclass(frozen=True, slots=True)
class SilverD24PlusGeometry:
    extreme_rays: tuple[SilverD24PlusExtremeRay, ...]
    insertion_profiles: tuple[SilverD24PlusInsertionProfile, ...]
    finite_vertex_count: int
    cusp_count: int
    hyperideal_ray_count: int
    irredundant: bool
    finite_volume: bool
    compact: bool


def _silver_projective_key(vector: Sequence[Silver]) -> Vector:
    pivot = next(value for value in vector if value)
    return tuple(value / pivot for value in vector)


def _silver_functional_value(
    functional: Sequence[Silver], vector: Sequence[Silver]
) -> Silver:
    return sum(
        (
            coefficient * coordinate
            for coefficient, coordinate in zip(functional, vector)
            if coefficient and coordinate
        ),
        start=SILVER_ZERO,
    )


def _silver_row_rank(rows: Sequence[Sequence[Silver]], width: int) -> int:
    work = [list(row) for row in rows]
    rank = 0
    for column in range(width):
        pivot = next((row for row in range(rank, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        for row in range(rank + 1, len(work)):
            if not work[row][column]:
                continue
            factor = work[row][column] / pivot_value
            for index in range(column, width):
                work[row][index] -= factor * work[rank][index]
        rank += 1
        if rank == len(work):
            break
    return rank


def _silver_kernel(rows: Sequence[Sequence[Silver]]) -> Vector | None:
    if not rows:
        return None
    width = len(rows[0])
    work = [list(row) for row in rows]
    pivots = []
    rank = 0
    for column in range(width):
        pivot = next((row for row in range(rank, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        value = work[rank][column]
        work[rank] = [entry / value for entry in work[rank]]
        for row in range(len(work)):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                work[row][index] - factor * work[rank][index]
                for index in range(width)
            ]
        pivots.append(column)
        rank += 1
    if width - rank != 1:
        return None
    free = next(index for index in range(width) if index not in pivots)
    vector = [SILVER_ZERO] * width
    vector[free] = SILVER_ONE
    for row, column in enumerate(pivots):
        vector[column] = -work[row][free]
    return tuple(vector)


def _d24_plus_functional(
    seed: SilverE8CubeSeed, wall: Sequence[Silver]
) -> Vector:
    gram = d24_plus_gram()
    return (
        -seed.timelike_coefficient * wall[0],
        *tuple(
            sum(
                (
                    Silver(gram[index][other]) * wall[1 + other]
                    for other in range(24)
                    if gram[index][other]
                ),
                start=SILVER_ZERO,
            )
            for index in range(24)
        ),
    )


@lru_cache(maxsize=None)
def analyze_silver_d24_plus_cone(
    search: SilverD24PlusSearch,
) -> SilverD24PlusGeometry:
    """Enumerate the exact projective rays of the current 24-cone.

    The first 25 walls form a simplicial cone.  Each remaining simultaneous
    wall is inserted by exact double description; adjacency is certified by
    rank 23 of the common active wall functionals.
    """

    ambient_rank = 25
    walls = search.walls
    if len(walls) < ambient_rank:
        raise ValueError("a pointed H24 cone needs at least 25 walls")
    functionals = tuple(_d24_plus_functional(search.seed, wall) for wall in walls)
    initial = tuple(range(ambient_rank))
    if _silver_row_rank(
        tuple(functionals[index] for index in initial), ambient_rank
    ) != ambient_rank:
        raise ValueError("the D24 simple walls and first cap are not independent")

    processed = list(initial)
    rays: dict[Vector, tuple[Vector, frozenset[int]]] = {}
    for omitted in initial:
        support = tuple(index for index in initial if index != omitted)
        ray = _silver_kernel(tuple(functionals[index] for index in support))
        if ray is None:
            raise AssertionError("the independent seed cone lost a ray")
        values = tuple(_silver_functional_value(functionals[index], ray) for index in initial)
        if all(value.sign(0) <= 0 for value in values):
            oriented = ray
        elif all(value.sign(0) >= 0 for value in values):
            oriented = tuple(-value for value in ray)
        else:
            raise AssertionError("the seed cone has inconsistent ray orientation")
        active = frozenset(
            index
            for index in initial
            if not _silver_functional_value(functionals[index], oriented)
        )
        rays[_silver_projective_key(oriented)] = (oriented, active)

    profiles = []
    rank_cache: dict[frozenset[int], int] = {}
    for new_facet in range(ambient_rank, len(walls)):
        functional = functionals[new_facet]
        negative = []
        zero = []
        positive = []
        for key, (ray, active) in rays.items():
            value = _silver_functional_value(functional, ray)
            destination = negative if value.sign(0) < 0 else positive if value.sign(0) > 0 else zero
            destination.append((key, ray, active, value))
        updated = {key: (ray, active) for key, ray, active, _value in negative}
        updated.update(
            {
                key: (ray, active | frozenset((new_facet,)))
                for key, ray, active, _value in zero
            }
        )
        for _negative_key, negative_ray, negative_active, negative_value in negative:
            for _positive_key, positive_ray, positive_active, positive_value in positive:
                common = negative_active & positive_active
                if len(common) < ambient_rank - 2:
                    continue
                rank = rank_cache.get(common)
                if rank is None:
                    rank = _silver_row_rank(
                        tuple(functionals[index] for index in common), ambient_rank
                    )
                    rank_cache[common] = rank
                if rank != ambient_rank - 2:
                    continue
                ray = tuple(
                    positive_value * negative_ray[index]
                    - negative_value * positive_ray[index]
                    for index in range(ambient_rank)
                )
                active = frozenset(
                    index
                    for index in (*processed, new_facet)
                    if not _silver_functional_value(functionals[index], ray)
                )
                updated[_silver_projective_key(ray)] = (ray, active)
        rays = updated
        processed.append(new_facet)
        profiles.append(
            SilverD24PlusInsertionProfile(
                new_facet,
                len(negative),
                len(zero),
                len(positive),
                len(rays),
            )
        )

    output = []
    for ray, active in rays.values():
        if _silver_row_rank(
            tuple(functionals[index] for index in active), ambient_rank
        ) != ambient_rank - 1:
            raise AssertionError("double description returned a non-extreme ray")
        norm = d24_plus_pair(search.seed, ray, ray)
        sign = norm.sign(0)
        output.append(
            SilverD24PlusExtremeRay(
                ray,
                tuple(sorted(active)),
                norm,
                "finite" if sign < 0 else "cusp" if sign == 0 else "hyperideal",
            )
        )
    ordered = tuple(sorted(output, key=lambda ray: _silver_projective_key(ray.vector)))
    finite = sum(ray.kind == "finite" for ray in ordered)
    cusps = sum(ray.kind == "cusp" for ray in ordered)
    hyperideal = sum(ray.kind == "hyperideal" for ray in ordered)
    irredundant = bool(ordered) and all(
        any(
            _silver_functional_value(functionals[facet], ray.vector).sign(0) < 0
            for ray in ordered
        )
        for facet in range(len(functionals))
    )
    finite_volume = irredundant and hyperideal == 0
    return SilverD24PlusGeometry(
        ordered,
        tuple(profiles),
        finite,
        cusps,
        hyperideal,
        irredundant,
        finite_volume,
        finite_volume and cusps == 0,
    )


__all__ = [
    "SilverD24PlusCandidate",
    "SilverD24PlusEnumeration",
    "SilverD24PlusPacket",
    "SilverD24PlusSearch",
    "SilverD24PlusExtremeRay",
    "SilverD24PlusGeometry",
    "SilverD24PlusInsertionProfile",
    "analyze_silver_d24_plus_cone",
    "d24_cartan_gram",
    "d24_height_zero_simple_roots",
    "d24_plus_gram",
    "d24_plus_pair",
    "enumerate_dominant_d24_plus_weights",
    "run_silver_d24_plus_first_shell",
]

"""Solver-ready global fusion constraints for golden Coxeter Gram matrices.

The input is a simple-polytope nerve (maximal vertices), finite labels in
``{2, 3, 5}``, and finite domains of dotted products in ``Z[phi]``.  Once the
labels internal to a finite vertex are fixed, the vertex Gram matrix is a
constant positive-definite matrix.  Its adjugate turns the rank-one Schur
condition into denominator-free equations over ``Z[phi]``.

For a vertex ``I``, complement ``J``, ``Delta=det(H_I)``, and outside
incidence columns ``v_j``, put

    D_j    = 2*Delta - v_j^t adj(H_I) v_j,
    R_jk   = Delta*h_jk - v_j^t adj(H_I) v_k.

The target rank and oriented-side conditions are exactly

    D_j(identity) < 0,       D_j(conjugate) > 0,
    R_jk**2 = D_j*D_k,       R_jk(identity) < 0.

Already before ``h_jk`` is chosen, rank one forces every nonzero ``D_j`` to
lie in one square class in ``F*/F*^2`` for ``F = Q(sqrt(5))``.

Thus only ``choose(|J|, 2)`` golden equations are needed.  For the 9/14
problem this is ten equations per vertex.  The module exposes both an exact
finite-domain propagator and a symbolic two-rational-component polynomial
system suitable for translation to SMT, SAT-modulo-tables, or a CAS.

This is a Gram feasibility compiler.  It does not enumerate all 9/14 nerves,
bound dotted-product domains, prove K3 effectivity, or certify Vinberg
completeness.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations
from math import isqrt
from typing import Mapping, Sequence

import sympy as sp

from coxeter_polytopes.exact.golden import Golden, PHI, ZERO, RationalLike
from coxeter_polytopes.exact.matrix import Inertia, determinant, inertia

from .k3 import k3_trace_square
from .golden_lanner import golden_lanner_blocks


Pair = tuple[int, int]
GoldenMatrix = tuple[tuple[Golden, ...], ...]


def _pair(left: int, right: int) -> Pair:
    if left == right:
        raise ValueError("a Gram pair needs two distinct facets")
    return (left, right) if left < right else (right, left)


def _golden_sort_key(value: Golden) -> tuple[Fraction, Fraction]:
    return value.a, value.b


def _rational_square_root(value: Fraction) -> Fraction | None:
    if value < 0:
        return None
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator * numerator != value.numerator:
        return None
    if denominator * denominator != value.denominator:
        return None
    return Fraction(numerator, denominator)


def golden_square_root(value: Golden | RationalLike) -> Golden | None:
    """Return an exact square root in ``Q(phi)``, if one exists.

    Write ``value = u + v*sqrt(5)`` and a possible root as
    ``p + q*sqrt(5)``.  Its field norm is one of the two rational square
    roots of ``Norm(value)``; this determines the rational candidates
    ``p**2`` and ``q**2``.  Testing their finitely many signs gives an exact
    square-class decision without numerical algebraic-number routines.
    """

    target = Golden.coerce(value)
    if not target:
        return ZERO
    norm_root = _rational_square_root(target.norm())
    if norm_root is None:
        return None
    rational_part = target.a + target.b / 2
    root_norms = (norm_root,) if norm_root == 0 else (norm_root, -norm_root)
    for root_norm in root_norms:
        p_abs = _rational_square_root((rational_part + root_norm) / 2)
        q_abs = _rational_square_root((rational_part - root_norm) / 10)
        if p_abs is None or q_abs is None:
            continue
        p_values = (p_abs,) if p_abs == 0 else (p_abs, -p_abs)
        q_values = (q_abs,) if q_abs == 0 else (q_abs, -q_abs)
        for p in p_values:
            for q in q_values:
                candidate = Golden(p - q, 2 * q)
                if candidate * candidate == target:
                    return candidate
    return None


def _coerce_domain(
    value: Golden | RationalLike | Sequence[Golden | RationalLike],
) -> tuple[Golden, ...]:
    if isinstance(value, (Golden, int, Fraction)):
        entries = (Golden.coerce(value),)
    else:
        entries = tuple(Golden.coerce(entry) for entry in value)
    if not entries:
        raise ValueError("a Gram-entry domain must be nonempty")
    return tuple(sorted(set(entries), key=_golden_sort_key))


def _normalize_vertices(
    vertices: Sequence[Sequence[int]],
    *,
    facet_count: int,
    dimension: int,
) -> tuple[tuple[int, ...], ...]:
    normalized: list[tuple[int, ...]] = []
    seen: set[tuple[int, ...]] = set()
    for raw_vertex in vertices:
        vertex = tuple(sorted(raw_vertex))
        if len(vertex) != dimension:
            raise ValueError(f"every vertex must contain exactly {dimension} facets")
        if len(vertex) != len(set(vertex)):
            raise ValueError("a vertex contains a repeated facet")
        if any(index < 0 or index >= facet_count for index in vertex):
            raise IndexError("a vertex facet index is out of range")
        if vertex not in seen:
            normalized.append(vertex)
            seen.add(vertex)
    if not normalized:
        raise ValueError("the nerve must contain at least one maximal vertex")
    return tuple(sorted(normalized))


def _intersecting_pairs(vertices: Sequence[Sequence[int]]) -> frozenset[Pair]:
    return frozenset(
        _pair(left, right)
        for vertex in vertices
        for left, right in combinations(vertex, 2)
    )


def minimal_missing_faces(
    vertices: Sequence[Sequence[int]],
    *,
    facet_count: int,
    max_order: int = 5,
) -> tuple[tuple[int, ...], ...]:
    """Return all minimal nonfaces through ``max_order`` from maximal faces."""

    maximal = tuple(frozenset(vertex) for vertex in vertices)

    def is_face(subset: tuple[int, ...]) -> bool:
        chosen = frozenset(subset)
        return any(chosen <= vertex for vertex in maximal)

    missing: list[tuple[int, ...]] = []
    for order in range(2, max_order + 1):
        for subset in combinations(range(facet_count), order):
            if is_face(subset):
                continue
            if all(is_face(subset[:i] + subset[i + 1 :]) for i in range(order)):
                missing.append(subset)
    return tuple(missing)


@dataclass(frozen=True, slots=True)
class GoldenPairDomain:
    """Finite domain of one off-diagonal field Gram entry."""

    pair: Pair
    kind: str
    values: tuple[Golden, ...]

    @property
    def is_singleton(self) -> bool:
        return len(self.values) == 1

    @property
    def value(self) -> Golden:
        if not self.is_singleton:
            raise ValueError("the pair domain is not resolved")
        return self.values[0]


@dataclass(frozen=True, slots=True)
class GoldenLannerPlacement:
    """One labelled embedding of an order-four/five Lanner block."""

    block_name: str
    facets: tuple[int, ...]
    requirements: tuple[tuple[Pair, Golden], ...]

    @property
    def requirement_map(self) -> dict[Pair, Golden]:
        return dict(self.requirements)


@dataclass(frozen=True, slots=True)
class TenWallSquareClassRelation:
    """Square-class relation between two vertex-plus-wall determinants."""

    facets: Pair
    determinant_ratio: Golden
    square_root: Golden | None

    @property
    def passes(self) -> bool:
        return self.square_root is not None


@dataclass(frozen=True, slots=True)
class VertexDeterminantSquareClassAudit:
    """Early square-class audit for individually resolved outside walls."""

    vertex: tuple[int, ...]
    resolved_determinants: tuple[tuple[int, Golden], ...]
    unresolved_facets: tuple[int, ...]
    relations: tuple[TenWallSquareClassRelation, ...]
    contradiction: str | None

    @property
    def passes(self) -> bool | None:
        if self.contradiction is not None:
            return False
        if self.unresolved_facets:
            return None
        return True


def golden_lanner_placement_options(
    facets: Sequence[int],
) -> tuple[GoldenLannerPlacement, ...]:
    """Enumerate all distinct labelled Lanner placements on a missing face.

    There are 54 placements on a fixed four-set and 240 on a fixed five-set.
    These finite tables are the first CSP branch forced by the missing-face
    combinatorics.
    """

    support = tuple(sorted(facets))
    if len(support) not in (4, 5) or len(support) != len(set(support)):
        raise ValueError("a golden Lanner placement needs four or five facets")
    options: dict[tuple[tuple[Pair, Golden], ...], GoldenLannerPlacement] = {}
    for block in golden_lanner_blocks():
        if block.order != len(support):
            continue
        assert block.name is not None
        for image in permutations(support):
            requirements = tuple(
                sorted(
                    (
                        _pair(image[i], image[j]),
                        block.field_gram[i][j],
                    )
                    for i, j in combinations(range(block.order), 2)
                )
            )
            options.setdefault(
                requirements,
                GoldenLannerPlacement(
                    block_name=block.name,
                    facets=support,
                    requirements=requirements,
                ),
            )
    return tuple(
        sorted(
            options.values(),
            key=lambda option: (
                option.block_name,
                tuple(
                    (pair, value.a, value.b) for pair, value in option.requirements
                ),
            ),
        )
    )


def _inverse_golden(matrix: GoldenMatrix) -> GoldenMatrix:
    size = len(matrix)
    work = [
        [*matrix[i], *(Golden(int(i == j)) for j in range(size))]
        for i in range(size)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]), None
        )
        if pivot is None:
            raise ZeroDivisionError("the vertex Gram matrix is singular")
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


def _bilinear(
    left: Sequence[Golden],
    matrix: Sequence[Sequence[Golden]],
    right: Sequence[Golden],
) -> Golden:
    return sum(
        (
            left[i] * matrix[i][j] * right[j]
            for i in range(len(left))
            for j in range(len(right))
        ),
        start=ZERO,
    )


@dataclass(frozen=True, slots=True)
class VertexFusionConstraint:
    """A fixed elliptic vertex anchor and its adjugate data."""

    vertex: tuple[int, ...]
    complement: tuple[int, ...]
    vertex_gram: GoldenMatrix
    determinant: Golden
    adjugate: GoldenMatrix

    def _resolved_row(
        self,
        outside: int,
        domains: Mapping[Pair, GoldenPairDomain],
    ) -> tuple[Golden, ...] | None:
        entries: list[Golden] = []
        for inside in self.vertex:
            domain = domains[_pair(inside, outside)]
            if not domain.is_singleton:
                return None
            entries.append(domain.value)
        return tuple(entries)

    def determinant_square_class_audit(
        self,
        domains: Mapping[Pair, GoldenPairDomain],
    ) -> VertexDeterminantSquareClassAudit:
        """Audit ten-wall determinant square classes before full fusion.

        The Schur numerator ``D_j`` is exactly
        ``det(G[I union {j}])``.  Rank one implies
        ``D_j / D_k = (R_jk / D_k)**2``, so all resolved nonzero
        determinants must occupy one class in ``F*/F*2``.  This check needs
        no value for the pair joining the two outside walls.
        """

        resolved: list[tuple[int, Golden]] = []
        unresolved: list[int] = []
        for outside in self.complement:
            row = self._resolved_row(outside, domains)
            if row is None:
                unresolved.append(outside)
                continue
            determinant_value = (
                2 * self.determinant - _bilinear(row, self.adjugate, row)
            )
            resolved.append((outside, determinant_value))

        sign_failure = next(
            (
                (outside, value)
                for outside, value in resolved
                if value.sign(0) >= 0 or value.sign(1) <= 0
            ),
            None,
        )
        if sign_failure is not None:
            outside, value = sign_failure
            return VertexDeterminantSquareClassAudit(
                vertex=self.vertex,
                resolved_determinants=tuple(resolved),
                unresolved_facets=tuple(unresolved),
                relations=(),
                contradiction=(
                    f"outside facet {outside} has ten-wall determinant {value}; "
                    "expected signs (-,+)"
                ),
            )

        relations: list[TenWallSquareClassRelation] = []
        for (left, left_value), (right, right_value) in combinations(resolved, 2):
            ratio = left_value / right_value
            witness = golden_square_root(ratio)
            relation = TenWallSquareClassRelation(
                facets=_pair(left, right),
                determinant_ratio=ratio,
                square_root=witness,
            )
            relations.append(relation)
            if not relation.passes:
                return VertexDeterminantSquareClassAudit(
                    vertex=self.vertex,
                    resolved_determinants=tuple(resolved),
                    unresolved_facets=tuple(unresolved),
                    relations=tuple(relations),
                    contradiction=(
                        f"ten-wall determinants for outside facets {left} and "
                        f"{right} have nonsquare ratio {ratio}"
                    ),
                )
        return VertexDeterminantSquareClassAudit(
            vertex=self.vertex,
            resolved_determinants=tuple(resolved),
            unresolved_facets=tuple(unresolved),
            relations=tuple(relations),
            contradiction=None,
        )

    def propagate(
        self,
        domains: Mapping[Pair, GoldenPairDomain],
    ) -> VertexFusionPropagation:
        """Apply the oriented rank-one rule when all outside rows resolve.

        The unsquared condition ``R_jk(identity) < 0`` is indispensable.  It
        chooses the common time orientation of all outside roots.  In the
        equivalent vertex-plus-anchor formulation, the coefficient of wall
        ``k`` along anchor ``j`` is ``R_jk / D_j``; the two negative signs
        below are therefore exactly the existing positive same-side test.
        """

        rows: dict[int, tuple[Golden, ...]] = {}
        unresolved: set[Pair] = set()
        for outside in self.complement:
            row = self._resolved_row(outside, domains)
            if row is None:
                unresolved.update(
                    _pair(inside, outside)
                    for inside in self.vertex
                    if not domains[_pair(inside, outside)].is_singleton
                )
            else:
                rows[outside] = row
        square_class_audit = self.determinant_square_class_audit(domains)
        if square_class_audit.contradiction is not None:
            return VertexFusionPropagation(
                vertex=self.vertex,
                ready=not unresolved,
                contradiction=square_class_audit.contradiction,
                unresolved_pairs=tuple(sorted(unresolved)),
                diagonal_numerators=square_class_audit.resolved_determinants,
                determinant_square_class_relations=square_class_audit.relations,
                residual_numerators=(),
                filtered_domains=(),
            )
        if unresolved:
            return VertexFusionPropagation(
                vertex=self.vertex,
                ready=False,
                contradiction=None,
                unresolved_pairs=tuple(sorted(unresolved)),
                diagonal_numerators=square_class_audit.resolved_determinants,
                determinant_square_class_relations=square_class_audit.relations,
                residual_numerators=(),
                filtered_domains=(),
            )

        diagonal = dict(square_class_audit.resolved_determinants)

        residuals: list[tuple[Pair, tuple[Golden, ...]]] = []
        filtered: list[tuple[Pair, tuple[Golden, ...]]] = []
        for left, right in combinations(self.complement, 2):
            pair = _pair(left, right)
            prediction = _bilinear(rows[left], self.adjugate, rows[right])
            allowed: list[Golden] = []
            pair_residuals: list[Golden] = []
            for value in domains[pair].values:
                residual = self.determinant * value - prediction
                pair_residuals.append(residual)
                if (
                    residual * residual == diagonal[left] * diagonal[right]
                    and residual.sign(0) < 0
                ):
                    allowed.append(value)
            residuals.append((pair, tuple(pair_residuals)))
            if not allowed:
                return VertexFusionPropagation(
                    vertex=self.vertex,
                    ready=True,
                    contradiction=(
                        f"pair {pair} has no value satisfying the oriented "
                        "rank-one Schur equation"
                    ),
                    unresolved_pairs=(),
                    diagonal_numerators=tuple(sorted(diagonal.items())),
                    determinant_square_class_relations=(
                        square_class_audit.relations
                    ),
                    residual_numerators=tuple(residuals),
                    filtered_domains=tuple(filtered),
                )
            filtered.append((pair, tuple(allowed)))
        return VertexFusionPropagation(
            vertex=self.vertex,
            ready=True,
            contradiction=None,
            unresolved_pairs=(),
            diagonal_numerators=tuple(sorted(diagonal.items())),
            determinant_square_class_relations=square_class_audit.relations,
            residual_numerators=tuple(residuals),
            filtered_domains=tuple(filtered),
        )

    def symbolic_system(
        self,
        domains: Mapping[Pair, GoldenPairDomain],
    ) -> VertexFusionPolynomialSystem:
        """Compile the denominator-free two-component polynomial system."""

        variables: dict[Pair, PhiDomainVariable] = {}

        def entry(pair: Pair) -> PhiExpression:
            domain = domains[pair]
            if domain.is_singleton:
                return PhiExpression.constant(domain.value)
            variable = variables.setdefault(pair, PhiDomainVariable.from_domain(domain))
            return PhiExpression(variable.a_symbol, variable.b_symbol)

        rows: dict[int, tuple[PhiExpression, ...]] = {
            outside: tuple(entry(_pair(inside, outside)) for inside in self.vertex)
            for outside in self.complement
        }
        delta = PhiExpression.constant(self.determinant)
        adjugate = tuple(
            tuple(PhiExpression.constant(value) for value in row)
            for row in self.adjugate
        )

        def symbolic_bilinear(
            left: Sequence[PhiExpression],
            right: Sequence[PhiExpression],
        ) -> PhiExpression:
            result = PhiExpression.constant(ZERO)
            for i in range(len(left)):
                for j in range(len(right)):
                    result += left[i] * adjugate[i][j] * right[j]
            return result.expanded()

        diagonal = {
            outside: (
                PhiExpression.constant(Golden(2)) * delta
                - symbolic_bilinear(row, row)
            ).expanded()
            for outside, row in rows.items()
        }
        equations: list[PhiPolynomialEquation] = []
        signs: list[PhiPolynomialSign] = []
        for outside, expression in diagonal.items():
            signs.append(
                PhiPolynomialSign(
                    name=f"v{self.vertex}_D_{outside}_identity_negative",
                    expression=expression,
                    embedding=0,
                    expected_sign=-1,
                )
            )
            signs.append(
                PhiPolynomialSign(
                    name=f"v{self.vertex}_D_{outside}_conjugate_positive",
                    expression=expression,
                    embedding=1,
                    expected_sign=1,
                )
            )
        for left, right in combinations(self.complement, 2):
            pair = _pair(left, right)
            residual = (
                delta * entry(pair) - symbolic_bilinear(rows[left], rows[right])
            ).expanded()
            equation = (
                residual * residual - diagonal[left] * diagonal[right]
            ).expanded()
            equations.append(
                PhiPolynomialEquation(
                    name=f"v{self.vertex}_rank1_{left}_{right}",
                    expression=equation,
                )
            )
            signs.append(
                PhiPolynomialSign(
                    name=f"v{self.vertex}_R_{left}_{right}_identity_negative",
                    expression=residual,
                    embedding=0,
                    expected_sign=-1,
                )
            )
        return VertexFusionPolynomialSystem(
            vertex=self.vertex,
            complement=self.complement,
            determinant=self.determinant,
            variables=tuple(variables[pair] for pair in sorted(variables)),
            equations=tuple(equations),
            signs=tuple(signs),
        )


@dataclass(frozen=True, slots=True)
class VertexFusionPropagation:
    vertex: tuple[int, ...]
    ready: bool
    contradiction: str | None
    unresolved_pairs: tuple[Pair, ...]
    diagonal_numerators: tuple[tuple[int, Golden], ...]
    determinant_square_class_relations: tuple[TenWallSquareClassRelation, ...]
    residual_numerators: tuple[tuple[Pair, tuple[Golden, ...]], ...]
    filtered_domains: tuple[tuple[Pair, tuple[Golden, ...]], ...]


def _sympy_rational(value: Fraction) -> sp.Rational:
    return sp.Rational(value.numerator, value.denominator)


@dataclass(frozen=True, slots=True)
class PhiExpression:
    """A symbolic element ``a+b*phi`` with rational-polynomial components."""

    a: sp.Expr
    b: sp.Expr

    @classmethod
    def constant(cls, value: Golden | RationalLike) -> PhiExpression:
        golden = Golden.coerce(value)
        return cls(_sympy_rational(golden.a), _sympy_rational(golden.b))

    def __add__(self, other: PhiExpression) -> PhiExpression:
        return PhiExpression(self.a + other.a, self.b + other.b)

    def __sub__(self, other: PhiExpression) -> PhiExpression:
        return PhiExpression(self.a - other.a, self.b - other.b)

    def __mul__(self, other: PhiExpression) -> PhiExpression:
        return PhiExpression(
            self.a * other.a + self.b * other.b,
            self.a * other.b + self.b * other.a + self.b * other.b,
        )

    def expanded(self) -> PhiExpression:
        return PhiExpression(sp.expand(self.a), sp.expand(self.b))

    def evaluate(self, substitutions: Mapping[sp.Symbol, sp.Rational]) -> Golden:
        a_value = sp.expand(self.a.subs(substitutions))
        b_value = sp.expand(self.b.subs(substitutions))
        if a_value.is_Rational is not True or b_value.is_Rational is not True:
            raise ValueError("a symbolic golden expression remains unresolved")
        return Golden(
            Fraction(int(a_value.p), int(a_value.q)),
            Fraction(int(b_value.p), int(b_value.q)),
        )


@dataclass(frozen=True, slots=True)
class PhiDomainVariable:
    pair: Pair
    a_symbol: sp.Symbol
    b_symbol: sp.Symbol
    allowed_coefficients: tuple[tuple[Fraction, Fraction], ...]

    @classmethod
    def from_domain(cls, domain: GoldenPairDomain) -> PhiDomainVariable:
        left, right = domain.pair
        return cls(
            pair=domain.pair,
            a_symbol=sp.Symbol(f"g_{left}_{right}_a"),
            b_symbol=sp.Symbol(f"g_{left}_{right}_b"),
            allowed_coefficients=tuple((value.a, value.b) for value in domain.values),
        )


@dataclass(frozen=True, slots=True)
class PhiPolynomialEquation:
    name: str
    expression: PhiExpression


@dataclass(frozen=True, slots=True)
class PhiPolynomialSign:
    name: str
    expression: PhiExpression
    embedding: int
    expected_sign: int


@dataclass(frozen=True, slots=True)
class FusionSymbolicEvaluation:
    equations_zero: bool
    signs_correct: bool

    @property
    def passes(self) -> bool:
        return self.equations_zero and self.signs_correct


@dataclass(frozen=True, slots=True)
class VertexFusionPolynomialSystem:
    """Two rational equations for every golden rank-one equation."""

    vertex: tuple[int, ...]
    complement: tuple[int, ...]
    determinant: Golden
    variables: tuple[PhiDomainVariable, ...]
    equations: tuple[PhiPolynomialEquation, ...]
    signs: tuple[PhiPolynomialSign, ...]

    @property
    def rational_equation_count(self) -> int:
        return 2 * len(self.equations)

    def evaluate(self, assignment: Mapping[Pair, Golden | RationalLike]) -> FusionSymbolicEvaluation:
        substitutions: dict[sp.Symbol, sp.Rational] = {}
        for variable in self.variables:
            if variable.pair not in assignment:
                raise KeyError(f"missing value for pair {variable.pair}")
            value = Golden.coerce(assignment[variable.pair])
            if (value.a, value.b) not in variable.allowed_coefficients:
                raise ValueError(f"value {value} is outside the domain of {variable.pair}")
            substitutions[variable.a_symbol] = _sympy_rational(value.a)
            substitutions[variable.b_symbol] = _sympy_rational(value.b)
        equations_zero = all(
            equation.expression.evaluate(substitutions) == ZERO
            for equation in self.equations
        )
        signs_correct = all(
            sign.expression.evaluate(substitutions).sign(sign.embedding)
            == sign.expected_sign
            for sign in self.signs
        )
        return FusionSymbolicEvaluation(equations_zero, signs_correct)


@dataclass(frozen=True, slots=True)
class LannerOverlapAudit:
    placements: tuple[str, str]
    support: tuple[int, ...]
    resolved: bool
    identity_inertia: Inertia | None
    conjugate_inertia: Inertia | None

    @property
    def passes(self) -> bool | None:
        if not self.resolved:
            return None
        assert self.identity_inertia is not None
        assert self.conjugate_inertia is not None
        return (
            self.identity_inertia.negative <= 1
            and self.conjugate_inertia.negative == 0
        )


def audit_lanner_overlaps(
    placements: Sequence[GoldenLannerPlacement],
    domains: Mapping[Pair, GoldenPairDomain],
) -> tuple[LannerOverlapAudit, ...]:
    """Apply principal-submatrix inertia to every selected seed pair."""

    audits: list[LannerOverlapAudit] = []
    for first, second in combinations(placements, 2):
        support = tuple(sorted(set(first.facets) | set(second.facets)))
        pairs = tuple(combinations(support, 2))
        if any(not domains[_pair(*pair)].is_singleton for pair in pairs):
            audits.append(
                LannerOverlapAudit(
                    placements=(first.block_name, second.block_name),
                    support=support,
                    resolved=False,
                    identity_inertia=None,
                    conjugate_inertia=None,
                )
            )
            continue
        gram = tuple(
            tuple(
                Golden(2)
                if i == j
                else domains[_pair(support[i], support[j])].value
                for j in range(len(support))
            )
            for i in range(len(support))
        )
        audits.append(
            LannerOverlapAudit(
                placements=(first.block_name, second.block_name),
                support=support,
                resolved=True,
                identity_inertia=inertia(gram, embedding=0),
                conjugate_inertia=inertia(gram, embedding=1),
            )
        )
    return tuple(audits)


@dataclass(frozen=True, slots=True)
class GoldenFusionPropagationResult:
    domains: tuple[GoldenPairDomain, ...]
    iterations: int
    pruned_value_count: int
    compiled_vertex_count: int
    ready_vertex_count: int
    unresolved_vertices: tuple[tuple[int, ...], ...]
    contradiction: str | None

    @property
    def passes(self) -> bool:
        return self.contradiction is None and not self.unresolved_vertices

    @property
    def domain_map(self) -> dict[Pair, GoldenPairDomain]:
        return {domain.pair: domain for domain in self.domains}


@dataclass(frozen=True, slots=True)
class GoldenFusionProgram:
    """Compiled finite-domain global fusion problem for one nerve."""

    facet_count: int
    dimension: int
    vertices: tuple[tuple[int, ...], ...]
    pair_domains: tuple[GoldenPairDomain, ...]
    lanner_placements: tuple[GoldenLannerPlacement, ...]
    initial_contradictions: tuple[str, ...]

    @property
    def domain_map(self) -> dict[Pair, GoldenPairDomain]:
        return {domain.pair: domain for domain in self.pair_domains}

    @property
    def missing_faces(self) -> tuple[tuple[int, ...], ...]:
        return minimal_missing_faces(
            self.vertices, facet_count=self.facet_count, max_order=5
        )

    def _constraint(
        self,
        vertex: tuple[int, ...],
        domains: Mapping[Pair, GoldenPairDomain],
    ) -> VertexFusionConstraint | None:
        internal = tuple(combinations(vertex, 2))
        if any(not domains[_pair(*pair)].is_singleton for pair in internal):
            return None
        gram = tuple(
            tuple(
                Golden(2)
                if i == j
                else domains[_pair(vertex[i], vertex[j])].value
                for j in range(self.dimension)
            )
            for i in range(self.dimension)
        )
        if inertia(gram, embedding=0) != Inertia(self.dimension, 0, 0):
            raise ValueError(f"vertex {vertex} is not elliptic at the identity place")
        if inertia(gram, embedding=1) != Inertia(self.dimension, 0, 0):
            raise ValueError(f"vertex {vertex} is not elliptic at the conjugate place")
        delta = determinant(gram)
        inverse = _inverse_golden(gram)
        adjugate = tuple(
            tuple(delta * inverse[i][j] for j in range(self.dimension))
            for i in range(self.dimension)
        )
        complement = tuple(
            index for index in range(self.facet_count) if index not in vertex
        )
        return VertexFusionConstraint(vertex, complement, gram, delta, adjugate)

    def constraints(
        self,
        domains: Mapping[Pair, GoldenPairDomain] | None = None,
    ) -> tuple[VertexFusionConstraint, ...]:
        current = self.domain_map if domains is None else domains
        constraints: list[VertexFusionConstraint] = []
        for vertex in self.vertices:
            constraint = self._constraint(vertex, current)
            if constraint is not None:
                constraints.append(constraint)
        return tuple(constraints)

    def symbolic_systems(self) -> tuple[VertexFusionPolynomialSystem, ...]:
        domains = self.domain_map
        return tuple(
            constraint.symbolic_system(domains)
            for constraint in self.constraints(domains)
        )

    def propagate(self) -> GoldenFusionPropagationResult:
        domains = self.domain_map
        if self.initial_contradictions:
            return GoldenFusionPropagationResult(
                domains=tuple(domains[pair] for pair in sorted(domains)),
                iterations=0,
                pruned_value_count=0,
                compiled_vertex_count=0,
                ready_vertex_count=0,
                unresolved_vertices=self.vertices,
                contradiction=self.initial_contradictions[0],
            )

        pruned = 0
        iterations = 0
        compiled_count = 0
        ready_count = 0
        while True:
            iterations += 1
            changed = False
            compiled_count = 0
            ready_count = 0
            try:
                constraints = self.constraints(domains)
            except ValueError as error:
                return GoldenFusionPropagationResult(
                    domains=tuple(domains[pair] for pair in sorted(domains)),
                    iterations=iterations,
                    pruned_value_count=pruned,
                    compiled_vertex_count=compiled_count,
                    ready_vertex_count=ready_count,
                    unresolved_vertices=(),
                    contradiction=str(error),
                )
            compiled_count = len(constraints)
            for constraint in constraints:
                result = constraint.propagate(domains)
                if result.ready:
                    ready_count += 1
                if result.contradiction is not None:
                    return GoldenFusionPropagationResult(
                        domains=tuple(domains[pair] for pair in sorted(domains)),
                        iterations=iterations,
                        pruned_value_count=pruned,
                        compiled_vertex_count=compiled_count,
                        ready_vertex_count=ready_count,
                        unresolved_vertices=(),
                        contradiction=f"vertex {constraint.vertex}: {result.contradiction}",
                    )
                if not result.ready:
                    continue
                for pair, allowed in result.filtered_domains:
                    old = domains[pair]
                    new_values = tuple(value for value in old.values if value in allowed)
                    if new_values != old.values:
                        pruned += len(old.values) - len(new_values)
                        domains[pair] = GoldenPairDomain(pair, old.kind, new_values)
                        changed = True

            overlaps = audit_lanner_overlaps(self.lanner_placements, domains)
            failed_overlap = next(
                (audit for audit in overlaps if audit.passes is False), None
            )
            if failed_overlap is not None:
                return GoldenFusionPropagationResult(
                    domains=tuple(domains[pair] for pair in sorted(domains)),
                    iterations=iterations,
                    pruned_value_count=pruned,
                    compiled_vertex_count=compiled_count,
                    ready_vertex_count=ready_count,
                    unresolved_vertices=(),
                    contradiction=(
                        f"Lanner placements {failed_overlap.placements} on "
                        f"{failed_overlap.support} have more than one identity "
                        "negative direction or a conjugate negative direction"
                    ),
                )
            if not changed:
                break

        unresolved = tuple(
            vertex
            for vertex in self.vertices
            if (
                (constraint := self._constraint(vertex, domains)) is None
                or not constraint.propagate(domains).ready
            )
        )
        return GoldenFusionPropagationResult(
            domains=tuple(domains[pair] for pair in sorted(domains)),
            iterations=iterations,
            pruned_value_count=pruned,
            compiled_vertex_count=compiled_count,
            ready_vertex_count=ready_count,
            unresolved_vertices=unresolved,
            contradiction=None,
        )


def _normalize_mapping_keys(
    mapping: Mapping[Pair, object],
) -> dict[Pair, object]:
    normalized: dict[Pair, object] = {}
    for raw_pair, value in mapping.items():
        if len(raw_pair) != 2:
            raise ValueError("a Gram-domain key must be a facet pair")
        pair = _pair(*raw_pair)
        if pair in normalized:
            raise ValueError(f"pair {pair} was supplied twice")
        normalized[pair] = value
    return normalized


def compile_golden_fusion_program(
    vertices: Sequence[Sequence[int]],
    finite_labels: Mapping[Pair, int | Sequence[int]],
    dotted_products: Mapping[
        Pair,
        Golden | RationalLike | Sequence[Golden | RationalLike],
    ],
    *,
    facet_count: int = 14,
    dimension: int = 9,
    lanner_placements: Sequence[GoldenLannerPlacement] = (),
    require_effective_trace: bool = True,
) -> GoldenFusionProgram:
    """Compile a nerve and finite Gram domains into fusion constraints."""

    normalized_vertices = _normalize_vertices(
        vertices, facet_count=facet_count, dimension=dimension
    )
    intersecting = _intersecting_pairs(normalized_vertices)
    all_pairs = frozenset(combinations(range(facet_count), 2))
    disjoint = all_pairs - intersecting
    finite_input = _normalize_mapping_keys(finite_labels)
    dotted_input = _normalize_mapping_keys(dotted_products)
    if set(finite_input) != set(intersecting):
        missing = set(intersecting) - set(finite_input)
        extra = set(finite_input) - set(intersecting)
        raise ValueError(f"finite-label pairs mismatch the nerve: missing={missing}, extra={extra}")
    if set(dotted_input) != set(disjoint):
        missing = set(disjoint) - set(dotted_input)
        extra = set(dotted_input) - set(disjoint)
        raise ValueError(f"dotted pairs mismatch the nerve: missing={missing}, extra={extra}")

    products = {2: ZERO, 3: Golden(-1), 5: -PHI}
    domains: dict[Pair, GoldenPairDomain] = {}
    for pair in sorted(intersecting):
        raw = finite_input[pair]
        labels = (raw,) if isinstance(raw, int) else tuple(raw)
        if not labels or any(label not in products for label in labels):
            raise ValueError(f"finite labels for {pair} must be a nonempty subset of 2,3,5")
        values = tuple(
            sorted({products[label] for label in labels}, key=_golden_sort_key)
        )
        domains[pair] = GoldenPairDomain(pair, "finite", values)
    for pair in sorted(disjoint):
        values = _coerce_domain(dotted_input[pair])
        for value in values:
            if (value + 2).sign(0) >= 0:
                raise ValueError(f"dotted product {value} at {pair} is not below -2")
            if (Golden(4) - value.conjugate() * value.conjugate()).sign(0) <= 0:
                raise ValueError(
                    f"dotted product {value} at {pair} violates conjugate positivity"
                )
            trace_value = k3_trace_square(value)
            if require_effective_trace and (
                trace_value.denominator != 1 or trace_value < 0
            ):
                raise ValueError(
                    f"dotted product {value} at {pair} has ineffective trace {trace_value}"
                )
        domains[pair] = GoldenPairDomain(pair, "dotted", values)

    contradictions: list[str] = []
    placements = tuple(lanner_placements)
    for placement in placements:
        if any(index < 0 or index >= facet_count for index in placement.facets):
            raise IndexError("a Lanner placement facet is out of range")
        for pair, required in placement.requirements:
            domain = domains[pair]
            if domain.kind != "finite":
                contradictions.append(
                    f"Lanner placement {placement.block_name} uses disjoint pair {pair}"
                )
                continue
            retained = tuple(value for value in domain.values if value == required)
            if not retained:
                contradictions.append(
                    f"Lanner placement {placement.block_name} conflicts at pair {pair}"
                )
            else:
                domains[pair] = GoldenPairDomain(pair, domain.kind, retained)

    return GoldenFusionProgram(
        facet_count=facet_count,
        dimension=dimension,
        vertices=normalized_vertices,
        pair_domains=tuple(domains[pair] for pair in sorted(domains)),
        lanner_placements=placements,
        initial_contradictions=tuple(contradictions),
    )


def fusion_program_from_field_gram(
    field_gram: Sequence[Sequence[Golden | RationalLike]],
    vertices: Sequence[Sequence[int]],
    *,
    pair_domain_overrides: Mapping[
        Pair,
        Golden | RationalLike | Sequence[Golden | RationalLike],
    ]
    | None = None,
    lanner_placements: Sequence[GoldenLannerPlacement] = (),
) -> GoldenFusionProgram:
    """Compile a completed field Gram matrix, optionally widening domains."""

    gram = tuple(tuple(Golden.coerce(entry) for entry in row) for row in field_gram)
    facet_count = len(gram)
    if any(len(row) != facet_count for row in gram):
        raise ValueError("the field Gram matrix must be square")
    if any(gram[index][index] != Golden(2) for index in range(facet_count)):
        raise ValueError("every field Gram diagonal entry must equal 2")
    if any(gram[i][j] != gram[j][i] for i in range(facet_count) for j in range(i)):
        raise ValueError("the field Gram matrix must be symmetric")
    raw_vertices = tuple(tuple(vertex) for vertex in vertices)
    if not raw_vertices:
        raise ValueError("the nerve must contain at least one maximal vertex")
    dimension = len(raw_vertices[0])
    normalized_vertices = _normalize_vertices(
        raw_vertices, facet_count=facet_count, dimension=dimension
    )
    intersecting = _intersecting_pairs(normalized_vertices)
    raw_overrides = _normalize_mapping_keys(pair_domain_overrides or {})
    all_pairs = set(combinations(range(facet_count), 2))
    if not set(raw_overrides) <= all_pairs:
        raise IndexError("a pair-domain override facet is out of range")
    overrides = {
        pair: _coerce_domain(value)
        for pair, value in raw_overrides.items()
    }
    inverse_labels = {ZERO: 2, Golden(-1): 3, -PHI: 5}
    finite: dict[Pair, int | tuple[int, ...]] = {}
    dotted: dict[Pair, Golden | tuple[Golden, ...]] = {}
    for pair in combinations(range(facet_count), 2):
        values = overrides.get(pair, (gram[pair[0]][pair[1]],))
        if pair in intersecting:
            try:
                labels = tuple(inverse_labels[value] for value in values)
            except KeyError as error:
                raise ValueError(
                    f"intersecting pair {pair} has an unrecognized product"
                ) from error
            finite[pair] = labels[0] if len(labels) == 1 else labels
        else:
            dotted[pair] = values[0] if len(values) == 1 else values
    return compile_golden_fusion_program(
        normalized_vertices,
        finite,
        dotted,
        facet_count=facet_count,
        dimension=dimension,
        lanner_placements=lanner_placements,
    )

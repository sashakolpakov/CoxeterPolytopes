"""K3 audits for the eight order-four/five golden Lanner seeds.

These are the five rank-four and three rank-five members relevant to missing
faces in the H9 fusion search.  The exhaustive equal-root class additionally
contains four rank-three triangles.

The trace lift of a golden Lanner block is an even Lorentzian lattice of
rank eight or ten.  This module separates three questions that are easy to
conflate:

* whether that seed embeds primitively in the K3 lattice;
* whether it has a proper even overlattice inside its own rational span; and
* whether it occurs primitively inside a rank-20, golden-action-stable
  Neron--Severi lattice.

The first question is decided by Nikulin's stable-range criterion.  The
second is an exact enumeration of isotropic subgroups of the seed
discriminant form, including invariance under multiplication by ``phi``.
For the third, explicit embeddings into the archived discriminant-11 and
discriminant-19 carriers provide integral gluing certificates.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import lcm
from typing import Iterable, Sequence

from sympy import Matrix as SympyMatrix
from sympy import ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.exact.matrix import Inertia, inertia
from coxeter_polytopes.exact.trace_lattice import integral_gram, trace_gram
from coxeter_polytopes.search.h11_vinberg import (
    h11_basis_gram,
    run_h11_vinberg_search,
)
from coxeter_polytopes.search.h19_vinberg import (
    h19_basis_gram,
    run_h19_vinberg_search,
)

from .golden_lanner import GoldenLannerAudit, golden_lanner_blocks
from .k3 import K3_TRACE_SCALAR


IntegerMatrix = tuple[tuple[int, ...], ...]
DiscriminantElement = tuple[int, ...]
GluingElement = tuple[DiscriminantElement, DiscriminantElement]


def _domain_matrix(matrix: SympyMatrix) -> DomainMatrix:
    return DomainMatrix(
        [
            [ZZ(int(matrix[i, j])) for j in range(matrix.cols)]
            for i in range(matrix.rows)
        ],
        matrix.shape,
        ZZ,
    )


def _matrix_tuple(matrix: SympyMatrix) -> IntegerMatrix:
    result: list[tuple[int, ...]] = []
    for i in range(matrix.rows):
        row: list[int] = []
        for j in range(matrix.cols):
            entry = matrix[i, j]
            if entry.is_Integer is not True:
                raise ArithmeticError("an expected integral matrix is nonintegral")
            row.append(int(entry))
        result.append(tuple(row))
    return tuple(result)


def _fraction(value: object) -> Fraction:
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    numerator, denominator = value.as_numer_denom()  # type: ignore[attr-defined]
    return Fraction(int(numerator), int(denominator))


def _mod_two(value: object) -> Fraction:
    rational = _fraction(value)
    return rational - 2 * (rational // 2)


def _smith_invariants(matrix: SympyMatrix) -> tuple[int, ...]:
    normal = smith_normal_form(matrix, domain=ZZ)
    return tuple(
        abs(int(normal[i, i]))
        for i in range(min(normal.shape))
        if abs(int(normal[i, i])) > 1
    )


def _interleaved_golden_action(rank: int) -> SympyMatrix:
    action = SympyMatrix.zeros(2 * rank)
    for index in range(rank):
        first = 2 * index
        action[first, first + 1] = 1
        action[first + 1, first] = 1
        action[first + 1, first + 1] = 1
    return action


class _DiscriminantModel:
    """Canonical Smith coordinates for one finite discriminant form."""

    def __init__(self, gram: SympyMatrix, action: SympyMatrix):
        if gram.det() == 0:
            raise ValueError("a discriminant form requires a nondegenerate lattice")
        normal, left, _ = smith_normal_decomp(_domain_matrix(gram))
        diagonal = normal.to_Matrix()
        self.gram = gram
        self.gram_inverse = gram.inv()
        self.left = left.to_Matrix()
        self.left_inverse = self.left.inv()
        self.diagonal = tuple(
            abs(int(diagonal[i, i])) for i in range(gram.rows)
        )
        self.active_indices = tuple(
            i for i, modulus in enumerate(self.diagonal) if modulus > 1
        )
        self.moduli = tuple(self.diagonal[i] for i in self.active_indices)
        self.action = action
        if action.T * gram != gram * action:
            raise ValueError("the proposed discriminant action is not self-adjoint")

    @property
    def zero(self) -> DiscriminantElement:
        return tuple(0 for _ in self.moduli)

    @property
    def elements(self) -> tuple[DiscriminantElement, ...]:
        return tuple(product(*(range(modulus) for modulus in self.moduli)))

    def _full_smith_vector(self, element: DiscriminantElement) -> SympyMatrix:
        if len(element) != len(self.moduli):
            raise ValueError("a discriminant element has the wrong length")
        vector = SympyMatrix.zeros(self.gram.rows, 1)
        for index, value, modulus in zip(
            self.active_indices, element, self.moduli
        ):
            vector[index] = value % modulus
        return vector

    def covector(self, element: DiscriminantElement) -> SympyMatrix:
        """Return ``z`` for the dual vector ``gram**-1 z``."""

        return self.left_inverse * self._full_smith_vector(element)

    def dual_vector(self, element: DiscriminantElement) -> SympyMatrix:
        return self.gram_inverse * self.covector(element)

    def canonical_from_lattice_coordinates(
        self, vector: SympyMatrix
    ) -> DiscriminantElement:
        smith = self.left * self.gram * vector
        result: list[int] = []
        for index, modulus in zip(self.active_indices, self.moduli):
            entry = smith[index]
            if entry.is_Integer is not True:
                raise ArithmeticError("a purported dual vector has nonintegral pairings")
            result.append(int(entry) % modulus)
        return tuple(result)

    def quadratic(self, element: DiscriminantElement) -> Fraction:
        covector = self.covector(element)
        return _mod_two((covector.T * self.gram_inverse * covector)[0])

    def apply_action(self, element: DiscriminantElement) -> DiscriminantElement:
        transformed_covector = self.action.T * self.covector(element)
        smith = self.left * transformed_covector
        return tuple(
            int(smith[index]) % modulus
            for index, modulus in zip(self.active_indices, self.moduli)
        )

    def add(
        self, left: DiscriminantElement, right: DiscriminantElement
    ) -> DiscriminantElement:
        return tuple(
            (a + b) % modulus
            for a, b, modulus in zip(left, right, self.moduli)
        )

    def span(
        self, generators: Iterable[DiscriminantElement]
    ) -> frozenset[DiscriminantElement]:
        generators = tuple(generators)
        subgroup = {self.zero}
        frontier = [self.zero]
        while frontier:
            current = frontier.pop()
            for generator in generators:
                value = self.add(current, generator)
                if value not in subgroup:
                    subgroup.add(value)
                    frontier.append(value)
        return frozenset(subgroup)


def _isotropic_subgroups(
    model: _DiscriminantModel,
) -> tuple[frozenset[DiscriminantElement], ...]:
    isotropic_elements = tuple(
        element for element in model.elements if model.quadratic(element) == 0
    )
    subgroups = {frozenset({model.zero})}
    frontier = list(subgroups)
    while frontier:
        subgroup = frontier.pop()
        for element in isotropic_elements:
            if element in subgroup:
                continue
            extension = model.span((*subgroup, element))
            if all(model.quadratic(value) == 0 for value in extension) and (
                extension not in subgroups
            ):
                subgroups.add(extension)
                frontier.append(extension)
    return tuple(sorted(subgroups, key=lambda group: (len(group), sorted(group))))


def _subgroup_generators(
    model: _DiscriminantModel,
    subgroup: frozenset[DiscriminantElement],
) -> tuple[DiscriminantElement, ...]:
    generators: list[DiscriminantElement] = []
    generated = model.span(())
    for element in sorted(subgroup):
        if element not in generated:
            generators.append(element)
            generated = model.span(generators)
    return tuple(generators)


def _overlattice_gram(
    model: _DiscriminantModel,
    subgroup: frozenset[DiscriminantElement],
) -> SympyMatrix:
    generators = _subgroup_generators(model, subgroup)
    dual_vectors = tuple(model.dual_vector(element) for element in generators)
    denominator = lcm(
        1,
        *(
            _fraction(entry).denominator
            for vector in dual_vectors
            for entry in vector
        ),
    )
    columns = tuple(
        SympyMatrix.eye(model.gram.rows).col(i)
        for i in range(model.gram.rows)
    )
    columns += dual_vectors
    scaled = SympyMatrix.hstack(
        *(column * denominator for column in columns)
    ).applyfunc(int)
    basis = hermite_normal_form(scaled) / denominator
    gram = basis.T * model.gram * basis
    if any(entry.is_Integer is not True for entry in gram):
        raise ArithmeticError("an isotropic subgroup produced a nonintegral overlattice")
    if any(int(gram[i, i]) % 2 for i in range(gram.rows)):
        raise ArithmeticError("an isotropic subgroup produced an odd overlattice")
    return gram


@dataclass(frozen=True, slots=True)
class EvenSeedOverlattice:
    """One exact even overlattice of a Lanner trace seed."""

    subgroup: tuple[DiscriminantElement, ...]
    index: int
    tau_stable: bool
    gram: IntegerMatrix
    determinant: int
    discriminant_group_invariants: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class GoldenLannerK3SeedAudit:
    """K3 and golden-saturation constraints for one Lanner trace seed."""

    name: str
    order: int
    trace_gram: IntegerMatrix
    trace_inertia: Inertia
    discriminant_moduli: tuple[int, ...]
    discriminant_quadratic_distribution: tuple[tuple[Fraction, int], ...]
    tau_basis_images: tuple[DiscriminantElement, ...]
    even_overlattices: tuple[EvenSeedOverlattice, ...]
    k3_complement_signature: Inertia
    k3_complement_rank: int
    nikulin_length_bound: int
    primitive_k3_embedding_exists: bool
    primitive_k3_embedding_unique: bool

    @property
    def discriminant_length(self) -> int:
        return len(self.discriminant_moduli)

    @property
    def nontrivial_tau_stable_overlattices(self) -> tuple[EvenSeedOverlattice, ...]:
        return tuple(
            overlattice
            for overlattice in self.even_overlattices
            if overlattice.index > 1 and overlattice.tau_stable
        )

    @property
    def primitive_in_every_tau_stable_completion(self) -> bool:
        return not self.nontrivial_tau_stable_overlattices

    @property
    def passes(self) -> bool:
        return (
            self.primitive_k3_embedding_exists
            and self.primitive_k3_embedding_unique
            and self.primitive_in_every_tau_stable_completion
        )


def audit_golden_lanner_k3_seed(
    block: GoldenLannerAudit,
) -> GoldenLannerK3SeedAudit:
    """Audit primitive embedding and every even self-overlattice of ``block``."""

    if block.name is None or not block.is_arithmetic_lanner:
        raise ValueError("the K3 seed audit requires a classified golden Lanner block")
    gram = SympyMatrix(block.restriction_trace_gram)
    action = _interleaved_golden_action(block.order)
    model = _DiscriminantModel(gram, action)
    overlattices: list[EvenSeedOverlattice] = []
    for subgroup in _isotropic_subgroups(model):
        overlattice_gram = _overlattice_gram(model, subgroup)
        index = len(subgroup)
        if abs(int(gram.det())) != abs(int(overlattice_gram.det())) * index * index:
            raise ArithmeticError("the overlattice determinant-index identity failed")
        overlattices.append(
            EvenSeedOverlattice(
                subgroup=tuple(sorted(subgroup)),
                index=index,
                tau_stable={model.apply_action(value) for value in subgroup}
                == set(subgroup),
                gram=_matrix_tuple(overlattice_gram),
                determinant=int(overlattice_gram.det()),
                discriminant_group_invariants=_smith_invariants(overlattice_gram),
            )
        )

    trace_rank = 2 * block.order
    complement_rank = 22 - trace_rank
    complement_signature = Inertia(2, 19 - (trace_rank - 1), 0)
    length_bound = complement_rank - 2
    strict_signature_room = (
        complement_signature.positive > 0 and complement_signature.negative > 0
    )
    stable_range = len(model.moduli) <= length_bound
    quadratic_distribution = Counter(
        model.quadratic(element) for element in model.elements
    )
    basis_images = tuple(
        model.apply_action(
            tuple(1 if i == j else 0 for i in range(len(model.moduli)))
        )
        for j in range(len(model.moduli))
    )
    return GoldenLannerK3SeedAudit(
        name=block.name,
        order=block.order,
        trace_gram=_matrix_tuple(gram),
        trace_inertia=inertia(block.restriction_trace_gram),
        discriminant_moduli=model.moduli,
        discriminant_quadratic_distribution=tuple(sorted(quadratic_distribution.items())),
        tau_basis_images=basis_images,
        even_overlattices=tuple(overlattices),
        k3_complement_signature=complement_signature,
        k3_complement_rank=complement_rank,
        nikulin_length_bound=length_bound,
        primitive_k3_embedding_exists=strict_signature_room and stable_range,
        primitive_k3_embedding_unique=strict_signature_room and stable_range,
    )


@lru_cache(maxsize=1)
def golden_lanner_k3_seed_audits() -> tuple[GoldenLannerK3SeedAudit, ...]:
    """Return the eight exact K3 seed audits in census order."""

    audits = tuple(audit_golden_lanner_k3_seed(block) for block in golden_lanner_blocks())
    if not all(audit.passes for audit in audits):
        raise AssertionError("a forced golden Lanner K3 seed audit failed")
    return audits


def _trace_generator_columns(
    roots: Sequence[Sequence[Golden]],
) -> SympyMatrix:
    columns: list[SympyMatrix] = []
    for root in roots:
        ordinary: list[int] = []
        translated: list[int] = []
        for entry in root:
            if not entry.is_integral:
                raise ValueError("a carrier seed root is not in Z[phi]")
            a = entry.a.numerator
            b = entry.b.numerator
            ordinary.extend((a, b))
            translated.extend((b, a + b))
        columns.extend((SympyMatrix(ordinary), SympyMatrix(translated)))
    return SympyMatrix.hstack(*columns)


def _rectangular_smith_diagonal(matrix: SympyMatrix) -> tuple[int, ...]:
    normal = smith_normal_form(matrix, domain=ZZ)
    return tuple(
        abs(int(normal[i, i]))
        for i in range(min(normal.shape))
        if normal[i, i] != 0
    )


def _combined_span(
    left: _DiscriminantModel,
    right: _DiscriminantModel,
    generators: Iterable[GluingElement],
) -> frozenset[GluingElement]:
    zero = (left.zero, right.zero)
    generators = tuple(generators)
    subgroup = {zero}
    frontier = [zero]
    while frontier:
        current_left, current_right = frontier.pop()
        for generator_left, generator_right in generators:
            value = (
                left.add(current_left, generator_left),
                right.add(current_right, generator_right),
            )
            if value not in subgroup:
                subgroup.add(value)
                frontier.append(value)
    return frozenset(subgroup)


@dataclass(frozen=True, slots=True)
class GoldenSeedCarrierCompletion:
    """A primitive seed and full discriminant gluing inside a rank-20 carrier."""

    carrier: str
    seed_name: str
    root_indices: tuple[int, ...]
    field_gram_verified: bool
    seed_trace_determinant: int
    seed_discriminant_group_invariants: tuple[int, ...]
    seed_tau_basis_images: tuple[DiscriminantElement, ...]
    seed_coordinate_smith_diagonal: tuple[int, ...]
    seed_primitive: bool
    complement_rank: int
    complement_inertia: Inertia
    complement_determinant: int
    complement_discriminant_group_invariants: tuple[int, ...]
    complement_tau_basis_images: tuple[DiscriminantElement, ...]
    complement_tau_polynomial_verified: bool
    complement_tau_self_adjoint_verified: bool
    direct_sum_index: int
    gluing_subgroup: tuple[GluingElement, ...]
    gluing_isotropic: bool
    gluing_tau_stable: bool
    gluing_projects_fully_to_seed: bool
    ambient_determinant: int
    ambient_discriminant_group_invariants: tuple[int, ...]
    ambient_tau_basis_images: tuple[DiscriminantElement, ...]

    @property
    def determinant_identity_verified(self) -> bool:
        return (
            abs(self.seed_trace_determinant * self.complement_determinant)
            == abs(self.ambient_determinant) * self.direct_sum_index**2
        )

    @property
    def passes(self) -> bool:
        return (
            self.field_gram_verified
            and self.seed_primitive
            and self.complement_inertia == Inertia(0, self.complement_rank, 0)
            and self.complement_tau_polynomial_verified
            and self.complement_tau_self_adjoint_verified
            and self.direct_sum_index == abs(self.seed_trace_determinant)
            and len(self.gluing_subgroup) == self.direct_sum_index
            and self.gluing_isotropic
            and self.gluing_tau_stable
            and self.gluing_projects_fully_to_seed
            and self.determinant_identity_verified
        )


def _carrier_completion(
    *,
    carrier: str,
    ambient_field_gram: Sequence[Sequence[Golden]],
    roots: Sequence[Sequence[Golden]],
    block: GoldenLannerAudit,
    root_indices: tuple[int, ...],
) -> GoldenSeedCarrierCompletion:
    selected_roots = tuple(roots[index] for index in root_indices)
    field_gram = tuple(
        tuple(
            sum(
                (
                    selected_roots[i][a]
                    * ambient_field_gram[a][b]
                    * selected_roots[j][b]
                    for a in range(10)
                    for b in range(10)
                ),
                start=Golden(0),
            )
            for j in range(block.order)
        )
        for i in range(block.order)
    )
    field_verified = field_gram == block.field_gram

    ambient_trace = SympyMatrix(
        integral_gram(trace_gram(ambient_field_gram, K3_TRACE_SCALAR))
    )
    ambient_action = _interleaved_golden_action(10)
    seed_columns = _trace_generator_columns(selected_roots)
    seed_gram = seed_columns.T * ambient_trace * seed_columns
    seed_action = _interleaved_golden_action(block.order)
    seed_model = _DiscriminantModel(seed_gram, seed_action)
    seed_smith = _rectangular_smith_diagonal(seed_columns)
    seed_primitive = all(value == 1 for value in seed_smith)

    pairing = seed_columns.T * ambient_trace
    normal, _, right = smith_normal_decomp(_domain_matrix(pairing))
    normal_matrix = normal.to_Matrix()
    pairing_rank = sum(
        normal_matrix[i, i] != 0 for i in range(min(normal_matrix.shape))
    )
    right_matrix = right.to_Matrix()
    complement_columns = right_matrix[:, pairing_rank:]
    if pairing * complement_columns != SympyMatrix.zeros(
        pairing.rows, complement_columns.cols
    ):
        raise ArithmeticError("the carrier complement basis is not orthogonal")
    complement_gram = complement_columns.T * ambient_trace * complement_columns

    transported_action = right_matrix.inv() * ambient_action * right_matrix
    if transported_action[:pairing_rank, pairing_rank:] != SympyMatrix.zeros(
        pairing_rank, complement_columns.cols
    ):
        raise ArithmeticError("the carrier complement is not tau-invariant")
    complement_action = transported_action[pairing_rank:, pairing_rank:]
    complement_model = _DiscriminantModel(complement_gram, complement_action)

    direct_sum = SympyMatrix.hstack(seed_columns, complement_columns)
    direct_sum_index = abs(int(direct_sum.det()))
    inverse_direct_sum = direct_sum.inv()
    gluing_generators: list[GluingElement] = []
    for index in range(ambient_trace.rows):
        coordinates = inverse_direct_sum[:, index]
        seed_coordinates = coordinates[: seed_columns.cols, :]
        complement_coordinates = coordinates[seed_columns.cols :, :]
        gluing_generators.append(
            (
                seed_model.canonical_from_lattice_coordinates(seed_coordinates),
                complement_model.canonical_from_lattice_coordinates(
                    complement_coordinates
                ),
            )
        )
    gluing = _combined_span(seed_model, complement_model, gluing_generators)
    gluing_isotropic = all(
        _mod_two(
            seed_model.quadratic(seed_element)
            + complement_model.quadratic(complement_element)
        )
        == 0
        for seed_element, complement_element in gluing
    )
    gluing_tau_stable = {
        (
            seed_model.apply_action(seed_element),
            complement_model.apply_action(complement_element),
        )
        for seed_element, complement_element in gluing
    } == set(gluing)
    seed_projection = {seed_element for seed_element, _ in gluing}
    complement_projection = {complement_element for _, complement_element in gluing}
    full_seed_projection = (
        seed_projection == set(seed_model.elements)
        and len(complement_projection) == len(gluing)
    )

    ambient_model = _DiscriminantModel(ambient_trace, ambient_action)
    seed_basis_images = tuple(
        seed_model.apply_action(
            tuple(1 if i == j else 0 for i in range(len(seed_model.moduli)))
        )
        for j in range(len(seed_model.moduli))
    )
    complement_basis_images = tuple(
        complement_model.apply_action(
            tuple(
                1 if i == j else 0
                for i in range(len(complement_model.moduli))
            )
        )
        for j in range(len(complement_model.moduli))
    )
    ambient_basis_images = tuple(
        ambient_model.apply_action(
            tuple(
                1 if i == j else 0
                for i in range(len(ambient_model.moduli))
            )
        )
        for j in range(len(ambient_model.moduli))
    )
    return GoldenSeedCarrierCompletion(
        carrier=carrier,
        seed_name=block.name or "unclassified",
        root_indices=root_indices,
        field_gram_verified=field_verified,
        seed_trace_determinant=int(seed_gram.det()),
        seed_discriminant_group_invariants=seed_model.moduli,
        seed_tau_basis_images=seed_basis_images,
        seed_coordinate_smith_diagonal=seed_smith,
        seed_primitive=seed_primitive,
        complement_rank=complement_columns.cols,
        complement_inertia=inertia(_matrix_tuple(complement_gram)),
        complement_determinant=int(complement_gram.det()),
        complement_discriminant_group_invariants=_smith_invariants(
            complement_gram
        ),
        complement_tau_basis_images=complement_basis_images,
        complement_tau_polynomial_verified=(
            complement_action * complement_action
            == complement_action + SympyMatrix.eye(complement_columns.cols)
        ),
        complement_tau_self_adjoint_verified=(
            complement_action.T * complement_gram
            == complement_gram * complement_action
        ),
        direct_sum_index=direct_sum_index,
        gluing_subgroup=tuple(sorted(gluing)),
        gluing_isotropic=gluing_isotropic,
        gluing_tau_stable=gluing_tau_stable,
        gluing_projects_fully_to_seed=full_seed_projection,
        ambient_determinant=int(ambient_trace.det()),
        ambient_discriminant_group_invariants=ambient_model.moduli,
        ambient_tau_basis_images=ambient_basis_images,
    )


_H19_ROOT_INDICES = {
    "A4(3-5-3)": (0, 9, 14, 13),
    "A4(5-3-5)": (0, 12, 8, 9),
    "D4(5,3,3)": (0, 12, 13, 9),
    "C4(one 5)": (1, 0, 13, 16),
    "C4(opposite 5s)": (11, 19, 13, 16),
    "A5(5-3-3-3)": (1, 0, 13, 9, 2),
    "A5(5-3-3-5)": (2, 4, 3, 16, 17),
    "D5(terminal 5)": (1, 13, 6, 2, 16),
}


@lru_cache(maxsize=1)
def h19_golden_seed_completions() -> tuple[GoldenSeedCarrierCompletion, ...]:
    """Embed all eight order-four/five seeds in the discriminant-19 carrier."""

    roots = run_h19_vinberg_search(5).roots
    completions = tuple(
        _carrier_completion(
            carrier="H19",
            ambient_field_gram=h19_basis_gram(),
            roots=roots,
            block=block,
            root_indices=_H19_ROOT_INDICES[block.name or ""],
        )
        for block in golden_lanner_blocks()
    )
    if not all(completion.passes for completion in completions):
        raise AssertionError("an H19 golden seed completion failed")
    return completions


@lru_cache(maxsize=1)
def h11_matching_seed_completion() -> GoldenSeedCarrierCompletion:
    """Embed the determinant-11 path seed in the H11 carrier."""

    block = next(
        block for block in golden_lanner_blocks() if block.name == "A4(3-5-3)"
    )
    completion = _carrier_completion(
        carrier="H11",
        ambient_field_gram=h11_basis_gram(),
        roots=run_h11_vinberg_search(6).roots,
        block=block,
        root_indices=(0, 11, 16, 10),
    )
    if not completion.passes:
        raise AssertionError("the H11 determinant-11 seed completion failed")
    return completion

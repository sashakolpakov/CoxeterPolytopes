"""Exact one-cross-edge shell for the three H9/14 multiwedge nerves.

This module deliberately studies only the first layer outside the frozen
sparse families of records 43, 44, and 47.  For every sparse symmetry-orbit
representative it changes exactly one of the 45 right-angled pairs joining
the isolated ``boundary(Delta_4)`` factor to the residual nine facets, to
one of the labels 3, 4, or 5.

The quotient is exact.  We generate the full automorphism group of the
minimal-nonface hypergraph, take the stabilizer of each labelled sparse
representative, and quotient its 45 cross pairs by that stabilizer.  Thus a
reported orbit is an orbit of *labelled one-edge extensions*, not merely an
orbit of unlabelled facet pairs.

The local test is also exact.  The sparse representative has already been
checked at all 90 rank-nine vertices and all four order-five Lanner
minimal nonfaces.  A cross edge belongs to no order-five minimal nonface, so
those four tests are unchanged; only rank-nine vertices containing both
endpoints can change, and every such vertex is tested with the complete
finite Coxeter classifier.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from typing import Iterable, Mapping, Sequence

import sympy as sp

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    Pair,
    audit_h9_14_coxeter_labels,
    is_finite_coxeter_labeling,
    iter_product_orthogonal_label_representatives,
    pair,
)
from coxeter_polytopes.combinatorics.h9_14_multiwedge import (
    RECORD43,
    RECORD44,
    RECORD47,
    MultiwedgeSpecification,
)
from coxeter_polytopes.combinatorics.h9_14_record43 import (
    record43_facet_complements,
)
from coxeter_polytopes.combinatorics.h9_14_record44 import (
    record44_facet_complements,
)
from coxeter_polytopes.combinatorics.h9_14_record47 import (
    record47_facet_complements,
)
from coxeter_polytopes.combinatorics.h9_14_sparse_coxeter_labels import (
    audit_sparse_h9_14_labels,
)


FACET_COUNT = 14
ONE_EDGE_LABELS = (3, 4, 5)

Permutation = tuple[int, ...]


def _permutation(*cycles: tuple[int, ...]) -> Permutation:
    image = list(range(FACET_COUNT))
    for cycle in cycles:
        if len(cycle) < 2:
            raise ValueError("a permutation cycle must have at least two entries")
        for source, target in zip(cycle, cycle[1:] + cycle[:1], strict=True):
            image[source] = target
    if len(set(image)) != FACET_COUNT:
        raise ValueError("cycles do not define a permutation")
    return tuple(image)


def _adjacent_transpositions(vertices: tuple[int, ...]) -> tuple[Permutation, ...]:
    return tuple(
        _permutation((vertices[index], vertices[index + 1]))
        for index in range(len(vertices) - 1)
    )


# Independently repeated here so that the one-edge certificate verifies its
# own group action.  Each generator is checked on all vertices and minimal
# nonfaces below before it is used.
_NERVE_GENERATORS: dict[int, tuple[Permutation, ...]] = {
    43: (
        *_adjacent_transpositions((3, 4, 5, 6, 7)),
        *_adjacent_transpositions((8, 9, 10)),
        _permutation((0, 13), (1, 12), (2, 11)),
    ),
    44: (
        *_adjacent_transpositions((2, 3, 4, 5, 6)),
        *_adjacent_transpositions((7, 8, 9, 10)),
        _permutation((0, 13), (1, 12)),
    ),
    47: (
        *_adjacent_transpositions((3, 4, 5, 6, 7)),
        *_adjacent_transpositions((8, 9, 10)),
        _permutation((0, 12), (11, 1)),
        _permutation((12, 13), (1, 2)),
    ),
}

_EXPECTED_GROUP_ORDERS = {43: 1_440, 44: 5_760, 47: 4_320}


@dataclass(frozen=True, slots=True)
class _Context:
    specification: MultiwedgeSpecification
    vertices: tuple[tuple[int, ...], ...]
    nonfaces: tuple[tuple[int, ...], ...]
    isolated: tuple[int, ...]
    residual: tuple[int, ...]
    cross_pairs: tuple[Pair, ...]
    group: tuple[Permutation, ...]


@dataclass(frozen=True, slots=True)
class _ProductContext:
    record_index: int
    nonfaces: tuple[tuple[int, ...], ...]
    vertices: tuple[tuple[int, ...], ...]
    factors: tuple[tuple[int, ...], ...]
    factor_vertex_sections: tuple[tuple[tuple[int, ...], ...], ...]
    cross_pairs: tuple[Pair, ...]
    full_group_order: int


def _compose(left: Permutation, right: Permutation) -> Permutation:
    return tuple(left[right[index]] for index in range(FACET_COUNT))


def _generated_group(generators: Sequence[Permutation]) -> tuple[Permutation, ...]:
    identity = tuple(range(FACET_COUNT))
    group = {identity}
    frontier = deque((identity,))
    while frontier:
        current = frontier.popleft()
        for generator in generators:
            image = _compose(generator, current)
            if image not in group:
                group.add(image)
                frontier.append(image)
    return tuple(sorted(group))


def _vertices_from_complements(complements: Iterable[int]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        sorted(
            tuple(index for index in range(FACET_COUNT) if not (mask & (1 << index)))
            for mask in complements
        )
    )


_SPECIFICATIONS = {43: RECORD43, 44: RECORD44, 47: RECORD47}
_COMPLEMENT_FUNCTIONS = {
    43: record43_facet_complements,
    44: record44_facet_complements,
    47: record47_facet_complements,
}


@lru_cache(maxsize=3)
def _context(record_index: int) -> _Context:
    specification = _SPECIFICATIONS[record_index]
    vertices = _vertices_from_complements(_COMPLEMENT_FUNCTIONS[record_index]())
    nonfaces = tuple(specification.missing_edges) + tuple(specification.missing_fives)
    isolated = tuple(specification.simplex_join_vertices)
    residual = tuple(sorted(set(range(FACET_COUNT)) - set(isolated)))
    cross_pairs = tuple(pair(left, right) for left in isolated for right in residual)
    if len(vertices) != 90 or len(cross_pairs) != 45:
        raise AssertionError("the archived 90-vertex join factorization changed")
    if any(
        {left, right} <= set(face)
        for left, right in cross_pairs
        for face in specification.missing_fives
    ):
        raise AssertionError("a declared cross pair lies in a Lanner support")

    generators = _NERVE_GENERATORS[record_index]
    vertex_set = {frozenset(vertex) for vertex in vertices}
    nonface_set = {frozenset(face) for face in nonfaces}
    for generator in generators:
        if {
            frozenset(generator[vertex] for vertex in face) for face in vertices
        } != vertex_set:
            raise AssertionError("a one-edge generator does not preserve vertices")
        if {
            frozenset(generator[vertex] for vertex in face) for face in nonfaces
        } != nonface_set:
            raise AssertionError("a one-edge generator does not preserve nonfaces")
    group = _generated_group(generators)
    if len(group) != _EXPECTED_GROUP_ORDERS[record_index]:
        raise AssertionError("the full nerve automorphism group order changed")
    return _Context(
        specification=specification,
        vertices=vertices,
        nonfaces=nonfaces,
        isolated=isolated,
        residual=residual,
        cross_pairs=tuple(sorted(cross_pairs)),
        group=group,
    )


def _label_code(labels: Mapping[Pair, int]) -> tuple[tuple[int, int, int], ...]:
    return tuple((left, right, labels[(left, right)]) for left, right in sorted(labels))


def _label_digest(labels: Mapping[Pair, int]) -> str:
    payload = json.dumps(_label_code(labels), separators=(",", ":")).encode()
    return sha256(payload).hexdigest()


def _preserves_labels(permutation: Permutation, labels: Mapping[Pair, int]) -> bool:
    return all(
        labels[pair(permutation[left], permutation[right])] == label
        for (left, right), label in labels.items()
    )


def _pair_orbits(
    cross_pairs: Sequence[Pair], stabilizer: Sequence[Permutation]
) -> tuple[tuple[Pair, tuple[Pair, ...]], ...]:
    remaining = set(cross_pairs)
    result: list[tuple[Pair, tuple[Pair, ...]]] = []
    while remaining:
        seed = min(remaining)
        orbit = tuple(
            sorted({pair(permutation[seed[0]], permutation[seed[1]]) for permutation in stabilizer})
        )
        if not set(orbit) <= set(cross_pairs):
            raise AssertionError("a labelled stabilizer left the cross-pair set")
        result.append((min(orbit), orbit))
        remaining.difference_update(orbit)
    if sum(len(orbit) for _representative, orbit in result) != len(cross_pairs):
        raise AssertionError("cross-pair stabilizer orbits do not partition 45 pairs")
    return tuple(result)


def _pair_orbits_from_generators(
    cross_pairs: Sequence[Pair], generators: Sequence[Permutation]
) -> tuple[tuple[Pair, tuple[Pair, ...]], ...]:
    cross_set = set(cross_pairs)
    remaining = set(cross_pairs)
    result: list[tuple[Pair, tuple[Pair, ...]]] = []
    while remaining:
        seed = min(remaining)
        orbit = {seed}
        frontier = deque((seed,))
        while frontier:
            current = frontier.popleft()
            for permutation in generators:
                image = pair(permutation[current[0]], permutation[current[1]])
                if image not in cross_set:
                    raise AssertionError("a product stabilizer left the cross-pair set")
                if image not in orbit:
                    orbit.add(image)
                    frontier.append(image)
        ordered = tuple(sorted(orbit))
        result.append((ordered[0], ordered))
        remaining.difference_update(orbit)
    if sum(len(orbit) for _edge, orbit in result) != len(cross_pairs):
        raise AssertionError("product cross-pair orbits do not partition the shell")
    return tuple(result)


def _changed_local_pass(
    context: _Context,
    baseline: Mapping[Pair, int],
    edge: Pair,
    label: int,
) -> tuple[bool, int, tuple[int, ...]]:
    if edge not in context.cross_pairs or baseline[edge] != 2:
        raise ValueError("the one-edge shell changes one right-angled cross pair")
    if label not in ONE_EDGE_LABELS:
        raise ValueError("the one-edge finite label must be 3, 4, or 5")
    changed = dict(baseline)
    changed[edge] = label
    affected = tuple(
        index
        for index, vertex in enumerate(context.vertices)
        if edge[0] in vertex and edge[1] in vertex
    )
    bad = tuple(
        index
        for index in affected
        if not is_finite_coxeter_labeling(context.vertices[index], changed)
    )
    return not bad, len(affected), bad


def _gram(support: tuple[int, ...], labels: Mapping[Pair, int]) -> sp.Matrix:
    matrix = 2 * sp.eye(len(support))
    for left_index, right_index in combinations(range(len(support)), 2):
        left, right = support[left_index], support[right_index]
        value = -2 * sp.cos(sp.pi / labels[pair(left, right)])
        matrix[left_index, right_index] = value
        matrix[right_index, left_index] = value
    return matrix


@lru_cache(maxsize=None)
def _schur_sign_from_codes(
    isolated_code: tuple[int, ...],
    residual_code: tuple[int, ...],
    isolated_position: int,
    residual_position: int,
    cross_label: int,
) -> tuple[int, str]:
    local_pairs = tuple(combinations(range(5), 2))
    isolated_labels = dict(zip(local_pairs, isolated_code, strict=True))
    residual_labels = dict(zip(local_pairs, residual_code, strict=True))
    isolated_gram = _gram(tuple(range(5)), isolated_labels)
    residual_gram = _gram(tuple(range(5)), residual_labels)
    isolated_det = sp.factor(
        isolated_gram.det(), extension=[sp.sqrt(2), sp.sqrt(5)]
    )
    residual_det = sp.factor(
        residual_gram.det(), extension=[sp.sqrt(2), sp.sqrt(5)]
    )
    isolated_minor = isolated_gram.minor_submatrix(
        isolated_position, isolated_position
    ).det()
    residual_minor = residual_gram.minor_submatrix(
        residual_position, residual_position
    ).det()
    if isolated_det.is_negative is not True or residual_det.is_negative is not True:
        raise AssertionError("the Schur test did not receive two Lanner blocks")
    if isolated_minor.is_positive is not True or residual_minor.is_positive is not True:
        raise AssertionError("a Lanner one-node deletion is not positive definite")
    inverse_diagonal_product = (
        isolated_minor * residual_minor / (isolated_det * residual_det)
    )
    coupling_square = (2 * sp.cos(sp.pi / cross_label)) ** 2
    schur_factor = sp.factor(
        1 - coupling_square * inverse_diagonal_product,
        extension=[sp.sqrt(2), sp.sqrt(5)],
    )
    if schur_factor.is_positive is True:
        sign = 1
    elif schur_factor.is_negative is True:
        sign = -1
    elif schur_factor.is_zero is True:
        sign = 0
    else:
        raise AssertionError("the exact algebraic Schur sign was undecided")
    return sign, sp.sstr(schur_factor)


def _induced_code(
    support: tuple[int, ...], labels: Mapping[Pair, int]
) -> tuple[int, ...]:
    return tuple(labels[pair(left, right)] for left, right in combinations(support, 2))


def _two_lanner_obstruction(
    context: _Context,
    baseline: Mapping[Pair, int],
    edge: Pair,
    label: int,
) -> dict[str, object]:
    isolated_endpoint = edge[0] if edge[0] in context.isolated else edge[1]
    residual_endpoint = edge[1] if edge[0] in context.isolated else edge[0]
    isolated = context.isolated
    isolated_code = _induced_code(isolated, baseline)
    witnesses = []
    for support in context.specification.missing_fives[1:]:
        if residual_endpoint not in support:
            return {
                "status": "negative_index_at_least_two",
                "method": "untouched_orthogonal_lanner_pair",
                "isolated_support": list(isolated),
                "residual_support": list(support),
            }
        sign, expression = _schur_sign_from_codes(
            isolated_code,
            _induced_code(support, baseline),
            isolated.index(isolated_endpoint),
            support.index(residual_endpoint),
            label,
        )
        witnesses.append(
            {
                "residual_support": list(support),
                "schur_factor_sign": sign,
                "schur_factor": expression,
            }
        )
        if sign > 0:
            return {
                "status": "negative_index_at_least_two",
                "method": "exact_two_lanner_schur_factor",
                "isolated_support": list(isolated),
                **witnesses[-1],
            }
    return {
        "status": "not_excluded_by_two_lanner_principal_test",
        "method": "exact_two_lanner_schur_factor",
        "isolated_support": list(isolated),
        "tests": witnesses,
    }


@lru_cache(maxsize=3)
def audit_multiwedge_one_cross_edge(record_index: int) -> dict[str, object]:
    """Return the exact one-cross-edge shell for one multiwedge record."""

    context = _context(record_index)
    sparse = audit_sparse_h9_14_labels(record_index)
    baselines = sparse.orbit_representative_label_maps()
    if len(baselines) != sparse.symmetry_orbit_count:
        raise AssertionError("the sparse representative export changed")

    total_labelled_trials = 0
    total_labelled_local_survivors = 0
    extension_orbit_count = 0
    local_survivor_orbit_count = 0
    metric_excluded_orbit_count = 0
    metric_unresolved_orbit_count = 0
    affected_histogram: Counter[int] = Counter()
    stabilizer_histogram: Counter[int] = Counter()
    pair_orbit_histogram: Counter[int] = Counter()
    survivor_records: list[dict[str, object]] = []
    baseline_records: list[dict[str, object]] = []

    for baseline_index, baseline in enumerate(baselines):
        baseline_audit = audit_h9_14_coxeter_labels(
            context.vertices, context.nonfaces, baseline
        )
        if not baseline_audit.passes_discrete_gate:
            raise AssertionError("a frozen sparse representative no longer passes")
        if any(baseline[edge] != 2 for edge in context.cross_pairs):
            raise AssertionError("a frozen sparse cross pair is not right-angled")

        stabilizer = tuple(
            permutation
            for permutation in context.group
            if _preserves_labels(permutation, baseline)
        )
        if not stabilizer or len(context.group) % len(stabilizer):
            raise AssertionError("invalid labelled stabilizer")
        baseline_orbit_size = len(context.group) // len(stabilizer)
        pair_orbits = _pair_orbits(context.cross_pairs, stabilizer)
        stabilizer_histogram[len(stabilizer)] += 1
        pair_orbit_histogram[len(pair_orbits)] += 1
        total_labelled_trials += baseline_orbit_size * len(context.cross_pairs) * 3
        extension_orbit_count += len(pair_orbits) * 3

        baseline_local_orbits = 0
        baseline_local_labelled = 0
        for edge, edge_orbit in pair_orbits:
            for label in ONE_EDGE_LABELS:
                passes, affected_count, bad_vertices = _changed_local_pass(
                    context, baseline, edge, label
                )
                affected_histogram[affected_count] += 1
                if not passes:
                    continue
                baseline_local_orbits += 1
                local_survivor_orbit_count += 1
                labelled_multiplicity = baseline_orbit_size * len(edge_orbit)
                baseline_local_labelled += labelled_multiplicity
                total_labelled_local_survivors += labelled_multiplicity
                metric = _two_lanner_obstruction(
                    context, baseline, edge, label
                )
                if metric["status"] == "negative_index_at_least_two":
                    metric_excluded_orbit_count += 1
                else:
                    metric_unresolved_orbit_count += 1
                survivor_records.append(
                    {
                        "baseline_orbit_index": baseline_index,
                        "baseline_label_sha256": _label_digest(baseline),
                        "edge": list(edge),
                        "label": label,
                        "edge_stabilizer_orbit_size": len(edge_orbit),
                        "full_labelled_orbit_size": labelled_multiplicity,
                        "affected_vertex_count": affected_count,
                        "bad_vertex_indices": list(bad_vertices),
                        "metric_principal_test": metric,
                    }
                )
        baseline_records.append(
            {
                "baseline_orbit_index": baseline_index,
                "baseline_label_sha256": _label_digest(baseline),
                "labelled_stabilizer_order": len(stabilizer),
                "baseline_full_nerve_orbit_size": baseline_orbit_size,
                "cross_pair_orbit_count": len(pair_orbits),
                "one_edge_extension_orbit_count": 3 * len(pair_orbits),
                "locally_admissible_extension_orbit_count": baseline_local_orbits,
                "locally_admissible_full_labelled_count": baseline_local_labelled,
            }
        )

    if total_labelled_trials != sparse.full_labelled_assignment_count * 45 * 3:
        raise AssertionError("full labelled one-edge trial count did not close")
    if metric_excluded_orbit_count + metric_unresolved_orbit_count != local_survivor_orbit_count:
        raise AssertionError("metric survivor partition did not close")

    return {
        "record_index": record_index,
        "graph6": context.specification.graph6,
        "scope": {
            "baseline_family": "frozen sparse Coxeter-label family",
            "changed_cross_edge_count": 1,
            "new_label_domain": list(ONE_EDGE_LABELS),
            "unrestricted_label_census": False,
            "dotted_entries_assigned": False,
        },
        "nerve": {
            "facet_count": FACET_COUNT,
            "rank_nine_vertex_count": len(context.vertices),
            "higher_lanner_nonface_count": len(context.specification.missing_fives),
            "isolated_factor": list(context.isolated),
            "residual_factor": list(context.residual),
            "right_cross_pair_count": len(context.cross_pairs),
            "full_automorphism_group_order": len(context.group),
        },
        "baseline": {
            "symmetry_orbit_count": len(baselines),
            "full_labelled_assignment_count": sparse.full_labelled_assignment_count,
            "stabilizer_order_histogram": {
                str(value): count for value, count in sorted(stabilizer_histogram.items())
            },
            "records": baseline_records,
        },
        "one_edge_shell": {
            "full_labelled_trial_count": total_labelled_trials,
            "full_nerve_symmetry_orbit_count": extension_orbit_count,
            "cross_pair_orbit_count_histogram_per_baseline": {
                str(value): count for value, count in sorted(pair_orbit_histogram.items())
            },
            "locally_admissible_full_labelled_count": total_labelled_local_survivors,
            "locally_admissible_full_nerve_symmetry_orbit_count": (
                local_survivor_orbit_count
            ),
            "affected_vertex_count_histogram_over_orbit_trials": {
                str(value): count for value, count in sorted(affected_histogram.items())
            },
            "all_unaffected_vertices_reused_from_frozen_exact_audit": True,
            "all_affected_vertices_checked": True,
            "all_higher_lanner_nonfaces_unchanged_and_checked": True,
        },
        "metric_principal_gate": {
            "negative_index_at_least_two_orbit_count": metric_excluded_orbit_count,
            "not_excluded_orbit_count": metric_unresolved_orbit_count,
            "method": (
                "Cauchy interlacing on either an untouched orthogonal pair of "
                "Lanner principal blocks or an exact one-edge two-Lanner "
                "Schur factor in Q(sqrt(2),sqrt(5))."
            ),
        },
        "survivor_orbit_representatives": survivor_records,
    }


_PRODUCT_NONFACES: dict[int, tuple[tuple[int, ...], ...]] = {
    30: (
        (0, 11),
        (1, 12),
        (2, 13),
        (3, 4, 5, 6),
        (7, 8, 9, 10),
    ),
    42: (
        (0, 11),
        (0, 12),
        (1, 12),
        (1, 13),
        (11, 13),
        (2, 3, 4, 5),
        (6, 7, 8, 9, 10),
    ),
}


@lru_cache(maxsize=2)
def _product_context(record_index: int) -> _ProductContext:
    nonfaces = _PRODUCT_NONFACES[record_index]
    vertices = tuple(
        subset
        for subset in combinations(range(FACET_COUNT), 9)
        if not any(set(face) <= set(subset) for face in nonfaces)
    )
    if record_index == 30:
        factors = ((0, 11), (1, 12), (2, 13), (3, 4, 5, 6), (7, 8, 9, 10))
        expected_vertices = 128
        full_group_order = 55_296
    elif record_index == 42:
        polygon = tuple(sorted({vertex for face in nonfaces[:5] for vertex in face}))
        factors = (polygon, (2, 3, 4, 5), (6, 7, 8, 9, 10))
        expected_vertices = 100
        full_group_order = 28_800
    else:
        raise ValueError("the product one-edge shell covers records 30 and 42")
    if len(vertices) != expected_vertices:
        raise AssertionError("the product rank-nine vertex census changed")
    if set().union(*map(set, factors)) != set(range(FACET_COUNT)) or any(
        set(left) & set(right) for left, right in combinations(factors, 2)
    ):
        raise AssertionError("the product factors no longer partition 14 facets")
    cross_pairs = tuple(
        sorted(
            pair(left, right)
            for first, second in combinations(factors, 2)
            for left in first
            for right in second
        )
    )
    expected_cross = {30: 76, 42: 65}[record_index]
    if len(cross_pairs) != expected_cross:
        raise AssertionError("the product cross-pair count changed")
    if any(
        {left, right} <= set(face)
        for left, right in cross_pairs
        for face in nonfaces
        if len(face) > 2
    ):
        raise AssertionError("a product cross edge lies inside a Lanner support")
    factor_vertex_sections = tuple(
        tuple(
            sorted(
                {
                    tuple(vertex for vertex in maximal if vertex in factor)
                    for maximal in vertices
                }
            )
        )
        for factor in factors
    )
    reconstructed = {
        tuple(sorted(vertex for section in sections for vertex in section))
        for sections in product(*factor_vertex_sections)
    }
    if reconstructed != set(vertices):
        raise AssertionError("the product vertices do not factor as declared")
    return _ProductContext(
        record_index=record_index,
        nonfaces=nonfaces,
        vertices=vertices,
        factors=factors,
        factor_vertex_sections=factor_vertex_sections,
        cross_pairs=cross_pairs,
        full_group_order=full_group_order,
    )


def _factor_automorphisms(
    factor: tuple[int, ...], labels: Mapping[Pair, int]
) -> tuple[Permutation, ...]:
    """All internal labelled automorphisms, dotted-vs-meeting included."""

    result = []
    internal_pairs = tuple(combinations(factor, 2))
    for image_vertices in permutations(factor):
        local_image = dict(zip(factor, image_vertices, strict=True))
        if all(
            labels.get(pair(left, right))
            == labels.get(pair(local_image[left], local_image[right]))
            for left, right in internal_pairs
        ):
            image = list(range(FACET_COUNT))
            for source, target in local_image.items():
                image[source] = target
            result.append(tuple(image))
    if not result:
        raise AssertionError("a labelled factor lost its identity automorphism")
    return tuple(sorted(result))


def _block_swap(left: tuple[int, ...], right: tuple[int, ...]) -> Permutation:
    if len(left) != len(right):
        raise ValueError("only equal factor blocks can be swapped")
    image = list(range(FACET_COUNT))
    for first, second in zip(left, right, strict=True):
        image[first], image[second] = second, first
    return tuple(image)


def _product_stabilizer_data(
    context: _ProductContext, labels: Mapping[Pair, int]
) -> tuple[int, tuple[Permutation, ...]]:
    factors = context.factors
    automorphisms = tuple(_factor_automorphisms(factor, labels) for factor in factors)
    generators: list[Permutation] = [
        permutation
        for factor_automorphisms in automorphisms
        for permutation in factor_automorphisms
    ]
    if context.record_index == 30:
        # The three dotted two-blocks are indistinguishable.  The two Lanner
        # blocks can be interchanged precisely when their canonical labelled
        # codes agree.
        generators.extend(
            (
                _block_swap(factors[0], factors[1]),
                _block_swap(factors[1], factors[2]),
            )
        )
        first_code = tuple(
            labels[pair(left, right)] for left, right in combinations(factors[3], 2)
        )
        second_code = tuple(
            labels[pair(left, right)] for left, right in combinations(factors[4], 2)
        )
        lanner_swap = first_code == second_code
        if lanner_swap:
            generators.append(_block_swap(factors[3], factors[4]))
        stabilizer_order = (
            48
            * len(automorphisms[3])
            * len(automorphisms[4])
            * (2 if lanner_swap else 1)
        )
    else:
        stabilizer_order = 1
        for group in automorphisms:
            stabilizer_order *= len(group)
    if context.full_group_order % stabilizer_order:
        raise AssertionError("a product labelled stabilizer has invalid order")
    if any(not _preserves_labels(generator, labels) for generator in generators):
        raise AssertionError("a product stabilizer generator changed a finite label")
    return stabilizer_order, tuple(generators)


def _changed_product_local_pass(
    context: _ProductContext,
    baseline: Mapping[Pair, int],
    edge: Pair,
    label: int,
) -> tuple[bool, int, tuple[int, ...]]:
    factor_indices = tuple(
        index
        for index, factor in enumerate(context.factors)
        if edge[0] in factor or edge[1] in factor
    )
    if len(factor_indices) != 2:
        raise AssertionError("a cross edge must join exactly two product factors")
    first_index, second_index = factor_indices
    first_factor, second_factor = context.factors[first_index], context.factors[second_index]
    first_endpoint = edge[0] if edge[0] in first_factor else edge[1]
    second_endpoint = edge[1] if edge[0] in first_factor else edge[0]
    first_code = tuple(
        baseline.get(pair(left, right), 0)
        for left, right in combinations(first_factor, 2)
    )
    second_code = tuple(
        baseline.get(pair(left, right), 0)
        for left, right in combinations(second_factor, 2)
    )
    return _cached_changed_product_local_pass(
        context.record_index,
        first_index,
        second_index,
        first_code,
        second_code,
        first_endpoint,
        second_endpoint,
        label,
    )


@lru_cache(maxsize=None)
def _cached_changed_product_local_pass(
    record_index: int,
    first_index: int,
    second_index: int,
    first_code: tuple[int, ...],
    second_code: tuple[int, ...],
    first_endpoint: int,
    second_endpoint: int,
    label: int,
) -> tuple[bool, int, tuple[int, ...]]:
    context = _product_context(record_index)
    first_factor, second_factor = context.factors[first_index], context.factors[second_index]
    first_sections = tuple(
        section
        for section in context.factor_vertex_sections[first_index]
        if first_endpoint in section
    )
    second_sections = tuple(
        section
        for section in context.factor_vertex_sections[second_index]
        if second_endpoint in section
    )
    changed: dict[Pair, int] = {
        pair(left, right): value
        for (left, right), value in zip(
            combinations(first_factor, 2), first_code, strict=True
        )
        if value
    }
    changed.update(
        {
            pair(left, right): value
            for (left, right), value in zip(
                combinations(second_factor, 2), second_code, strict=True
            )
            if value
        }
    )
    for left in first_factor:
        for right in second_factor:
            changed[pair(left, right)] = 2
    changed[pair(first_endpoint, second_endpoint)] = label
    bad_sections = tuple(
        tuple(sorted((*first, *second)))
        for first in first_sections
        for second in second_sections
        if not is_finite_coxeter_labeling((*first, *second), changed)
    )
    other_multiplicity = 1
    for index, sections in enumerate(context.factor_vertex_sections):
        if index not in (first_index, second_index):
            other_multiplicity *= len(sections)
    affected_count = len(first_sections) * len(second_sections) * other_multiplicity
    # ``bad_sections`` are the complete distinct affected factor components.
    # Their Cartesian products with the unchanged orthogonal factors are bad
    # exactly when the section itself is bad.
    return not bad_sections, affected_count, tuple(range(len(bad_sections)))


@lru_cache(maxsize=2)
def audit_product_one_cross_edge(record_index: int) -> dict[str, object]:
    """Exact labelled-stabilizer quotient for a product family's one-edge shell."""

    context = _product_context(record_index)
    baselines = tuple(
        iter_product_orthogonal_label_representatives(context.nonfaces)
    )
    expected_baseline_orbits = {30: 45, 42: 6_120}[record_index]
    expected_full_baselines = {30: 11_664, 42: 46_448_640}[record_index]
    if len(baselines) != expected_baseline_orbits:
        raise AssertionError("the frozen product representative count changed")

    full_baseline_count = 0
    total_labelled_trials = 0
    total_labelled_local_survivors = 0
    extension_orbit_count = 0
    local_survivor_orbit_count = 0
    stabilizer_histogram: Counter[int] = Counter()
    pair_orbit_histogram: Counter[int] = Counter()
    affected_histogram: Counter[int] = Counter()
    survivors: list[dict[str, object]] = []
    baseline_records: list[dict[str, object]] = []

    for baseline_index, representative in enumerate(baselines):
        baseline = representative.label_map()
        audit = audit_h9_14_coxeter_labels(
            context.vertices, context.nonfaces, baseline
        )
        if not audit.passes_discrete_gate:
            raise AssertionError("a frozen product representative no longer passes")
        if any(baseline[edge] != 2 for edge in context.cross_pairs):
            raise AssertionError("a frozen product cross pair is not right-angled")
        stabilizer_order, generators = _product_stabilizer_data(context, baseline)
        baseline_orbit_size = context.full_group_order // stabilizer_order
        full_baseline_count += baseline_orbit_size
        pair_orbits = _pair_orbits_from_generators(context.cross_pairs, generators)
        stabilizer_histogram[stabilizer_order] += 1
        pair_orbit_histogram[len(pair_orbits)] += 1
        total_labelled_trials += baseline_orbit_size * len(context.cross_pairs) * 3
        extension_orbit_count += len(pair_orbits) * 3
        baseline_local_orbits = 0
        baseline_local_labelled = 0
        for edge, edge_orbit in pair_orbits:
            for label in ONE_EDGE_LABELS:
                passes, affected_count, bad = _changed_product_local_pass(
                    context, baseline, edge, label
                )
                affected_histogram[affected_count] += 1
                if not passes:
                    continue
                baseline_local_orbits += 1
                local_survivor_orbit_count += 1
                multiplicity = baseline_orbit_size * len(edge_orbit)
                baseline_local_labelled += multiplicity
                total_labelled_local_survivors += multiplicity
                survivors.append(
                    {
                        "baseline_orbit_index": baseline_index,
                        "baseline_label_sha256": _label_digest(baseline),
                        "baseline_parameters": {
                            "lanner_type_indices": list(representative.lanner_type_indices),
                            "cycle_color_code": list(representative.cycle_color_code),
                        },
                        "edge": list(edge),
                        "label": label,
                        "edge_stabilizer_orbit_size": len(edge_orbit),
                        "full_labelled_orbit_size": multiplicity,
                        "affected_vertex_count": affected_count,
                        "bad_vertex_indices": list(bad),
                        "metric_principal_test": {
                            "status": "negative_index_at_least_two",
                            "method": "product_component_interlacing",
                            "record_30_lower_bound": 4 if record_index == 30 else None,
                            "record_42_lower_bound": 2 if record_index == 42 else None,
                        },
                    }
                )
        baseline_records.append(
            {
                "baseline_orbit_index": baseline_index,
                "baseline_label_sha256": _label_digest(baseline),
                "labelled_stabilizer_order": stabilizer_order,
                "baseline_full_nerve_orbit_size": baseline_orbit_size,
                "cross_pair_orbit_count": len(pair_orbits),
                "one_edge_extension_orbit_count": 3 * len(pair_orbits),
                "locally_admissible_extension_orbit_count": baseline_local_orbits,
                "locally_admissible_full_labelled_count": baseline_local_labelled,
            }
        )

    if full_baseline_count != expected_full_baselines:
        raise AssertionError("product stabilizer weights did not recover the sparse census")
    if total_labelled_trials != expected_full_baselines * len(context.cross_pairs) * 3:
        raise AssertionError("product full labelled one-edge count did not close")
    return {
        "record_index": record_index,
        "scope": {
            "baseline_family": "frozen factor-orthogonal Coxeter-label family",
            "changed_cross_edge_count": 1,
            "new_label_domain": list(ONE_EDGE_LABELS),
            "unrestricted_label_census": False,
            "dotted_entries_assigned": False,
        },
        "nerve": {
            "facet_count": FACET_COUNT,
            "rank_nine_vertex_count": len(context.vertices),
            "higher_lanner_nonface_count": sum(len(face) > 2 for face in context.nonfaces),
            "factor_blocks": [list(factor) for factor in context.factors],
            "right_cross_pair_count": len(context.cross_pairs),
            "full_automorphism_group_order": context.full_group_order,
        },
        "baseline": {
            "symmetry_orbit_count": len(baselines),
            "full_labelled_assignment_count": full_baseline_count,
            "stabilizer_order_histogram": {
                str(value): count for value, count in sorted(stabilizer_histogram.items())
            },
            "records": baseline_records,
        },
        "one_edge_shell": {
            "full_labelled_trial_count": total_labelled_trials,
            "full_nerve_symmetry_orbit_count": extension_orbit_count,
            "cross_pair_orbit_count_histogram_per_baseline": {
                str(value): count for value, count in sorted(pair_orbit_histogram.items())
            },
            "locally_admissible_full_labelled_count": total_labelled_local_survivors,
            "locally_admissible_full_nerve_symmetry_orbit_count": local_survivor_orbit_count,
            "affected_vertex_count_histogram_over_orbit_trials": {
                str(value): count for value, count in sorted(affected_histogram.items())
            },
            "all_unaffected_vertices_reused_from_frozen_exact_audit": True,
            "all_affected_vertices_checked": True,
            "all_higher_lanner_nonfaces_unchanged_and_checked": True,
        },
        "metric_principal_gate": {
            "negative_index_at_least_two_orbit_count": local_survivor_orbit_count,
            "not_excluded_orbit_count": 0,
            "method": (
                "A single cross edge couples at most two product factors. "
                "Record 30 retains at least three orthogonal negative factors "
                "and the coupled piece contains a negative principal block, "
                "so n_- >= 4.  In record 42 either one Lanner remains "
                "orthogonal to a coupled piece containing the other, or the "
                "two Lanners are coupled while an orthogonal dotted 2x2 "
                "polygon principal block remains; hence n_- >= 2."
            ),
            "record_negative_index_lower_bound": 4 if record_index == 30 else 2,
        },
        "survivor_orbit_representatives": survivors,
    }


def audit_all_multiwedge_one_cross_edges() -> tuple[dict[str, object], ...]:
    return tuple(audit_multiwedge_one_cross_edge(index) for index in (43, 44, 47))


def audit_all_one_cross_edges() -> tuple[dict[str, object], ...]:
    """Return the exact one-edge shell on all five frozen ordinary families."""

    return (
        audit_product_one_cross_edge(30),
        audit_product_one_cross_edge(42),
        *audit_all_multiwedge_one_cross_edges(),
    )


__all__ = [
    "ONE_EDGE_LABELS",
    "audit_all_multiwedge_one_cross_edges",
    "audit_all_one_cross_edges",
    "audit_multiwedge_one_cross_edge",
    "audit_product_one_cross_edge",
]

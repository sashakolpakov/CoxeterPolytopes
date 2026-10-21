"""Exact second-layer filters for the H9/14 candidate nerves.

The input is one of the missing-edge graphs surviving
``combinatorics.missing_edges``.  A facet of the putative simplicial
9-polytope has nine of the fourteen vertices.  Its five-vertex complement
must therefore be a vertex cover of the missing-edge graph.

This module does *not* enumerate complete nerves.  It proves two further
necessary conditions:

* the generalized lower-bound theorem forces ``g4 > 0``; and
* every actual facet complement lies in a deterministic ridge-support core.

The latter follows from the closed-pseudomanifold equation saying that, for
each six-set C, the number of selected five-set complements B contained in C
is either zero or two.  Zero is allowed: the complementary eight-set need
not be a ridge.  Iteratively deleting a candidate B that cannot have a
partner in one of its nine six-set supersets never deletes an actual facet.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations
from math import comb
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.missing_edges import (
    TARGET_POLYTOPE_DIMENSION,
    TARGET_VERTEX_COUNT,
    MissingEdgeAudit,
    analyze_missing_edge_graph,
    decode_graph6,
)


FACET_COMPLEMENT_SIZE = TARGET_VERTEX_COUNT - TARGET_POLYTOPE_DIMENSION
GALE_RANK = TARGET_VERTEX_COUNT - (TARGET_POLYTOPE_DIMENSION + 1)

EXPECTED_INPUT_GRAPH6_SHA256 = (
    "c53a9d5466b1a7ff2fda61dc2d3542d2fca0d3f2ecd3ca6ad551dd502fc77062"
)
EXPECTED_SOURCE_MANIFEST_SHA256 = (
    "55ac9ef133781dd9e6701d181538ea5f00d7a36d29a7d31711520c6b70e63ff3"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "1bc70a1b5bd0cc237f41d7d851fb90274567f4819bdd5bbcd73d033c91bd5f98"
)
EXPECTED_INPUT_COUNTS: dict[int, int] = {
    2: 2,
    3: 5,
    4: 11,
    5: 23,
    6: 50,
    7: 89,
    8: 118,
    9: 146,
}
EXPECTED_STRICT_GLBT_COUNTS: dict[int, int] = {
    2: 2,
    3: 5,
    4: 11,
    5: 21,
    6: 42,
    7: 60,
    8: 63,
    9: 59,
}
EXPECTED_RIDGE_SUPPORT_COUNTS: dict[int, int] = {
    2: 2,
    3: 5,
    4: 11,
    5: 19,
    6: 35,
    7: 46,
    8: 44,
    9: 35,
}
EXPECTED_JOINT_SURVIVOR_COUNTS: dict[int, int] = {
    2: 2,
    3: 5,
    4: 11,
    5: 17,
    6: 29,
    7: 33,
    8: 28,
    9: 17,
}


def _mask(vertices: Iterable[int]) -> int:
    value = 0
    for vertex in vertices:
        if not 0 <= vertex < TARGET_VERTEX_COUNT:
            raise ValueError("vertex outside the H9/14 ground set")
        value |= 1 << vertex
    return value


def _vertices(mask: int) -> tuple[int, ...]:
    return tuple(
        vertex
        for vertex in range(TARGET_VERTEX_COUNT)
        if mask & (1 << vertex)
    )


def _edge_masks(graph6: str) -> tuple[int, ...]:
    vertex_count, edges = decode_graph6(graph6)
    if vertex_count != TARGET_VERTEX_COUNT:
        raise ValueError("ridge-support filter is specific to 14 vertices")
    return tuple((1 << first) | (1 << second) for first, second in edges)


def candidate_facet_complement_masks(graph6: str) -> tuple[int, ...]:
    """Return every possible five-vertex complement of a nerve facet.

    If B is the complement of a nine-vertex facet, then B intersects every
    missing edge.  Conversely, every five-vertex cover is retained as a
    *candidate*; missing four- and five-faces may rule it out later.
    """

    edge_masks = _edge_masks(graph6)
    candidates = []
    for vertices in combinations(range(TARGET_VERTEX_COUNT), FACET_COMPLEMENT_SIZE):
        mask = _mask(vertices)
        if all(mask & edge_mask for edge_mask in edge_masks):
            candidates.append(mask)
    return tuple(sorted(candidates))


def ridge_support_core(
    candidates: Iterable[int],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Return the greatest iteratively ridge-supported candidate set.

    A selected complement B lies in nine six-sets ``B union {v}``.  In each,
    a closed pseudomanifold needs exactly one other selected complement.
    Thus B is impossible if any such six-set contains fewer than two current
    candidates.  Repeating this monotone deletion reaches a unique greatest
    fixed point.  The second return value lists deletion counts by round.
    """

    full_mask = (1 << TARGET_VERTEX_COUNT) - 1
    active = set(candidates)
    if any(
        mask < 0
        or mask & ~full_mask
        or mask.bit_count() != FACET_COMPLEMENT_SIZE
        for mask in active
    ):
        raise ValueError("facet complements must be distinct five-set masks")

    removed_by_round: list[int] = []
    while True:
        six_set_counts: Counter[int] = Counter()
        for complement in active:
            outside = full_mask ^ complement
            while outside:
                bit = outside & -outside
                six_set_counts[complement | bit] += 1
                outside ^= bit

        doomed: set[int] = set()
        for complement in active:
            outside = full_mask ^ complement
            while outside:
                bit = outside & -outside
                if six_set_counts[complement | bit] < 2:
                    doomed.add(complement)
                    break
                outside ^= bit
        if not doomed:
            break
        active.difference_update(doomed)
        removed_by_round.append(len(doomed))

    return tuple(sorted(active)), tuple(removed_by_round)


def independent_set_count(graph6: str, size: int) -> int:
    """Count independent sets of the requested size exactly."""

    if not 0 <= size <= TARGET_VERTEX_COUNT:
        raise ValueError("invalid independent-set size")
    edge_masks = _edge_masks(graph6)
    return sum(
        1
        for vertices in combinations(range(TARGET_VERTEX_COUNT), size)
        if not any(_mask(vertices) & edge == edge for edge in edge_masks)
    )


def first_uncovered_independent_triple(
    graph6: str, facet_complement_masks: Sequence[int]
) -> tuple[int, int, int] | None:
    """Return the first target triple not covered by any candidate facet.

    A triple T lies in the facet complementary to B exactly when ``T & B``
    is empty.  Missing edges are already nonfaces.  The target no-missing-
    triples theorem requires every independent triple to be covered.
    """

    edge_masks = _edge_masks(graph6)
    for vertices in combinations(range(TARGET_VERTEX_COUNT), 3):
        triple = _mask(vertices)
        if any(triple & edge == edge for edge in edge_masks):
            continue
        if not any(triple & complement == 0 for complement in facet_complement_masks):
            return vertices
    return None


def dehn_sommerville_f_vector(
    missing_edge_audit: MissingEdgeAudit, missing_four_count: int
) -> tuple[int, ...]:
    """Return ``(f0,...,f8)`` forced by DS for a proposed four-face count."""

    if not 0 <= missing_four_count <= missing_edge_audit.independent_four_set_count:
        raise ValueError("invalid number of missing four-faces")
    h4 = (
        missing_edge_audit.h3
        + missing_edge_audit.g4_if_no_missing_fours
        - missing_four_count
    )
    h = (
        1,
        missing_edge_audit.h1,
        missing_edge_audit.h2,
        missing_edge_audit.h3,
        h4,
        h4,
        missing_edge_audit.h3,
        missing_edge_audit.h2,
        missing_edge_audit.h1,
        1,
    )
    dimension = TARGET_POLYTOPE_DIMENSION
    return tuple(
        sum(comb(dimension - index, dimension - size) * h[index] for index in range(size + 1))
        for size in range(1, dimension + 1)
    )


def _strict_missing_four_interval(
    audit: MissingEdgeAudit,
) -> tuple[int, int] | None:
    """Apply the generalized lower-bound obstruction ``g4 != 0``.

    Murai--Nevo gives a 3-stacked triangulation when g4=0.  Since the target
    has 14>10 vertices, that triangulation is nontrivial and has an
    inclusion-minimal interior face.  No interior face has dimension at most
    9-4=5, so this missing boundary face has at least seven vertices,
    contradicting the Coxeter/Lanner bound five.
    """

    lower = audit.missing_four_count_min
    upper = min(
        audit.missing_four_count_max,
        audit.g4_if_no_missing_fours - 1,
    )
    return (lower, upper) if lower <= upper else None


def _ds_schedule(
    audit: MissingEdgeAudit, graph6: str, interval: tuple[int, int] | None
) -> tuple[dict[str, object], ...]:
    if interval is None:
        return ()
    independent_counts = {
        size: independent_set_count(graph6, size)
        for size in range(4, TARGET_POLYTOPE_DIMENSION + 1)
    }
    independent_five_count = independent_counts[5]
    schedule: list[dict[str, object]] = []
    for missing_four_count in range(interval[0], interval[1] + 1):
        f_vector = dehn_sommerville_f_vector(audit, missing_four_count)
        schedule.append(
            {
                "missing_four_count": missing_four_count,
                "g4": audit.g4_if_no_missing_fours - missing_four_count,
                "f_vector": list(f_vector),
                # If c(A4) is the number of independent five-sets containing
                # a chosen member of A4, then |A5| equals this offset-c(A4).
                "missing_five_count_offset": independent_five_count - f_vector[4],
                "within_independent_set_bounds": all(
                    0 <= f_vector[size - 1] <= independent_counts[size]
                    for size in range(4, TARGET_POLYTOPE_DIMENSION + 1)
                ),
            }
        )
    return tuple(schedule)


@dataclass(frozen=True)
class NerveSupportAudit:
    """One exact second-layer audit record."""

    graph6: str
    edge_count: int
    gale_rank: int
    strict_missing_four_count_min: int
    strict_missing_four_count_max: int
    passes_strict_generalized_lower_bound: bool
    independent_five_set_count: int
    dehn_sommerville_schedule: tuple[dict[str, object], ...]
    passes_dehn_sommerville_independent_bounds: bool
    candidate_facet_complement_count: int
    ridge_support_core_masks: tuple[int, ...]
    ridge_pruning_removed_by_round: tuple[int, ...]
    uncovered_independent_triple: tuple[int, int, int] | None
    passes_ridge_support_coverage: bool
    survives_second_layer: bool

    def to_record(self) -> dict[str, object]:
        if not self.passes_strict_generalized_lower_bound:
            disposition = "reject_g4_zero_generalized_lower_bound"
        elif not self.passes_dehn_sommerville_independent_bounds:
            disposition = "reject_dehn_sommerville_independent_bounds"
        elif not self.passes_ridge_support_coverage:
            disposition = "reject_ridge_support_uncovered_triple"
        else:
            disposition = "survives_second_layer"
        return {
            "graph6": self.graph6,
            "edges": self.edge_count,
            "gale_rank": self.gale_rank,
            "strict_missing_four_count_min": self.strict_missing_four_count_min,
            "strict_missing_four_count_max": self.strict_missing_four_count_max,
            "passes_strict_generalized_lower_bound": self.passes_strict_generalized_lower_bound,
            "independent_five_set_count": self.independent_five_set_count,
            "dehn_sommerville_schedule": list(self.dehn_sommerville_schedule),
            "passes_dehn_sommerville_independent_bounds": self.passes_dehn_sommerville_independent_bounds,
            "candidate_facet_complement_count": self.candidate_facet_complement_count,
            "ridge_support_core_count": len(self.ridge_support_core_masks),
            "ridge_support_core_masks": list(self.ridge_support_core_masks),
            "ridge_pruning_removed_by_round": list(self.ridge_pruning_removed_by_round),
            "uncovered_independent_triple": (
                list(self.uncovered_independent_triple)
                if self.uncovered_independent_triple is not None
                else None
            ),
            "passes_ridge_support_coverage": self.passes_ridge_support_coverage,
            "survives_second_layer": self.survives_second_layer,
            "disposition": disposition,
        }


def analyze_nerve_support(graph6: str) -> NerveSupportAudit:
    """Run the complete strict-GLBT and ridge-support audit for one graph."""

    missing_edge_audit = analyze_missing_edge_graph(graph6)
    if not missing_edge_audit.survives_graph_prefilter:
        raise ValueError("second layer accepts only first-layer survivors")

    interval = _strict_missing_four_interval(missing_edge_audit)
    schedule = _ds_schedule(missing_edge_audit, graph6, interval)
    # The missing-four count is not fixed at graph level.  A graph passes
    # this necessary DS layer when at least one admissible count has a forced
    # f-vector within the independent-set capacities of H.
    passes_ds_bounds = bool(schedule) and any(
        bool(entry["within_independent_set_bounds"]) for entry in schedule
    )

    candidates = candidate_facet_complement_masks(graph6)
    core, removed_by_round = ridge_support_core(candidates)
    uncovered_triple = first_uncovered_independent_triple(graph6, core)
    passes_ridge_support = uncovered_triple is None
    passes_glbt = interval is not None

    return NerveSupportAudit(
        graph6=graph6,
        edge_count=missing_edge_audit.edge_count,
        gale_rank=GALE_RANK,
        strict_missing_four_count_min=interval[0] if interval else -1,
        strict_missing_four_count_max=interval[1] if interval else -1,
        passes_strict_generalized_lower_bound=passes_glbt,
        independent_five_set_count=independent_set_count(graph6, 5),
        dehn_sommerville_schedule=schedule,
        passes_dehn_sommerville_independent_bounds=passes_ds_bounds,
        candidate_facet_complement_count=len(candidates),
        ridge_support_core_masks=core,
        ridge_pruning_removed_by_round=removed_by_round,
        uncovered_independent_triple=uncovered_triple,
        passes_ridge_support_coverage=passes_ridge_support,
        survives_second_layer=(passes_glbt and passes_ds_bounds and passes_ridge_support),
    )


@dataclass(frozen=True)
class FacetCompletionAudit:
    """Direct audit of a concrete selected facet-complement collection."""

    facet_count: int
    ridge_count: int
    all_ridge_degrees_two: bool
    dual_graph_connected: bool
    covers_all_independent_triples: bool
    minimal_nonfaces_by_size: dict[int, tuple[int, ...]]


def audit_facet_completion(
    graph6: str, facet_complement_masks: Iterable[int]
) -> FacetCompletionAudit:
    """Verify incidence and minimal nonfaces of a concrete pure 8-complex.

    The supplied distinct five-sets are complements of its nine-vertex
    facets.  This routine does not prove geometric realizability.
    """

    selected = tuple(sorted(set(facet_complement_masks)))
    if not selected:
        raise ValueError("a facet completion must contain a facet")
    candidates = set(candidate_facet_complement_masks(graph6))
    if any(mask not in candidates for mask in selected):
        raise ValueError("selected complement is not a five-vertex cover")

    full_mask = (1 << TARGET_VERTEX_COUNT) - 1
    ridge_groups: dict[int, list[int]] = {}
    for complement in selected:
        outside = full_mask ^ complement
        while outside:
            bit = outside & -outside
            ridge_groups.setdefault(complement | bit, []).append(complement)
            outside ^= bit
    all_ridge_degrees_two = all(len(group) == 2 for group in ridge_groups.values())

    adjacency = {complement: set() for complement in selected}
    for group in ridge_groups.values():
        if len(group) == 2:
            first, second = group
            adjacency[first].add(second)
            adjacency[second].add(first)
    seen = {selected[0]}
    queue = deque([selected[0]])
    while queue:
        current = queue.popleft()
        for neighbor in adjacency[current] - seen:
            seen.add(neighbor)
            queue.append(neighbor)

    uncovered = first_uncovered_independent_triple(graph6, selected)
    is_face = [
        any(subset & complement == 0 for complement in selected)
        for subset in range(1 << TARGET_VERTEX_COUNT)
    ]
    minimal_nonfaces: dict[int, list[int]] = {}
    for subset in range(1, 1 << TARGET_VERTEX_COUNT):
        if is_face[subset]:
            continue
        vertices = _vertices(subset)
        if all(is_face[subset ^ (1 << vertex)] for vertex in vertices):
            minimal_nonfaces.setdefault(len(vertices), []).append(subset)

    return FacetCompletionAudit(
        facet_count=len(selected),
        ridge_count=len(ridge_groups),
        all_ridge_degrees_two=all_ridge_degrees_two,
        dual_graph_connected=len(seen) == len(selected),
        covers_all_independent_triples=uncovered is None,
        minimal_nonfaces_by_size={
            size: tuple(values) for size, values in sorted(minimal_nonfaces.items())
        },
    )


def _edge_count_map(counter: Counter[int]) -> dict[str, int]:
    return {str(edge_count): counter.get(edge_count, 0) for edge_count in range(2, 10)}


def summarize_support_audits(
    audits: Sequence[NerveSupportAudit],
) -> dict[str, object]:
    input_counts = Counter(audit.edge_count for audit in audits)
    glbt_counts = Counter(
        audit.edge_count
        for audit in audits
        if audit.passes_strict_generalized_lower_bound
    )
    ds_counts = Counter(
        audit.edge_count
        for audit in audits
        if audit.passes_strict_generalized_lower_bound
        and audit.passes_dehn_sommerville_independent_bounds
    )
    ridge_counts = Counter(
        audit.edge_count for audit in audits if audit.passes_ridge_support_coverage
    )
    survivors = Counter(
        audit.edge_count for audit in audits if audit.survives_second_layer
    )
    exclusive_glbt = Counter(
        audit.edge_count
        for audit in audits
        if not audit.passes_strict_generalized_lower_bound
    )
    exclusive_ridge = Counter(
        audit.edge_count
        for audit in audits
        if audit.passes_strict_generalized_lower_bound
        and audit.passes_dehn_sommerville_independent_bounds
        and not audit.passes_ridge_support_coverage
    )
    return {
        "input_total": len(audits),
        "input_by_edges": _edge_count_map(input_counts),
        "strict_glbt_survivors_total": sum(glbt_counts.values()),
        "strict_glbt_survivors_by_edges": _edge_count_map(glbt_counts),
        "dehn_sommerville_bound_survivors_total": sum(ds_counts.values()),
        "dehn_sommerville_bound_survivors_by_edges": _edge_count_map(ds_counts),
        "ridge_support_coverage_survivors_total": sum(ridge_counts.values()),
        "ridge_support_coverage_survivors_by_edges": _edge_count_map(ridge_counts),
        "exclusive_glbt_rejections_total": sum(exclusive_glbt.values()),
        "exclusive_glbt_rejections_by_edges": _edge_count_map(exclusive_glbt),
        "exclusive_ridge_rejections_total": sum(exclusive_ridge.values()),
        "exclusive_ridge_rejections_by_edges": _edge_count_map(exclusive_ridge),
        "second_layer_survivors_total": sum(survivors.values()),
        "second_layer_survivors_by_edges": _edge_count_map(survivors),
    }


def _graph6_bytes(values: Sequence[str]) -> bytes:
    return "".join(f"{value}\n" for value in values).encode("ascii")


def _certificate_bytes(audits: Sequence[NerveSupportAudit]) -> bytes:
    lines = []
    for audit in audits:
        triple = (
            "-"
            if audit.uncovered_independent_triple is None
            else ",".join(map(str, audit.uncovered_independent_triple))
        )
        core = ",".join(map(str, audit.ridge_support_core_masks))
        lines.append(
            f"{audit.graph6}|{audit.strict_missing_four_count_min}:"
            f"{audit.strict_missing_four_count_max}|{core}|{triple}\n"
        )
    return "".join(lines).encode("ascii")


def _independent_macaulay_raise(value: int, degree: int) -> int:
    """Small verifier-local Macaulay implementation."""

    remaining = value
    upper: int | None = None
    raised = 0
    for lower in range(degree, 0, -1):
        top = lower - 1
        while (
            (upper is None or top + 1 < upper)
            and comb(top + 1, lower) <= remaining
        ):
            top += 1
        remaining -= comb(top, lower)
        raised += comb(top + 1, lower + 1)
        upper = top
    if remaining:
        raise AssertionError("verifier failed to expand a Macaulay integer")
    return raised


def _verify_critical_record_with_frozensets(
    record: Mapping[str, object], index: int
) -> None:
    """Independently check GLBT and ridge certificates using set objects.

    The production analyzer represents subsets by integer bitmasks.  This
    verifier deliberately regenerates five-covers and the deletion fixed
    point with ``frozenset`` unions/intersections, then derives the strict
    four-face interval directly from the low face-number formulas.
    """

    graph6 = record.get("graph6")
    if not isinstance(graph6, str):
        raise ValueError(f"record {index} lacks graph6")
    vertex_count, edge_tuples = decode_graph6(graph6)
    ordered_vertices = tuple(range(vertex_count))
    vertices = frozenset(ordered_vertices)
    edges = tuple(frozenset(edge) for edge in edge_tuples)

    degrees = Counter(vertex for edge in edges for vertex in edge)
    wedge_count = sum(comb(degrees[vertex], 2) for vertex in vertices)
    triangle_count = sum(
        1
        for triple in combinations(ordered_vertices, 3)
        if all(frozenset(pair) in edges for pair in combinations(triple, 2))
    )
    edge_count = len(edges)
    f1 = comb(vertex_count, 2) - edge_count
    f2 = (
        comb(vertex_count, 3)
        - (vertex_count - 2) * edge_count
        + wedge_count
        - triangle_count
    )
    dimension = TARGET_POLYTOPE_DIMENSION
    h1 = vertex_count - dimension
    h2 = f1 - (dimension - 1) * vertex_count + comb(dimension, 2)
    h3 = (
        f2
        - (dimension - 2) * f1
        + comb(dimension - 1, 2) * vertex_count
        - comb(dimension, 3)
    )
    g3 = h3 - h2
    independent_fours = sum(
        1
        for four in combinations(ordered_vertices, 4)
        if not any(edge <= frozenset(four) for edge in edges)
    )
    h4_if_no_missing_fours = (
        independent_fours
        - (dimension - 3) * f2
        + comb(dimension - 2, 2) * f1
        - comb(dimension - 1, 3) * vertex_count
        + comb(dimension, 4)
    )
    beta = h4_if_no_missing_fours - h3
    lower = max(0, beta - _independent_macaulay_raise(g3, 3))
    upper = min(independent_fours, beta, beta - 1)
    passes_glbt = lower <= upper
    expected_lower = lower if passes_glbt else -1
    expected_upper = upper if passes_glbt else -1
    if (
        record.get("strict_missing_four_count_min") != expected_lower
        or record.get("strict_missing_four_count_max") != expected_upper
        or record.get("passes_strict_generalized_lower_bound") != passes_glbt
    ):
        raise ValueError(f"record {index} fails independent GLBT verification")

    active = {
        frozenset(five)
        for five in combinations(ordered_vertices, FACET_COMPLEMENT_SIZE)
        if all(frozenset(five) & edge for edge in edges)
    }
    if record.get("candidate_facet_complement_count") != len(active):
        raise ValueError(f"record {index} has incorrect candidate-cover count")
    removed_by_round: list[int] = []
    while True:
        six_counts = Counter(
            complement | {vertex}
            for complement in active
            for vertex in vertices - complement
        )
        doomed = {
            complement
            for complement in active
            if any(
                six_counts[complement | {vertex}] < 2
                for vertex in vertices - complement
            )
        }
        if not doomed:
            break
        active.difference_update(doomed)
        removed_by_round.append(len(doomed))

    core_masks = sorted(
        sum(1 << vertex for vertex in complement) for complement in active
    )
    if record.get("ridge_support_core_masks") != core_masks:
        raise ValueError(f"record {index} fails independent ridge-core verification")
    if record.get("ridge_support_core_count") != len(core_masks):
        raise ValueError(f"record {index} has incorrect ridge-core count")
    if record.get("ridge_pruning_removed_by_round") != removed_by_round:
        raise ValueError(f"record {index} has incorrect ridge-deletion rounds")

    uncovered: list[int] | None = None
    for triple_tuple in combinations(ordered_vertices, 3):
        triple = frozenset(triple_tuple)
        if any(edge <= triple for edge in edges):
            continue
        if not any(complement.isdisjoint(triple) for complement in active):
            uncovered = list(triple_tuple)
            break
    if record.get("uncovered_independent_triple") != uncovered:
        raise ValueError(f"record {index} has incorrect uncovered triple")
    if record.get("passes_ridge_support_coverage") != (uncovered is None):
        raise ValueError(f"record {index} has incorrect ridge-support disposition")


def build_support_manifest(
    graph6_values: Iterable[str], source: Mapping[str, object]
) -> dict[str, object]:
    values = [value.strip() for value in graph6_values if value.strip()]
    audits = [analyze_nerve_support(value) for value in values]
    manifest = {
        "schema_version": "1.0.0",
        "object": "second-layer graph shadows, not completed nerves",
        "target": {
            "simple_polytope_dimension": TARGET_POLYTOPE_DIMENSION,
            "facet_count": TARGET_VERTEX_COUNT,
            "polar_dimension": TARGET_POLYTOPE_DIMENSION,
            "polar_vertex_count": TARGET_VERTEX_COUNT,
            "gale_rank": GALE_RANK,
        },
        "coverage_claim": (
            "Every polytopal H9/14 Coxeter nerve surviving the first layer has "
            "positive g4 and all its facet complements inside the archived "
            "ridge-support core. Rejected records cannot be target nerves."
        ),
        "nonclaims": [
            "A surviving record is not a globally selected pseudomanifold.",
            "No homology-sphere or PL-sphere certificate is supplied.",
            "No oriented-matroid or polytopality certificate is supplied.",
            "No Coxeter labeling or Gram matrix is supplied.",
        ],
        "source": dict(source),
        "input_graph6_sha256": sha256(_graph6_bytes(values)).hexdigest(),
        "certificate_sha256": sha256(_certificate_bytes(audits)).hexdigest(),
        "counts": summarize_support_audits(audits),
        "records": [audit.to_record() for audit in audits],
    }
    verify_support_manifest(manifest)
    return manifest


def verify_support_manifest(manifest: Mapping[str, object]) -> dict[str, object]:
    """Recompute every second-layer record and all frozen coverage counts."""

    if manifest.get("schema_version") != "1.0.0":
        raise ValueError("unsupported nerve-support manifest schema")
    if manifest.get("object") != "second-layer graph shadows, not completed nerves":
        raise ValueError("manifest does not identify its intermediate scope")
    target = manifest.get("target")
    expected_target = {
        "simple_polytope_dimension": TARGET_POLYTOPE_DIMENSION,
        "facet_count": TARGET_VERTEX_COUNT,
        "polar_dimension": TARGET_POLYTOPE_DIMENSION,
        "polar_vertex_count": TARGET_VERTEX_COUNT,
        "gale_rank": GALE_RANK,
    }
    if not isinstance(target, Mapping) or dict(target) != expected_target:
        raise ValueError("manifest has incorrect H9/14 target metadata")

    raw_records = manifest.get("records")
    if not isinstance(raw_records, list):
        raise ValueError("manifest records must be a list")
    graph6_values: list[str] = []
    audits: list[NerveSupportAudit] = []
    for index, raw_record in enumerate(raw_records):
        if not isinstance(raw_record, Mapping):
            raise ValueError(f"record {index} is not an object")
        graph6 = raw_record.get("graph6")
        if not isinstance(graph6, str):
            raise ValueError(f"record {index} lacks graph6")
        audit = analyze_nerve_support(graph6)
        if dict(raw_record) != audit.to_record():
            raise ValueError(f"record {index} invariants do not verify")
        _verify_critical_record_with_frozensets(raw_record, index)
        graph6_values.append(graph6)
        audits.append(audit)

    if len(set(graph6_values)) != len(graph6_values):
        raise ValueError("manifest repeats a graph")
    graph_hash = sha256(_graph6_bytes(graph6_values)).hexdigest()
    if graph_hash != EXPECTED_INPUT_GRAPH6_SHA256:
        raise ValueError("input graph list does not match the frozen first layer")
    if manifest.get("input_graph6_sha256") != graph_hash:
        raise ValueError("manifest input graph hash mismatch")
    source = manifest.get("source")
    if not isinstance(source, Mapping):
        raise ValueError("manifest lacks first-layer provenance")
    if source.get("sha256") != EXPECTED_SOURCE_MANIFEST_SHA256:
        raise ValueError("first-layer manifest hash is not the frozen source")
    certificate_hash = sha256(_certificate_bytes(audits)).hexdigest()
    if certificate_hash != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("critical certificates do not match the frozen archive")
    if manifest.get("certificate_sha256") != certificate_hash:
        raise ValueError("manifest certificate hash mismatch")

    summary = summarize_support_audits(audits)
    if manifest.get("counts") != summary:
        raise ValueError("manifest summary does not match its records")
    frozen_vectors = (
        ("input_by_edges", EXPECTED_INPUT_COUNTS),
        ("strict_glbt_survivors_by_edges", EXPECTED_STRICT_GLBT_COUNTS),
        ("dehn_sommerville_bound_survivors_by_edges", EXPECTED_STRICT_GLBT_COUNTS),
        ("ridge_support_coverage_survivors_by_edges", EXPECTED_RIDGE_SUPPORT_COUNTS),
        ("second_layer_survivors_by_edges", EXPECTED_JOINT_SURVIVOR_COUNTS),
    )
    for key, expected in frozen_vectors:
        observed = {int(edge_count): count for edge_count, count in summary[key].items()}
        if observed != expected:
            raise ValueError(f"unexpected frozen count vector for {key}")
    return summary

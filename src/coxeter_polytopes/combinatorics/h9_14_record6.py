"""Exact incidence no-go for H9/14 graph-shadow record 6.

The missing-edge graph is a triangle on three ``special`` vertices and
eleven isolates.  Facet complements split into three two-graphs ``Y_s`` on
the isolates and one ordinary edge set ``Z``.  This module classifies every
possible special-vertex link and exhausts the ten numerical ways in which
three admissible links can satisfy Dehn--Sommerville facet counting.

The conclusion concerns ordinary simplicial-sphere incidence already: no
compact hyperbolic Coxeter realization can have this graph shadow.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, combinations_with_replacement
from math import comb
from typing import Iterable, Sequence

from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


TARGET_RECORD_INDEX = 6
TARGET_GRAPH6 = "M???????????_?_?_"
TARGET_EDGES = ((0, 12), (0, 13), (12, 13))

LOCAL_ORDER = 11
LOCAL_VERTICES = tuple(range(LOCAL_ORDER))
LOCAL_FULL_MASK = (1 << LOCAL_ORDER) - 1
LOCAL_EDGES = tuple(combinations(LOCAL_VERTICES, 2))
LOCAL_TRIPLES = tuple(combinations(LOCAL_VERTICES, 3))
TRIPLE_INDEX = {triple: index for index, triple in enumerate(LOCAL_TRIPLES)}

EXPECTED_TRIANGLE_FREE_GRAPH_COUNT = 12_172
EXPECTED_RIDGE_GRAPH_COUNT = 308
EXPECTED_LINK_GRAPH_COUNT = 42
EXPECTED_FULL_LINK_CLASS_COUNT = 12

# The exact output of the necessary link audit on ``geng -q -t 10``.
EXPECTED_LINK_GRAPH6 = (
    "I???FB_w?", "I???FB_}?", "I???FB_FW", "I???FBo{?",
    "I???FBo}?", "I???FBw}?", "I???F?\\Fo", "I???Fbo{?",
    "I???Fbo}?", "I???FboBw", "I???Frw}?", "I??EEBo{?",
    "I??EEBo}?", "I??EEBoBw", "I??EEBw}?", "I??EFB_}?",
    "I??EFBo{?", "I??EFBo}?", "I??EFBoBw", "I??EFBw}?",
    "I??EFbo}?", "I??EFboBw", "I??EFrw}?", "I??E@cypW",
    "I??E@syN_", "I??E@sypW", "I??FFB_}?", "I??FFBo{?",
    "I??FFBo}?", "I??FFBw}?", "I??FF?^Fo", "I??Ffbo{?",
    "I??Ffbo}?", "I??FfboBw", "I?AEFBo}?", "I?AEFbo}?",
    "I?AEFboBw", "I?AFF?^Fo", "I?AFfbo}?", "I?AFfboBw",
    "I?BEFbo}?", "I?BfFBw}?",
)

EXPECTED_LINK_INVARIANT_HISTOGRAM = (
    ((45, 2, 1), 2),
    ((48, 1, 1), 2),
    ((49, 1, 2), 10),
    ((50, 1, 3), 4),
    ((52, 0, 2), 6),
    ((53, 0, 3), 12),
    ((54, 0, 4), 5),
    ((55, 0, 5), 1),
)

ADMISSIBLE_LINK_SIGNATURES = tuple(
    key for key, _multiplicity in EXPECTED_LINK_INVARIANT_HISTOGRAM
)

# (three link signatures, isolate-only A4 count, Z-edge count).
EXPECTED_NUMERICAL_BRANCHES = (
    (((45, 2, 1), (45, 2, 1), (45, 2, 1)), 0, 3),
    (((45, 2, 1), (45, 2, 1), (45, 2, 1)), 1, 1),
    (((45, 2, 1), (45, 2, 1), (48, 1, 1)), 0, 2),
    (((45, 2, 1), (45, 2, 1), (48, 1, 1)), 1, 0),
    (((45, 2, 1), (45, 2, 1), (49, 1, 2)), 0, 1),
    (((45, 2, 1), (45, 2, 1), (50, 1, 3)), 0, 0),
    (((45, 2, 1), (45, 2, 1), (52, 0, 2)), 0, 0),
    (((45, 2, 1), (48, 1, 1), (48, 1, 1)), 0, 1),
    (((45, 2, 1), (48, 1, 1), (49, 1, 2)), 0, 0),
    (((48, 1, 1), (48, 1, 1), (48, 1, 1)), 0, 0),
)


def vertex_mask(vertices: Iterable[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


TRIPLE_VERTEX_MASKS = tuple(vertex_mask(triple) for triple in LOCAL_TRIPLES)
FOUR_TRIPLE_INDICES = tuple(
    tuple(
        index
        for index, triple_mask in enumerate(TRIPLE_VERTEX_MASKS)
        if triple_mask & vertex_mask(four) == triple_mask
    )
    for four in combinations(LOCAL_VERTICES, 4)
)
DISJOINT_TRIPLES = tuple(
    sum(
        1 << index
        for index, triple_mask in enumerate(TRIPLE_VERTEX_MASKS)
        if not subset & triple_mask
    )
    for subset in range(1 << LOCAL_ORDER)
)
PAIR_VERTEX_MASKS = tuple(vertex_mask(edge) for edge in LOCAL_EDGES)
TRIPLE_ADJACENCY = tuple(
    sum(
        1 << other
        for other, other_mask in enumerate(TRIPLE_VERTEX_MASKS)
        if other != index and (triple_mask | other_mask).bit_count() == 4
    )
    for index, triple_mask in enumerate(TRIPLE_VERTEX_MASKS)
)
EDGE_COBoundaries = tuple(
    sum(
        1 << index
        for index, triple in enumerate(LOCAL_TRIPLES)
        if first in triple and second in triple
    )
    for first, second in LOCAL_EDGES
)


def _macaulay_raise(value: int, degree: int) -> int:
    remaining = value
    next_top: int | None = None
    result = 0
    for lower in range(degree, 0, -1):
        top = lower - 1
        while (
            (next_top is None or top + 1 < next_top)
            and comb(top + 1, lower) <= remaining
        ):
            top += 1
        remaining -= comb(top, lower)
        result += comb(top + 1, lower + 1)
        next_top = top
    if remaining:
        raise AssertionError("incomplete Macaulay expansion")
    return result


def graph_two_graph(graph6: str) -> int:
    """Return ``delta X`` on vertices 0,...,10, with 0 the switch base."""

    order, raw_edges = decode_graph6(graph6)
    if order != 10:
        raise ValueError("a record-6 link switching graph has order ten")
    edges = {
        tuple(sorted((first + 1, second + 1)))
        for first, second in raw_edges
    }
    return sum(
        1 << index
        for index, triple in enumerate(LOCAL_TRIPLES)
        if sum(tuple(sorted(edge)) in edges for edge in combinations(triple, 2)) % 2
    )


def _switch_graph_edges(graph6: str) -> tuple[tuple[int, int], ...]:
    order, edges = decode_graph6(graph6)
    if order != 10:
        raise ValueError("the switching graph must have order ten")
    return edges


def ridge_graph_admissible(graph6: str) -> bool:
    """Check triangle-free and induced-2K2-free directly."""

    edges = set(_switch_graph_edges(graph6))
    for triple in combinations(range(10), 3):
        if all(tuple(sorted(edge)) in edges for edge in combinations(triple, 2)):
            return False
    for four in combinations(range(10), 4):
        induced = [edge for edge in edges if set(edge) <= set(four)]
        if len(induced) == 2 and set(induced[0]).isdisjoint(induced[1]):
            return False
    return True


def link_face_table(two_graph: int) -> tuple[bool, ...]:
    return tuple(
        subset == 0 or bool(two_graph & DISJOINT_TRIPLES[subset])
        for subset in range(1 << LOCAL_ORDER)
    )


def link_minimal_nonfaces(two_graph: int) -> tuple[int, ...]:
    table = link_face_table(two_graph)
    return tuple(
        subset
        for subset in range(1, 1 << LOCAL_ORDER)
        if not table[subset]
        and all(
            table[subset ^ (1 << vertex)]
            for vertex in LOCAL_VERTICES
            if subset & (1 << vertex)
        )
    )


def _dual_connected(two_graph: int) -> bool:
    selected = [index for index in range(len(LOCAL_TRIPLES)) if two_graph >> index & 1]
    if not selected:
        return False
    seen = 1 << selected[0]
    stack = [selected[0]]
    while stack:
        index = stack.pop()
        fresh = TRIPLE_ADJACENCY[index] & two_graph & ~seen
        while fresh:
            bit = fresh & -fresh
            other = bit.bit_length() - 1
            seen |= bit
            stack.append(other)
            fresh ^= bit
    return seen.bit_count() == len(selected)


@lru_cache(maxsize=None)
def admissible_link_signature(two_graph: int) -> tuple[int, int, int] | None:
    """Apply every necessary special-link test used by the finite census.

    The result is ``(facet_count, missing_triangle_count, g4)``.  The input
    is already a two-graph, but the stronger 0/2 ridge rule is checked.
    """

    if any(
        sum(bool(two_graph >> index & 1) for index in indices) not in (0, 2)
        for indices in FOUR_TRIPLE_INDICES
    ):
        return None
    if any(not two_graph & DISJOINT_TRIPLES[pair] for pair in PAIR_VERTEX_MASKS):
        return None
    table = link_face_table(two_graph)
    nonfaces = []
    for subset in range(1, 1 << LOCAL_ORDER):
        if table[subset]:
            continue
        if all(
            table[subset ^ (1 << vertex)]
            for vertex in LOCAL_VERTICES
            if subset & (1 << vertex)
        ):
            if subset.bit_count() > 5:
                return None
            nonfaces.append(subset)
    if not nonfaces or not _dual_connected(two_graph):
        return None
    missing_triangles = sum(value.bit_count() == 3 for value in nonfaces)
    four_faces = sum(
        table[subset]
        for subset in range(1 << LOCAL_ORDER)
        if subset.bit_count() == 4
    )
    g4 = four_faces - 325 + 6 * missing_triangles
    g3 = 4 - missing_triangles
    if not (
        0 <= missing_triangles <= 3
        and 1 <= g4 <= _macaulay_raise(g3, 3)
    ):
        return None
    return two_graph.bit_count(), missing_triangles, g4


@dataclass(frozen=True)
class LinkCensusAudit:
    input_count: int
    ridge_graph_count: int
    survivor_graph6: tuple[str, ...]
    invariant_histogram: tuple[tuple[tuple[int, int, int], int], ...]

    @property
    def verified(self) -> bool:
        return (
            self.input_count == EXPECTED_TRIANGLE_FREE_GRAPH_COUNT
            and self.ridge_graph_count == EXPECTED_RIDGE_GRAPH_COUNT
            and self.survivor_graph6 == EXPECTED_LINK_GRAPH6
            and self.invariant_histogram == EXPECTED_LINK_INVARIANT_HISTOGRAM
        )

    def to_record(self) -> dict[str, object]:
        return {
            "input_triangle_free_unlabelled_graphs": self.input_count,
            "triangle_and_induced_2k2_free_types": self.ridge_graph_count,
            "audited_fixed_base_link_types": len(self.survivor_graph6),
            "survivor_graph6": list(self.survivor_graph6),
            "invariant_histogram": [
                {
                    "facet_count": key[0],
                    "missing_triangle_count": key[1],
                    "g4": key[2],
                    "multiplicity": count,
                }
                for key, count in self.invariant_histogram
            ],
        }


def audit_link_census(graph6_words: Sequence[str]) -> LinkCensusAudit:
    survivors: list[str] = []
    signatures: list[tuple[int, int, int]] = []
    ridge_count = 0
    for word in graph6_words:
        if not ridge_graph_admissible(word):
            continue
        ridge_count += 1
        signature = admissible_link_signature(graph_two_graph(word))
        if signature is not None:
            survivors.append(word)
            signatures.append(signature)
    return LinkCensusAudit(
        input_count=len(graph6_words),
        ridge_graph_count=ridge_count,
        survivor_graph6=tuple(survivors),
        invariant_histogram=tuple(sorted(Counter(signatures).items())),
    )


def partition_two_graph(blocks: Sequence[Iterable[int]]) -> int:
    block_sets = tuple(frozenset(block) for block in blocks)
    return sum(
        1 << index
        for index, triple in enumerate(LOCAL_TRIPLES)
        if all(any(vertex in block for vertex in triple) for block in block_sets)
    )


def e_link_orbit() -> tuple[int, ...]:
    """All 5,775 labelled 3+4+4 transversal two-graphs."""

    result = []
    for small in combinations(LOCAL_VERTICES, 3):
        remaining = tuple(vertex for vertex in LOCAL_VERTICES if vertex not in small)
        for first in combinations(remaining, 4):
            second = tuple(vertex for vertex in remaining if vertex not in first)
            if first <= second:
                result.append(partition_two_graph((small, first, second)))
    if len(result) != 5_775 or len(set(result)) != 5_775:
        raise AssertionError("E-link labelled orbit count changed")
    return tuple(result)


def h_link_orbit() -> tuple[int, ...]:
    """All 4,620 labelled 3+3+5 transversal two-graphs."""

    result = []
    for large in combinations(LOCAL_VERTICES, 5):
        remaining = tuple(vertex for vertex in LOCAL_VERTICES if vertex not in large)
        for first in combinations(remaining, 3):
            second = tuple(vertex for vertex in remaining if vertex not in first)
            if first <= second:
                result.append(partition_two_graph((first, second, large)))
    if len(result) != 4_620 or len(set(result)) != 4_620:
        raise AssertionError("H-link labelled orbit count changed")
    return tuple(result)


E_CANONICAL = partition_two_graph((range(3), range(3, 7), range(7, 11)))
H_CANONICAL = partition_two_graph((range(3), range(3, 6), range(6, 11)))


def coboundary(edge_indices: Iterable[int]) -> int:
    result = 0
    for index in edge_indices:
        result ^= EDGE_COBoundaries[index]
    return result


def derive_numerical_branches() -> tuple[
    tuple[tuple[tuple[int, int, int], ...], int, int], ...
]:
    """Derive, rather than assume, the ten facet-count branches.

    For a link signature ``(y,q,r)``, ``y`` is its facet count and ``q``
    counts global A4s through that special vertex.  If ``a0`` is the number
    of isolate-only A4s, global Dehn--Sommerville gives

    ``z = 150 - 2*(sum(q)+a0) - sum(y)``.

    We enumerate every unordered signature triple and every nonnegative
    ``a0`` for which this number of all-special complement edges is
    nonnegative.
    """

    result = []
    for signatures in combinations_with_replacement(ADMISSIBLE_LINK_SIGNATURES, 3):
        total_y = sum(signature[0] for signature in signatures)
        total_q = sum(signature[1] for signature in signatures)
        a0 = 0
        while True:
            z_size = 150 - 2 * (total_q + a0) - total_y
            if z_size < 0:
                break
            result.append((tuple(sorted(signatures)), a0, z_size))
            a0 += 1
    return tuple(sorted(result))


def missing_isolate_faces(
    ys: Sequence[int], z_edge_indices: Sequence[int], size: int
) -> tuple[tuple[int, ...], ...]:
    """Return isolate-only nonfaces of one fixed size.

    When ``size==3`` these violate the no-missing-triples theorem.  Once the
    size-three result is empty, the size-four result consists precisely of
    the isolate-only minimal four-nonfaces.
    """

    union = 0
    for two_graph in ys:
        union |= two_graph
    z_edges = tuple(LOCAL_EDGES[index] for index in z_edge_indices)
    result = []
    for face in combinations(LOCAL_VERTICES, size):
        face_mask = vertex_mask(face)
        y_cover = bool(union & DISJOINT_TRIPLES[face_mask])
        z_cover = any(set(face).isdisjoint(edge) for edge in z_edges)
        if not y_cover and not z_cover:
            result.append(face)
    return tuple(result)


@dataclass(frozen=True)
class BranchAudit:
    name: str
    link_types: str
    missing_four_count: int
    isolate_missing_four_count: int
    z_size: int
    target_signature: tuple[int, int, int]
    facet_count_candidates: int
    link_compatible_candidates: int
    independent_triple_compatible_candidates: int
    expected_a0_compatible_candidates: int
    missing_triple_histogram: tuple[tuple[int, int], ...]
    isolate_four_histogram: tuple[tuple[int, int], ...]
    first_missing_triple_witness: tuple[int, ...] | None
    first_isolate_four_witness: tuple[int, ...] | None

    @property
    def excluded(self) -> bool:
        return self.expected_a0_compatible_candidates == 0

    def to_record(self) -> dict[str, object]:
        return {
            "name": self.name,
            "link_types": self.link_types,
            "missing_four_count": self.missing_four_count,
            "isolate_missing_four_count": self.isolate_missing_four_count,
            "z_size": self.z_size,
            "target_signature": list(self.target_signature),
            "facet_count_candidates": self.facet_count_candidates,
            "link_compatible_candidates": self.link_compatible_candidates,
            "independent_triple_compatible_candidates": (
                self.independent_triple_compatible_candidates
            ),
            "expected_a0_compatible_candidates": (
                self.expected_a0_compatible_candidates
            ),
            "missing_triple_histogram": {
                str(key): value for key, value in self.missing_triple_histogram
            },
            "isolate_four_histogram": {
                str(key): value for key, value in self.isolate_four_histogram
            },
            "first_missing_triple_witness": (
                list(self.first_missing_triple_witness)
                if self.first_missing_triple_witness is not None
                else None
            ),
            "first_isolate_four_witness": (
                list(self.first_isolate_four_witness)
                if self.first_isolate_four_witness is not None
                else None
            ),
            "excluded": self.excluded,
        }


def _scan_derived_branch(
    *,
    name: str,
    link_types: str,
    missing_four_count: int,
    expected_a0: int,
    base: int,
    orbit: Sequence[int],
    z_supports: Sequence[tuple[int, ...]],
    target_signature: tuple[int, int, int],
) -> BranchAudit:
    facet_count_candidates = 0
    link_candidates: list[tuple[int, int, tuple[int, ...]]] = []
    for second in orbit:
        for z_edges in z_supports:
            third = base ^ second ^ coboundary(z_edges)
            if third.bit_count() != target_signature[0]:
                continue
            facet_count_candidates += 1
            if admissible_link_signature(third) == target_signature:
                link_candidates.append((second, third, z_edges))

    histogram: Counter[int] = Counter()
    four_histogram: Counter[int] = Counter()
    independent = []
    first_missing_triple: tuple[int, ...] | None = None
    first_isolate_four: tuple[int, ...] | None = None
    expected = 0
    for second, third, z_edges in link_candidates:
        missing_triples = missing_isolate_faces((base, second, third), z_edges, 3)
        histogram[len(missing_triples)] += 1
        if missing_triples:
            if first_missing_triple is None:
                first_missing_triple = missing_triples[0]
            continue
        independent.append((second, third, z_edges))
        missing_fours = missing_isolate_faces((base, second, third), z_edges, 4)
        four_histogram[len(missing_fours)] += 1
        if len(missing_fours) == expected_a0:
            expected += 1
        elif first_isolate_four is None and missing_fours:
            first_isolate_four = missing_fours[0]
    return BranchAudit(
        name=name,
        link_types=link_types,
        missing_four_count=missing_four_count,
        isolate_missing_four_count=expected_a0,
        z_size=len(z_supports[0]) if z_supports else 0,
        target_signature=target_signature,
        facet_count_candidates=facet_count_candidates,
        link_compatible_candidates=len(link_candidates),
        independent_triple_compatible_candidates=len(independent),
        expected_a0_compatible_candidates=expected,
        missing_triple_histogram=tuple(sorted(histogram.items())),
        isolate_four_histogram=tuple(sorted(four_histogram.items())),
        first_missing_triple_witness=first_missing_triple,
        first_isolate_four_witness=first_isolate_four,
    )


def _hhh_z3_branch(h_orbit: Sequence[int]) -> BranchAudit:
    delta_three = {
        coboundary(indices): indices
        for indices in combinations(range(len(LOCAL_EDGES)), 3)
    }
    if len(delta_three) != comb(55, 3):
        raise AssertionError("three-edge coboundaries unexpectedly collided")
    unordered = []
    for first_index, first in enumerate(h_orbit):
        base_xor = H_CANONICAL ^ first
        for second_index in range(first_index, len(h_orbit)):
            z_edges = delta_three.get(base_xor ^ h_orbit[second_index])
            if z_edges is not None:
                unordered.append((first_index, second_index, z_edges))
    ordered_count = sum(
        1 if first == second else 2 for first, second, _ in unordered
    )
    histogram: Counter[int] = Counter()
    first_missing_triple = None
    first_isolate_four = None
    four_histogram: Counter[int] = Counter()
    independent_ordered = 0
    expected_ordered = 0
    for first_index, second_index, z_edges in unordered:
        multiplicity = 1 if first_index == second_index else 2
        ys = (H_CANONICAL, h_orbit[first_index], h_orbit[second_index])
        missing_triples = missing_isolate_faces(ys, z_edges, 3)
        histogram[len(missing_triples)] += multiplicity
        if missing_triples:
            if first_missing_triple is None:
                first_missing_triple = missing_triples[0]
            continue
        independent_ordered += multiplicity
        missing_fours = missing_isolate_faces(ys, z_edges, 4)
        four_histogram[len(missing_fours)] += multiplicity
        if not missing_fours:
            expected_ordered += multiplicity
        elif first_isolate_four is None:
            first_isolate_four = missing_fours[0]
    return BranchAudit(
        name="HHH_z3_a0_0",
        link_types="H,H,H",
        missing_four_count=6,
        isolate_missing_four_count=0,
        z_size=3,
        target_signature=(45, 2, 1),
        facet_count_candidates=ordered_count,
        link_compatible_candidates=ordered_count,
        independent_triple_compatible_candidates=independent_ordered,
        expected_a0_compatible_candidates=expected_ordered,
        missing_triple_histogram=tuple(sorted(histogram.items())),
        isolate_four_histogram=tuple(sorted(four_histogram.items())),
        first_missing_triple_witness=first_missing_triple,
        first_isolate_four_witness=first_isolate_four,
    )


def audit_ten_branches() -> tuple[BranchAudit, ...]:
    """Exhaust every relative labelled alignment in the ten DS branches.

    Common isolate relabelling sends one E or H link to the displayed
    canonical partition.  ``e_link_orbit``/``h_link_orbit`` then enumerate
    *every* labelled copy of the second link, and the Z lists enumerate every
    labelled edge support.  Hence these are full relative-orbit searches,
    not samples of canonical graph labels.
    """

    if derive_numerical_branches() != EXPECTED_NUMERICAL_BRANCHES:
        raise AssertionError("record-6 numerical branch derivation changed")
    e_orbit = e_link_orbit()
    h_orbit = h_link_orbit()
    e_set = set(e_orbit)
    h_set = set(h_orbit)
    e_survivors = {
        graph_two_graph(word)
        for word in EXPECTED_LINK_GRAPH6
        if admissible_link_signature(graph_two_graph(word)) == (48, 1, 1)
    }
    h_survivors = {
        graph_two_graph(word)
        for word in EXPECTED_LINK_GRAPH6
        if admissible_link_signature(graph_two_graph(word)) == (45, 2, 1)
    }
    if len(e_survivors) != 2 or not e_survivors <= e_set:
        raise AssertionError("E-link partition orbit is not exhaustive")
    if len(h_survivors) != 2 or not h_survivors <= h_set:
        raise AssertionError("H-link partition orbit is not exhaustive")
    z0 = ((),)
    z1 = tuple((index,) for index in range(55))
    z2 = tuple(combinations(range(55), 2))
    branches = (
        _scan_derived_branch(
            name="EEE_z0_a0_0", link_types="E,E,E", missing_four_count=3,
            expected_a0=0, base=E_CANONICAL, orbit=e_orbit, z_supports=z0,
            target_signature=(48, 1, 1),
        ),
        _scan_derived_branch(
            name="AHH_z0_a0_0", link_types="H,H,A", missing_four_count=4,
            expected_a0=0, base=H_CANONICAL, orbit=h_orbit, z_supports=z0,
            target_signature=(52, 0, 2),
        ),
        _scan_derived_branch(
            name="EEH_z1_a0_0", link_types="H,E,E", missing_four_count=4,
            expected_a0=0, base=H_CANONICAL, orbit=e_orbit, z_supports=z1,
            target_signature=(48, 1, 1),
        ),
        _scan_derived_branch(
            name="EFH_z0_a0_0", link_types="H,E,F", missing_four_count=4,
            expected_a0=0, base=H_CANONICAL, orbit=e_orbit, z_supports=z0,
            target_signature=(49, 1, 2),
        ),
        _scan_derived_branch(
            name="EHH_z2_a0_0", link_types="H,H,E", missing_four_count=5,
            expected_a0=0, base=H_CANONICAL, orbit=h_orbit, z_supports=z2,
            target_signature=(48, 1, 1),
        ),
        _scan_derived_branch(
            name="FHH_z1_a0_0", link_types="H,H,F", missing_four_count=5,
            expected_a0=0, base=H_CANONICAL, orbit=h_orbit, z_supports=z1,
            target_signature=(49, 1, 2),
        ),
        _scan_derived_branch(
            name="GHH_z0_a0_0", link_types="H,H,G", missing_four_count=5,
            expected_a0=0, base=H_CANONICAL, orbit=h_orbit, z_supports=z0,
            target_signature=(50, 1, 3),
        ),
        _hhh_z3_branch(h_orbit),
        _scan_derived_branch(
            name="EHH_z0_a0_1", link_types="H,H,E", missing_four_count=6,
            expected_a0=1, base=H_CANONICAL, orbit=h_orbit, z_supports=z0,
            target_signature=(48, 1, 1),
        ),
        _scan_derived_branch(
            name="HHH_z1_a0_1", link_types="H,H,H", missing_four_count=7,
            expected_a0=1, base=H_CANONICAL, orbit=h_orbit, z_supports=z1,
            target_signature=(45, 2, 1),
        ),
    )
    if not all(branch.excluded for branch in branches):
        raise AssertionError("record-6 branch census found an unexcluded case")
    return branches


def verify_target_graph() -> bool:
    return decode_graph6(TARGET_GRAPH6) == (14, TARGET_EDGES)

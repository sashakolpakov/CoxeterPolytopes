"""Exact whole-shadow reduction for H9/14 second-layer record 17.

The missing-edge graph has edges ``0-12, 0-13, 1-12, 2-13, 12-13``
and nine isolates.  Blocker duality and the complete ridge equations force
two independent ``K5,5`` facet families.  The remaining equations are an
exact finite problem consisting of a triangle-free graph on the nine
isolates and two five-subsets.  This module exhausts that problem without a
SAT/SMT backend.

The eight surviving incidence-complex orbits are only *candidates*: passing
the closed-pseudomanifold and mod-two homology-link audits does not prove that
they are PL, convex, or nerves of hyperbolic Coxeter polytopes.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.local_link import analyze_local_link_graph
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


TARGET_RECORD_INDEX = 17
TARGET_GRAPH6 = "M???????????o?g?_"
TARGET_VERTEX_COUNT = 14
TARGET_DIMENSION = 9
FULL_MASK = (1 << TARGET_VERTEX_COUNT) - 1
MAX_LANNER_ORDER = 5

SPECIAL_VERTICES = (12, 13)
WORK_VERTICES = tuple(range(12))
ACTIVE_VERTICES = (0, 1, 2)
ISOLATED_VERTICES = tuple(range(3, 12))
LOCAL_ISOLATES = tuple(range(9))
P_BASE = (0, 1, 13)
P_GROUND = (2,) + ISOLATED_VERTICES
Q_BASE = (0, 2, 12)
Q_GROUND = (1,) + ISOLATED_VERTICES

TARGET_EDGES = ((0, 12), (1, 12), (0, 13), (2, 13), (12, 13))
EDGE_MASKS = tuple(sum(1 << vertex for vertex in edge) for edge in TARGET_EDGES)
MINIMAL_GRAPH_COVERS = ((12, 13), P_BASE, Q_BASE)

EXPECTED_STREAM_COUNT = 1_897
EXPECTED_ADMISSIBLE_X_COUNT = 156
EXPECTED_COUPLED_STATE_COUNT = 11_545
EXPECTED_STAGE_COUNTS = {
    "target_y_count": {0: 1_434, 1: 1_376, 2: 1_455, 3: 1_869, 4: 521},
    "independent_triples": {0: 654, 1: 584, 2: 159, 3: 136, 4: 21},
    "exact_a4": {0: 36, 1: 82, 2: 76, 3: 53, 4: 19},
    "exact_a5_and_support": {0: 36, 1: 82, 2: 76, 3: 53, 4: 19},
    "full_blocker": {0: 1, 1: 17, 2: 3, 3: 0, 4: 0},
    "full_face_vector": {0: 1, 1: 17, 2: 3, 3: 0, 4: 0},
}

UNIQUE_X_GRAPH6 = "H??EEBo"
UNIQUE_P_LOCAL_MASK = sum(1 << vertex for vertex in (4, 5, 6, 7, 8))
UNIQUE_Q_LOCAL_MASK = UNIQUE_P_LOCAL_MASK
UNIQUE_A5_MASKS = (121, 122, 124, 481, 1921, 2078, 3590, 3968)
EXPECTED_RESIDUAL_ORBITS = (
    (1, "H???FB_", 31, 6),
    (1, "H???FB_", 440, 4),
    (2, "H???F?^", 248, 1),
    (2, "H???Fbo", 496, 1),
    (0, "H??EEBo", 496, 1),
    (1, "H??EEBw", 31, 1),
    (1, "H??FFB_", 31, 6),
    (2, "H??FF?^", 248, 1),
)


def mask(vertices: Iterable[int]) -> int:
    result = 0
    for vertex in vertices:
        result |= 1 << vertex
    return result


def vertices(value: int, order: int = TARGET_VERTEX_COUNT) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(order) if value & (1 << vertex))


def _is_independent(value: int) -> bool:
    return not any(value & edge == edge for edge in EDGE_MASKS)


def _minimal_vertex_covers() -> tuple[int, ...]:
    result = []
    for value in range(1 << TARGET_VERTEX_COUNT):
        if not all(value & edge for edge in EDGE_MASKS):
            continue
        if all(
            not all((value ^ bit) & edge for edge in EDGE_MASKS)
            for bit in (1 << vertex for vertex in vertices(value))
        ):
            result.append(value)
    return tuple(result)


@dataclass(frozen=True)
class ShadowAudit:
    decoded_edges: tuple[tuple[int, int], ...]
    minimal_vertex_cover_masks: tuple[int, ...]
    independent_counts_3_4_5: tuple[int, int, int]

    @property
    def verified(self) -> bool:
        return (
            self.decoded_edges == TARGET_EDGES
            and self.minimal_vertex_cover_masks
            == tuple(sorted(mask(cover) for cover in MINIMAL_GRAPH_COVERS))
            and self.independent_counts_3_4_5 == (310, 735, 1_212)
        )

    def to_record(self) -> dict[str, object]:
        return {
            "decoded_edges": [list(edge) for edge in self.decoded_edges],
            "minimal_vertex_covers": [
                list(vertices(value)) for value in self.minimal_vertex_cover_masks
            ],
            "independent_counts_3_4_5": list(self.independent_counts_3_4_5),
            "verified": self.verified,
        }


def audit_shadow() -> ShadowAudit:
    order, decoded = decode_graph6(TARGET_GRAPH6)
    if order != TARGET_VERTEX_COUNT:
        raise AssertionError("record-17 graph order changed")
    audit = ShadowAudit(
        decoded_edges=decoded,
        minimal_vertex_cover_masks=_minimal_vertex_covers(),
        independent_counts_3_4_5=tuple(
            sum(
                _is_independent(mask(chosen))
                for chosen in combinations(range(TARGET_VERTEX_COUNT), size)
            )
            for size in (3, 4, 5)
        ),
    )
    if not audit.verified:
        raise AssertionError("record-17 graph-shadow audit changed")
    return audit


@dataclass(frozen=True)
class SupportAudit:
    special_missing_four_intervals: tuple[tuple[int, int], tuple[int, int]]
    four_on_work_count: int
    allowed_four_count: int
    five_on_work_count: int
    every_forbidden_four_extends_at_a_special: bool

    @property
    def verified(self) -> bool:
        return (
            self.special_missing_four_intervals == ((0, 0), (0, 0))
            and self.four_on_work_count == 495
            and self.allowed_four_count == 201
            and self.five_on_work_count == 792
            and self.every_forbidden_four_extends_at_a_special
        )

    def to_record(self) -> dict[str, object]:
        return {
            "special_missing_four_intervals": [
                list(interval) for interval in self.special_missing_four_intervals
            ],
            "four_on_work_count": self.four_on_work_count,
            "allowed_four_count": self.allowed_four_count,
            "five_on_work_count": self.five_on_work_count,
            "every_forbidden_four_extends_at_a_special": (
                self.every_forbidden_four_extends_at_a_special
            ),
            "support_conclusion": (
                "A4 lies on W and contains 0 or both 1,2; extra A5 lies on W"
            ),
            "verified": self.verified,
        }


def _allowed_a4(value: int) -> bool:
    if value & ~mask(WORK_VERTICES):
        return False
    return bool(value & 1) or value & mask((1, 2)) == mask((1, 2))


def audit_support() -> SupportAudit:
    local = analyze_local_link_graph(TARGET_GRAPH6)
    intervals = tuple(
        local.vertex_figures[vertex].missing_four_incidence_interval
        for vertex in SPECIAL_VERTICES
    )
    work_fours = tuple(mask(chosen) for chosen in combinations(WORK_VERTICES, 4))
    forbidden = tuple(value for value in work_fours if not _allowed_a4(value))
    every_extends = all(
        (
            (value & mask((0, 1)))
            or _is_independent(value | (1 << 12))
        )
        and (
            (value & mask((0, 2)))
            or _is_independent(value | (1 << 13))
        )
        for value in forbidden
    )
    audit = SupportAudit(
        special_missing_four_intervals=intervals,  # type: ignore[arg-type]
        four_on_work_count=len(work_fours),
        allowed_four_count=sum(_allowed_a4(value) for value in work_fours),
        five_on_work_count=sum(1 for _ in combinations(WORK_VERTICES, 5)),
        every_forbidden_four_extends_at_a_special=every_extends,
    )
    if not audit.verified:
        raise AssertionError("record-17 support audit changed")
    return audit


def _graph_cover_five_sets() -> tuple[int, ...]:
    return tuple(
        value
        for chosen in combinations(range(TARGET_VERTEX_COUNT), 5)
        if all((value := mask(chosen)) & edge for edge in EDGE_MASKS)
    )


@dataclass(frozen=True)
class FacetTypeAudit:
    graph_cover_count: int
    y_count: int
    p_count: int
    q_count: int
    every_cover_classified_once: bool
    p_and_q_blocker_part_sizes: tuple[int, int]

    @property
    def verified(self) -> bool:
        return (
            (self.graph_cover_count, self.y_count, self.p_count, self.q_count)
            == (310, 220, 45, 45)
            and self.every_cover_classified_once
            and self.p_and_q_blocker_part_sizes == (5, 5)
        )

    def to_record(self) -> dict[str, object]:
        return {
            "graph_cover_count": self.graph_cover_count,
            "facet_complement_types": {
                "y_{12,13,T}": self.y_count,
                "p_{0,1,13,S}": self.p_count,
                "q_{0,2,12,S}": self.q_count,
            },
            "every_cover_classified_once": self.every_cover_classified_once,
            "blocker_cut_conclusion": (
                "each p/q graph is K5,5; its two blocker sides are A5 rows"
            ),
            "p_and_q_blocker_part_sizes": list(self.p_and_q_blocker_part_sizes),
            "verified": self.verified,
        }


def audit_facet_types() -> FacetTypeAudit:
    covers = _graph_cover_five_sets()
    special = mask(SPECIAL_VERTICES)
    p_base = mask(P_BASE)
    q_base = mask(Q_BASE)
    p_ground = mask(P_GROUND)
    q_ground = mask(Q_GROUND)
    counts = Counter()
    classified = True
    for value in covers:
        types = []
        if value & special == special and (value ^ special).bit_count() == 3:
            types.append("y")
        if value & p_base == p_base and (value ^ p_base) & ~p_ground == 0:
            types.append("p")
        if value & q_base == q_base and (value ^ q_base) & ~q_ground == 0:
            types.append("q")
        classified &= len(types) == 1
        counts.update(types)
    audit = FacetTypeAudit(
        graph_cover_count=len(covers),
        y_count=counts["y"],
        p_count=counts["p"],
        q_count=counts["q"],
        every_cover_classified_once=classified,
        p_and_q_blocker_part_sizes=(5, 5),
    )
    if not audit.verified:
        raise AssertionError("record-17 facet-type audit changed")
    return audit


def _symbolic_terms() -> tuple[tuple[int, tuple[str, int]], ...]:
    terms = []
    special = mask(SPECIAL_VERTICES)
    for triple in combinations(WORK_VERTICES, 3):
        terms.append((special | mask(triple), ("y", mask(triple))))
    p_base = mask(P_BASE)
    for pair in combinations(P_GROUND, 2):
        terms.append((p_base | mask(pair), ("p", mask(pair))))
    q_base = mask(Q_BASE)
    for pair in combinations(Q_GROUND, 2):
        terms.append((q_base | mask(pair), ("q", mask(pair))))
    return tuple(terms)


@dataclass(frozen=True)
class RidgeSystemAudit:
    total_six_sets: int
    nonempty_groups: int
    both_special_groups: int
    p_triangle_groups: int
    q_triangle_groups: int
    y_incidences: int
    p_incidences: int
    q_incidences: int
    signature_histogram: tuple[tuple[tuple[int, ...], int], ...]

    @property
    def verified(self) -> bool:
        return (
            self.total_six_sets == 3_003
            and self.nonempty_groups == 735
            and self.both_special_groups == 495
            and self.p_triangle_groups == 120
            and self.q_triangle_groups == 120
            and (self.y_incidences, self.p_incidences, self.q_incidences)
            == (1_980, 405, 405)
        )

    def to_record(self) -> dict[str, object]:
        return {
            "total_six_sets": self.total_six_sets,
            "nonempty_groups": self.nonempty_groups,
            "both_special_groups": self.both_special_groups,
            "p_triangle_groups": self.p_triangle_groups,
            "q_triangle_groups": self.q_triangle_groups,
            "incidences": {
                "y": self.y_incidences,
                "p": self.p_incidences,
                "q": self.q_incidences,
            },
            "signature_histogram": [
                {"coefficients": list(signature), "count": count}
                for signature, count in self.signature_histogram
            ],
            "verified": self.verified,
        }


def audit_ridge_system() -> RidgeSystemAudit:
    terms = _symbolic_terms()
    special = mask(SPECIAL_VERTICES)
    p_base = mask(P_BASE)
    q_base = mask(Q_BASE)
    kinds = Counter()
    incidences = Counter()
    signatures = Counter()
    nonempty = 0
    for chosen in combinations(range(TARGET_VERTEX_COUNT), 6):
        six = mask(chosen)
        equation = Counter(key for facet, key in terms if facet & six == facet)
        if not equation:
            continue
        nonempty += 1
        signatures[tuple(sorted(equation.values()))] += 1
        for (kind, _), coefficient in equation.items():
            incidences[kind] += coefficient
        if six & special == special:
            kinds["both"] += 1
        elif six & (1 << 13) and six & p_base == p_base:
            kinds["p"] += 1
        elif six & (1 << 12) and six & q_base == q_base:
            kinds["q"] += 1
        else:
            raise AssertionError("unclassified nonempty record-17 ridge group")
    audit = RidgeSystemAudit(
        total_six_sets=sum(1 for _ in combinations(range(TARGET_VERTEX_COUNT), 6)),
        nonempty_groups=nonempty,
        both_special_groups=kinds["both"],
        p_triangle_groups=kinds["p"],
        q_triangle_groups=kinds["q"],
        y_incidences=incidences["y"],
        p_incidences=incidences["p"],
        q_incidences=incidences["q"],
        signature_histogram=tuple(sorted(signatures.items())),
    )
    if not audit.verified:
        raise AssertionError("record-17 complete ridge system changed")
    return audit


def _adjacency(order: int, edges: Iterable[tuple[int, int]]) -> tuple[int, ...]:
    result = [0] * order
    for first, second in edges:
        result[first] |= 1 << second
        result[second] |= 1 << first
    return tuple(result)


def _triangle_free(adjacency: Sequence[int]) -> bool:
    return all(
        not (
            adjacency[first] & (1 << second)
            and adjacency[first] & (1 << third)
            and adjacency[second] & (1 << third)
        )
        for first, second, third in combinations(range(len(adjacency)), 3)
    )


def _induced_2k2_free(adjacency: Sequence[int]) -> bool:
    for chosen in combinations(range(len(adjacency)), 4):
        chosen_mask = mask(chosen)
        degrees = tuple(
            (adjacency[vertex] & chosen_mask).bit_count() for vertex in chosen
        )
        if sorted(degrees) == [1, 1, 1, 1]:
            return False
    return True


def _vertex_signature(adjacency: Sequence[int], vertex: int) -> tuple[object, ...]:
    degrees = tuple(value.bit_count() for value in adjacency)
    return (
        degrees[vertex],
        tuple(
            sorted(
                degrees[neighbor]
                for neighbor in range(len(adjacency))
                if adjacency[vertex] & (1 << neighbor)
            )
        ),
    )


def _isomorphic(first: Sequence[int], second: Sequence[int]) -> bool:
    if len(first) != len(second):
        return False
    order = len(first)
    first_sig = tuple(_vertex_signature(first, v) for v in range(order))
    second_sig = tuple(_vertex_signature(second, v) for v in range(order))
    if sorted(first_sig) != sorted(second_sig):
        return False
    mapping: dict[int, int] = {}
    used: set[int] = set()

    def search() -> bool:
        if len(mapping) == order:
            return True
        options = []
        for vertex in range(order):
            if vertex in mapping:
                continue
            images = [
                image
                for image in range(order)
                if image not in used
                and first_sig[vertex] == second_sig[image]
                and all(
                    bool(first[vertex] & (1 << other))
                    == bool(second[image] & (1 << mapped))
                    for other, mapped in mapping.items()
                )
            ]
            if not images:
                return False
            options.append((len(images), vertex, images))
        _, vertex, images = min(options)
        for image in images:
            mapping[vertex] = image
            used.add(image)
            if search():
                return True
            used.remove(image)
            del mapping[vertex]
        return False

    return search()


def _pairwise_nonisomorphic(graphs: Sequence[Sequence[int]]) -> bool:
    groups: dict[tuple[object, ...], list[Sequence[int]]] = {}
    for adjacency in graphs:
        signature = (
            sum(value.bit_count() for value in adjacency) // 2,
            tuple(sorted(value.bit_count() for value in adjacency)),
            tuple(sorted(_vertex_signature(adjacency, v) for v in range(len(adjacency)))),
        )
        for previous in groups.setdefault(signature, []):
            if _isomorphic(previous, adjacency):
                return False
        groups[signature].append(adjacency)
    return True


def _canonical_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


@dataclass(frozen=True)
class XArchiveAudit:
    archive_certificate_sha256: str
    raw_stream_sha256: str
    stream_count: int
    triangle_free_count: int
    induced_2k2_free_count: int
    pairwise_nonisomorphic: bool
    admissible_graphs: tuple[tuple[str, tuple[int, ...]], ...]

    @property
    def verified(self) -> bool:
        return (
            self.stream_count == EXPECTED_STREAM_COUNT
            and self.triangle_free_count == EXPECTED_STREAM_COUNT
            and self.induced_2k2_free_count == EXPECTED_ADMISSIBLE_X_COUNT
            and self.pairwise_nonisomorphic
            and len(self.admissible_graphs) == EXPECTED_ADMISSIBLE_X_COUNT
        )

    def to_record(self) -> dict[str, object]:
        return {
            "archive_certificate_sha256": self.archive_certificate_sha256,
            "raw_stream_sha256": self.raw_stream_sha256,
            "stream_count": self.stream_count,
            "triangle_free_count": self.triangle_free_count,
            "induced_2k2_free_count": self.induced_2k2_free_count,
            "pairwise_nonisomorphic": self.pairwise_nonisomorphic,
            "verified": self.verified,
        }


def audit_x_archive(path: Path) -> XArchiveAudit:
    document = json.loads(path.read_text(encoding="utf-8"))
    digest = _canonical_digest(document)
    if document.get("certificate_sha256") != digest:
        raise ValueError("record-17 X archive canonical digest changed")
    words = document.get("graph6")
    if not isinstance(words, list) or not all(isinstance(word, str) for word in words):
        raise ValueError("record-17 X archive is malformed")
    raw_digest = sha256(("".join(word + "\n" for word in words)).encode()).hexdigest()
    if document.get("raw_stream_sha256") != raw_digest or len(set(words)) != len(words):
        raise ValueError("record-17 X raw stream changed")
    all_graphs = []
    admissible = []
    triangle_count = 0
    for word in words:
        order, edges = decode_graph6(word)
        if order != len(LOCAL_ISOLATES):
            raise ValueError("record-17 X stream has mixed orders")
        adjacency = _adjacency(order, edges)
        all_graphs.append(adjacency)
        if not _triangle_free(adjacency):
            continue
        triangle_count += 1
        if _induced_2k2_free(adjacency):
            admissible.append((word, adjacency))
    audit = XArchiveAudit(
        archive_certificate_sha256=digest,
        raw_stream_sha256=raw_digest,
        stream_count=len(words),
        triangle_free_count=triangle_count,
        induced_2k2_free_count=len(admissible),
        pairwise_nonisomorphic=_pairwise_nonisomorphic(all_graphs),
        admissible_graphs=tuple(admissible),
    )
    if not audit.verified:
        raise AssertionError("record-17 X archive audit changed")
    return audit


def _cut_edge(first: int, second: int, cut_mask: int) -> bool:
    return bool((cut_mask >> first) & 1) != bool((cut_mask >> second) & 1)


def _switched_triangle_free(adjacency: Sequence[int], cut_mask: int) -> bool:
    return all(
        not all(
            bool(adjacency[first] & (1 << second))
            != _cut_edge(first, second, cut_mask)
            for first, second in combinations(triple, 2)
        )
        for triple in combinations(LOCAL_ISOLATES, 3)
    )


def _difference_constraint(
    adjacency: Sequence[int], p_mask: int, q_mask: int
) -> bool:
    difference = p_mask ^ q_mask
    return all(
        bool(adjacency[first] & (1 << second))
        == _cut_edge(first, second, p_mask)
        for first, second in combinations(vertices(difference, 9), 2)
    )


def _selected_complements(
    adjacency: Sequence[int], p_mask: int, q_mask: int
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...], tuple[int, ...]]:
    y_triples = []
    for first, second in combinations(LOCAL_ISOLATES, 2):
        global_pair = (first + 3, second + 3)
        x = bool(adjacency[first] & (1 << second))
        if x:
            y_triples.append(mask((0,) + global_pair))
        if x != _cut_edge(first, second, p_mask):
            y_triples.append(mask((1,) + global_pair))
        if x != _cut_edge(first, second, q_mask):
            y_triples.append(mask((2,) + global_pair))
    for local in LOCAL_ISOLATES:
        if bool(p_mask & (1 << local)) != bool(q_mask & (1 << local)):
            y_triples.append(mask((1, 2, local + 3)))
    for triple in combinations(LOCAL_ISOLATES, 3):
        parity = sum(
            bool(adjacency[first] & (1 << second))
            for first, second in combinations(triple, 2)
        ) % 2
        if parity:
            y_triples.append(mask(vertex + 3 for vertex in triple))

    p_pairs = []
    p_colours = {2: False} | {
        local + 3: bool(p_mask & (1 << local)) for local in LOCAL_ISOLATES
    }
    for first, second in combinations(P_GROUND, 2):
        if p_colours[first] != p_colours[second]:
            p_pairs.append(mask((first, second)))
    q_pairs = []
    q_colours = {1: False} | {
        local + 3: bool(q_mask & (1 << local)) for local in LOCAL_ISOLATES
    }
    for first, second in combinations(Q_GROUND, 2):
        if q_colours[first] != q_colours[second]:
            q_pairs.append(mask((first, second)))

    y = tuple(sorted(set(y_triples)))
    p = tuple(sorted(p_pairs))
    q = tuple(sorted(q_pairs))
    complements = tuple(
        sorted(
            [mask(SPECIAL_VERTICES) | triple for triple in y]
            + [mask(P_BASE) | pair for pair in p]
            + [mask(Q_BASE) | pair for pair in q]
        )
    )
    return y, p, q, complements


def _all_ridges_two(complements: Sequence[int]) -> bool:
    groups = Counter()
    for complement in complements:
        outside = FULL_MASK ^ complement
        while outside:
            bit = outside & -outside
            groups[complement | bit] += 1
            outside ^= bit
    return bool(groups) and set(groups.values()) == {2}


def _is_face(face: int, complements: Sequence[int]) -> bool:
    return any(not face & complement for complement in complements)


def _minimal_a4s(complements: Sequence[int]) -> tuple[int, ...]:
    return tuple(
        value
        for chosen in combinations(range(TARGET_VERTEX_COUNT), 4)
        if _is_independent(value := mask(chosen))
        and not _is_face(value, complements)
    )


def _minimal_a5s(
    complements: Sequence[int], a4_masks: Sequence[int]
) -> tuple[int, ...]:
    return tuple(
        value
        for chosen in combinations(range(TARGET_VERTEX_COUNT), 5)
        if _is_independent(value := mask(chosen))
        and not any(value & four == four for four in a4_masks)
        and not _is_face(value, complements)
    )


def face_table_from_complements(complements: Sequence[int]) -> tuple[bool, ...]:
    table = [False] * (1 << TARGET_VERTEX_COUNT)
    table[0] = True
    for complement in complements:
        facet = FULL_MASK ^ complement
        subset = facet
        while True:
            table[subset] = True
            if subset == 0:
                break
            subset = (subset - 1) & facet
    return tuple(table)


def minimal_nonfaces_from_table(
    face_table: Sequence[bool],
) -> dict[int, tuple[int, ...]]:
    result: dict[int, list[int]] = {}
    for value in range(1, 1 << TARGET_VERTEX_COUNT):
        if face_table[value]:
            continue
        bits = value
        minimal = True
        while bits:
            bit = bits & -bits
            if not face_table[value ^ bit]:
                minimal = False
                break
            bits ^= bit
        if minimal:
            result.setdefault(value.bit_count(), []).append(value)
    return {size: tuple(values) for size, values in sorted(result.items())}


def _face_vector(face_table: Sequence[bool]) -> tuple[int, ...]:
    counts = Counter(
        value.bit_count() for value, is_face in enumerate(face_table) if is_face
    )
    return tuple(counts[size] for size in range(1, TARGET_DIMENSION + 1))


@dataclass(frozen=True)
class CandidateWitness:
    missing_four_count: int
    x_graph6: str
    p_local_mask: int
    q_local_mask: int
    selected_facet_complement_masks: tuple[int, ...]
    minimal_nonfaces_by_size: dict[int, tuple[int, ...]]
    face_vector: tuple[int, ...]

    def to_record(self) -> dict[str, object]:
        return {
            "missing_four_count": self.missing_four_count,
            "x_graph6": self.x_graph6,
            "p_local_vertices": list(vertices(self.p_local_mask, 9)),
            "q_local_vertices": list(vertices(self.q_local_mask, 9)),
            "p_global_isolates": [v + 3 for v in vertices(self.p_local_mask, 9)],
            "q_global_isolates": [v + 3 for v in vertices(self.q_local_mask, 9)],
            "selected_facet_complement_masks": list(
                self.selected_facet_complement_masks
            ),
            "selected_facet_complement_sha256": sha256(
                (",".join(map(str, self.selected_facet_complement_masks))).encode()
            ).hexdigest(),
            "facet_count": len(self.selected_facet_complement_masks),
            "minimal_nonfaces_by_size": {
                str(size): list(values)
                for size, values in self.minimal_nonfaces_by_size.items()
            },
            "minimal_nonface_vertices_by_size": {
                str(size): [list(vertices(value)) for value in values]
                for size, values in self.minimal_nonfaces_by_size.items()
            },
            "face_vector": list(self.face_vector),
        }


def _coloured_isomorphic(
    adjacency: Sequence[int],
    first_p: int,
    first_q: int,
    second_p: int,
    second_q: int,
) -> bool:
    """Test cut-pair equivalence under ``Aut(X)`` and branch interchange."""

    order = len(adjacency)
    degrees = tuple(value.bit_count() for value in adjacency)

    def colour(value_p: int, value_q: int) -> tuple[int, ...]:
        return tuple(
            int(bool(value_p & (1 << vertex)))
            + 2 * int(bool(value_q & (1 << vertex)))
            for vertex in range(order)
        )

    first_colours = colour(first_p, first_q)

    def fixed_branch(second_colours: Sequence[int]) -> bool:
        def signature(colours: Sequence[int], vertex: int) -> tuple[object, ...]:
            return (
                colours[vertex],
                degrees[vertex],
                tuple(
                    sorted(
                        (colours[neighbor], degrees[neighbor])
                        for neighbor in range(order)
                        if adjacency[vertex] & (1 << neighbor)
                    )
                ),
            )

        first_signatures = tuple(
            signature(first_colours, vertex) for vertex in range(order)
        )
        second_signatures = tuple(
            signature(second_colours, vertex) for vertex in range(order)
        )
        if sorted(first_signatures) != sorted(second_signatures):
            return False
        mapping: dict[int, int] = {}
        used: set[int] = set()

        def search() -> bool:
            if len(mapping) == order:
                return True
            options = []
            for vertex in range(order):
                if vertex in mapping:
                    continue
                images = [
                    image
                    for image in range(order)
                    if image not in used
                    and first_signatures[vertex] == second_signatures[image]
                    and all(
                        bool(adjacency[vertex] & (1 << other))
                        == bool(adjacency[image] & (1 << mapped))
                        for other, mapped in mapping.items()
                    )
                ]
                if not images:
                    return False
                options.append((len(images), vertex, images))
            _, vertex, images = min(options)
            for image in images:
                mapping[vertex] = image
                used.add(image)
                if search():
                    return True
                used.remove(image)
                del mapping[vertex]
            return False

        return search()

    return fixed_branch(colour(second_p, second_q)) or fixed_branch(
        colour(second_q, second_p)
    )


@dataclass(frozen=True)
class ResidualOrbit:
    representative: CandidateWitness
    representation_count: int

    def to_record(self) -> dict[str, object]:
        return {
            "missing_four_count": self.representative.missing_four_count,
            "x_graph6": self.representative.x_graph6,
            "p_local_vertices": list(
                vertices(self.representative.p_local_mask, 9)
            ),
            "q_local_vertices": list(
                vertices(self.representative.q_local_mask, 9)
            ),
            "representation_count": self.representation_count,
            "facet_count": len(
                self.representative.selected_facet_complement_masks
            ),
            "minimal_nonface_histogram": {
                str(size): len(values)
                for size, values in self.representative.minimal_nonfaces_by_size.items()
            },
            "selected_facet_complement_sha256": sha256(
                (
                    ",".join(
                        map(
                            str,
                            self.representative.selected_facet_complement_masks,
                        )
                    )
                ).encode()
            ).hexdigest(),
        }


def _residual_orbits(
    witnesses: Sequence[CandidateWitness], archive: XArchiveAudit
) -> tuple[ResidualOrbit, ...]:
    adjacency_by_word = dict(archive.admissible_graphs)
    groups: list[list[CandidateWitness]] = []
    for witness in witnesses:
        adjacency = adjacency_by_word[witness.x_graph6]
        matching = None
        for group in groups:
            representative = group[0]
            if representative.x_graph6 != witness.x_graph6:
                continue
            if _coloured_isomorphic(
                adjacency,
                witness.p_local_mask,
                witness.q_local_mask,
                representative.p_local_mask,
                representative.q_local_mask,
            ):
                matching = group
                break
        if matching is None:
            groups.append([witness])
        else:
            matching.append(witness)
    return tuple(
        ResidualOrbit(group[0], len(group)) for group in groups
    )


@dataclass(frozen=True)
class EnumerationAudit:
    coupled_state_count: int
    every_coupled_state_passes_all_ridge_equations: bool
    stage_counts: dict[str, dict[int, int]]
    witnesses: tuple[CandidateWitness, ...]
    residual_orbits: tuple[ResidualOrbit, ...]

    @property
    def verified(self) -> bool:
        branch_counts = Counter(
            witness.missing_four_count for witness in self.witnesses
        )
        a0 = tuple(
            witness for witness in self.witnesses if witness.missing_four_count == 0
        )
        orbit_signatures = tuple(
            (
                orbit.representative.missing_four_count,
                orbit.representative.x_graph6,
                orbit.representative.p_local_mask,
                orbit.representation_count,
            )
            for orbit in self.residual_orbits
        )
        return (
            self.coupled_state_count == EXPECTED_COUPLED_STATE_COUNT
            and self.every_coupled_state_passes_all_ridge_equations
            and self.stage_counts == EXPECTED_STAGE_COUNTS
            and branch_counts == {0: 1, 1: 17, 2: 3}
            and all(
                witness.p_local_mask == witness.q_local_mask
                for witness in self.witnesses
            )
            and len(a0) == 1
            and a0[0].x_graph6 == UNIQUE_X_GRAPH6
            and a0[0].p_local_mask == UNIQUE_P_LOCAL_MASK
            and a0[0].q_local_mask == UNIQUE_Q_LOCAL_MASK
            and a0[0].minimal_nonfaces_by_size.get(5)
            == UNIQUE_A5_MASKS
            and orbit_signatures == EXPECTED_RESIDUAL_ORBITS
        )

    def to_record(self) -> dict[str, object]:
        return {
            "coupled_state_count": self.coupled_state_count,
            "every_coupled_state_passes_all_ridge_equations": (
                self.every_coupled_state_passes_all_ridge_equations
            ),
            "stage_counts": {
                stage: {str(a): count for a, count in counts.items()}
                for stage, counts in self.stage_counts.items()
            },
            "witnesses": [witness.to_record() for witness in self.witnesses],
            "residual_orbit_count": len(self.residual_orbits),
            "residual_orbits": [
                orbit.to_record() for orbit in self.residual_orbits
            ],
            "verified": self.verified,
        }


def enumerate_whole_shadow(
    schedules: Sequence[Mapping[str, object]], archive: XArchiveAudit
) -> EnumerationAudit:
    parsed = {
        int(schedule["missing_four_count"]): {
            "missing_five_count": int(schedule["missing_five_count_offset"]),
            "face_vector": tuple(int(value) for value in schedule["f_vector"]),
        }
        for schedule in schedules
    }
    if set(parsed) != set(range(7)):
        raise ValueError("record-17 schedules changed")
    cuts = tuple(mask(chosen) for chosen in combinations(LOCAL_ISOLATES, 5))
    independent_triples = tuple(
        mask(chosen)
        for chosen in combinations(range(TARGET_VERTEX_COUNT), 3)
        if _is_independent(mask(chosen))
    )
    stage = {name: Counter() for name in EXPECTED_STAGE_COUNTS}
    coupled_count = 0
    all_ridges = True
    witnesses = []
    work_mask = mask(WORK_VERTICES)

    for word, adjacency in archive.admissible_graphs:
        valid_cuts = tuple(
            cut for cut in cuts if _switched_triangle_free(adjacency, cut)
        )
        for p_mask in valid_cuts:
            for q_mask in valid_cuts:
                if not _difference_constraint(adjacency, p_mask, q_mask):
                    continue
                coupled_count += 1
                y, p_edges, q_edges, complements = _selected_complements(
                    adjacency, p_mask, q_mask
                )
                if len(p_edges) != 25 or len(q_edges) != 25:
                    raise AssertionError("record-17 K5,5 facet count changed")
                ridge_ok = _all_ridges_two(complements)
                all_ridges &= ridge_ok
                if not ridge_ok:
                    raise AssertionError("coupled characterization missed a ridge")
                y_count = len(y)
                if y_count < 52 or y_count > 60 or y_count % 2:
                    continue
                a = (60 - y_count) // 2
                stage["target_y_count"][a] += 1
                if not all(_is_face(triple, complements) for triple in independent_triples):
                    continue
                stage["independent_triples"][a] += 1
                a4_masks = _minimal_a4s(complements)
                if len(a4_masks) != a:
                    continue
                stage["exact_a4"][a] += 1
                a5_masks = _minimal_a5s(complements, a4_masks)
                missing_five_count = sum(
                    _is_independent(value := mask(chosen))
                    and not _is_face(value, complements)
                    for chosen in combinations(range(TARGET_VERTEX_COUNT), 5)
                )
                if (
                    missing_five_count != parsed[a]["missing_five_count"]
                    or any(not _allowed_a4(value) for value in a4_masks)
                    or any(value & ~work_mask for value in a5_masks)
                ):
                    continue
                stage["exact_a5_and_support"][a] += 1
                face_table = face_table_from_complements(complements)
                minimal = minimal_nonfaces_from_table(face_table)
                if (
                    set(minimal) - {2, 4, 5}
                    or set(minimal.get(2, ())) != set(EDGE_MASKS)
                    or set(minimal.get(4, ())) != set(a4_masks)
                    or set(minimal.get(5, ())) != set(a5_masks)
                ):
                    continue
                stage["full_blocker"][a] += 1
                vector = _face_vector(face_table)
                if vector != parsed[a]["face_vector"]:
                    continue
                stage["full_face_vector"][a] += 1
                witnesses.append(
                    CandidateWitness(
                        missing_four_count=a,
                        x_graph6=word,
                        p_local_mask=p_mask,
                        q_local_mask=q_mask,
                        selected_facet_complement_masks=complements,
                        minimal_nonfaces_by_size=minimal,
                        face_vector=vector,
                    )
                )
    normalized_stage = {
        name: {a: counts[a] for a in range(5)} for name, counts in stage.items()
    }
    audit = EnumerationAudit(
        coupled_state_count=coupled_count,
        every_coupled_state_passes_all_ridge_equations=all_ridges,
        stage_counts=normalized_stage,
        witnesses=tuple(witnesses),
        residual_orbits=_residual_orbits(witnesses, archive),
    )
    if not audit.verified:
        raise AssertionError(
            "record-17 whole-shadow census changed: "
            f"coupled={coupled_count}, stages={normalized_stage}, "
            f"witnesses={len(witnesses)}"
        )
    return audit


@dataclass(frozen=True)
class Record17Audit:
    shadow: ShadowAudit
    support: SupportAudit
    facet_types: FacetTypeAudit
    ridges: RidgeSystemAudit
    x_archive: XArchiveAudit
    enumeration: EnumerationAudit

    @property
    def verified(self) -> bool:
        return (
            self.shadow.verified
            and self.support.verified
            and self.facet_types.verified
            and self.ridges.verified
            and self.x_archive.verified
            and self.enumeration.verified
        )

    def to_record(self) -> dict[str, object]:
        return {
            "shadow": self.shadow.to_record(),
            "support": self.support.to_record(),
            "facet_types": self.facet_types.to_record(),
            "ridges": self.ridges.to_record(),
            "x_archive": self.x_archive.to_record(),
            "enumeration": self.enumeration.to_record(),
            "conclusion": (
                "a=1,2,3,4 are excluded and a=0 has one incidence candidate; "
                "a=5,6 are excluded separately by strict vertex-link g4"
            ),
            "verified": self.verified,
        }


def audit_record17(
    schedules: Sequence[Mapping[str, object]], x_archive_path: Path
) -> Record17Audit:
    x_archive = audit_x_archive(x_archive_path)
    audit = Record17Audit(
        shadow=audit_shadow(),
        support=audit_support(),
        facet_types=audit_facet_types(),
        ridges=audit_ridge_system(),
        x_archive=x_archive,
        enumeration=enumerate_whole_shadow(schedules, x_archive),
    )
    if not audit.verified:
        raise AssertionError("record-17 audit changed")
    return audit

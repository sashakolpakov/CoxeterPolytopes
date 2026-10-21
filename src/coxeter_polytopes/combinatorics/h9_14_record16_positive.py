"""Positive-``a`` completion obstruction for H9/14 record 16.

This module continues the analytic ridge reduction in
``h9_14_record16``.  The same ``y/z`` facet types remain valid when missing
four-faces are present because every admissible four-face lies on
``W={0,...,11}`` and meets ``{0,1}``.  A coarse degree bound eliminates the
``a=1,2`` schedules.  Only five small graph invariants remain for ``a=3,4``;
an archived nauty ``geng -t`` stream on at most eight vertices is replayed
independently to finish them.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record16 import (
    ACTIVE_VERTICES,
    EDGE_MASKS,
    ISOLATED_VERTICES,
    SPECIAL_VERTICES,
    TARGET_GRAPH6,
    TARGET_RECORD_INDEX,
    TARGET_VERTEX_COUNT,
    WORK_VERTICES,
    mask,
    vertices,
)
from coxeter_polytopes.combinatorics.local_link import (
    analyze_local_link_graph,
)
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


EXPECTED_STREAM_COUNTS = {6: 38, 7: 107, 8: 410}
EXPECTED_ADMISSIBLE_COUNTS = {6: 22, 7: 41, 8: 82}
TARGET_X_INVARIANTS = ((8, 72), (7, 62), (8, 70), (7, 60), (6, 54))


def _is_independent(candidate: int) -> bool:
    return not any(candidate & edge == edge for edge in EDGE_MASKS)


@dataclass(frozen=True)
class PositiveSupportAudit:
    independent_four_count: int
    four_on_work_vertices: int
    four_meeting_active: int
    isolate_only_four_count: int
    special_vertex_four_incidence_upper: tuple[int, int]
    isolate_only_four_has_special_extensions: bool
    independent_five_count: int
    five_on_work_vertices: int
    five_meeting_special: int

    @property
    def verified(self) -> bool:
        return (
            self.independent_four_count == 735
            and self.four_on_work_vertices == 495
            and self.four_meeting_active == 285
            and self.isolate_only_four_count == 210
            and self.special_vertex_four_incidence_upper == (0, 0)
            and self.isolate_only_four_has_special_extensions
            and self.independent_five_count == 1_212
            and self.five_on_work_vertices == 792
            and self.five_meeting_special == 420
        )

    def to_record(self) -> dict[str, object]:
        return {
            "independent_four_count": self.independent_four_count,
            "four_on_work_vertices": self.four_on_work_vertices,
            "four_meeting_active": self.four_meeting_active,
            "isolate_only_four_count": self.isolate_only_four_count,
            "special_vertex_four_incidence_upper": list(
                self.special_vertex_four_incidence_upper
            ),
            "isolate_only_four_has_special_extensions": (
                self.isolate_only_four_has_special_extensions
            ),
            "independent_five_count": self.independent_five_count,
            "five_on_work_vertices": self.five_on_work_vertices,
            "five_meeting_special": self.five_meeting_special,
            "conclusion": (
                "every admissible A4 lies on W and meets {0,1}; every "
                "additional minimal A5 lies on W"
            ),
            "verified": self.verified,
        }


def audit_positive_support() -> PositiveSupportAudit:
    """Replay the strict-link support reduction for ``A4`` and extra ``A5``."""

    independent_fours = tuple(
        mask(chosen)
        for chosen in combinations(range(TARGET_VERTEX_COUNT), 4)
        if _is_independent(mask(chosen))
    )
    independent_fives = tuple(
        mask(chosen)
        for chosen in combinations(range(TARGET_VERTEX_COUNT), 5)
        if _is_independent(mask(chosen))
    )
    work_mask = mask(WORK_VERTICES)
    active_mask = mask(ACTIVE_VERTICES)
    special_mask = mask(SPECIAL_VERTICES)
    on_work_fours = tuple(
        face for face in independent_fours if not face & ~work_mask
    )
    isolate_only_fours = tuple(
        face for face in on_work_fours if not face & active_mask
    )
    special_extensions = all(
        all(_is_independent(face | (1 << special)) for special in SPECIAL_VERTICES)
        for face in isolate_only_fours
    )
    local = analyze_local_link_graph(TARGET_GRAPH6)
    special_upper = tuple(
        local.vertex_figures[vertex].missing_four_incidence_interval[1]
        for vertex in SPECIAL_VERTICES
    )
    audit = PositiveSupportAudit(
        independent_four_count=len(independent_fours),
        four_on_work_vertices=len(on_work_fours),
        four_meeting_active=sum(bool(face & active_mask) for face in on_work_fours),
        isolate_only_four_count=len(isolate_only_fours),
        special_vertex_four_incidence_upper=special_upper,
        isolate_only_four_has_special_extensions=special_extensions,
        independent_five_count=len(independent_fives),
        five_on_work_vertices=sum(
            not face & ~work_mask for face in independent_fives
        ),
        five_meeting_special=sum(bool(face & special_mask) for face in independent_fives),
    )
    if not audit.verified:
        raise AssertionError("record-16 positive support audit changed")
    return audit


@dataclass(frozen=True)
class ScheduleCutCase:
    missing_four_count: int
    target_facet_count: int
    smaller_cut_side: int
    larger_cut_side: int
    z_edge_count: int
    required_y_count: int
    cauchy_y_upper_bound: int

    @property
    def survives_coarse_bound(self) -> bool:
        return self.required_y_count <= self.cauchy_y_upper_bound

    def to_record(self) -> dict[str, int | bool]:
        return {
            "missing_four_count": self.missing_four_count,
            "target_facet_count": self.target_facet_count,
            "smaller_cut_side": self.smaller_cut_side,
            "larger_cut_side": self.larger_cut_side,
            "z_edge_count": self.z_edge_count,
            "required_y_count": self.required_y_count,
            "cauchy_y_upper_bound": self.cauchy_y_upper_bound,
            "survives_coarse_bound": self.survives_coarse_bound,
        }


def positive_cut_cases(
    schedules: Sequence[Mapping[str, object]],
) -> tuple[ScheduleCutCase, ...]:
    parsed = {
        int(schedule["missing_four_count"]): int(schedule["f_vector"][-1])
        for schedule in schedules
    }
    result = []
    for a in range(1, 5):
        target = parsed[a]
        for smaller in range(6):
            larger = 10 - smaller
            z_edges = smaller * larger
            result.append(
                ScheduleCutCase(
                    missing_four_count=a,
                    target_facet_count=target,
                    smaller_cut_side=smaller,
                    larger_cut_side=larger,
                    z_edge_count=z_edges,
                    required_y_count=target - 2 * z_edges,
                    cauchy_y_upper_bound=9 * larger,
                )
            )
    expected_survivors = {
        (3, 2, 8, 72),
        (3, 3, 7, 62),
        (4, 2, 8, 70),
        (4, 3, 7, 60),
        (4, 4, 6, 54),
    }
    observed = {
        (
            case.missing_four_count,
            case.smaller_cut_side,
            case.larger_cut_side,
            case.required_y_count,
        )
        for case in result
        if case.survives_coarse_bound
    }
    if observed != expected_survivors:
        raise AssertionError("record-16 residual cut cases changed")
    return tuple(result)


def _adjacency(order: int, edges: Iterable[tuple[int, int]]) -> tuple[int, ...]:
    result = [0] * order
    for first, second in edges:
        result[first] |= 1 << second
        result[second] |= 1 << first
    return tuple(result)


def _triangle_free(adjacency: Sequence[int]) -> bool:
    order = len(adjacency)
    return all(
        not (
            adjacency[first] & (1 << second)
            and adjacency[first] & (1 << third)
            and adjacency[second] & (1 << third)
        )
        for first, second, third in combinations(range(order), 3)
    )


def _induced_2k2_free(adjacency: Sequence[int]) -> bool:
    for chosen in combinations(range(len(adjacency)), 4):
        degrees = []
        edge_count = 0
        chosen_mask = mask(chosen)
        for vertex in chosen:
            degree = (adjacency[vertex] & chosen_mask).bit_count()
            degrees.append(degree)
            edge_count += degree
        if edge_count // 2 == 2 and sorted(degrees) == [1, 1, 1, 1]:
            return False
    return True


def _vertex_signature(adjacency: Sequence[int], vertex: int) -> tuple[object, ...]:
    degrees = tuple(value.bit_count() for value in adjacency)
    neighbor_degrees = tuple(
        sorted(
            degrees[neighbor]
            for neighbor in range(len(adjacency))
            if adjacency[vertex] & (1 << neighbor)
        )
    )
    return degrees[vertex], neighbor_degrees


def _isomorphic(first: Sequence[int], second: Sequence[int]) -> bool:
    """Exact small-graph isomorphism test used to audit the geng quotient."""

    if len(first) != len(second):
        return False
    order = len(first)
    first_signatures = tuple(
        _vertex_signature(first, vertex) for vertex in range(order)
    )
    second_signatures = tuple(
        _vertex_signature(second, vertex) for vertex in range(order)
    )
    if sorted(first_signatures) != sorted(second_signatures):
        return False
    mapping: dict[int, int] = {}
    used: set[int] = set()

    def search() -> bool:
        if len(mapping) == order:
            return True
        candidates_by_vertex = []
        for vertex in range(order):
            if vertex in mapping:
                continue
            candidates = []
            for image in range(order):
                if image in used or first_signatures[vertex] != second_signatures[image]:
                    continue
                if all(
                    bool(first[vertex] & (1 << other))
                    == bool(second[image] & (1 << mapped))
                    for other, mapped in mapping.items()
                ):
                    candidates.append(image)
            if not candidates:
                return False
            candidates_by_vertex.append((len(candidates), vertex, candidates))
        _, vertex, candidates = min(candidates_by_vertex)
        for image in candidates:
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


def _x_y_triples(
    order: int, edges: Sequence[tuple[int, int]]
) -> tuple[int, ...]:
    x_edges = {
        tuple(sorted((first + 2, second + 2))) for first, second in edges
    }
    result = []
    active = set(ACTIVE_VERTICES)
    for triple in combinations(WORK_VERTICES, 3):
        active_count = len(set(triple) & active)
        isolates = tuple(vertex for vertex in triple if vertex not in active)
        if active_count == 2:
            selected = False
        elif active_count == 1:
            selected = tuple(sorted(isolates)) in x_edges
        else:
            selected = (
                sum(tuple(sorted(pair)) in x_edges for pair in combinations(isolates, 2))
                % 2
                == 1
            )
        if selected:
            result.append(mask(triple))
    return tuple(result)


def _a4_blockers(y_triples: Sequence[int]) -> tuple[int, ...]:
    active_mask = mask(ACTIVE_VERTICES)
    return tuple(
        candidate
        for chosen in combinations(WORK_VERTICES, 4)
        if (candidate := mask(chosen)) & active_mask
        and all(candidate & triple for triple in y_triples)
    )


def _complete_bipartite_parts(
    adjacency: Sequence[int],
) -> tuple[tuple[int, ...], tuple[int, ...]] | None:
    """Return the two parts exactly when the graph is complete bipartite."""

    order = len(adjacency)
    if order == 0:
        return None
    colour: dict[int, int] = {0: 0}
    queue = [0]
    while queue:
        vertex = queue.pop()
        for neighbor in range(order):
            if not adjacency[vertex] & (1 << neighbor):
                continue
            expected = 1 - colour[vertex]
            if neighbor in colour:
                if colour[neighbor] != expected:
                    return None
            else:
                colour[neighbor] = expected
                queue.append(neighbor)
    if len(colour) != order:
        return None
    parts = tuple(
        tuple(vertex for vertex in range(order) if colour[vertex] == value)
        for value in (0, 1)
    )
    if not all(
        bool(adjacency[first] & (1 << second)) == (colour[first] != colour[second])
        for first, second in combinations(range(order), 2)
    ):
        return None
    return parts


@dataclass(frozen=True)
class XGraphAudit:
    graph6: str
    order: int
    edge_count: int
    degree_sequence: tuple[int, ...]
    y_count: int
    a4_blocker_masks: tuple[int, ...]
    complete_bipartite_3_3: bool

    def to_record(self) -> dict[str, object]:
        return {
            "graph6": self.graph6,
            "order": self.order,
            "edge_count": self.edge_count,
            "degree_sequence": list(self.degree_sequence),
            "y_count": self.y_count,
            "a4_blocker_count": len(self.a4_blocker_masks),
            "a4_blocker_masks": list(self.a4_blocker_masks),
            "a4_blocker_vertices": [
                list(vertices(value)) for value in self.a4_blocker_masks
            ],
            "complete_bipartite_3_3": self.complete_bipartite_3_3,
        }


def _audit_x_graph(word: str) -> tuple[XGraphAudit, tuple[int, ...]]:
    order, edges = decode_graph6(word)
    adjacency = _adjacency(order, edges)
    if not _triangle_free(adjacency):
        raise ValueError("archived geng -t graph contains a triangle")
    degrees = tuple(sorted((value.bit_count() for value in adjacency), reverse=True))
    edge_count = len(edges)
    y_count = 12 * edge_count - sum(degree * degree for degree in degrees)
    y_triples = _x_y_triples(order, edges)
    if len(y_triples) != y_count:
        raise AssertionError("X degree formula disagrees with direct y replay")
    blockers = _a4_blockers(y_triples)
    bipartition = _complete_bipartite_parts(adjacency)
    complete_3_3 = (
        bipartition is not None
        and tuple(sorted(map(len, bipartition))) == (3, 3)
    )
    return (
        XGraphAudit(
            graph6=word,
            order=order,
            edge_count=edge_count,
            degree_sequence=degrees,
            y_count=y_count,
            a4_blocker_masks=blockers,
            complete_bipartite_3_3=complete_3_3,
        ),
        adjacency,
    )


def _canonical_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _stream_digest(words: Sequence[str]) -> str:
    return sha256(("".join(word + "\n" for word in words)).encode()).hexdigest()


@dataclass(frozen=True)
class XArchiveAudit:
    archive_certificate_sha256: str
    stream_counts: Mapping[int, int]
    admissible_counts: Mapping[int, int]
    pairwise_nonisomorphic: Mapping[int, bool]
    target_graphs: Mapping[tuple[int, int], tuple[XGraphAudit, ...]]

    @property
    def verified(self) -> bool:
        expected_targets = {
            (8, 72): 0,
            (7, 62): 1,
            (8, 70): 1,
            (7, 60): 2,
            (6, 54): 1,
        }
        return (
            dict(self.stream_counts) == EXPECTED_STREAM_COUNTS
            and dict(self.admissible_counts) == EXPECTED_ADMISSIBLE_COUNTS
            and all(self.pairwise_nonisomorphic.values())
            and {
                key: len(self.target_graphs.get(key, ()))
                for key in expected_targets
            }
            == expected_targets
        )

    def to_record(self) -> dict[str, object]:
        return {
            "archive_certificate_sha256": self.archive_certificate_sha256,
            "stream_counts": {
                str(key): value for key, value in sorted(self.stream_counts.items())
            },
            "admissible_triangle_and_2k2_free_counts": {
                str(key): value
                for key, value in sorted(self.admissible_counts.items())
            },
            "pairwise_nonisomorphic": {
                str(key): value
                for key, value in sorted(self.pairwise_nonisomorphic.items())
            },
            "target_graphs": {
                f"n{order}_y{y_count}": [
                    graph.to_record() for graph in self.target_graphs[(order, y_count)]
                ]
                for order, y_count in TARGET_X_INVARIANTS
            },
            "verified": self.verified,
        }


def audit_x_graph_archive(path: str | Path) -> XArchiveAudit:
    path = Path(path)
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema_version") != 1 or document.get("kind") != (
        "h9_14_record16_triangle_free_x_graph_streams"
    ):
        raise ValueError("wrong record-16 X archive metadata")
    if document.get("certificate_sha256") != _canonical_digest(document):
        raise ValueError("record-16 X archive digest changed")
    streams = document.get("streams")
    if not isinstance(streams, Mapping):
        raise ValueError("record-16 X archive lacks streams")
    counts: dict[int, int] = {}
    admissible_counts: dict[int, int] = {}
    nonisomorphic: dict[int, bool] = {}
    targets: dict[tuple[int, int], list[XGraphAudit]] = {
        key: [] for key in TARGET_X_INVARIANTS
    }
    for order, expected_count in EXPECTED_STREAM_COUNTS.items():
        stream = streams.get(str(order))
        if not isinstance(stream, Mapping):
            raise ValueError(f"record-16 X archive lacks order {order}")
        words = stream.get("graph6")
        if not isinstance(words, list) or not all(isinstance(word, str) for word in words):
            raise ValueError("record-16 X archive has malformed graph6 words")
        if len(words) != expected_count or len(set(words)) != expected_count:
            raise ValueError("record-16 X archive count or distinctness changed")
        if stream.get("raw_stream_sha256") != _stream_digest(words):
            raise ValueError("record-16 X raw stream digest changed")
        all_adjacencies = []
        retained = 0
        for word in words:
            graph, adjacency = _audit_x_graph(word)
            if graph.order != order:
                raise ValueError("record-16 X archive mixes graph orders")
            all_adjacencies.append(adjacency)
            if not _induced_2k2_free(adjacency):
                continue
            retained += 1
            key = (order, graph.y_count)
            if key in targets:
                targets[key].append(graph)
        counts[order] = len(words)
        admissible_counts[order] = retained
        nonisomorphic[order] = _pairwise_nonisomorphic(all_adjacencies)
    audit = XArchiveAudit(
        archive_certificate_sha256=str(document["certificate_sha256"]),
        stream_counts=counts,
        admissible_counts=admissible_counts,
        pairwise_nonisomorphic=nonisomorphic,
        target_graphs={key: tuple(value) for key, value in targets.items()},
    )
    if not audit.verified:
        raise AssertionError("record-16 X archive audit changed")
    return audit


@dataclass(frozen=True)
class PositiveBranchAudit:
    support: PositiveSupportAudit
    cut_cases: tuple[ScheduleCutCase, ...]
    x_archive: XArchiveAudit
    forced_a4_masks: tuple[int, ...]
    forced_five_shadow_size: int
    forced_additional_minimal_a5_count: int
    final_a5_audit: "FinalA5Audit"
    full_blocker_audit: "FullBlockerAudit"
    ridge_z_cut_part_sizes: tuple[int, int]

    @property
    def verified(self) -> bool:
        residual = [case for case in self.cut_cases if case.survives_coarse_bound]
        target = self.x_archive.target_graphs
        return (
            self.support.verified
            and len(residual) == 5
            and len(target[(8, 72)]) == 0
            and len(target[(7, 62)]) == 1
            and len(target[(7, 62)][0].a4_blocker_masks) == 0
            and len(target[(8, 70)]) == 1
            and len(target[(8, 70)][0].a4_blocker_masks) == 0
            and len(target[(7, 60)]) == 2
            and all(len(graph.a4_blocker_masks) == 2 for graph in target[(7, 60)])
            and len(target[(6, 54)]) == 1
            and target[(6, 54)][0].complete_bipartite_3_3
            and self.forced_a4_masks == target[(6, 54)][0].a4_blocker_masks
            and len(self.forced_a4_masks) == 4
            and self.forced_five_shadow_size == 30
            and self.forced_additional_minimal_a5_count == 2
            and self.final_a5_audit.verified
            and self.full_blocker_audit.verified
            and tuple(
                sorted(dict(self.full_blocker_audit.minimal_transversals_by_size)[4])
            )
            == tuple(sorted(self.forced_a4_masks))
            and self.ridge_z_cut_part_sizes == (4, 6)
        )

    def to_record(self) -> dict[str, object]:
        residual = [
            case.to_record() for case in self.cut_cases if case.survives_coarse_bound
        ]
        return {
            "support_reduction": self.support.to_record(),
            "cut_cases": [case.to_record() for case in self.cut_cases],
            "coarse_survivors": residual,
            "x_graph_archive_audit": self.x_archive.to_record(),
            "forced_final_case": {
                "missing_four_count": 4,
                "z_cut_part_sizes": list(self.ridge_z_cut_part_sizes),
                "x_graph": "K3,3 on the six-vertex cut side",
                "a4_masks": list(self.forced_a4_masks),
                "a4_vertices": [
                    list(vertices(value)) for value in self.forced_a4_masks
                ],
                "forced_five_shadow_size": self.forced_five_shadow_size,
                "schedule_missing_five_offset": 32,
                "additional_minimal_a5_count": (
                    self.forced_additional_minimal_a5_count
                ),
                "extra_a5_exhaustive_audit": self.final_a5_audit.to_record(),
                "direct_full_blocker_replay": (
                    self.full_blocker_audit.to_record()
                ),
                "contradiction": "K5,5 is not the ridge-forced K4,6",
            },
            "conclusion": (
                "branches a=1,2,3,4 are impossible; together with the "
                "a=0 and strict-g4 a=5,6 certificates, record 16 is excluded"
            ),
            "verified": self.verified,
        }


@dataclass(frozen=True)
class FinalA5Audit:
    candidate_minimal_a5_count: int
    unordered_two_a5_families: int
    tau_z_at_least_two_families: int
    all_tau_valid_rows_avoid_active: bool
    all_tau_valid_intersections_empty: bool
    all_tau_valid_unions_are_isolates: bool
    z_edge_count_histogram: tuple[tuple[int, int], ...]
    z_part_size_histogram: tuple[tuple[tuple[int, int], int], ...]
    example_rows: tuple[int, int]

    @property
    def verified(self) -> bool:
        return (
            self.candidate_minimal_a5_count == 762
            and self.unordered_two_a5_families == 289_941
            and self.tau_z_at_least_two_families == 126
            and self.all_tau_valid_rows_avoid_active
            and self.all_tau_valid_intersections_empty
            and self.all_tau_valid_unions_are_isolates
            and self.z_edge_count_histogram == ((25, 126),)
            and self.z_part_size_histogram == (((5, 5), 126),)
        )

    def to_record(self) -> dict[str, object]:
        return {
            "candidate_minimal_a5_count": self.candidate_minimal_a5_count,
            "unordered_two_a5_families": self.unordered_two_a5_families,
            "tau_z_at_least_two_families": self.tau_z_at_least_two_families,
            "all_tau_valid_rows_avoid_active": (
                self.all_tau_valid_rows_avoid_active
            ),
            "all_tau_valid_intersections_empty": (
                self.all_tau_valid_intersections_empty
            ),
            "all_tau_valid_unions_are_isolates": (
                self.all_tau_valid_unions_are_isolates
            ),
            "z_edge_count_histogram": {
                str(edge_count): multiplicity
                for edge_count, multiplicity in self.z_edge_count_histogram
            },
            "z_part_size_histogram": [
                {"part_sizes": list(parts), "multiplicity": multiplicity}
                for parts, multiplicity in self.z_part_size_histogram
            ],
            "example_rows": [list(vertices(row)) for row in self.example_rows],
            "verified": self.verified,
        }


def audit_final_extra_a5s(
    forced_a4: Sequence[int], shadow: set[int]
) -> FinalA5Audit:
    """Enumerate the final case's two additional minimal five-faces."""

    active_mask = mask(ACTIVE_VERTICES)
    isolate_mask = mask(ISOLATED_VERTICES)
    candidates = tuple(
        candidate
        for chosen in combinations(WORK_VERTICES, 5)
        if (candidate := mask(chosen)) not in shadow
        and all(candidate & face != face for face in forced_a4)
    )
    valid: list[tuple[int, int]] = []
    edge_histogram: Counter[int] = Counter()
    part_histogram: Counter[tuple[int, int]] = Counter()
    all_avoid_active = True
    all_empty = True
    all_union = True
    for first, second in combinations(candidates, 2):
        z_rows = tuple(
            row for row in (first, second) if not row & active_mask
        )
        # For a nonempty clutter, tau(Z)>=2 is exactly the absence of a
        # common vertex.  With fewer than two rows the intersection is
        # necessarily nonempty, so this also forces both rows into Z.
        if not z_rows:
            continue
        intersection = z_rows[0]
        for row in z_rows[1:]:
            intersection &= row
        if intersection:
            continue
        valid.append((first, second))
        all_avoid_active &= len(z_rows) == 2
        all_empty &= first & second == 0
        all_union &= first | second == isolate_mask
        z_edges = tuple(
            mask(pair)
            for pair in combinations(ISOLATED_VERTICES, 2)
            if all(mask(pair) & row for row in z_rows)
        )
        edge_histogram[len(z_edges)] += 1
        adjacency = {vertex: set() for vertex in ISOLATED_VERTICES}
        for pair_mask in z_edges:
            first_vertex, second_vertex = vertices(pair_mask)
            adjacency[first_vertex].add(second_vertex)
            adjacency[second_vertex].add(first_vertex)
        seed = ISOLATED_VERTICES[0]
        first_part = {seed}
        second_part = set(adjacency[seed])
        first_part.update(set(ISOLATED_VERTICES) - second_part - {seed})
        if not all(
            (second_vertex in adjacency[first_vertex])
            == ((first_vertex in first_part) != (second_vertex in first_part))
            for first_vertex, second_vertex in combinations(ISOLATED_VERTICES, 2)
        ):
            raise AssertionError("final z blocker is not complete bipartite")
        part_histogram[tuple(sorted((len(first_part), len(second_part))))] += 1
    audit = FinalA5Audit(
        candidate_minimal_a5_count=len(candidates),
        unordered_two_a5_families=(
            len(candidates) * (len(candidates) - 1) // 2
        ),
        tau_z_at_least_two_families=len(valid),
        all_tau_valid_rows_avoid_active=all_avoid_active,
        all_tau_valid_intersections_empty=all_empty,
        all_tau_valid_unions_are_isolates=all_union,
        z_edge_count_histogram=tuple(sorted(edge_histogram.items())),
        z_part_size_histogram=tuple(sorted(part_histogram.items())),
        example_rows=valid[0] if valid else (),
    )
    if not audit.verified:
        raise AssertionError("record-16 final A5 enumeration changed")
    return audit


@dataclass(frozen=True)
class FullBlockerAudit:
    facet_complement_count: int
    minimal_transversals_by_size: tuple[tuple[int, tuple[int, ...]], ...]

    @property
    def histogram(self) -> tuple[tuple[int, int], ...]:
        return tuple(
            (size, len(values))
            for size, values in self.minimal_transversals_by_size
        )

    @property
    def verified(self) -> bool:
        by_size = dict(self.minimal_transversals_by_size)
        return (
            self.facet_complement_count == 102
            and self.histogram == ((2, 5), (4, 4), (5, 2), (6, 2), (7, 2))
            and tuple(sorted(by_size[2]))
            == tuple(sorted(mask(edge) for edge in (
                (0, 12),
                (1, 12),
                (0, 13),
                (1, 13),
                (12, 13),
            )))
            and tuple(sorted(by_size[5]))
            == tuple(
                sorted(
                    (
                        mask((8, 9, 10, 11, 12)),
                        mask((8, 9, 10, 11, 13)),
                    )
                )
            )
            and bool(by_size.get(6))
            and bool(by_size.get(7))
        )

    def to_record(self) -> dict[str, object]:
        return {
            "facet_complement_count": self.facet_complement_count,
            "minimal_transversal_histogram": {
                str(size): len(values)
                for size, values in self.minimal_transversals_by_size
            },
            "minimal_transversal_masks": {
                str(size): list(values)
                for size, values in self.minimal_transversals_by_size
            },
            "minimal_transversal_vertices": {
                str(size): [list(vertices(value)) for value in values]
                for size, values in self.minimal_transversals_by_size
            },
            "special_a5_support_violation": [
                [8, 9, 10, 11, 12],
                [8, 9, 10, 11, 13],
            ],
            "lanner_order_bound_violations": {
                str(size): len(values)
                for size, values in self.minimal_transversals_by_size
                if size > 5
            },
            "verified": self.verified,
        }


def audit_final_full_blocker(final_graph: XGraphAudit) -> FullBlockerAudit:
    """Dualize the final case's 102 facets back to all minimal nonfaces."""

    order, edges = decode_graph6(final_graph.graph6)
    if order != 6:
        raise ValueError("final record-16 X graph must have order six")
    y_triples = _x_y_triples(order, edges)
    special_mask = mask(SPECIAL_VERTICES)
    active_mask = mask(ACTIVE_VERTICES)
    large_side = tuple(range(2, 8))
    small_side = tuple(range(8, 12))
    facets = [special_mask | triple for triple in y_triples]
    for first in large_side:
        for second in small_side:
            pair = (1 << first) | (1 << second)
            for special in SPECIAL_VERTICES:
                facets.append(active_mask | pair | (1 << special))
    facets = sorted(set(facets))
    minimal = []
    for candidate in range(1, 1 << TARGET_VERTEX_COUNT):
        if not all(candidate & facet for facet in facets):
            continue
        if all(
            not all((candidate ^ (1 << vertex)) & facet for facet in facets)
            for vertex in vertices(candidate)
        ):
            minimal.append(candidate)
    by_size = tuple(
        (size, tuple(value for value in minimal if value.bit_count() == size))
        for size in sorted({value.bit_count() for value in minimal})
    )
    audit = FullBlockerAudit(
        facet_complement_count=len(facets),
        minimal_transversals_by_size=by_size,
    )
    if not audit.verified:
        raise AssertionError("record-16 final full blocker replay changed")
    return audit


def audit_record16_positive_branches(
    schedules: Sequence[Mapping[str, object]],
    x_archive_path: str | Path,
) -> PositiveBranchAudit:
    support = audit_positive_support()
    cases = positive_cut_cases(schedules)
    archive = audit_x_graph_archive(x_archive_path)
    final_graph = archive.target_graphs[(6, 54)][0]
    forced_a4 = final_graph.a4_blocker_masks
    shadow = {
        face | (1 << vertex)
        for face in forced_a4
        for vertex in WORK_VERTICES
        if not face & (1 << vertex)
    }
    schedule_by_a = {
        int(schedule["missing_four_count"]): schedule for schedule in schedules
    }
    offset = int(schedule_by_a[4]["missing_five_count_offset"])
    final_a5_audit = audit_final_extra_a5s(forced_a4, shadow)
    full_blocker_audit = audit_final_full_blocker(final_graph)
    audit = PositiveBranchAudit(
        support=support,
        cut_cases=cases,
        x_archive=archive,
        forced_a4_masks=forced_a4,
        forced_five_shadow_size=len(shadow),
        forced_additional_minimal_a5_count=offset - len(shadow),
        final_a5_audit=final_a5_audit,
        full_blocker_audit=full_blocker_audit,
        ridge_z_cut_part_sizes=(4, 6),
    )
    if not audit.verified:
        raise AssertionError("record-16 positive branch certificate changed")
    return audit

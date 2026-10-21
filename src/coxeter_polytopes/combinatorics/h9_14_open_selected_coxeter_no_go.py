"""Compact-Coxeter label obstructions for seven selected H9/14 nerves.

Records 0, 1, 3, 4, 5, 7, and 9 each have a separately certified ordinary
convex realization.  This module uses only their exact minimal-nonface
clutters and rank-nine facets.  Every order-four/five minimal nonface is
joined with the complete labelled classical Lanner table, and every induced
restriction on an actual rank-nine vertex is required to be elliptic.

Pairs not fixed by a joined Lanner support are temporarily right angled.  In
normalized Gram form ``G = I-A``, adding any finite non-right label increases
the symmetric nonnegative matrix ``A`` entrywise.  Perron--Frobenius
monotonicity therefore makes every sparse non-ellipticity rejection lossless.

The result is selected-nerve only.  It neither closes another completion of
the same graph shadow nor addresses noncompact finite-volume chambers.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
from importlib import import_module
from itertools import combinations
import json
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    Pair,
    is_finite_coxeter_labeling,
    lanner_placement_codes,
    pair,
)


TARGET_VERTEX_COUNT = 14
TARGET_DIMENSION = 9
FULL_MASK = (1 << TARGET_VERTEX_COUNT) - 1
TARGET_RECORDS = (0, 1, 3, 4, 5, 7, 9)
LabelState = tuple[tuple[Pair, int], ...]


def _state_digest(states: Sequence[LabelState]) -> str:
    payload = [
        [[left, right, label] for (left, right), label in state]
        for state in sorted(states)
    ]
    return sha256(
        json.dumps(payload, separators=(",", ":"), ensure_ascii=True).encode(
            "utf-8"
        )
    ).hexdigest()


@lru_cache(maxsize=None)
def _placement_options(support: tuple[int, ...]) -> tuple[LabelState, ...]:
    support = tuple(sorted(support))
    edges = tuple(combinations(support, 2))
    return tuple(
        tuple(zip(edges, code, strict=True))
        for code in lanner_placement_codes(len(support))
    )


@lru_cache(maxsize=262_144)
def _finite_code(order: int, code: tuple[int, ...]) -> bool:
    labels = dict(zip(combinations(range(order), 2), code, strict=True))
    return is_finite_coxeter_labeling(range(order), labels)


def _restriction_is_finite(
    restriction: tuple[int, ...], labels: Mapping[Pair, int]
) -> bool:
    code = tuple(
        labels[pair(left, right)]
        for left, right in combinations(restriction, 2)
    )
    return _finite_code(len(restriction), code)


def _next_support_index(
    supports: Sequence[tuple[int, ...]],
    remaining: Iterable[int],
    cluster_union: Iterable[int],
) -> int:
    """Frozen outcome-independent overlap ordering used by both replays."""

    choices = tuple(sorted(remaining))
    if not choices:
        raise ValueError("no Lanner support remains")
    union = set(cluster_union)
    if not union:
        return choices[0]
    return max(
        choices,
        key=lambda index: (
            len(set(supports[index]) & union),
            -len(set(supports[index]) - union),
            -index,
        ),
    )


@dataclass(frozen=True, slots=True)
class SelectedNerve:
    record_index: int
    missing_edges: tuple[tuple[int, int], ...]
    supports: tuple[tuple[int, ...], ...]
    rank_nine_vertices: tuple[tuple[int, ...], ...]

    @property
    def support_histogram(self) -> tuple[tuple[int, int], ...]:
        return tuple(sorted(Counter(map(len, self.supports)).items()))


def load_selected_nerve(record_index: int) -> SelectedNerve:
    if record_index not in TARGET_RECORDS:
        raise ValueError(f"record must lie in {TARGET_RECORDS}")
    module = import_module(
        f"coxeter_polytopes.combinatorics.h9_14_record{record_index}"
    )
    missing_edges = tuple(tuple(face) for face in module.TARGET_MISSING_EDGES)
    missing_fours = tuple(
        sorted(tuple(face) for face in getattr(module, "TARGET_MISSING_FOURS", ()))
    )
    missing_fives = tuple(sorted(tuple(face) for face in module.TARGET_MISSING_FIVES))
    complements = tuple(sorted(module.target_facet_complements()))
    vertices = tuple(
        sorted(
            tuple(
                vertex
                for vertex in range(TARGET_VERTEX_COUNT)
                if (FULL_MASK ^ complement) >> vertex & 1
            )
            for complement in complements
        )
    )
    if any(len(vertex) != TARGET_DIMENSION for vertex in vertices):
        raise AssertionError("a selected facet is not rank nine")
    nonfaces = tuple(map(frozenset, (*missing_edges, *missing_fours, *missing_fives)))
    reconstructed = tuple(
        face
        for face in combinations(range(TARGET_VERTEX_COUNT), TARGET_DIMENSION)
        if not any(nonface <= set(face) for nonface in nonfaces)
    )
    if vertices != reconstructed:
        raise AssertionError("the minimal nonfaces do not recover the selected facets")
    return SelectedNerve(
        record_index=record_index,
        missing_edges=missing_edges,
        supports=missing_fours + missing_fives,
        rank_nine_vertices=vertices,
    )


@dataclass(frozen=True, slots=True)
class LabelStage:
    support_index: int
    support: tuple[int, ...]
    cluster_union: tuple[int, ...]
    restriction_size_histogram: tuple[tuple[int, int], ...]
    lanner_placement_count: int
    raw_extension_count: int
    compatible_extension_count: int
    passing_state_count: int
    passing_state_sha256: str

    @property
    def count_signature(self) -> tuple[int, int, int]:
        return (
            self.raw_extension_count,
            self.compatible_extension_count,
            self.passing_state_count,
        )

    def to_record(self) -> dict[str, object]:
        return {
            "support_index": self.support_index,
            "support": list(self.support),
            "support_order": len(self.support),
            "cumulative_cluster_union": list(self.cluster_union),
            "actual_vertex_restriction_size_histogram": {
                str(size): count for size, count in self.restriction_size_histogram
            },
            "complete_lanner_placement_count": self.lanner_placement_count,
            "raw_extension_count": self.raw_extension_count,
            "compatible_lanner_extension_count": self.compatible_extension_count,
            "elliptic_actual_vertex_survivor_count": self.passing_state_count,
            "passing_state_sha256": self.passing_state_sha256,
        }


def _extend_stage(
    nerve: SelectedNerve,
    states: Sequence[LabelState],
    previous_union: Iterable[int],
    support_index: int,
) -> tuple[tuple[LabelState, ...], LabelStage]:
    support = nerve.supports[support_index]
    cluster_union = tuple(sorted(set(previous_union) | set(support)))
    cluster_set = set(cluster_union)
    restrictions = tuple(
        sorted(
            {
                tuple(vertex for vertex in face if vertex in cluster_set)
                for face in nerve.rank_nine_vertices
            },
            key=lambda value: (len(value), value),
        )
    )
    if any(
        set(dotted) <= set(restriction)
        for dotted in nerve.missing_edges
        for restriction in restrictions
    ):
        raise AssertionError("an actual vertex restriction contains a dotted pair")

    options = _placement_options(support)
    passing: set[LabelState] = set()
    compatible = 0
    for state in states:
        old = dict(state)
        for option in options:
            if not all(old.get(edge, label) == label for edge, label in option):
                continue
            compatible += 1
            merged = old | dict(option)
            labels = {
                pair(left, right): 2
                for left, right in combinations(cluster_union, 2)
            }
            labels.update(merged)
            if all(
                _restriction_is_finite(restriction, labels)
                for restriction in restrictions
            ):
                passing.add(tuple(sorted(merged.items())))
    passing_states = tuple(sorted(passing))
    stage = LabelStage(
        support_index=support_index,
        support=support,
        cluster_union=cluster_union,
        restriction_size_histogram=tuple(
            sorted(Counter(map(len, restrictions)).items())
        ),
        lanner_placement_count=len(options),
        raw_extension_count=len(states) * len(options),
        compatible_extension_count=compatible,
        passing_state_count=len(passing_states),
        passing_state_sha256=_state_digest(passing_states),
    )
    return passing_states, stage


@dataclass(frozen=True, slots=True)
class SelectedNerveLabelAudit:
    nerve: SelectedNerve
    stages: tuple[LabelStage, ...]

    @property
    def compact_coxeter_no_go(self) -> bool:
        return bool(self.stages) and self.stages[-1].passing_state_count == 0

    @property
    def verified(self) -> bool:
        remaining = set(range(len(self.nerve.supports)))
        union: tuple[int, ...] = ()
        for stage in self.stages:
            expected = _next_support_index(self.nerve.supports, remaining, union)
            if (
                stage.support_index != expected
                or stage.support != self.nerve.supports[expected]
                or stage.lanner_placement_count
                != len(lanner_placement_codes(len(stage.support)))
            ):
                return False
            remaining.remove(expected)
            union = stage.cluster_union
        return self.compact_coxeter_no_go

    def to_record(self) -> dict[str, object]:
        return {
            "record_index": self.nerve.record_index,
            "rank_nine_vertex_count": len(self.nerve.rank_nine_vertices),
            "minimal_nonface_histogram": {
                "2": len(self.nerve.missing_edges),
                **{
                    str(size): count
                    for size, count in self.nerve.support_histogram
                },
            },
            "decisive_support_indices": [
                stage.support_index for stage in self.stages
            ],
            "incremental_lanner_and_actual_vertex_gate": [
                stage.to_record() for stage in self.stages
            ],
            "selected_nerve_compact_coxeter_no_go": self.compact_coxeter_no_go,
            "verified": self.verified,
        }


def audit_selected_nerve(record_index: int) -> SelectedNerveLabelAudit:
    nerve = load_selected_nerve(record_index)
    remaining = set(range(len(nerve.supports)))
    states: tuple[LabelState, ...] = ((),)
    union: tuple[int, ...] = ()
    stages: list[LabelStage] = []
    while states and remaining:
        support_index = _next_support_index(nerve.supports, remaining, union)
        states, stage = _extend_stage(
            nerve, states, union, support_index
        )
        stages.append(stage)
        remaining.remove(support_index)
        union = stage.cluster_union
    audit = SelectedNerveLabelAudit(nerve=nerve, stages=tuple(stages))
    if not audit.verified:
        raise AssertionError(
            f"record {record_index} did not reach an exact selected-nerve no-go"
        )
    return audit


@dataclass(frozen=True, slots=True)
class AllSelectedNervesAudit:
    audits: tuple[SelectedNerveLabelAudit, ...]

    @property
    def all_seven_selected_nerves_excluded(self) -> bool:
        return (
            tuple(audit.nerve.record_index for audit in self.audits)
            == TARGET_RECORDS
            and all(audit.compact_coxeter_no_go for audit in self.audits)
        )

    @property
    def verified(self) -> bool:
        return self.all_seven_selected_nerves_excluded and all(
            audit.verified for audit in self.audits
        )

    def to_record(self) -> dict[str, object]:
        return {
            "record_indices": list(TARGET_RECORDS),
            "selected_nerves_excluded": len(self.audits),
            "finite_label_alphabet": [2, 3, 4, 5],
            "complete_lanner_tables": {
                "order_four_placements": len(lanner_placement_codes(4)),
                "order_five_placements": len(lanner_placement_codes(5)),
            },
            "record_gates": [audit.to_record() for audit in self.audits],
            "all_seven_selected_nerves_excluded": (
                self.all_seven_selected_nerves_excluded
            ),
            "verified": self.verified,
        }


def audit_all_open_selected_nerves() -> AllSelectedNervesAudit:
    audit = AllSelectedNervesAudit(
        audits=tuple(audit_selected_nerve(record) for record in TARGET_RECORDS)
    )
    if not audit.verified:
        raise AssertionError("the seven selected-nerve no-gos did not verify")
    return audit


__all__ = [
    "AllSelectedNervesAudit",
    "LabelStage",
    "SelectedNerve",
    "SelectedNerveLabelAudit",
    "TARGET_RECORDS",
    "audit_all_open_selected_nerves",
    "audit_selected_nerve",
    "load_selected_nerve",
]

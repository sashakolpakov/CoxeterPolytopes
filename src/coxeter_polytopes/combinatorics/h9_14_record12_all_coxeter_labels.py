"""Whole-shadow compact-Coxeter obstruction for H9/14 record 12.

The corrected record-12 incidence census contains exactly 371 residual
complexes.  This module reconstructs every one from its two archived balanced
branch matrices, recovers every minimal nonface and every rank-nine face, and
then applies a deterministic compact Coxeter label gate.

For each selected order-four/five minimal nonface the complete labelled
Lanner table is joined to the current states.  After every join, every
restriction induced by an actual rank-nine face must be elliptic.  Pairs not
yet fixed by a selected Lanner support are temporarily right angled.  This is
proof-safe: in normalized Gram form ``G=I-A``, adding any finite non-right
edge increases the symmetric nonnegative matrix ``A`` entrywise, so
Perron--Frobenius monotonicity cannot repair a non-elliptic restriction.

Burcroff's global compact-H9/14 reduction (EJC 120 (2024), Lemma 10.10)
excludes multi-multiple edges.  Thus ``2,3,4,5`` is the complete finite label
alphabet in this target; dotted pairs remain continuous but never occur in
an actual vertex.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    Pair,
    is_finite_coxeter_labeling,
    lanner_placement_codes,
    pair,
)
from coxeter_polytopes.combinatorics.h9_14_record12_residual import (
    EXPECTED_FULL_ORBIT_COUNT,
    EXPECTED_SCHEDULE_ORBIT_COUNTS,
    FULL_MASK,
    TARGET_DIMENSION,
    TARGET_EDGE_MASKS,
    TARGET_VERTEX_COUNT,
    face_table_from_complements,
    minimal_nonfaces_from_table,
    selected_complements_from_matrices,
    verify_archive,
    vertices,
)


EXPECTED_INPUT_FILE_SHA256 = (
    "faca367738bba58afc0ec78a43e9d4aa4b02721802c5d9e77b7d155c8e7faa05"
)
EXPECTED_INPUT_CERTIFICATE_SHA256 = (
    "50ca012334e7a4973a4c29bde4226d4460bc8ade2d644a3b9dbf56b1bf92e8f4"
)


LabelState = tuple[tuple[Pair, int], ...]


def _canonical_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(
            payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("utf-8")
    ).hexdigest()


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


def _records_digest(records: Sequence[Mapping[str, object]]) -> str:
    return sha256(
        json.dumps(
            list(records), sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("utf-8")
    ).hexdigest()


def _rank_nine_vertices(
    minimal_nonfaces: Iterable[Iterable[int]],
) -> tuple[tuple[int, ...], ...]:
    forbidden = tuple(frozenset(face) for face in minimal_nonfaces)
    return tuple(
        facets
        for facets in combinations(range(TARGET_VERTEX_COUNT), TARGET_DIMENSION)
        if not any(face <= set(facets) for face in forbidden)
    )


@dataclass(frozen=True, slots=True)
class Record12ResidualNerve:
    orbit_index: int
    x_matrix: int
    z_matrix: int
    missing_four_count: int
    relative_embedding_multiplicity: int
    archived_facet_count: int
    selected_complement_sha256: str
    minimal_nonfaces_by_size: tuple[
        tuple[int, tuple[tuple[int, ...], ...]], ...
    ]
    rank_nine_vertices: tuple[tuple[int, ...], ...]

    @property
    def higher_minimal_nonfaces(self) -> tuple[tuple[int, ...], ...]:
        by_size = dict(self.minimal_nonfaces_by_size)
        return by_size.get(4, ()) + by_size.get(5, ())

    @property
    def dotted_pairs(self) -> tuple[Pair, ...]:
        return tuple(pair(*face) for face in dict(self.minimal_nonfaces_by_size)[2])

    @property
    def signature(self) -> tuple[object, ...]:
        return (
            self.orbit_index,
            self.x_matrix,
            self.z_matrix,
            self.missing_four_count,
            self.relative_embedding_multiplicity,
            self.archived_facet_count,
            self.selected_complement_sha256,
        )


def load_record12_residual_nerves(
    archive_path: Path,
) -> tuple[Record12ResidualNerve, ...]:
    """Load and reconstruct all 371 frozen residual complexes exactly."""

    raw = archive_path.read_bytes()
    if sha256(raw).hexdigest() != EXPECTED_INPUT_FILE_SHA256:
        raise ValueError("the frozen record-12 residual archive changed")
    document = json.loads(raw)
    if document.get("certificate_sha256") != EXPECTED_INPUT_CERTIFICATE_SHA256:
        raise ValueError("the record-12 residual certificate digest changed")
    if _canonical_digest(document) != EXPECTED_INPUT_CERTIFICATE_SHA256:
        raise ValueError("the record-12 residual canonical digest does not replay")
    summary = verify_archive(document)
    if (
        summary["residual_orbits"] != EXPECTED_FULL_ORBIT_COUNT
        or summary["schedule_orbit_counts"] != EXPECTED_SCHEDULE_ORBIT_COUNTS
    ):
        raise ValueError("the exact record-12 residual census changed")

    nerves: list[Record12ResidualNerve] = []
    for orbit_index, orbit in enumerate(document["residual_orbits"]):
        x_matrix = int(orbit["x_matrix"])
        z_matrix = int(orbit["z_matrix"])
        complements = selected_complements_from_matrices(x_matrix, z_matrix)
        table = face_table_from_complements(complements)
        nonface_masks = minimal_nonfaces_from_table(table)
        by_size = tuple(
            (
                size,
                tuple(vertices(value) for value in rows),
            )
            for size, rows in nonface_masks.items()
        )
        all_nonfaces = tuple(
            face for _size, faces in by_size for face in faces
        )
        reconstructed_vertices = _rank_nine_vertices(all_nonfaces)
        direct_vertices = tuple(
            sorted(vertices(FULL_MASK ^ complement) for complement in complements)
        )
        if reconstructed_vertices != direct_vertices:
            raise AssertionError("minimal nonfaces do not recover every rank-nine face")

        nerve = Record12ResidualNerve(
            orbit_index=orbit_index,
            x_matrix=x_matrix,
            z_matrix=z_matrix,
            missing_four_count=int(orbit["missing_four_count"]),
            relative_embedding_multiplicity=int(
                orbit["relative_embedding_multiplicity"]
            ),
            archived_facet_count=int(orbit["facet_count"]),
            selected_complement_sha256=str(
                orbit["selected_complement_sha256"]
            ),
            minimal_nonfaces_by_size=by_size,
            rank_nine_vertices=direct_vertices,
        )
        histogram = {size: len(faces) for size, faces in by_size}
        archived_histogram = {
            int(size): int(count)
            for size, count in orbit["minimal_nonface_histogram"].items()
        }
        complement_digest = sha256(
            b"".join(value.to_bytes(2, "little") for value in complements)
        ).hexdigest()
        if (
            set(histogram) - {2, 4, 5}
            or histogram != archived_histogram
            or set(nonface_masks[2]) != set(TARGET_EDGE_MASKS)
            or histogram.get(4, 0) != nerve.missing_four_count
            or len(direct_vertices) != nerve.archived_facet_count
            or complement_digest != nerve.selected_complement_sha256
        ):
            raise AssertionError(f"record-12 residual orbit {orbit_index} changed")
        nerves.append(nerve)
    if len(nerves) != EXPECTED_FULL_ORBIT_COUNT:
        raise AssertionError("the record-12 residual orbit count changed")
    return tuple(nerves)


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


@dataclass(frozen=True, slots=True)
class Record12LabelStage:
    support_index: int
    support: tuple[int, ...]
    lanner_placement_count: int
    cluster_union: tuple[int, ...]
    restriction_size_histogram: tuple[tuple[int, int], ...]
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
            "complete_lanner_placement_count": self.lanner_placement_count,
            "cumulative_cluster_union": list(self.cluster_union),
            "actual_vertex_restriction_size_histogram": {
                str(size): count for size, count in self.restriction_size_histogram
            },
            "raw_extension_count": self.raw_extension_count,
            "compatible_lanner_extension_count": self.compatible_extension_count,
            "elliptic_actual_vertex_survivor_count": self.passing_state_count,
            "passing_state_sha256": self.passing_state_sha256,
        }


def _extend_stage(
    nerve: Record12ResidualNerve,
    states: Sequence[LabelState],
    previous_union: Iterable[int],
    support_index: int,
    support: tuple[int, ...],
) -> tuple[tuple[LabelState, ...], Record12LabelStage]:
    cluster_union = tuple(sorted(set(previous_union) | set(support)))
    cluster_set = set(cluster_union)
    restrictions = tuple(
        sorted(
            {
                tuple(facet for facet in vertex if facet in cluster_set)
                for vertex in nerve.rank_nine_vertices
            },
            key=lambda value: (len(value), value),
        )
    )
    if any(
        set(dotted) <= set(restriction)
        for dotted in nerve.dotted_pairs
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
    stage = Record12LabelStage(
        support_index=support_index,
        support=support,
        lanner_placement_count=len(options),
        cluster_union=cluster_union,
        restriction_size_histogram=tuple(
            sorted(Counter(map(len, restrictions)).items())
        ),
        raw_extension_count=len(states) * len(options),
        compatible_extension_count=compatible,
        passing_state_count=len(passing_states),
        passing_state_sha256=_state_digest(passing_states),
    )
    return passing_states, stage


def _next_support_index(
    supports: Sequence[tuple[int, ...]],
    remaining: Iterable[int],
    cluster_union: Iterable[int],
) -> int:
    """Choose the next support by one frozen, state/outcome-independent rule.

    Start with support zero.  Afterwards maximize overlap with the cumulative
    cluster, then minimize the number of new vertices, then use the smallest
    archived support index.  The rule affects speed only: every chosen support
    is joined with its complete placement table.
    """

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
class Record12OrbitLabelAudit:
    nerve: Record12ResidualNerve
    stages: tuple[Record12LabelStage, ...]

    @property
    def decisive_support_indices(self) -> tuple[int, ...]:
        return tuple(stage.support_index for stage in self.stages)

    @property
    def compact_coxeter_no_go(self) -> bool:
        return bool(self.stages) and self.stages[-1].passing_state_count == 0

    @property
    def verified(self) -> bool:
        supports = self.nerve.higher_minimal_nonfaces
        remaining = set(range(len(supports)))
        union: tuple[int, ...] = ()
        for stage in self.stages:
            expected = _next_support_index(supports, remaining, union)
            if (
                stage.support_index != expected
                or stage.support != supports[expected]
                or stage.lanner_placement_count
                != len(lanner_placement_codes(len(stage.support)))
            ):
                return False
            remaining.remove(expected)
            union = stage.cluster_union
        return self.compact_coxeter_no_go

    def to_record(self) -> dict[str, object]:
        return {
            "orbit_index": self.nerve.orbit_index,
            "orbit_signature": {
                "x_matrix": self.nerve.x_matrix,
                "z_matrix": self.nerve.z_matrix,
                "missing_four_count": self.nerve.missing_four_count,
                "relative_embedding_multiplicity": (
                    self.nerve.relative_embedding_multiplicity
                ),
                "facet_count": self.nerve.archived_facet_count,
                "selected_complement_sha256": (
                    self.nerve.selected_complement_sha256
                ),
            },
            "rank_nine_vertex_count": len(self.nerve.rank_nine_vertices),
            "minimal_nonface_histogram": {
                str(size): len(faces)
                for size, faces in self.nerve.minimal_nonfaces_by_size
            },
            "decisive_support_indices": list(self.decisive_support_indices),
            "incremental_lanner_and_actual_vertex_gate": [
                stage.to_record() for stage in self.stages
            ],
            "compact_coxeter_no_go": self.compact_coxeter_no_go,
            "verified": self.verified,
        }


def audit_record12_orbit_labels(
    nerve: Record12ResidualNerve,
) -> Record12OrbitLabelAudit:
    supports = nerve.higher_minimal_nonfaces
    remaining = set(range(len(supports)))
    states: tuple[LabelState, ...] = ((),)
    cluster_union: tuple[int, ...] = ()
    stages: list[Record12LabelStage] = []
    while states and remaining:
        support_index = _next_support_index(supports, remaining, cluster_union)
        states, stage = _extend_stage(
            nerve,
            states,
            cluster_union,
            support_index,
            supports[support_index],
        )
        stages.append(stage)
        remaining.remove(support_index)
        cluster_union = stage.cluster_union
    audit = Record12OrbitLabelAudit(nerve=nerve, stages=tuple(stages))
    if not audit.verified:
        raise AssertionError(
            f"record-12 orbit {nerve.orbit_index} did not reach an exact no-go"
        )
    return audit


@dataclass(frozen=True, slots=True)
class Record12AllOrbitLabelAudit:
    input_file_sha256: str
    input_certificate_sha256: str
    orbits: tuple[Record12OrbitLabelAudit, ...]

    @property
    def whole_shadow_compact_coxeter_no_go(self) -> bool:
        return len(self.orbits) == EXPECTED_FULL_ORBIT_COUNT and all(
            orbit.compact_coxeter_no_go for orbit in self.orbits
        )

    @property
    def support_count_histogram(self) -> tuple[tuple[int, int], ...]:
        return tuple(sorted(Counter(len(orbit.stages) for orbit in self.orbits).items()))

    @property
    def schedule_exclusion_histogram(self) -> tuple[tuple[int, int], ...]:
        return tuple(
            sorted(
                Counter(
                    orbit.nerve.missing_four_count for orbit in self.orbits
                ).items()
            )
        )

    @property
    def orbit_transcript_sha256(self) -> str:
        return _records_digest([orbit.to_record() for orbit in self.orbits])

    @property
    def verified(self) -> bool:
        return (
            self.input_file_sha256 == EXPECTED_INPUT_FILE_SHA256
            and self.input_certificate_sha256
            == EXPECTED_INPUT_CERTIFICATE_SHA256
            and tuple(dict(self.schedule_exclusion_histogram).items())
            == tuple(EXPECTED_SCHEDULE_ORBIT_COUNTS.items())
            and all(orbit.verified for orbit in self.orbits)
            and self.whole_shadow_compact_coxeter_no_go
        )

    def to_record(self) -> dict[str, object]:
        return {
            "input_residual_certificate": {
                "file_sha256": self.input_file_sha256,
                "certificate_sha256": self.input_certificate_sha256,
                "exact_residual_orbit_count": len(self.orbits),
                "schedule_orbit_counts": {
                    str(key): value for key, value in self.schedule_exclusion_histogram
                },
            },
            "finite_label_alphabet": [2, 3, 4, 5],
            "complete_lanner_tables": {
                "order_four_placements": len(lanner_placement_codes(4)),
                "order_five_placements": len(lanner_placement_codes(5)),
            },
            "deterministic_support_rule": (
                "start with archived support 0; then maximize overlap with "
                "the cumulative cluster, minimize new vertices, and break "
                "ties by smallest archived support index"
            ),
            "orbit_gates": [orbit.to_record() for orbit in self.orbits],
            "support_count_histogram": {
                str(count): multiplicity
                for count, multiplicity in self.support_count_histogram
            },
            "orbit_transcript_sha256": self.orbit_transcript_sha256,
            "orbits_excluded": sum(
                orbit.compact_coxeter_no_go for orbit in self.orbits
            ),
            "whole_shadow_compact_coxeter_no_go": (
                self.whole_shadow_compact_coxeter_no_go
            ),
            "verified": self.verified,
        }


@lru_cache(maxsize=1)
def audit_record12_all_coxeter_labels(
    archive_path: Path,
) -> Record12AllOrbitLabelAudit:
    nerves = load_record12_residual_nerves(archive_path)
    audit = Record12AllOrbitLabelAudit(
        input_file_sha256=sha256(archive_path.read_bytes()).hexdigest(),
        input_certificate_sha256=EXPECTED_INPUT_CERTIFICATE_SHA256,
        orbits=tuple(audit_record12_orbit_labels(nerve) for nerve in nerves),
    )
    if not audit.verified:
        raise AssertionError("record-12 all-orbit compact-Coxeter gate changed")
    return audit


__all__ = [
    "EXPECTED_INPUT_CERTIFICATE_SHA256",
    "EXPECTED_INPUT_FILE_SHA256",
    "Record12AllOrbitLabelAudit",
    "Record12LabelStage",
    "Record12OrbitLabelAudit",
    "Record12ResidualNerve",
    "audit_record12_all_coxeter_labels",
    "audit_record12_orbit_labels",
    "load_record12_residual_nerves",
]

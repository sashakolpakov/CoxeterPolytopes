"""Compact-Coxeter obstruction for every residual record-17 nerve.

The exact whole-shadow certificate for H9/14 record 17 leaves eight
incidence-complex orbits.  This module treats *all eight*: it extracts each
minimal-nonface clutter from that frozen certificate, reconstructs every
rank-nine face, and incrementally joins complete labelled Lanner tables on a
small decisive collection of order-four/five minimal nonfaces.

After each join, every restriction induced by an actual rank-nine face is
required to be an elliptic finite Coxeter diagram.  Pairs not yet forced by a
Lanner support are put at right angle.  This pruning is lossless for a
nonexistence proof: in normalized Gram form ``G = I-A``, changing a right
angle to any finite non-right angle increases the symmetric nonnegative
matrix ``A`` entrywise, so Perron--Frobenius monotonicity cannot turn a
non-elliptic principal diagram into an elliptic one.

Five orbits die after three supports and three after four.  Consequently no
one of the eight residual incidence complexes is the nerve of a compact
hyperbolic Coxeter polytope.  This is a whole-shadow compact theorem.  It does
not assert that the seven non-selected complexes are PL or polytopal, and it
does not exclude noncompact finite-volume chambers with different ideal-face
data.
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


FACET_COUNT = 14
DIMENSION = 9
EXPECTED_INPUT_FILE_SHA256 = (
    "db613bd1b3f6606f7b140f4183f5fe2897db6ef8292d72670e38523f7c397d94"
)
EXPECTED_INPUT_CERTIFICATE_SHA256 = (
    "b98e58b390cbf8d9b40e08da16da29f78a2f04805bb0b13d7f2efb623b9afb51"
)

# Supports are indexed after concatenating order-four and order-five minimal
# nonfaces in the order stored by the exact whole-shadow certificate.
DECISIVE_SUPPORT_INDICES = (
    (1, 2, 0),
    (2, 3, 4, 0),
    (2, 3, 0),
    (2, 3, 4),
    (0, 1, 2),
    (3, 4, 5, 0),
    (4, 5, 6, 0),
    (1, 2, 0),
)

EXPECTED_ORBIT_SIGNATURES = (
    (1, "H???FB_", (0, 1, 2, 3, 4), 6, 108),
    (1, "H???FB_", (3, 4, 5, 7, 8), 4, 108),
    (2, "H???F?^", (3, 4, 5, 6, 7), 1, 106),
    (2, "H???Fbo", (4, 5, 6, 7, 8), 1, 106),
    (0, "H??EEBo", (4, 5, 6, 7, 8), 1, 110),
    (1, "H??EEBw", (0, 1, 2, 3, 4), 1, 108),
    (1, "H??FFB_", (0, 1, 2, 3, 4), 6, 108),
    (2, "H??FF?^", (3, 4, 5, 6, 7), 1, 106),
)

# (raw extensions, compatible extensions, elliptic survivors) at every
# cumulative support stage.  These counts are independently replayed by the
# hostile verifier; they also prevent a weakened or reordered search from
# silently replacing the frozen proof.
EXPECTED_STAGE_COUNTS = (
    ((420, 420, 420), (176_400, 11_634, 24), (2_592, 432, 0)),
    (
        (420, 420, 420),
        (176_400, 1_332, 96),
        (40_320, 384, 18),
        (1_944, 1_944, 0),
    ),
    ((420, 420, 420), (176_400, 1_332, 96), (10_368, 10_368, 0)),
    ((420, 420, 420), (176_400, 1_332, 96), (40_320, 384, 0)),
    ((420, 420, 420), (176_400, 1_332, 96), (40_320, 384, 0)),
    (
        (420, 420, 420),
        (176_400, 11_634, 24),
        (10_080, 88, 8),
        (864, 32, 0),
    ),
    (
        (420, 420, 420),
        (176_400, 1_332, 96),
        (40_320, 384, 18),
        (1_944, 828, 0),
    ),
    ((108, 108, 108), (45_360, 1_542, 84), (9_072, 9_072, 0)),
)

EXPECTED_STAGE_DIGESTS = (
    (
        "20fec51dea65e9dc5594f7d63f8dec7e81506b7fee8a65752f815d769542b8d1",
        "f2739a77eb5ff316065b3698a2129fa0010c0544866b09817fe58e58b108133c",
        "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945",
    ),
    (
        "edd92955baeaa3c2cd81d7af99ef744334dde3813c2dff15e6ce68ace226081c",
        "1e1dfec7c2c5e50dc6004c589cd3b11b302850bd97dd2c9ff85c194b87fc9695",
        "de3c06f3354244de4b1fdc0065eebe44e25e5ba234fa5f9998b35bb14eb87526",
        "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945",
    ),
    (
        "b12aa301f700f0bb64097a85cd8849b3148d398a13099151fe581decdf1a2f9c",
        "ab8244d7ba0c74ec2157441f833f22e9a11552b367ff1e058a59ad1805dfd113",
        "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945",
    ),
    (
        "6653a78154113ed28fa8257d896c11e11f809876790ce5d7dba23a73f5082cc7",
        "8bad258bca95f7a4b51e88253ffac9d49a8b95807bf18462a6fcf1092142907f",
        "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945",
    ),
    (
        "6653a78154113ed28fa8257d896c11e11f809876790ce5d7dba23a73f5082cc7",
        "8bad258bca95f7a4b51e88253ffac9d49a8b95807bf18462a6fcf1092142907f",
        "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945",
    ),
    (
        "fb8d75842c42190ec9f8988edf38c2b873cbf1e015a371aca54c7e25c5cb6173",
        "93e2e8fce4eee6c988cea4e41edf6e3d5227504900a4555e96e9a7741d1de500",
        "f99ebfdec3e49744bc3ac6a62dd6c7d6f60de6f0d2d2542c2af34a0f43bc6c2f",
        "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945",
    ),
    (
        "bcdfeb8bb8016c4f00286d86c5971d41503bab9509f45ca9d56f5cdbf2bd2041",
        "fdba768376ba6464f80b61b59ebe6906478f1d0df6c0c11cb2d25b8912834dbb",
        "83f794045d9a0c420af0c1f78eedce07fb2777cbb785616190a4be406681db6c",
        "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945",
    ),
    (
        "5702c97e3ec3c514110ce9280c2db69e887c436a4b74e397b2d3c8d163b0cc03",
        "a833067ac762cb683a75d0f55916b304d88354e0692c17229f2e4ed1c8269b97",
        "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945",
    ),
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


@dataclass(frozen=True, slots=True)
class Record17ResidualNerve:
    orbit_index: int
    missing_four_count: int
    x_graph6: str
    p_local_vertices: tuple[int, ...]
    representation_count: int
    archived_facet_count: int
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
            self.missing_four_count,
            self.x_graph6,
            self.p_local_vertices,
            self.representation_count,
            self.archived_facet_count,
        )


def _rank_nine_vertices(
    minimal_nonfaces: Iterable[Iterable[int]],
) -> tuple[tuple[int, ...], ...]:
    forbidden = tuple(frozenset(face) for face in minimal_nonfaces)
    return tuple(
        facets
        for facets in combinations(range(FACET_COUNT), DIMENSION)
        if not any(face <= set(facets) for face in forbidden)
    )


def load_record17_residual_nerves(
    certificate_path: Path,
) -> tuple[Record17ResidualNerve, ...]:
    """Load and strictly validate the eight frozen whole-shadow orbits."""

    raw = certificate_path.read_bytes()
    if sha256(raw).hexdigest() != EXPECTED_INPUT_FILE_SHA256:
        raise ValueError("the frozen record-17 whole-shadow input file changed")
    document = json.loads(raw)
    if document.get("kind") != "h9_14_record17_exact_whole_shadow_residual":
        raise ValueError("unexpected record-17 whole-shadow certificate kind")
    if document.get("certificate_sha256") != EXPECTED_INPUT_CERTIFICATE_SHA256:
        raise ValueError("unexpected record-17 whole-shadow certificate digest")
    if _canonical_digest(document) != EXPECTED_INPUT_CERTIFICATE_SHA256:
        raise ValueError("record-17 whole-shadow canonical digest does not replay")

    enumeration = document["whole_shadow_audit"]["enumeration"]
    if enumeration.get("residual_orbit_count") != 8:
        raise ValueError("the record-17 residual orbit count changed")
    witnesses = enumeration["witnesses"]
    nerves: list[Record17ResidualNerve] = []
    for orbit_index, orbit in enumerate(enumeration["residual_orbits"]):
        signature = (
            int(orbit["missing_four_count"]),
            str(orbit["x_graph6"]),
            tuple(int(value) for value in orbit["p_local_vertices"]),
        )
        matches = [
            witness
            for witness in witnesses
            if (
                int(witness["missing_four_count"]),
                str(witness["x_graph6"]),
                tuple(int(value) for value in witness["p_local_vertices"]),
            )
            == signature
        ]
        if len(matches) != 1:
            raise ValueError("a residual orbit does not select one archived witness")
        witness = matches[0]
        if witness["p_local_vertices"] != witness["q_local_vertices"]:
            raise ValueError("the frozen record-17 residual cut equality changed")
        by_size = tuple(
            (
                int(size),
                tuple(tuple(int(vertex) for vertex in face) for face in faces),
            )
            for size, faces in sorted(
                witness["minimal_nonface_vertices_by_size"].items(),
                key=lambda item: int(item[0]),
            )
        )
        all_nonfaces = tuple(
            face for _size, faces in by_size for face in faces
        )
        vertices = _rank_nine_vertices(all_nonfaces)
        nerve = Record17ResidualNerve(
            orbit_index=orbit_index,
            missing_four_count=signature[0],
            x_graph6=signature[1],
            p_local_vertices=signature[2],
            representation_count=int(orbit["representation_count"]),
            archived_facet_count=int(orbit["facet_count"]),
            minimal_nonfaces_by_size=by_size,
            rank_nine_vertices=vertices,
        )
        histogram = {size: len(faces) for size, faces in by_size}
        archived_histogram = {
            int(size): int(count)
            for size, count in orbit["minimal_nonface_histogram"].items()
        }
        if (
            nerve.signature != EXPECTED_ORBIT_SIGNATURES[orbit_index]
            or len(vertices) != nerve.archived_facet_count
            or histogram != archived_histogram
            or set(histogram) - {2, 4, 5}
            or histogram.get(2) != 5
            or histogram.get(4, 0) != nerve.missing_four_count
        ):
            raise AssertionError("a frozen record-17 residual nerve changed")
        nerves.append(nerve)
    return tuple(nerves)


@lru_cache(maxsize=None)
def _placement_options(support: tuple[int, ...]) -> tuple[LabelState, ...]:
    support = tuple(sorted(support))
    edges = tuple(combinations(support, 2))
    return tuple(
        tuple(zip(edges, code, strict=True))
        for code in lanner_placement_codes(len(support))
    )


@lru_cache(maxsize=None)
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
class Record17AllLabelStage:
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
                str(size): count
                for size, count in self.restriction_size_histogram
            },
            "raw_extension_count": self.raw_extension_count,
            "compatible_lanner_extension_count": self.compatible_extension_count,
            "elliptic_actual_vertex_survivor_count": self.passing_state_count,
            "passing_state_sha256": self.passing_state_sha256,
        }


def _extend_stage(
    nerve: Record17ResidualNerve,
    states: Sequence[LabelState],
    previous_union: Iterable[int],
    support_index: int,
    support: tuple[int, ...],
) -> tuple[tuple[LabelState, ...], Record17AllLabelStage]:
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
    stage = Record17AllLabelStage(
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


@dataclass(frozen=True, slots=True)
class Record17OrbitLabelAudit:
    nerve: Record17ResidualNerve
    decisive_support_indices: tuple[int, ...]
    stages: tuple[Record17AllLabelStage, ...]

    @property
    def compact_coxeter_no_go(self) -> bool:
        return bool(self.stages) and self.stages[-1].passing_state_count == 0

    @property
    def verified(self) -> bool:
        index = self.nerve.orbit_index
        return (
            self.decisive_support_indices == DECISIVE_SUPPORT_INDICES[index]
            and tuple(stage.count_signature for stage in self.stages)
            == EXPECTED_STAGE_COUNTS[index]
            and tuple(stage.passing_state_sha256 for stage in self.stages)
            == EXPECTED_STAGE_DIGESTS[index]
            and self.compact_coxeter_no_go
        )

    def to_record(self) -> dict[str, object]:
        by_size = dict(self.nerve.minimal_nonfaces_by_size)
        return {
            "orbit_index": self.nerve.orbit_index,
            "orbit_signature": {
                "missing_four_count": self.nerve.missing_four_count,
                "x_graph6": self.nerve.x_graph6,
                "p_equals_q_local_vertices": list(self.nerve.p_local_vertices),
                "labelled_representation_count": self.nerve.representation_count,
            },
            "rank_nine_vertex_count": len(self.nerve.rank_nine_vertices),
            "minimal_nonfaces_by_size": {
                str(size): [list(face) for face in faces]
                for size, faces in self.nerve.minimal_nonfaces_by_size
            },
            "minimal_nonface_histogram": {
                str(size): len(faces)
                for size, faces in self.nerve.minimal_nonfaces_by_size
            },
            "dotted_pair_count": len(by_size[2]),
            "decisive_support_indices": list(self.decisive_support_indices),
            "decisive_supports": [
                list(self.nerve.higher_minimal_nonfaces[index])
                for index in self.decisive_support_indices
            ],
            "incremental_lanner_and_actual_vertex_gate": [
                stage.to_record() for stage in self.stages
            ],
            "compact_coxeter_no_go": self.compact_coxeter_no_go,
            "verified": self.verified,
        }


def audit_record17_orbit_labels(
    nerve: Record17ResidualNerve,
) -> Record17OrbitLabelAudit:
    support_indices = DECISIVE_SUPPORT_INDICES[nerve.orbit_index]
    supports = nerve.higher_minimal_nonfaces
    states: tuple[LabelState, ...] = ((),)
    cluster_union: tuple[int, ...] = ()
    stages: list[Record17AllLabelStage] = []
    for support_index in support_indices:
        states, stage = _extend_stage(
            nerve,
            states,
            cluster_union,
            support_index,
            supports[support_index],
        )
        stages.append(stage)
        cluster_union = stage.cluster_union
        if not states:
            break
    audit = Record17OrbitLabelAudit(
        nerve=nerve,
        decisive_support_indices=support_indices,
        stages=tuple(stages),
    )
    if not audit.verified:
        raise AssertionError(
            f"record-17 orbit {nerve.orbit_index} label gate changed"
        )
    return audit


@dataclass(frozen=True, slots=True)
class Record17AllOrbitLabelAudit:
    input_file_sha256: str
    input_certificate_sha256: str
    orbits: tuple[Record17OrbitLabelAudit, ...]

    @property
    def whole_shadow_compact_coxeter_no_go(self) -> bool:
        return len(self.orbits) == 8 and all(
            orbit.compact_coxeter_no_go for orbit in self.orbits
        )

    @property
    def verified(self) -> bool:
        return (
            self.input_file_sha256 == EXPECTED_INPUT_FILE_SHA256
            and self.input_certificate_sha256
            == EXPECTED_INPUT_CERTIFICATE_SHA256
            and tuple(orbit.nerve.signature for orbit in self.orbits)
            == EXPECTED_ORBIT_SIGNATURES
            and all(orbit.verified for orbit in self.orbits)
            and self.whole_shadow_compact_coxeter_no_go
        )

    def to_record(self) -> dict[str, object]:
        return {
            "input_whole_shadow_certificate": {
                "file_sha256": self.input_file_sha256,
                "certificate_sha256": self.input_certificate_sha256,
                "exact_residual_orbit_count": len(self.orbits),
            },
            "complete_lanner_tables": {
                "order_four_placements": len(lanner_placement_codes(4)),
                "order_five_placements": len(lanner_placement_codes(5)),
            },
            "orbit_gates": [orbit.to_record() for orbit in self.orbits],
            "orbits_excluded": sum(
                orbit.compact_coxeter_no_go for orbit in self.orbits
            ),
            "whole_shadow_compact_coxeter_no_go": (
                self.whole_shadow_compact_coxeter_no_go
            ),
            "verified": self.verified,
        }


@lru_cache(maxsize=1)
def audit_record17_all_coxeter_labels(
    certificate_path: Path,
) -> Record17AllOrbitLabelAudit:
    nerves = load_record17_residual_nerves(certificate_path)
    audit = Record17AllOrbitLabelAudit(
        input_file_sha256=sha256(certificate_path.read_bytes()).hexdigest(),
        input_certificate_sha256=EXPECTED_INPUT_CERTIFICATE_SHA256,
        orbits=tuple(audit_record17_orbit_labels(nerve) for nerve in nerves),
    )
    if not audit.verified:
        raise AssertionError("record-17 all-orbit compact-Coxeter gate changed")
    return audit


__all__ = [
    "DECISIVE_SUPPORT_INDICES",
    "EXPECTED_ORBIT_SIGNATURES",
    "EXPECTED_STAGE_COUNTS",
    "EXPECTED_STAGE_DIGESTS",
    "Record17AllLabelStage",
    "Record17AllOrbitLabelAudit",
    "Record17OrbitLabelAudit",
    "Record17ResidualNerve",
    "audit_record17_all_coxeter_labels",
    "audit_record17_orbit_labels",
    "load_record17_residual_nerves",
]

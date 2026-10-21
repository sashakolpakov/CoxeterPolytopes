"""Exact minimal-transversal obstruction for the forced-A4-empty stratum.

For a simplicial complex on fourteen vertices, complements identify maximal
faces with inclusion-minimal transversals of its minimal-nonface hypergraph.
A pure eight-dimensional nerve therefore has only minimal transversals of
size five.

In the H9/14 graph shadows whose strict Dehn--Sommerville interval forces A4
to be empty, the minimal-nonface hypergraph is exactly the missing edges
together with the Dehn--Sommerville-forced family A5. This module enumerates
every possible A5 family and archives one explicit wrong-sized minimal
transversal for every choice.

The verifier does not use Z3 or trust the generator's search. It reconstructs
the complete A5 choice space and checks directly that every archived
transversal hits every hyperedge, is inclusion-minimal, and has size other
than five.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.local_link import analyze_local_link_graph
from coxeter_polytopes.combinatorics.missing_edges import (
    TARGET_VERTEX_COUNT,
    analyze_missing_edge_graph,
    decode_graph6,
    macaulay_raise,
)


FULL_MASK = (1 << TARGET_VERTEX_COUNT) - 1
TARGET_TRANSVERSAL_SIZE = TARGET_VERTEX_COUNT - 9

EXPECTED_FORCED_GRAPH_COUNTS: dict[int, int] = {
    4: 1,
    6: 2,
    7: 5,
    8: 11,
    9: 17,
}
EXPECTED_FORCED_GRAPH_TOTAL = 36
EXPECTED_MISSING_FIVE_COUNT_HISTOGRAM: dict[int, int] = {0: 1, 1: 35}
EXPECTED_CHOICE_TOTAL = 30_934
EXPECTED_BAD_TRANSVERSAL_SIZE_HISTOGRAM: dict[int, int] = {
    3: 3_192,
    4: 27_742,
}
EXPECTED_LOCAL_LINK_INTERSECTION = 1
EXPECTED_COMBINED_SURVIVORS = 29

# Filled after generation and then frozen.
EXPECTED_GRAPH6_TRANSCRIPT_SHA256 = (
    "75592ebd001da892ecc226da16aa259569ee47e9bb876e5837b3601c97f32e93"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "4abf18addd6301b15edf7b9ef99a257aa61c94f160c524127ab23f0c828af587"
)


def _mask(vertices: Iterable[int]) -> int:
    result = 0
    for vertex in vertices:
        result |= 1 << vertex
    return result


def _edge_masks(graph6: str) -> tuple[int, ...]:
    vertex_count, edges = decode_graph6(graph6)
    if vertex_count != TARGET_VERTEX_COUNT:
        raise ValueError("forced-a=0 filter is specific to 14 vertices")
    return tuple(sorted((1 << first) | (1 << second) for first, second in edges))


def independent_masks(graph6: str, size: int) -> tuple[int, ...]:
    """Return all independent masks of the requested size."""

    if not 0 <= size <= TARGET_VERTEX_COUNT:
        raise ValueError("invalid independent-set size")
    edges = _edge_masks(graph6)
    result = []
    for vertices in combinations(range(TARGET_VERTEX_COUNT), size):
        mask = _mask(vertices)
        if not any(mask & edge == edge for edge in edges):
            result.append(mask)
    return tuple(result)


def strict_missing_four_interval(graph6: str) -> tuple[int, int] | None:
    """Recompute the strict global A4 interval directly from H."""

    audit = analyze_missing_edge_graph(graph6)
    if not audit.survives_graph_prefilter:
        return None
    lower = max(
        0,
        audit.g4_if_no_missing_fours - macaulay_raise(audit.g3, 3),
    )
    upper = min(
        audit.independent_four_set_count,
        audit.g4_if_no_missing_fours - 1,
    )
    return (lower, upper) if lower <= upper else None


def forced_missing_five_count(graph6: str) -> tuple[int, int, int]:
    """Return (b, I5, f4) under the forced choice A4=empty."""

    edge_count = len(_edge_masks(graph6))
    independent_three_count = len(independent_masks(graph6, 3))
    independent_four_count = len(independent_masks(graph6, 4))
    independent_five_count = len(independent_masks(graph6, 5))
    f1 = comb(TARGET_VERTEX_COUNT, 2) - edge_count
    f2 = independent_three_count
    f3 = independent_four_count
    f4 = 6 * f3 - 21 * f2 + 56 * f1 - 126 * TARGET_VERTEX_COUNT + 252
    return independent_five_count - f4, independent_five_count, f4


def is_transversal(mask: int, hyperedges: Sequence[int]) -> bool:
    """Whether mask meets every nonempty hyperedge."""

    return (
        0 <= mask <= FULL_MASK
        and all(edge != 0 and mask & edge for edge in hyperedges)
    )


def is_minimal_transversal(mask: int, hyperedges: Sequence[int]) -> bool:
    """Check inclusion-minimality directly, one deleted vertex at a time."""

    if not is_transversal(mask, hyperedges):
        return False
    remaining = mask
    while remaining:
        bit = remaining & -remaining
        if is_transversal(mask ^ bit, hyperedges):
            return False
        remaining ^= bit
    return True


def minimal_vertex_covers(graph6: str) -> tuple[int, ...]:
    """Enumerate every inclusion-minimal cover of the sparse graph H."""

    edges = _edge_masks(graph6)
    return tuple(
        mask
        for mask in range(1 << TARGET_VERTEX_COUNT)
        if is_minimal_transversal(mask, edges)
    )


def first_wrong_sized_transversal(
    graph6: str,
    missing_five_mask: int,
    minimal_covers: Sequence[int] | None = None,
) -> int | None:
    """Find a deterministic wrong-sized transversal for one A5 choice.

    A minimal transversal of E(H) union {A} either is a minimal vertex cover
    of H that already meets A, or is obtained by adding one member of A to a
    minimal vertex cover disjoint from A. Testing those candidates is
    complete, although the public verifier only validates the witness.
    """

    if missing_five_mask:
        if missing_five_mask.bit_count() != 5:
            raise ValueError("A5 member must be a five-set")
        hyperedges = _edge_masks(graph6) + (missing_five_mask,)
    else:
        hyperedges = _edge_masks(graph6)
    covers = (
        tuple(minimal_covers)
        if minimal_covers is not None
        else minimal_vertex_covers(graph6)
    )
    for cover in covers:
        if (not missing_five_mask or cover & missing_five_mask) and (
            cover.bit_count() != TARGET_TRANSVERSAL_SIZE
        ):
            return cover
    if missing_five_mask:
        for cover in covers:
            if cover & missing_five_mask:
                continue
            remaining = missing_five_mask
            while remaining:
                bit = remaining & -remaining
                candidate = cover | bit
                if (
                    candidate.bit_count() != TARGET_TRANSVERSAL_SIZE
                    and is_minimal_transversal(candidate, hyperedges)
                ):
                    return candidate
                remaining ^= bit
    return None


@dataclass(frozen=True)
class ForcedA0GraphCertificate:
    """Complete choice-by-choice certificate for one forced-A4-empty graph."""

    graph6: str
    edge_count: int
    independent_five_count: int
    forced_missing_five_count: int
    forced_f4: int
    choice_certificates: tuple[tuple[int, int], ...]

    def to_record(self) -> dict[str, object]:
        return {
            "graph6": self.graph6,
            "edges": self.edge_count,
            "strict_missing_four_interval": [0, 0],
            "independent_five_count": self.independent_five_count,
            "forced_missing_five_count": self.forced_missing_five_count,
            "forced_f4": self.forced_f4,
            "choice_count": len(self.choice_certificates),
            "choice_certificates": [
                [missing_five, transversal]
                for missing_five, transversal in self.choice_certificates
            ],
        }


def build_graph_certificate(graph6: str) -> ForcedA0GraphCertificate:
    """Generate every explicit transversal witness for one forced graph."""

    if strict_missing_four_interval(graph6) != (0, 0):
        raise ValueError("graph does not force A4 to be empty")
    b, independent_five_count, f4 = forced_missing_five_count(graph6)
    if b not in (0, 1):
        raise ValueError("forced-a=0 certificate supports only b in {0,1}")
    choices = independent_masks(graph6, 5) if b == 1 else (0,)
    covers = minimal_vertex_covers(graph6)
    certificates: list[tuple[int, int]] = []
    for missing_five in choices:
        bad = first_wrong_sized_transversal(graph6, missing_five, covers)
        if bad is None:
            raise ValueError("choice has no wrong-sized minimal transversal")
        certificates.append((missing_five, bad))
    return ForcedA0GraphCertificate(
        graph6=graph6,
        edge_count=len(_edge_masks(graph6)),
        independent_five_count=independent_five_count,
        forced_missing_five_count=b,
        forced_f4=f4,
        choice_certificates=tuple(certificates),
    )


def _canonical_digest(payload: Mapping[str, object]) -> str:
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return sha256(canonical.encode("utf-8")).hexdigest()


def build_forced_a0_manifest(
    second_layer_graph6_values: Sequence[str],
    source: Mapping[str, object],
) -> dict[str, object]:
    """Build the exhaustive forced-A4-empty transversal archive."""

    forced_graphs = tuple(
        graph6
        for graph6 in second_layer_graph6_values
        if strict_missing_four_interval(graph6) == (0, 0)
    )
    certificates = tuple(build_graph_certificate(graph6) for graph6 in forced_graphs)
    edge_counts = Counter(certificate.edge_count for certificate in certificates)
    missing_five_counts = Counter(
        certificate.forced_missing_five_count for certificate in certificates
    )
    transversal_sizes = Counter(
        transversal.bit_count()
        for certificate in certificates
        for _, transversal in certificate.choice_certificates
    )
    forced_set = set(forced_graphs)
    local_link_survivors = {
        graph6
        for graph6 in second_layer_graph6_values
        if analyze_local_link_graph(graph6).passes_edge_and_triple_feasibility
    }
    summary = {
        "second_layer_input_total": len(second_layer_graph6_values),
        "forced_a0_graphs_by_edges": {
            str(edge_count): edge_counts[edge_count]
            for edge_count in sorted(edge_counts)
        },
        "forced_a0_graphs_total": len(certificates),
        "forced_missing_five_count_histogram": {
            str(count): missing_five_counts[count]
            for count in sorted(missing_five_counts)
        },
        "a5_choice_total": sum(
            len(certificate.choice_certificates) for certificate in certificates
        ),
        "bad_transversal_size_histogram": {
            str(size): transversal_sizes[size] for size in sorted(transversal_sizes)
        },
        "local_link_survivors_total": len(local_link_survivors),
        "forced_a0_local_link_intersection": len(forced_set & local_link_survivors),
        "combined_local_link_survivors": len(local_link_survivors - forced_set),
    }
    payload: dict[str, object] = {
        "schema_version": "1.0",
        "kind": "h9_14_forced_a0_minimal_transversals",
        "source": dict(source),
        "summary": summary,
        "records": [certificate.to_record() for certificate in certificates],
    }
    payload["certificate_sha256"] = _canonical_digest(payload)
    return payload


def _verify_choice_certificate(
    graph6: str,
    forced_b: int,
    pair: object,
) -> tuple[int, int]:
    if (
        not isinstance(pair, list)
        or len(pair) != 2
        or not all(isinstance(value, int) for value in pair)
    ):
        raise ValueError("invalid forced-a=0 choice certificate")
    missing_five, transversal = pair
    if missing_five < 0 or missing_five > FULL_MASK:
        raise ValueError("missing-five mask outside the ground set")
    if forced_b == 0:
        if missing_five != 0:
            raise ValueError("b=0 choice must use the zero sentinel")
    elif missing_five.bit_count() != 5:
        raise ValueError("b=1 choice must contain one five-set")
    hyperedges = _edge_masks(graph6)
    if missing_five:
        hyperedges += (missing_five,)
    if transversal.bit_count() == TARGET_TRANSVERSAL_SIZE:
        raise ValueError("archived transversal has target size five")
    if not is_transversal(transversal, hyperedges):
        raise ValueError("archived mask does not hit every minimal nonface")
    if not is_minimal_transversal(transversal, hyperedges):
        raise ValueError("archived transversal is not inclusion-minimal")
    return missing_five, transversal


def verify_forced_a0_manifest(
    manifest: Mapping[str, object],
    second_layer_graph6_values: Sequence[str] | None = None,
) -> dict[str, object]:
    """Verify coverage and every wrong-sized transversal directly."""

    if manifest.get("schema_version") != "1.0" or manifest.get("kind") != (
        "h9_14_forced_a0_minimal_transversals"
    ):
        raise ValueError("unexpected forced-a=0 manifest metadata")
    records = manifest.get("records")
    if not isinstance(records, list):
        raise ValueError("forced-a=0 records must be a list")

    expected_forced_graphs: tuple[str, ...] | None = None
    if second_layer_graph6_values is not None:
        expected_forced_graphs = tuple(
            graph6
            for graph6 in second_layer_graph6_values
            if strict_missing_four_interval(graph6) == (0, 0)
        )

    graph6_values: list[str] = []
    edge_counts: Counter[int] = Counter()
    missing_five_counts: Counter[int] = Counter()
    transversal_sizes: Counter[int] = Counter()
    choice_total = 0
    for record in records:
        if not isinstance(record, dict) or not isinstance(record.get("graph6"), str):
            raise ValueError("invalid forced-a=0 graph record")
        graph6 = record["graph6"]
        graph6_values.append(graph6)
        if strict_missing_four_interval(graph6) != (0, 0):
            raise ValueError("recorded graph does not force A4=empty")
        b, independent_five_count, f4 = forced_missing_five_count(graph6)
        if b not in (0, 1):
            raise ValueError("forced missing-five count is outside {0,1}")
        edge_count = len(_edge_masks(graph6))
        if (
            record.get("edges") != edge_count
            or record.get("strict_missing_four_interval") != [0, 0]
            or record.get("independent_five_count") != independent_five_count
            or record.get("forced_missing_five_count") != b
            or record.get("forced_f4") != f4
        ):
            raise ValueError("forced-a=0 graph invariants do not verify")

        raw_choices = record.get("choice_certificates")
        if not isinstance(raw_choices, list):
            raise ValueError("choice certificates must be a list")
        expected_choices = independent_masks(graph6, 5) if b == 1 else (0,)
        observed_choices: list[int] = []
        for pair in raw_choices:
            missing_five, transversal = _verify_choice_certificate(graph6, b, pair)
            observed_choices.append(missing_five)
            transversal_sizes[transversal.bit_count()] += 1
        if tuple(observed_choices) != expected_choices:
            raise ValueError("archived A5 choice space is not exhaustive")
        if record.get("choice_count") != len(expected_choices):
            raise ValueError("forced-a=0 choice count does not verify")
        edge_counts[edge_count] += 1
        missing_five_counts[b] += 1
        choice_total += len(expected_choices)

    if expected_forced_graphs is not None and tuple(graph6_values) != (
        expected_forced_graphs
    ):
        raise ValueError("forced-a=0 graph stratum is not exhaustive")
    graph_transcript = "\n".join(graph6_values) + "\n"
    transcript_digest = sha256(graph_transcript.encode("utf-8")).hexdigest()
    if EXPECTED_GRAPH6_TRANSCRIPT_SHA256 and transcript_digest != (
        EXPECTED_GRAPH6_TRANSCRIPT_SHA256
    ):
        raise ValueError("unexpected forced-a=0 graph transcript")

    if second_layer_graph6_values is None:
        local_intersection = EXPECTED_LOCAL_LINK_INTERSECTION
        combined_survivors = EXPECTED_COMBINED_SURVIVORS
        local_total = 30
        second_layer_total = 142
    else:
        local_link_survivors = {
            graph6
            for graph6 in second_layer_graph6_values
            if analyze_local_link_graph(graph6).passes_edge_and_triple_feasibility
        }
        forced_set = set(graph6_values)
        local_intersection = len(forced_set & local_link_survivors)
        combined_survivors = len(local_link_survivors - forced_set)
        local_total = len(local_link_survivors)
        second_layer_total = len(second_layer_graph6_values)

    computed_summary = {
        "second_layer_input_total": second_layer_total,
        "forced_a0_graphs_by_edges": {
            str(edge_count): edge_counts[edge_count] for edge_count in sorted(edge_counts)
        },
        "forced_a0_graphs_total": len(records),
        "forced_missing_five_count_histogram": {
            str(count): missing_five_counts[count]
            for count in sorted(missing_five_counts)
        },
        "a5_choice_total": choice_total,
        "bad_transversal_size_histogram": {
            str(size): transversal_sizes[size] for size in sorted(transversal_sizes)
        },
        "local_link_survivors_total": local_total,
        "forced_a0_local_link_intersection": local_intersection,
        "combined_local_link_survivors": combined_survivors,
    }
    if computed_summary != manifest.get("summary"):
        raise ValueError("forced-a=0 summary does not verify")

    expected_edge_counts = {
        str(key): value for key, value in EXPECTED_FORCED_GRAPH_COUNTS.items()
    }
    expected_b_histogram = {
        str(key): value
        for key, value in EXPECTED_MISSING_FIVE_COUNT_HISTOGRAM.items()
    }
    expected_size_histogram = {
        str(key): value
        for key, value in EXPECTED_BAD_TRANSVERSAL_SIZE_HISTOGRAM.items()
    }
    if computed_summary["forced_a0_graphs_by_edges"] != expected_edge_counts:
        raise ValueError("unexpected forced-a=0 graph count vector")
    if computed_summary["forced_a0_graphs_total"] != EXPECTED_FORCED_GRAPH_TOTAL:
        raise ValueError("unexpected forced-a=0 graph total")
    if (
        computed_summary["forced_missing_five_count_histogram"]
        != expected_b_histogram
    ):
        raise ValueError("unexpected forced missing-five histogram")
    if computed_summary["a5_choice_total"] != EXPECTED_CHOICE_TOTAL:
        raise ValueError("unexpected A5 choice total")
    if computed_summary["bad_transversal_size_histogram"] != expected_size_histogram:
        raise ValueError("unexpected bad-transversal size histogram")
    if local_intersection != EXPECTED_LOCAL_LINK_INTERSECTION:
        raise ValueError("unexpected local-link intersection")
    if combined_survivors != EXPECTED_COMBINED_SURVIVORS:
        raise ValueError("unexpected combined survivor total")

    payload = dict(manifest)
    archived_digest = payload.pop("certificate_sha256", None)
    computed_digest = _canonical_digest(payload)
    if archived_digest != computed_digest:
        raise ValueError("forced-a=0 certificate digest does not verify")
    if EXPECTED_CERTIFICATE_SHA256 and archived_digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("unexpected forced-a=0 certificate digest")
    return computed_summary

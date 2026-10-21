"""Exact missing-edge prefilter for the 9-dimensional, 14-facet case.

Let ``K`` be the nerve of a hypothetical compact Coxeter 9-polytope with
14 facets, equivalently the boundary of its simplicial polar.  This module
audits the graph ``H`` whose edges are the missing two-faces of ``K``.

The published reduction that ``K`` has no missing triples makes its first
three face numbers functions of ``H`` alone.  The g-theorem then gives a
small, exact graph-level prefilter.  A surviving graph is *not* a nerve: the
missing faces of orders four and five, the sphere condition, polytopality,
and all Coxeter/Gram constraints remain to be supplied.

The implementation deliberately has no graph-library dependency.  In
particular, the archived nauty output can be checked using pure Python.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations
from math import comb
from typing import Iterable, Mapping, Sequence


TARGET_VERTEX_COUNT = 14
TARGET_POLYTOPE_DIMENSION = 9
MIN_MISSING_EDGES = 2
MAX_MISSING_EDGES = 10
MAX_COXETER_MISSING_FACE_ORDER = 5

# Independent nauty/geng counts for unlabeled simple graphs on 14 vertices
# in the sparse range forced by the argument below.  They are frozen in the
# verifier so that editing both a manifest's records and its summary cannot
# conceal a dropped graph.
EXPECTED_UNLABELED_COUNTS: dict[int, int] = {
    2: 2,
    3: 5,
    4: 11,
    5: 26,
    6: 68,
    7: 177,
    8: 495,
    9: 1464,
    10: 4541,
}

# Digest of the exact ordered graph6 transcript produced by the archived
# nauty/geng run.  Unlike the generator executable digest, this is portable:
# the pure-Python verifier freezes it so that coordinated edits to records and
# the manifest's two self-reported transcript hashes are still detected.
EXPECTED_GRAPH6_STDOUT_SHA256 = (
    "780c577cc7344d270d4551f1ba6db7834e0eef97e96a7e35da1c8abe8728d64f"
)


def macaulay_raise(value: int, degree: int) -> int:
    """Return Macaulay's ``value^(<degree>)`` exactly.

    If

    ``value = C(a_degree, degree) + ... + C(a_j, j)``

    is its unique Macaulay expansion, the result is the sum of
    ``C(a_i + 1, i + 1)``.  Zero summands are harmless and make the greedy
    implementation uniform at ``value == 0``.
    """

    if value < 0:
        raise ValueError("Macaulay raising is defined here only for value >= 0")
    if degree < 1:
        raise ValueError("Macaulay degree must be positive")

    remaining = value
    next_top: int | None = None
    raised = 0
    for lower in range(degree, 0, -1):
        top = lower - 1
        while (
            (next_top is None or top + 1 < next_top)
            and comb(top + 1, lower) <= remaining
        ):
            top += 1
        remaining -= comb(top, lower)
        raised += comb(top + 1, lower + 1)
        next_top = top
    if remaining != 0:  # defensive: the greedy Macaulay expansion is exact
        raise AssertionError("incomplete Macaulay expansion")
    return raised


def decode_graph6(value: str) -> tuple[int, tuple[tuple[int, int], ...]]:
    """Decode a small graph6 string into its vertex count and edge tuple.

    The target has only 14 vertices, so accepting the one-byte graph6 order
    field (orders 0 through 62) keeps this verifier small and auditable.
    """

    if value.startswith(">>graph6<<"):
        value = value[len(">>graph6<<") :]
    if not value:
        raise ValueError("empty graph6 string")
    try:
        codes = [ord(character) - 63 for character in value]
    except TypeError as exc:  # pragma: no cover - protects a public API
        raise ValueError("graph6 input must be a string") from exc
    if any(code < 0 or code > 63 for code in codes):
        raise ValueError("graph6 characters must lie in ASCII range 63..126")
    if codes[0] == 63:
        raise ValueError("extended graph6 orders are outside this verifier's scope")

    vertex_count = codes[0]
    edge_bit_count = comb(vertex_count, 2)
    payload_length = (edge_bit_count + 5) // 6
    if len(codes) != payload_length + 1:
        raise ValueError(
            f"graph6 payload has length {len(codes) - 1}, expected {payload_length}"
        )

    bits: list[int] = []
    for code in codes[1:]:
        bits.extend((code >> shift) & 1 for shift in range(5, -1, -1))
    if any(bits[edge_bit_count:]):
        raise ValueError("nonzero graph6 padding bits")

    edges: list[tuple[int, int]] = []
    bit_index = 0
    for upper in range(1, vertex_count):
        for lower in range(upper):
            if bits[bit_index]:
                edges.append((lower, upper))
            bit_index += 1
    return vertex_count, tuple(edges)


def encode_graph6(
    vertex_count: int, edges: Iterable[tuple[int, int]]
) -> str:
    """Encode a simple graph of order at most 62 in graph6 format."""

    if not 0 <= vertex_count <= 62:
        raise ValueError("only graph6 orders 0..62 are supported")
    normalized: set[tuple[int, int]] = set()
    for first, second in edges:
        if first == second:
            raise ValueError("graph6 encodes simple graphs, not loops")
        lower, upper = sorted((first, second))
        if lower < 0 or upper >= vertex_count:
            raise ValueError("edge endpoint outside graph order")
        normalized.add((lower, upper))

    bits = [
        int((lower, upper) in normalized)
        for upper in range(1, vertex_count)
        for lower in range(upper)
    ]
    bits.extend([0] * ((-len(bits)) % 6))
    payload = []
    for start in range(0, len(bits), 6):
        code = 0
        for bit in bits[start : start + 6]:
            code = (code << 1) | bit
        payload.append(chr(code + 63))
    return chr(vertex_count + 63) + "".join(payload)


@dataclass(frozen=True)
class MissingEdgeAudit:
    """Exact low-skeleton invariants of one missing-edge graph."""

    graph6: str
    vertex_count: int
    edge_count: int
    wedge_count: int
    triangle_count: int
    f1: int
    f2: int
    h1: int
    h2: int
    h3: int
    g2: int
    g3: int
    macaulay_g3_bound: int
    passes_g3_macaulay: bool
    independent_four_set_count: int
    g4_if_no_missing_fours: int
    macaulay_g4_bound: int
    missing_four_count_min: int
    missing_four_count_max: int
    admits_missing_four_count: bool
    passes_g_vector_prefilter: bool
    stacked_obstruction: bool
    survives_graph_prefilter: bool

    def to_record(self) -> dict[str, object]:
        """Return the canonical JSON record used in the archive."""

        if not self.passes_g3_macaulay:
            disposition = "reject_g3_macaulay"
        elif not self.admits_missing_four_count:
            disposition = "reject_g4_macaulay_feasibility"
        elif self.stacked_obstruction:
            disposition = "reject_g2_zero_stacked_missing_9_face"
        else:
            disposition = "survives_graph_prefilter"
        return {
            "graph6": self.graph6,
            "edges": self.edge_count,
            "wedges": self.wedge_count,
            "triangles": self.triangle_count,
            "f1": self.f1,
            "f2": self.f2,
            "h1": self.h1,
            "h2": self.h2,
            "h3": self.h3,
            "g2": self.g2,
            "g3": self.g3,
            "macaulay_g3_bound": self.macaulay_g3_bound,
            "passes_g3_macaulay": self.passes_g3_macaulay,
            "independent_four_set_count": self.independent_four_set_count,
            "g4_if_no_missing_fours": self.g4_if_no_missing_fours,
            "macaulay_g4_bound": self.macaulay_g4_bound,
            "missing_four_count_min": self.missing_four_count_min,
            "missing_four_count_max": self.missing_four_count_max,
            "admits_missing_four_count": self.admits_missing_four_count,
            "passes_g_vector_prefilter": self.passes_g_vector_prefilter,
            "stacked_obstruction": self.stacked_obstruction,
            "survives_graph_prefilter": self.survives_graph_prefilter,
            "disposition": disposition,
        }


def _graph_statistics(
    vertex_count: int, edges: Sequence[tuple[int, int]]
) -> tuple[int, int]:
    adjacency = [set() for _ in range(vertex_count)]
    for first, second in edges:
        adjacency[first].add(second)
        adjacency[second].add(first)
    wedge_count = sum(comb(len(neighbors), 2) for neighbors in adjacency)
    triangle_count = sum(
        1
        for first in range(vertex_count)
        for second in adjacency[first]
        if first < second
        for third in adjacency[first] & adjacency[second]
        if second < third
    )
    return wedge_count, triangle_count


def _independent_four_set_count(
    vertex_count: int, edges: Sequence[tuple[int, int]]
) -> int:
    edge_set = set(edges)
    return sum(
        1
        for vertices in combinations(range(vertex_count), 4)
        if all(
            (first, second) not in edge_set
            for first, second in combinations(vertices, 2)
        )
    )


def analyze_missing_edge_graph(graph6: str) -> MissingEdgeAudit:
    """Audit one graph under the no-missing-triples target assumptions."""

    vertex_count, edges = decode_graph6(graph6)
    if vertex_count != TARGET_VERTEX_COUNT:
        raise ValueError(
            f"target missing-edge graph must have {TARGET_VERTEX_COUNT} vertices"
        )

    edge_count = len(edges)
    wedge_count, triangle_count = _graph_statistics(vertex_count, edges)
    dimension = TARGET_POLYTOPE_DIMENSION

    # Inclusion-exclusion over triples containing an edge of H.  With no
    # missing triples, precisely the triples containing no H-edge are faces.
    bad_triples = (
        (vertex_count - 2) * edge_count - wedge_count + triangle_count
    )
    f1 = comb(vertex_count, 2) - edge_count
    f2 = comb(vertex_count, 3) - bad_triples

    h1 = vertex_count - dimension
    h2 = f1 - (dimension - 1) * vertex_count + comb(dimension, 2)
    h3 = (
        f2
        - (dimension - 2) * f1
        + comb(dimension - 1, 2) * vertex_count
        - comb(dimension, 3)
    )
    g2 = h2 - h1
    g3 = h3 - h2
    macaulay_g3_bound = macaulay_raise(g2, 2) if g2 >= 0 else -1
    passes_g3_macaulay = g2 >= 0 and 0 <= g3 <= macaulay_g3_bound

    # Let a be the number of missing faces of order four.  An independent
    # four-set in H has all its edges, and the no-missing-triples assumption
    # gives all of its triangles.  Therefore it is either a 3-face or exactly
    # a missing four-face, and f3 = I4(H)-a.  Consequently g4 is affine in a.
    independent_four_set_count = _independent_four_set_count(vertex_count, edges)
    h4_if_no_missing_fours = (
        independent_four_set_count
        - (dimension - 3) * f2
        + comb(dimension - 2, 2) * f1
        - comb(dimension - 1, 3) * vertex_count
        + comb(dimension, 4)
    )
    g4_if_no_missing_fours = h4_if_no_missing_fours - h3
    if passes_g3_macaulay:
        macaulay_g4_bound = macaulay_raise(g3, 3)
        missing_four_count_min = max(
            0, g4_if_no_missing_fours - macaulay_g4_bound
        )
        missing_four_count_max = min(
            independent_four_set_count, g4_if_no_missing_fours
        )
        admits_missing_four_count = (
            missing_four_count_min <= missing_four_count_max
        )
    else:
        macaulay_g4_bound = -1
        missing_four_count_min = -1
        missing_four_count_max = -1
        admits_missing_four_count = False
    passes_g_vector_prefilter = (
        passes_g3_macaulay and admits_missing_four_count
    )

    # If e=10 then g2=0.  The equality case of the lower-bound theorem
    # makes the polar a stacked 9-polytope.  Since it has 14>10 vertices,
    # at least one stacking occurred and its removed facet remains a missing
    # face of order 9.  Coxeter missing faces are Lanner diagrams of order at
    # most 5, so this entire edge stratum is impossible.
    stacked_obstruction = edge_count == MAX_MISSING_EDGES
    survives_graph_prefilter = passes_g_vector_prefilter and not stacked_obstruction

    return MissingEdgeAudit(
        graph6=graph6,
        vertex_count=vertex_count,
        edge_count=edge_count,
        wedge_count=wedge_count,
        triangle_count=triangle_count,
        f1=f1,
        f2=f2,
        h1=h1,
        h2=h2,
        h3=h3,
        g2=g2,
        g3=g3,
        macaulay_g3_bound=macaulay_g3_bound,
        passes_g3_macaulay=passes_g3_macaulay,
        independent_four_set_count=independent_four_set_count,
        g4_if_no_missing_fours=g4_if_no_missing_fours,
        macaulay_g4_bound=macaulay_g4_bound,
        missing_four_count_min=missing_four_count_min,
        missing_four_count_max=missing_four_count_max,
        admits_missing_four_count=admits_missing_four_count,
        passes_g_vector_prefilter=passes_g_vector_prefilter,
        stacked_obstruction=stacked_obstruction,
        survives_graph_prefilter=survives_graph_prefilter,
    )


def _edge_count_map(counter: Counter[int]) -> dict[str, int]:
    return {
        str(edge_count): counter.get(edge_count, 0)
        for edge_count in range(MIN_MISSING_EDGES, MAX_MISSING_EDGES + 1)
    }


def summarize_audits(audits: Sequence[MissingEdgeAudit]) -> dict[str, object]:
    generated = Counter(audit.edge_count for audit in audits)
    g3_macaulay = Counter(
        audit.edge_count for audit in audits if audit.passes_g3_macaulay
    )
    g4_feasible = Counter(
        audit.edge_count for audit in audits if audit.passes_g_vector_prefilter
    )
    stacked = Counter(
        audit.edge_count
        for audit in audits
        if audit.passes_g_vector_prefilter and audit.stacked_obstruction
    )
    survivors = Counter(
        audit.edge_count
        for audit in audits
        if audit.survives_graph_prefilter
    )
    return {
        "generated_total": len(audits),
        "generated_by_edges": _edge_count_map(generated),
        "g3_macaulay_survivors_total": sum(g3_macaulay.values()),
        "g3_macaulay_survivors_by_edges": _edge_count_map(g3_macaulay),
        "g4_feasible_survivors_total": sum(g4_feasible.values()),
        "g4_feasible_survivors_by_edges": _edge_count_map(g4_feasible),
        "stacked_rejections_after_g4_total": sum(stacked.values()),
        "stacked_rejections_after_g4_by_edges": _edge_count_map(stacked),
        "graph_prefilter_survivors_total": sum(survivors.values()),
        "graph_prefilter_survivors_by_edges": _edge_count_map(survivors),
    }


def _stdout_bytes(graph6_values: Sequence[str]) -> bytes:
    return "".join(f"{value}\n" for value in graph6_values).encode("ascii")


def build_manifest(
    graph6_values: Iterable[str], generator: Mapping[str, object]
) -> dict[str, object]:
    """Build the deterministic graph-level coverage manifest."""

    values = [value.strip() for value in graph6_values if value.strip()]
    audits = [analyze_missing_edge_graph(value) for value in values]
    manifest = {
        "schema_version": "1.1.0",
        "object": "missing-edge graphs, not complete nerves",
        "target": {
            "simple_polytope_dimension": TARGET_POLYTOPE_DIMENSION,
            "facet_count": TARGET_VERTEX_COUNT,
            "polar": "simplicial 9-polytope with 14 vertices",
            "missing_edge_graph_order": TARGET_VERTEX_COUNT,
            "assumed_minimum_missing_edges": MIN_MISSING_EDGES,
            "assumed_no_missing_faces_of_order_3": True,
            "assumed_maximum_missing_face_order": (
                MAX_COXETER_MISSING_FACE_ORDER
            ),
        },
        "coverage_claim": (
            "Conditional on the stated published reductions, every target nerve's "
            "missing-edge graph occurs in this unlabeled graph list. Surviving "
            "records are graph-level candidates only."
        ),
        "nonclaims": [
            "A surviving graph need not extend to a simplicial sphere or polytope nerve.",
            "Missing faces of orders 4 and 5 have not been enumerated.",
            "No Coxeter labeling or Gram matrix is supplied by this manifest.",
        ],
        "generator": dict(generator),
        "graph6_stdout_sha256": sha256(_stdout_bytes(values)).hexdigest(),
        "counts": summarize_audits(audits),
        "records": [audit.to_record() for audit in audits],
    }
    verify_manifest(manifest)
    return manifest


def verify_manifest(manifest: Mapping[str, object]) -> dict[str, object]:
    """Independently verify every archived graph and summary in pure Python.

    This verifies the graph6 transcript, all low-skeleton arithmetic, every
    recorded disposition, and the frozen nauty count vector.  It does not
    re-prove nauty's canonical-generation algorithm.
    """

    if manifest.get("schema_version") != "1.1.0":
        raise ValueError("unsupported missing-edge manifest schema")
    if manifest.get("object") != "missing-edge graphs, not complete nerves":
        raise ValueError("manifest does not identify its graph-only scope")

    raw_records = manifest.get("records")
    if not isinstance(raw_records, list):
        raise ValueError("manifest records must be a list")

    graph6_values: list[str] = []
    audits: list[MissingEdgeAudit] = []
    for index, raw_record in enumerate(raw_records):
        if not isinstance(raw_record, Mapping):
            raise ValueError(f"record {index} is not an object")
        graph6_value = raw_record.get("graph6")
        if not isinstance(graph6_value, str):
            raise ValueError(f"record {index} has no graph6 string")
        audit = analyze_missing_edge_graph(graph6_value)
        if dict(raw_record) != audit.to_record():
            raise ValueError(f"record {index} invariants do not verify")
        if not MIN_MISSING_EDGES <= audit.edge_count <= MAX_MISSING_EDGES:
            raise ValueError(f"record {index} lies outside generated edge range")
        graph6_values.append(graph6_value)
        audits.append(audit)

    if len(set(graph6_values)) != len(graph6_values):
        raise ValueError("manifest repeats a graph6 record")
    transcript_hash = sha256(_stdout_bytes(graph6_values)).hexdigest()
    if manifest.get("graph6_stdout_sha256") != transcript_hash:
        raise ValueError("graph6 transcript hash mismatch")
    if transcript_hash != EXPECTED_GRAPH6_STDOUT_SHA256:
        raise ValueError("graph6 transcript does not match the frozen archive")

    summary = summarize_audits(audits)
    if manifest.get("counts") != summary:
        raise ValueError("manifest summary does not match its records")
    generated_counts = {
        int(edge_count): count
        for edge_count, count in summary["generated_by_edges"].items()
    }
    if generated_counts != EXPECTED_UNLABELED_COUNTS:
        raise ValueError("unlabeled graph count vector is incomplete or unexpected")

    generator = manifest.get("generator")
    if not isinstance(generator, Mapping):
        raise ValueError("manifest lacks generator provenance")
    if generator.get("command") != ["geng", "-q", "14", "2:10"]:
        raise ValueError("unexpected geng command")
    if generator.get("stdout_sha256") != transcript_hash:
        raise ValueError("generator stdout hash does not match archived transcript")
    for required in ("version", "executable", "executable_sha256"):
        if not isinstance(generator.get(required), str) or not generator[required]:
            raise ValueError(f"generator provenance lacks {required}")

    return summary

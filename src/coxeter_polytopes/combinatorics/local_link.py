"""Polytopal-link obstructions for the H9/14 nerve search.

Let ``K`` be the boundary of the putative simplicial polar, and let ``H``
be its graph of missing edges.  Every link in ``K`` is again the boundary of
a simplicial polytope.  Moreover, a minimal nonface of a link is obtained by
deleting the linked face from a minimal nonface of ``K``.  The Coxeter/Lanner
bound therefore says that link minimal nonfaces also have order at most five.

Combining this observation with the equality case of the generalized lower
bound theorem makes several entries of every local g-vector *strictly*
positive.  The vertex-figure conditions through g3 depend only on ``H`` and
on the incidence numbers of the still-unknown missing four-faces.  This file
implements those exact formulas without a SAT or graph-library dependency.

The module is a necessary filter only.  Passing it does not construct a
nerve, a sphere, an oriented matroid, or a Coxeter Gram matrix.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations
import json
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.missing_edges import (
    TARGET_VERTEX_COUNT,
    decode_graph6,
    macaulay_raise,
)


EXPECTED_INPUT_COUNTS: dict[int, int] = {
    2: 2,
    3: 5,
    4: 11,
    5: 17,
    6: 29,
    7: 33,
    8: 28,
    9: 17,
}

EXPECTED_VERTEX_G1_COUNTS: dict[int, int] = {
    2: 2,
    3: 5,
    4: 10,
    5: 15,
    6: 22,
    7: 14,
    8: 4,
    9: 0,
}

EXPECTED_VERTEX_G2_COUNTS: dict[int, int] = {
    2: 2,
    3: 5,
    4: 9,
    5: 10,
    6: 9,
    7: 3,
    8: 1,
    9: 0,
}

EXPECTED_VERTEX_G3_COUNTS: dict[int, int] = {
    2: 2,
    3: 5,
    4: 8,
    5: 9,
    6: 8,
    7: 3,
    8: 1,
    9: 0,
}

EXPECTED_LOWER_FACE_COUNTS: dict[int, int] = {
    2: 2,
    3: 5,
    4: 7,
    5: 8,
    6: 6,
    7: 2,
    8: 0,
    9: 0,
}

EXPECTED_GRAPH6_TRANSCRIPT_SHA256 = (
    "30c8ba87c5669b8ec714343d6648617ebc47820ce381cb571aa67f00c3236868"
)

EXPECTED_SURVIVOR_TRANSCRIPT_SHA256 = (
    "9ffa3e2f998cb6349e34adfa3df8c616c89aee257dbb474c415093249cd29968"
)

EXPECTED_CERTIFICATE_SHA256 = (
    "946137ef8fe8f38b7f9f7df46261a0ad617060d286d28ab7e853e7b245a75b74"
)


def _edge_set(graph6: str) -> tuple[int, set[tuple[int, int]]]:
    vertex_count, edges = decode_graph6(graph6)
    if vertex_count != TARGET_VERTEX_COUNT:
        raise ValueError("local-link filter is specific to 14 vertices")
    return vertex_count, {tuple(sorted(edge)) for edge in edges}


def _is_independent(
    vertices: Sequence[int], edges: set[tuple[int, int]]
) -> bool:
    return all(tuple(sorted(pair)) not in edges for pair in combinations(vertices, 2))


def _strict_interval(
    *,
    previous_g: int,
    next_g_if_no_missing_fours: int,
    candidate_missing_four_count: int,
    macaulay_degree: int,
) -> tuple[int, int] | None:
    """Return the allowed incidence interval after a strict local g-test.

    The actual next g-entry is ``next_g_if_no_missing_fours - a``, where
    ``a`` is the number of selected missing four-faces incident to the linked
    face. The inherited Lanner bound excludes the simplex figure itself and
    makes the next g-entry strictly positive. The g-theorem supplies the
    Macaulay upper bound.
    """

    beta = next_g_if_no_missing_fours
    upper_g = macaulay_raise(previous_g, macaulay_degree)
    lower = max(0, beta - upper_g)
    upper = beta - 1
    lower = max(0, lower)
    upper = min(candidate_missing_four_count, upper)
    return (lower, upper) if lower <= upper else None


@dataclass(frozen=True)
class VertexFigureAudit:
    """Low g-vector data for the 8-dimensional figure at one vertex."""

    vertex: int
    link_vertex_count: int
    inherited_missing_edge_count: int
    candidate_link_triangle_count: int
    g1: int
    g2: int
    g3_if_no_missing_fours: int
    passes_strict_g1: bool
    passes_strict_g2: bool
    missing_four_incidence_interval: tuple[int, int] | None

    @property
    def admits_strict_g3(self) -> bool:
        return self.missing_four_incidence_interval is not None

    @property
    def passes_graph_filter(self) -> bool:
        return self.passes_strict_g2 and self.admits_strict_g3

    def to_record(self) -> dict[str, object]:
        return {
            "vertex": self.vertex,
            "link_vertex_count": self.link_vertex_count,
            "inherited_missing_edge_count": self.inherited_missing_edge_count,
            "candidate_link_triangle_count": self.candidate_link_triangle_count,
            "g1": self.g1,
            "g2": self.g2,
            "g3_if_no_missing_fours": self.g3_if_no_missing_fours,
            "passes_strict_g1": self.passes_strict_g1,
            "passes_strict_g2": self.passes_strict_g2,
            "missing_four_incidence_interval": (
                list(self.missing_four_incidence_interval)
                if self.missing_four_incidence_interval is not None
                else None
            ),
            "admits_strict_g3": self.admits_strict_g3,
        }


@dataclass(frozen=True)
class EdgeFigureConstraint:
    """The g2 constraint in the 7-dimensional figure of one edge."""

    edge: tuple[int, int]
    link_vertex_count: int
    candidate_link_edge_count: int
    g1: int
    g2_if_no_missing_fours: int
    missing_four_incidence_interval: tuple[int, int] | None

    @property
    def feasible(self) -> bool:
        return self.missing_four_incidence_interval is not None


@dataclass(frozen=True)
class TripleFigureConstraint:
    """Vertex-count constraint in the 6-dimensional figure of a triangle."""

    triple: tuple[int, int, int]
    candidate_link_vertex_count: int
    missing_four_incidence_max: int

    @property
    def feasible(self) -> bool:
        return self.missing_four_incidence_max >= 0


@dataclass(frozen=True)
class LocalLinkGraphAudit:
    """All graph-determined local-link constraints for one shadow."""

    graph6: str
    edge_count: int
    vertex_figures: tuple[VertexFigureAudit, ...]
    edge_figures: tuple[EdgeFigureConstraint, ...]
    triple_figures: tuple[TripleFigureConstraint, ...]

    @property
    def passes_all_vertex_g1(self) -> bool:
        return all(figure.passes_strict_g1 for figure in self.vertex_figures)

    @property
    def passes_all_vertex_g2(self) -> bool:
        return self.passes_all_vertex_g1 and all(
            figure.passes_strict_g2 for figure in self.vertex_figures
        )

    @property
    def passes_all_vertex_g3_feasibility(self) -> bool:
        return self.passes_all_vertex_g2 and all(
            figure.admits_strict_g3 for figure in self.vertex_figures
        )

    @property
    def passes_edge_and_triple_feasibility(self) -> bool:
        return self.passes_all_vertex_g3_feasibility and all(
            figure.feasible for figure in self.edge_figures
        ) and all(figure.feasible for figure in self.triple_figures)

    def to_record(self) -> dict[str, object]:
        if not self.passes_all_vertex_g1:
            disposition = "reject_vertex_figure_strict_g1"
        elif not self.passes_all_vertex_g2:
            disposition = "reject_vertex_figure_strict_g2"
        elif not self.passes_all_vertex_g3_feasibility:
            disposition = "reject_vertex_figure_strict_g3"
        elif not self.passes_edge_and_triple_feasibility:
            disposition = "reject_lower_face_figure_feasibility"
        else:
            disposition = "survives_local_link_graph_filter"
        return {
            "graph6": self.graph6,
            "edges": self.edge_count,
            "vertex_figures": [figure.to_record() for figure in self.vertex_figures],
            "passes_all_vertex_g1": self.passes_all_vertex_g1,
            "passes_all_vertex_g2": self.passes_all_vertex_g2,
            "passes_all_vertex_g3_feasibility": (
                self.passes_all_vertex_g3_feasibility
            ),
            "passes_edge_and_triple_feasibility": (
                self.passes_edge_and_triple_feasibility
            ),
            "disposition": disposition,
        }


def vertex_figure_audits(graph6: str) -> tuple[VertexFigureAudit, ...]:
    """Compute exact vertex-figure g2 and missing-four/g3 intervals."""

    vertex_count, edges = _edge_set(graph6)
    result: list[VertexFigureAudit] = []
    for vertex in range(vertex_count):
        link_vertices = tuple(
            other
            for other in range(vertex_count)
            if other != vertex and tuple(sorted((vertex, other))) not in edges
        )
        n = len(link_vertices)
        inherited_missing_edges = sum(
            tuple(sorted(pair)) in edges for pair in combinations(link_vertices, 2)
        )
        f1 = len(tuple(combinations(link_vertices, 2))) - inherited_missing_edges
        candidate_f2 = sum(
            _is_independent(triple, edges)
            for triple in combinations(link_vertices, 3)
        )

        # The link is the boundary of a simplicial 8-polytope.
        g1 = n - 9
        g2 = f1 - 8 * n + 36
        beta = candidate_f2 - 7 * f1 + 28 * n - 84
        passes_g1 = g1 > 0
        passes_g2 = passes_g1 and 1 <= g2 <= macaulay_raise(g1, 1)

        interval = None
        if passes_g2:
            interval = _strict_interval(
                previous_g=g2,
                next_g_if_no_missing_fours=beta,
                candidate_missing_four_count=candidate_f2,
                macaulay_degree=2,
            )

        result.append(
            VertexFigureAudit(
                vertex=vertex,
                link_vertex_count=n,
                inherited_missing_edge_count=inherited_missing_edges,
                candidate_link_triangle_count=candidate_f2,
                g1=g1,
                g2=g2,
                g3_if_no_missing_fours=beta,
                passes_strict_g1=passes_g1,
                passes_strict_g2=passes_g2,
                missing_four_incidence_interval=interval,
            )
        )
    return tuple(result)


def edge_figure_constraints(graph6: str) -> tuple[EdgeFigureConstraint, ...]:
    """Return linear missing-four incidence constraints for every face edge."""

    vertex_count, edges = _edge_set(graph6)
    result: list[EdgeFigureConstraint] = []
    for edge in combinations(range(vertex_count), 2):
        if edge in edges:
            continue
        link_vertices = tuple(
            other
            for other in range(vertex_count)
            if other not in edge
            and all(tuple(sorted((other, endpoint))) not in edges for endpoint in edge)
        )
        n = len(link_vertices)
        candidate_f1 = sum(
            tuple(sorted(pair)) not in edges
            for pair in combinations(link_vertices, 2)
        )
        g1 = n - 8
        beta = candidate_f1 - 7 * n + 28
        interval = None
        if g1 > 0:
            interval = _strict_interval(
                previous_g=g1,
                next_g_if_no_missing_fours=beta,
                candidate_missing_four_count=candidate_f1,
                macaulay_degree=1,
            )
        result.append(
            EdgeFigureConstraint(
                edge=edge,
                link_vertex_count=n,
                candidate_link_edge_count=candidate_f1,
                g1=g1,
                g2_if_no_missing_fours=beta,
                missing_four_incidence_interval=interval,
            )
        )
    return tuple(result)


def triple_figure_constraints(graph6: str) -> tuple[TripleFigureConstraint, ...]:
    """Return the minimum-vertex constraint for every triangular face.

    An independent triple is a face by the published no-missing-triples
    reduction. Its face figure has dimension six. A seven-vertex figure would
    be a simplex whose boundary has a forbidden missing face of order seven,
    so it must have at least eight vertices. Each missing four-face containing
    the triple removes exactly one candidate link vertex.
    """

    vertex_count, edges = _edge_set(graph6)
    result: list[TripleFigureConstraint] = []
    for triple in combinations(range(vertex_count), 3):
        if not _is_independent(triple, edges):
            continue
        candidate_vertices = sum(
            all(tuple(sorted((other, endpoint))) not in edges for endpoint in triple)
            for other in range(vertex_count)
            if other not in triple
        )
        result.append(
            TripleFigureConstraint(
                triple=triple,
                candidate_link_vertex_count=candidate_vertices,
                missing_four_incidence_max=candidate_vertices - 8,
            )
        )
    return tuple(result)


def analyze_local_link_graph(graph6: str) -> LocalLinkGraphAudit:
    """Run every graph-determined local-link obstruction."""

    _, edges = _edge_set(graph6)
    return LocalLinkGraphAudit(
        graph6=graph6,
        edge_count=len(edges),
        vertex_figures=vertex_figure_audits(graph6),
        edge_figures=edge_figure_constraints(graph6),
        triple_figures=triple_figure_constraints(graph6),
    )


def _count_by_edges(
    audits: Sequence[LocalLinkGraphAudit], predicate: str | None = None
) -> dict[str, int]:
    # String keys make the in-memory object identical before and after its
    # JSON round trip; JSON object keys are strings by definition.
    result = {str(edge_count): 0 for edge_count in EXPECTED_INPUT_COUNTS}
    for audit in audits:
        if predicate is None or bool(getattr(audit, predicate)):
            result[str(audit.edge_count)] += 1
    return result


def build_local_link_manifest(
    graph6_values: Sequence[str], source: Mapping[str, object]
) -> dict[str, object]:
    """Build a compact exact archive of the vertex-figure reduction."""

    audits = tuple(analyze_local_link_graph(value) for value in graph6_values)
    records = [audit.to_record() for audit in audits]
    summary = {
        "input_by_edges": _count_by_edges(audits),
        "input_total": len(audits),
        "vertex_g1_survivors_by_edges": _count_by_edges(
            audits, "passes_all_vertex_g1"
        ),
        "vertex_g1_survivors_total": sum(
            audit.passes_all_vertex_g1 for audit in audits
        ),
        "vertex_g2_survivors_by_edges": _count_by_edges(
            audits, "passes_all_vertex_g2"
        ),
        "vertex_g2_survivors_total": sum(
            audit.passes_all_vertex_g2 for audit in audits
        ),
        "vertex_g3_survivors_by_edges": _count_by_edges(
            audits, "passes_all_vertex_g3_feasibility"
        ),
        "vertex_g3_survivors_total": sum(
            audit.passes_all_vertex_g3_feasibility for audit in audits
        ),
        "lower_face_feasibility_survivors_total": sum(
            audit.passes_edge_and_triple_feasibility for audit in audits
        ),
        "lower_face_feasibility_survivors_by_edges": _count_by_edges(
            audits, "passes_edge_and_triple_feasibility"
        ),
    }
    payload = {
        "schema_version": "1.0",
        "kind": "h9_14_local_polytopal_link_filter",
        "source": dict(source),
        "summary": summary,
        "records": records,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["certificate_sha256"] = sha256(canonical.encode("utf-8")).hexdigest()
    return payload


def verify_local_link_manifest(manifest: Mapping[str, object]) -> dict[str, object]:
    """Independently recompute every archived graph and frozen count."""

    if manifest.get("schema_version") != "1.0" or manifest.get("kind") != (
        "h9_14_local_polytopal_link_filter"
    ):
        raise ValueError("unexpected local-link manifest metadata")
    records = manifest.get("records")
    if not isinstance(records, list):
        raise ValueError("local-link records must be a list")
    graph6_values = []
    for record in records:
        if not isinstance(record, dict) or not isinstance(record.get("graph6"), str):
            raise ValueError("invalid local-link record")
        graph6_values.append(record["graph6"])
    rebuilt = build_local_link_manifest(graph6_values, manifest.get("source", {}))
    if rebuilt["records"] != records or rebuilt["summary"] != manifest.get("summary"):
        raise ValueError("local-link invariants do not verify")
    if rebuilt["certificate_sha256"] != manifest.get("certificate_sha256"):
        raise ValueError("local-link certificate digest does not verify")
    graph_transcript = "\n".join(graph6_values) + "\n"
    if sha256(graph_transcript.encode("utf-8")).hexdigest() != (
        EXPECTED_GRAPH6_TRANSCRIPT_SHA256
    ):
        raise ValueError("unexpected local-link graph transcript")
    survivor_transcript = "\n".join(
        graph6
        for graph6 in graph6_values
        if analyze_local_link_graph(graph6).passes_edge_and_triple_feasibility
    ) + "\n"
    if sha256(survivor_transcript.encode("utf-8")).hexdigest() != (
        EXPECTED_SURVIVOR_TRANSCRIPT_SHA256
    ):
        raise ValueError("unexpected local-link survivor transcript")
    if EXPECTED_CERTIFICATE_SHA256 and rebuilt["certificate_sha256"] != (
        EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("unexpected local-link certificate digest")
    summary = rebuilt["summary"]
    expected_input = {str(key): value for key, value in EXPECTED_INPUT_COUNTS.items()}
    expected_g1 = {
        str(key): value for key, value in EXPECTED_VERTEX_G1_COUNTS.items()
    }
    expected_g2 = {
        str(key): value for key, value in EXPECTED_VERTEX_G2_COUNTS.items()
    }
    expected_g3 = {
        str(key): value for key, value in EXPECTED_VERTEX_G3_COUNTS.items()
    }
    expected_lower = {
        str(key): value for key, value in EXPECTED_LOWER_FACE_COUNTS.items()
    }
    if summary["input_by_edges"] != expected_input:
        raise ValueError("unexpected local-link input count vector")
    if summary["vertex_g1_survivors_by_edges"] != expected_g1:
        raise ValueError("unexpected vertex-g1 survivor count vector")
    if summary["vertex_g2_survivors_by_edges"] != expected_g2:
        raise ValueError("unexpected vertex-g2 survivor count vector")
    if summary["vertex_g3_survivors_by_edges"] != expected_g3:
        raise ValueError("unexpected vertex-g3 survivor count vector")
    if summary["lower_face_feasibility_survivors_by_edges"] != expected_lower:
        raise ValueError("unexpected lower-face survivor count vector")
    return summary

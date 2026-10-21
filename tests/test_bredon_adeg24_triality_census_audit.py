"""Focused audit of the Bredon ADEG integral-triality precompiler."""

from __future__ import annotations

from collections import Counter, deque
import hashlib
import json
from pathlib import Path

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from scripts.certify_bredon_adeg24_triality_census import (
    DIAGRAMS,
    exact_quadratic_inertia,
    normalized_symmetric_gram,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/bredon_adeg24_triality_census.json"
SOURCE_TRANSCRIPTION_SHA256 = (
    "f4a540278fd465199b3eeaf5f873c24496adc7c8cd0b215d3af891471ece24ad"
)
ARCHIVE_CANONICAL_SHA256 = (
    "00ce6afa559a4fd84217b6e38cfa57867747da19964bbe0c1e03e94b1fc14722"
)
FAILING_IDS = {
    "ADEG-N2-02-triangle-633",
    "ADEG-N2-04-triangle-666",
    "ADEG-N3-06-square-one-6",
}
NON_STABLE_PACKAGES = {
    ("ADEG-N5-02-path-triangle", 1),
    ("ADEG-N6-01-diamond-tail", 0),
    ("ADEG-N7-01-pentagon-tail", 0),
    ("ADEG-N7-02-path-star", 1),
    ("ADEG-N7-03-three-armed-path", 1),
    ("ADEG-N7-04-triangle-branch-path", 1),
    ("ADEG-N7-05-two-triangles-branch", 0),
    ("ADEG-N9-01-branched-path", 1),
    ("ADEG-N9-02-P-star", 0),
    ("ADEG-N9-02-P-star", 1),
    ("ADEG-N11-01-branched-path", 0),
}


def _archive() -> dict:
    return json.loads(ARCHIVE.read_text(encoding="utf-8"))


def _source_transcription_digest() -> str:
    rows = [
        {
            "id": diagram.id,
            "dimension": diagram.dimension,
            "vertices": diagram.vertices,
            "m3": [list(edge) for edge in diagram.m3],
            "m6": [list(edge) for edge in diagram.m6],
            "source_order": diagram.source_order,
        }
        for diagram in DIAGRAMS
    ]
    encoded = json.dumps(rows, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _independent_triality(diagram) -> tuple[bool, tuple[frozenset[int], ...]]:
    """Contract m3 edges and two-colour the m6 quotient, independently."""

    parent = list(range(diagram.vertices))

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    def union(left: int, right: int) -> None:
        left, right = find(left), find(right)
        if left != right:
            parent[right] = left

    for edge in diagram.m3:
        union(*edge)
    components = {find(vertex) for vertex in range(diagram.vertices)}
    quotient = {component: set() for component in components}
    for left, right in diagram.m6:
        left_component, right_component = find(left), find(right)
        if left_component == right_component:
            return False, ()
        quotient[left_component].add(right_component)
        quotient[right_component].add(left_component)

    colours: dict[int, int] = {}
    for start in sorted(components):
        if start in colours:
            continue
        colours[start] = 0
        queue = deque([start])
        while queue:
            current = queue.popleft()
            for neighbour in quotient[current]:
                expected = 1 - colours[current]
                if neighbour in colours and colours[neighbour] != expected:
                    return False, ()
                if neighbour not in colours:
                    colours[neighbour] = expected
                    queue.append(neighbour)

    colour_classes = tuple(
        frozenset(
            vertex
            for vertex in range(diagram.vertices)
            if colours[find(vertex)] == colour
        )
        for colour in (0, 1)
    )
    return True, colour_classes


def _independent_folded_gram(diagram, long_vertices: frozenset[int]) -> sp.Matrix:
    gram = sp.zeros(diagram.vertices)
    for vertex in range(diagram.vertices):
        gram[vertex, vertex] = 6 if vertex in long_vertices else 2
    for left, right in diagram.m3:
        assert (left in long_vertices) == (right in long_vertices)
        value = -3 if left in long_vertices else -1
        gram[left, right] = gram[right, left] = value
    for left, right in diagram.m6:
        assert (left in long_vertices) != (right in long_vertices)
        gram[left, right] = gram[right, left] = -3
    return gram


def _independent_unfold(diagram, long_vertices: frozenset[int]) -> sp.Matrix:
    labels = [
        (vertex, phase)
        for vertex in range(diagram.vertices)
        for phase in (range(3) if vertex in long_vertices else range(1))
    ]
    position = {label: index for index, label in enumerate(labels)}
    gram = 2 * sp.eye(len(labels))
    for left, right in diagram.m3:
        phases = range(3) if left in long_vertices else range(1)
        for phase in phases:
            i, j = position[left, phase], position[right, phase]
            gram[i, j] = gram[j, i] = -1
    for left, right in diagram.m6:
        long_vertex = left if left in long_vertices else right
        short_vertex = right if left in long_vertices else left
        for phase in range(3):
            i, j = position[long_vertex, phase], position[short_vertex, 0]
            gram[i, j] = gram[j, i] = -1
    return gram


def _nonzero_smith_diagonal(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(smith[index, index]))
        for index in range(min(smith.shape))
        if smith[index, index] != 0
    ]


def test_archive_digest_scope_and_source_inventory() -> None:
    payload = _archive()
    claimed = payload.pop("canonical_sha256")
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    assert hashlib.sha256(encoded).hexdigest() == claimed == ARCHIVE_CANONICAL_SHA256
    assert _source_transcription_digest() == SOURCE_TRANSCRIPTION_SHA256
    assert [diagram.source_order for diagram in DIAGRAMS] == list(range(1, 25))
    assert Counter(diagram.dimension for diagram in DIAGRAMS) == {
        2: 4,
        3: 7,
        5: 4,
        6: 1,
        7: 5,
        9: 2,
        11: 1,
    }
    assert "precompiler" in payload["scope"]
    assert "not a complete K3-nef return theorem" in payload["scope"]


def test_all_24_normalized_grams_have_exact_source_rank_and_inertia() -> None:
    payload = _archive()
    archived_rows = {row["id"]: row for row in payload["rows"]}
    for diagram in DIAGRAMS:
        gram = normalized_symmetric_gram(diagram)
        expected = (
            diagram.dimension,
            1,
            diagram.vertices - diagram.dimension - 1,
        )
        assert gram.rank() == diagram.dimension + 1
        assert exact_quadratic_inertia(gram) == expected
        assert archived_rows[diagram.id]["normalized_symmetric_rank_inertia"] == [
            diagram.dimension + 1,
            list(expected),
        ]


def test_independent_bipartite_quotient_gives_21_rows_and_42_orientations() -> None:
    payload = _archive()
    archived_rows = {row["id"]: row for row in payload["rows"]}
    passing = []
    orientation_count = 0
    for diagram in DIAGRAMS:
        passes, colour_classes = _independent_triality(diagram)
        row = archived_rows[diagram.id]
        assert row["triality_criterion"]["passes"] is passes
        assert passes is (diagram.id not in FAILING_IDS)
        if not passes:
            assert not row["orientations"]
            continue
        passing.append(diagram.id)
        archived_long_sets = {
            frozenset(orientation["long_vertices"])
            for orientation in row["orientations"]
        }
        assert archived_long_sets == set(colour_classes)
        assert colour_classes[0] == frozenset(range(diagram.vertices)) - colour_classes[1]
        assert len(archived_long_sets) == 2
        orientation_count += 2
    assert len(passing) == 21
    assert orientation_count == 42
    assert payload["passing_ids"] == passing
    assert set(payload["failing_ids"]) == FAILING_IDS


def test_all_42_folded_and_unfolded_packages_replay_exactly() -> None:
    payload = _archive()
    diagram_by_id = {diagram.id: diagram for diagram in DIAGRAMS}
    for row in payload["rows"]:
        diagram = diagram_by_id[row["id"]]
        normalized = normalized_symmetric_gram(diagram)
        for orientation in row["orientations"]:
            long_vertices = frozenset(orientation["long_vertices"])
            folded = _independent_folded_gram(diagram, long_vertices)
            assert folded == sp.Matrix(orientation["folded_integral_gram"])
            scale = sp.diag(
                *[
                    sp.sqrt(3) if vertex in long_vertices else 1
                    for vertex in range(diagram.vertices)
                ]
            )
            assert folded == scale * normalized * scale
            assert orientation["folded_rank_inertia"] == [
                folded.rank(),
                list(exact_inertia(folded)),
            ]

            unfolded = _independent_unfold(diagram, long_vertices)
            assert unfolded == sp.Matrix(orientation["unfolded_gram"])
            assert orientation["unfolded_root_count"] == unfolded.rows
            assert orientation["unfolded_rank_inertia"] == [
                unfolded.rank(),
                list(exact_inertia(unfolded)),
            ]


def test_root_quotients_and_31_stable_range_candidates_replay_from_snf() -> None:
    payload = _archive()
    stable_packages = set()
    non_stable_packages = set()
    for row in payload["rows"]:
        for orientation in row["orientations"]:
            package = (row["id"], orientation["orientation"])
            unfolded = sp.Matrix(orientation["unfolded_gram"])
            smith = _nonzero_smith_diagonal(unfolded)
            lattice = orientation["root_generated_lattice"]
            assert smith == lattice["smith_diagonal"]
            assert abs(int(sp.prod(smith))) == abs(lattice["determinant"])
            length = sum(value > 1 for value in smith)
            assert length == lattice["discriminant_length"]

            # The preceding test checks these archived inertias independently.
            # Here the nonzero Smith factors give the rank without repeating
            # the comparatively expensive congruence elimination.
            rank = len(smith)
            inertia = tuple(orientation["unfolded_rank_inertia"][1])
            lorentzian = inertia == (rank - 1, 1, unfolded.rows - rank)
            stable = lorentzian and rank <= 19 and 22 - rank >= length + 2
            assert orientation["unfolded_lorentzian"] is lorentzian
            assert orientation["nikulin_stable_range_primitive_k3_embedding"] is stable
            (stable_packages if stable else non_stable_packages).add(package)
            assert orientation["strong_same_P_k3_return_certified"] is False

    assert len(stable_packages) == 31
    assert non_stable_packages == NON_STABLE_PACKAGES
    assert stable_packages == {
        tuple(package) for package in payload["stable_range_packages"]
    }
    assert payload["counts"] == {
        "published_diagrams": 24,
        "integral_triality_pass": 21,
        "integral_triality_fail": 3,
        "passing_dimension_at_least_5": 13,
        "orientation_packages": 42,
        "stable_range_k3_embedding_packages": 31,
        "strong_same_P_k3_returns": 0,
    }

import json
from itertools import combinations, product
from pathlib import Path

import pytest

from coxeter_polytopes.combinatorics.facet_completion import (
    FULL_VERTEX_MASK,
    audit_global_facet_selection,
    explicit_product_facet_selection,
    forced_small_transversal_obstruction,
    generated_face_table,
    minimal_nonfaces_from_table,
    solve_facet_selection,
)
from coxeter_polytopes.combinatorics.missing_edges import encode_graph6
from coxeter_polytopes.combinatorics.nerve_support import analyze_nerve_support


ROOT = Path(__file__).resolve().parents[1]
SUPPORT_MANIFEST = ROOT / "data" / "enumerations" / "h9_14_nerve_support.json"
REFINEMENT_MANIFEST = ROOT / "data" / "enumerations" / "h9_14_facet_refinement.json"
SCHEDULE_SPLIT_MANIFEST = ROOT / "data" / "enumerations" / "h9_14_schedule_split.json"


def _mask(vertices) -> int:
    return sum(1 << vertex for vertex in vertices)


def _matching_product_fixture() -> tuple[str, tuple[int, ...]]:
    missing_edges = {(0, 1), (2, 3), (4, 5)}
    simplex_blocks = ((6, 7, 8, 9), (10, 11, 12, 13))
    facets = []
    for binary_choices in product((0, 1), (2, 3), (4, 5)):
        for omitted_first in simplex_blocks[0]:
            for omitted_second in simplex_blocks[1]:
                facet = set(binary_choices)
                facet.update(set(simplex_blocks[0]) - {omitted_first})
                facet.update(set(simplex_blocks[1]) - {omitted_second})
                facets.append(_mask(facet))
    complements = tuple(sorted(FULL_VERTEX_MASK ^ facet for facet in facets))
    return encode_graph6(14, missing_edges), complements


def _pentagon_product_fixture() -> tuple[str, tuple[int, ...]]:
    pentagon_edges = {(0, 1), (1, 2), (2, 3), (3, 4), (0, 4)}
    missing_edges = {
        pair
        for pair in combinations(range(5), 2)
        if pair not in pentagon_edges
    }
    tetrahedron = (5, 6, 7, 8)
    four_simplex = (9, 10, 11, 12, 13)
    facets = []
    for pentagon_edge in pentagon_edges:
        for omitted_first in tetrahedron:
            for omitted_second in four_simplex:
                facet = set(pentagon_edge)
                facet.update(set(tetrahedron) - {omitted_first})
                facet.update(set(four_simplex) - {omitted_second})
                facets.append(_mask(facet))
    complements = tuple(sorted(FULL_VERTEX_MASK ^ facet for facet in facets))
    return encode_graph6(14, missing_edges), complements


@pytest.mark.parametrize(
    ("fixture", "expected_missing_sizes"),
    [
        (_matching_product_fixture, {2: 3, 4: 2}),
        (_pentagon_product_fixture, {2: 5, 4: 1, 5: 1}),
    ],
)
def test_exact_global_audit_accepts_polytopal_product_fixtures(
    fixture, expected_missing_sizes: dict[int, int]
) -> None:
    graph6, selected = fixture()
    assert explicit_product_facet_selection(graph6) == selected
    support = analyze_nerve_support(graph6)
    audit = audit_global_facet_selection(
        graph6,
        support.ridge_support_core_masks,
        support.dehn_sommerville_schedule,
        selected,
    )
    assert audit.passes_selection_model is True
    assert audit.passes_connected_euler_audit is True
    assert audit.passes_connected_homology_audit is True
    assert audit.first_strict_local_g_failure is None
    assert audit.first_local_macaulay_failure is None
    assert audit.passes_connected_homology_audit is True
    assert {
        size: len(values)
        for size, values in audit.minimal_nonfaces_by_size.items()
    } == expected_missing_sizes


def test_large_minimal_nonface_is_detected_exactly() -> None:
    # Nerve of Delta_5 x Delta_1^4: a join of a 5-simplex boundary and four
    # 0-spheres.  It is a genuine simplicial 8-sphere, but its six-element
    # missing face violates the Lanner order bound.
    large_block = tuple(range(6))
    pairs = ((6, 7), (8, 9), (10, 11), (12, 13))
    facets = []
    for omitted in large_block:
        for choices in product(*pairs):
            facet = (set(large_block) - {omitted}) | set(choices)
            facets.append(_mask(facet))
    complements = tuple(sorted(FULL_VERTEX_MASK ^ facet for facet in facets))
    nonfaces = minimal_nonfaces_from_table(generated_face_table(complements))
    assert {size: len(values) for size, values in nonfaces.items()} == {
        2: 4,
        6: 1,
    }


def test_small_archived_shadow_is_exactly_unsatisfiable() -> None:
    pytest.importorskip("z3")
    manifest = json.loads(SUPPORT_MANIFEST.read_text(encoding="utf-8"))
    record = min(
        (entry for entry in manifest["records"] if entry["survives_second_layer"]),
        key=lambda entry: entry["ridge_support_core_count"],
    )
    result = solve_facet_selection(
        record["graph6"],
        record["ridge_support_core_masks"],
        record["dehn_sommerville_schedule"],
        timeout_ms=10_000,
    )
    assert result.status == "unsat"
    assert result.selected_facet_complement_masks == ()


def test_forced_empty_higher_nonfaces_give_small_transversal_certificate() -> None:
    manifest = json.loads(SUPPORT_MANIFEST.read_text(encoding="utf-8"))
    survivors = [
        entry for entry in manifest["records"] if entry["survives_second_layer"]
    ]
    record = survivors[134]
    certificate = forced_small_transversal_obstruction(
        record["graph6"], record["dehn_sommerville_schedule"]
    )
    assert certificate == {
        "missing_four_count": 0,
        "missing_five_count": 0,
        "missing_edge_count": 4,
        "minimal_vertex_cover_upper_bound": 4,
        "required_facet_complement_size": 5,
    }


def test_refinement_archive_has_unambiguous_theorem_partition() -> None:
    support_manifest = json.loads(SUPPORT_MANIFEST.read_text(encoding="utf-8"))
    survivors = [
        entry
        for entry in support_manifest["records"]
        if entry["survives_second_layer"]
    ]
    refinement = json.loads(REFINEMENT_MANIFEST.read_text(encoding="utf-8"))
    assert refinement["certificate_sha256"] == (
        "c976ac262bb78dc4020210e3eab8ed662d5719de0d058f1a94089c60c7f4ae25"
    )
    assert refinement["counts"] == {
        "input": 142,
        "certified_unsat": 113,
        "certified_explicit_witness": 2,
        "theorem_undecided": 27,
        "bounded_timeout_unknown": 25,
        "backend_unsat_pending_independent_certificate": 2,
    }
    for index in (30, 42):
        record = refinement["records"][index]
        assert record["theorem_status"] == "certified_explicit_witness"
        source = survivors[index]
        audit = audit_global_facet_selection(
            source["graph6"],
            source["ridge_support_core_masks"],
            source["dehn_sommerville_schedule"],
            record["selected_facet_complement_masks"],
        )
        assert audit.passes_connected_homology_audit is True
        assert record["witness_audit"] == audit.to_record()
    assert refinement["records"][134]["theorem_status"] == "certified_unsat"
    assert refinement["records"][134]["small_transversal_certificate"] == (
        forced_small_transversal_obstruction(
            survivors[134]["graph6"],
            survivors[134]["dehn_sommerville_schedule"],
        )
    )


def test_schedule_split_archive_never_promotes_backend_unsat() -> None:
    manifest = json.loads(SCHEDULE_SPLIT_MANIFEST.read_text(encoding="utf-8"))
    assert manifest["certificate_sha256"] == (
        "2d4e6aa60912a9060c62a07a8d6997dda00fd4eec82d60e850c41fab523f2e1c"
    )
    records = {record["index"]: record for record in manifest["records"]}
    assert records[92]["aggregate_backend_status"] == "unsat"
    assert records[93]["aggregate_backend_status"] == "unknown"
    assert all(record["theorem_status"] == "undecided" for record in records.values())

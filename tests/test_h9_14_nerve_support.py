import json
from itertools import combinations, product
from pathlib import Path

import pytest

from coxeter_polytopes.combinatorics.missing_edges import (
    analyze_missing_edge_graph,
    encode_graph6,
)
from coxeter_polytopes.combinatorics.nerve_support import (
    EXPECTED_CERTIFICATE_SHA256,
    EXPECTED_INPUT_COUNTS,
    EXPECTED_JOINT_SURVIVOR_COUNTS,
    EXPECTED_RIDGE_SUPPORT_COUNTS,
    EXPECTED_STRICT_GLBT_COUNTS,
    GALE_RANK,
    analyze_nerve_support,
    audit_facet_completion,
    dehn_sommerville_f_vector,
    verify_support_manifest,
)


ROOT = Path(__file__).resolve().parents[1]
FIRST_LAYER = ROOT / "data" / "enumerations" / "h9_14_missing_edges.json"
MANIFEST = ROOT / "data" / "enumerations" / "h9_14_nerve_support.json"
FULL_MASK = (1 << 14) - 1


def _mask(vertices: set[int] | tuple[int, ...]) -> int:
    return sum(1 << vertex for vertex in vertices)


def _matching_product_fixture() -> tuple[str, tuple[int, ...]]:
    # Nerve of Delta_1^3 x Delta_3^2: join of three 0-spheres and two
    # tetrahedral boundary spheres.
    missing_edges = {(0, 1), (2, 3), (4, 5)}
    simplex_blocks = ((6, 7, 8, 9), (10, 11, 12, 13))
    facets = []
    for binary_choices in product(*((block[0], block[1]) for block in ((0, 1), (2, 3), (4, 5)))):
        for omitted_first in simplex_blocks[0]:
            for omitted_second in simplex_blocks[1]:
                facet = set(binary_choices)
                facet.update(set(simplex_blocks[0]) - {omitted_first})
                facet.update(set(simplex_blocks[1]) - {omitted_second})
                facets.append(_mask(facet))
    complements = tuple(sorted(FULL_MASK ^ facet for facet in facets))
    return encode_graph6(14, missing_edges), complements


def _pentagon_product_fixture() -> tuple[str, tuple[int, ...]]:
    # Nerve of C5 x Delta_3 x Delta_4.
    pentagon_edges = ((0, 1), (1, 2), (2, 3), (3, 4), (0, 4))
    missing_edges = {
        pair for pair in combinations(range(5), 2) if pair not in pentagon_edges
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
    complements = tuple(sorted(FULL_MASK ^ facet for facet in facets))
    return encode_graph6(14, missing_edges), complements


def _independent_five_sets(graph6: str) -> tuple[int, ...]:
    from coxeter_polytopes.combinatorics.missing_edges import decode_graph6

    _, edges = decode_graph6(graph6)
    edge_masks = tuple((1 << first) | (1 << second) for first, second in edges)
    return tuple(
        _mask(vertices)
        for vertices in combinations(range(14), 5)
        if not any(_mask(vertices) & edge == edge for edge in edge_masks)
    )


def test_gale_rank_and_strict_generalized_lower_bound_obstruction() -> None:
    assert GALE_RANK == 14 - (9 + 1) == 4
    first_layer = json.loads(FIRST_LAYER.read_text(encoding="utf-8"))
    forced_g4_zero = next(
        record
        for record in first_layer["records"]
        if record["survives_graph_prefilter"]
        and record["missing_four_count_min"]
        == record["missing_four_count_max"]
        == record["g4_if_no_missing_fours"]
    )
    audit = analyze_nerve_support(forced_g4_zero["graph6"])
    assert audit.passes_strict_generalized_lower_bound is False
    assert audit.survives_second_layer is False


@pytest.mark.parametrize(
    ("fixture", "expected_facets", "expected_ridges", "missing_sizes", "actual_a", "actual_b", "expected_f"),
    [
        (
            _matching_product_fixture,
            128,
            576,
            {2: 3, 4: 2},
            2,
            0,
            (14, 88, 328, 804, 1352, 1568, 1216, 576, 128),
        ),
        (
            _pentagon_product_fixture,
            100,
            450,
            {2: 5, 4: 1, 5: 1},
            1,
            1,
            (14, 86, 309, 725, 1165, 1295, 970, 450, 100),
        ),
    ],
)
def test_polytopal_product_fixtures_survive_and_realize_ds_schedule(
    fixture,
    expected_facets: int,
    expected_ridges: int,
    missing_sizes: dict[int, int],
    actual_a: int,
    actual_b: int,
    expected_f: tuple[int, ...],
) -> None:
    graph6, complements = fixture()
    support = analyze_nerve_support(graph6)
    assert support.survives_second_layer is True
    assert all(complement in support.ridge_support_core_masks for complement in complements)

    completion = audit_facet_completion(graph6, complements)
    assert completion.facet_count == expected_facets
    assert completion.ridge_count == expected_ridges
    assert completion.all_ridge_degrees_two is True
    assert completion.dual_graph_connected is True
    assert completion.covers_all_independent_triples is True
    assert {
        size: len(nonfaces)
        for size, nonfaces in completion.minimal_nonfaces_by_size.items()
    } == missing_sizes

    first = analyze_missing_edge_graph(graph6)
    assert dehn_sommerville_f_vector(first, actual_a) == expected_f
    schedule = next(
        entry
        for entry in support.dehn_sommerville_schedule
        if entry["missing_four_count"] == actual_a
    )
    missing_fours = completion.minimal_nonfaces_by_size[4]
    covered_by_a4 = sum(
        any(five_set & missing_four == missing_four for missing_four in missing_fours)
        for five_set in _independent_five_sets(graph6)
    )
    assert schedule["missing_five_count_offset"] - covered_by_a4 == actual_b


def test_six_sets_with_no_selected_complement_are_not_asserted_ridges() -> None:
    graph6, complements = _matching_product_fixture()
    selected = set(complements)
    empty_six_sets = [
        vertices
        for vertices in combinations(range(14), 6)
        if not any(_mask(set(vertices) - {removed}) in selected for removed in vertices)
    ]
    assert empty_six_sets
    assert audit_facet_completion(graph6, complements).all_ridge_degrees_two is True


def test_archived_nerve_support_manifest_verifies() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    assert manifest["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    summary = verify_support_manifest(manifest)
    assert summary["input_total"] == sum(EXPECTED_INPUT_COUNTS.values()) == 444
    assert summary["strict_glbt_survivors_total"] == sum(
        EXPECTED_STRICT_GLBT_COUNTS.values()
    ) == 263
    assert summary["ridge_support_coverage_survivors_total"] == sum(
        EXPECTED_RIDGE_SUPPORT_COUNTS.values()
    ) == 197
    assert summary["second_layer_survivors_total"] == sum(
        EXPECTED_JOINT_SURVIVOR_COUNTS.values()
    ) == 142


def test_nerve_support_manifest_tampering_is_detected() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    manifest["records"][0]["ridge_support_core_masks"].pop()
    with pytest.raises(ValueError, match="invariants do not verify"):
        verify_support_manifest(manifest)

    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    manifest["target"]["facet_count"] = 99
    with pytest.raises(ValueError, match="target metadata"):
        verify_support_manifest(manifest)

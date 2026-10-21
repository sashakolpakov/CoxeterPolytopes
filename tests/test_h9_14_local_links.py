import json
from pathlib import Path

import pytest

from coxeter_polytopes.combinatorics.local_link import (
    EXPECTED_CERTIFICATE_SHA256,
    EXPECTED_INPUT_COUNTS,
    EXPECTED_LOWER_FACE_COUNTS,
    EXPECTED_VERTEX_G1_COUNTS,
    EXPECTED_VERTEX_G2_COUNTS,
    EXPECTED_VERTEX_G3_COUNTS,
    analyze_local_link_graph,
    edge_figure_constraints,
    triple_figure_constraints,
    verify_local_link_manifest,
)
from coxeter_polytopes.combinatorics.missing_edges import encode_graph6


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "enumerations" / "h9_14_local_links.json"


def _missing_matching_graph() -> str:
    return encode_graph6(14, {(0, 1), (2, 3), (4, 5)})


def _missing_pentagon_graph() -> str:
    cycle_edges = {(0, 1), (1, 2), (2, 3), (3, 4), (0, 4)}
    diagonals = {
        pair
        for first in range(5)
        for second in range(first + 1, 5)
        if (pair := (first, second)) not in cycle_edges
    }
    return encode_graph6(14, diagonals)


@pytest.mark.parametrize("graph6", [_missing_matching_graph(), _missing_pentagon_graph()])
def test_polytopal_product_shadows_pass_local_link_filter(graph6: str) -> None:
    audit = analyze_local_link_graph(graph6)
    assert audit.passes_all_vertex_g2 is True
    assert audit.passes_all_vertex_g3_feasibility is True
    assert audit.passes_edge_and_triple_feasibility is True


def test_matching_fixture_missing_fours_obey_all_linear_intervals() -> None:
    graph6 = _missing_matching_graph()
    missing_fours = ({6, 7, 8, 9}, {10, 11, 12, 13})
    audit = analyze_local_link_graph(graph6)
    for figure in audit.vertex_figures:
        incidence = sum(figure.vertex in missing for missing in missing_fours)
        assert figure.missing_four_incidence_interval is not None
        lower, upper = figure.missing_four_incidence_interval
        assert lower <= incidence <= upper
    for figure in edge_figure_constraints(graph6):
        incidence = sum(set(figure.edge) <= missing for missing in missing_fours)
        assert figure.missing_four_incidence_interval is not None
        lower, upper = figure.missing_four_incidence_interval
        assert lower <= incidence <= upper
    for figure in triple_figure_constraints(graph6):
        incidence = sum(set(figure.triple) <= missing for missing in missing_fours)
        assert incidence <= figure.missing_four_incidence_max


def test_simplex_vertex_figure_is_rejected_by_inherited_missing_face_bound() -> None:
    # At vertex 0 the four missing edges leave exactly nine link vertices,
    # giving the boundary of an 8-simplex. Its own missing face has order
    # nine, so it violates the inherited Lanner bound five.
    graph6 = encode_graph6(14, {(0, 1), (0, 2), (0, 3), (0, 4)})
    audit = analyze_local_link_graph(graph6)
    figure = audit.vertex_figures[0]
    assert figure.g1 == 0
    assert figure.passes_strict_g1 is False
    assert audit.passes_all_vertex_g1 is False
    assert audit.passes_all_vertex_g2 is False


def test_archived_local_link_manifest_verifies_and_freezes_counts() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    assert manifest["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    summary = verify_local_link_manifest(manifest)
    assert summary["input_total"] == sum(EXPECTED_INPUT_COUNTS.values()) == 142
    assert summary["vertex_g1_survivors_total"] == sum(
        EXPECTED_VERTEX_G1_COUNTS.values()
    ) == 72
    assert summary["vertex_g2_survivors_total"] == sum(
        EXPECTED_VERTEX_G2_COUNTS.values()
    ) == 39
    assert summary["vertex_g3_survivors_total"] == sum(
        EXPECTED_VERTEX_G3_COUNTS.values()
    ) == 36
    assert summary["lower_face_feasibility_survivors_total"] == sum(
        EXPECTED_LOWER_FACE_COUNTS.values()
    ) == 30


def test_local_link_manifest_tampering_is_detected() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    manifest["records"][0]["vertex_figures"][0]["g2"] += 1
    with pytest.raises(ValueError, match="invariants do not verify"):
        verify_local_link_manifest(manifest)

    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    manifest["records"].pop()
    with pytest.raises(ValueError, match="invariants|graph transcript|input count"):
        verify_local_link_manifest(manifest)

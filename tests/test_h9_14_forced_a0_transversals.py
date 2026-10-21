import json
from pathlib import Path

import pytest

from coxeter_polytopes.combinatorics.forced_a0_transversals import (
    EXPECTED_BAD_TRANSVERSAL_SIZE_HISTOGRAM,
    EXPECTED_CERTIFICATE_SHA256,
    EXPECTED_CHOICE_TOTAL,
    EXPECTED_COMBINED_SURVIVORS,
    EXPECTED_FORCED_GRAPH_COUNTS,
    EXPECTED_FORCED_GRAPH_TOTAL,
    EXPECTED_LOCAL_LINK_INTERSECTION,
    EXPECTED_MISSING_FIVE_COUNT_HISTOGRAM,
    forced_missing_five_count,
    independent_masks,
    is_minimal_transversal,
    is_transversal,
    strict_missing_four_interval,
    verify_forced_a0_manifest,
)
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "data" / "enumerations" / "h9_14_nerve_support.json"
MANIFEST = (
    ROOT / "data" / "enumerations" / "h9_14_forced_a0_transversals.json"
)


def _second_layer_graphs() -> tuple[str, ...]:
    source = json.loads(SUPPORT.read_text(encoding="utf-8"))
    return tuple(
        record["graph6"]
        for record in source["records"]
        if record["survives_second_layer"]
    )


def test_complements_exchange_maximal_faces_and_minimal_transversals() -> None:
    # Five disjoint minimal nonfaces model a pure complex whose maximal-face
    # complements choose exactly one point from every block.
    hyperedges = (
        (1 << 0) | (1 << 1),
        (1 << 2) | (1 << 3),
        (1 << 4) | (1 << 5),
        sum(1 << vertex for vertex in (6, 7, 8, 9)),
        sum(1 << vertex for vertex in (10, 11, 12, 13)),
    )
    transversal = sum(1 << vertex for vertex in (0, 2, 4, 6, 10))
    assert is_transversal(transversal, hyperedges)
    assert is_minimal_transversal(transversal, hyperedges)
    assert transversal.bit_count() == 5
    assert not is_transversal(transversal ^ (1 << 0), hyperedges)
    assert not is_minimal_transversal(transversal | (1 << 1), hyperedges)


def test_forced_a0_stratum_and_a5_choice_counts_are_rederived() -> None:
    forced = [
        graph6
        for graph6 in _second_layer_graphs()
        if strict_missing_four_interval(graph6) == (0, 0)
    ]
    assert len(forced) == EXPECTED_FORCED_GRAPH_TOTAL == 36
    edge_histogram: dict[int, int] = {}
    b_histogram: dict[int, int] = {}
    choice_total = 0
    for graph6 in forced:
        _, edges = decode_graph6(graph6)
        edge_count = len(edges)
        edge_histogram[edge_count] = edge_histogram.get(edge_count, 0) + 1
        b, _, _ = forced_missing_five_count(graph6)
        b_histogram[b] = b_histogram.get(b, 0) + 1
        choice_total += len(independent_masks(graph6, 5)) if b == 1 else 1
    assert edge_histogram == EXPECTED_FORCED_GRAPH_COUNTS
    assert b_histogram == EXPECTED_MISSING_FIVE_COUNT_HISTOGRAM
    assert choice_total == EXPECTED_CHOICE_TOTAL == 30_934


def test_archived_forced_a0_manifest_verifies_every_choice() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    assert manifest["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    summary = verify_forced_a0_manifest(manifest, _second_layer_graphs())
    assert summary["forced_a0_graphs_total"] == 36
    assert summary["a5_choice_total"] == 30_934
    assert summary["bad_transversal_size_histogram"] == {
        str(size): count
        for size, count in EXPECTED_BAD_TRANSVERSAL_SIZE_HISTOGRAM.items()
    }
    assert summary["forced_a0_local_link_intersection"] == (
        EXPECTED_LOCAL_LINK_INTERSECTION
    ) == 1
    assert summary["combined_local_link_survivors"] == (
        EXPECTED_COMBINED_SURVIVORS
    ) == 29


def test_unique_b0_graph_and_all_bad_witnesses_have_size_three_or_four() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    b0_records = [
        record for record in manifest["records"] if record["forced_missing_five_count"] == 0
    ]
    assert len(b0_records) == 1
    assert b0_records[0]["graph6"] == "M???????C?G?G?C??"
    assert b0_records[0]["choice_certificates"][0][0] == 0
    sizes = {
        transversal.bit_count()
        for record in manifest["records"]
        for _, transversal in record["choice_certificates"]
    }
    assert sizes == {3, 4}


def test_forced_a0_manifest_tampering_is_detected() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    manifest["records"][0]["choice_certificates"].pop()
    with pytest.raises(ValueError, match="choice space is not exhaustive"):
        verify_forced_a0_manifest(manifest, _second_layer_graphs())

    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    manifest["records"][0]["choice_certificates"][0][1] = 31
    with pytest.raises(ValueError, match="target size five"):
        verify_forced_a0_manifest(manifest, _second_layer_graphs())

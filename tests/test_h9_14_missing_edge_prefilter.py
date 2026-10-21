from hashlib import sha256
import json
from pathlib import Path

import pytest

from coxeter_polytopes.combinatorics.missing_edges import (
    EXPECTED_GRAPH6_STDOUT_SHA256,
    EXPECTED_UNLABELED_COUNTS,
    analyze_missing_edge_graph,
    decode_graph6,
    encode_graph6,
    macaulay_raise,
    verify_manifest,
)


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "enumerations" / "h9_14_missing_edges.json"


def test_graph6_round_trip_and_macaulay_values() -> None:
    edges = {(0, 1), (0, 13), (6, 9), (12, 13)}
    encoded = encode_graph6(14, edges)
    vertex_count, decoded = decode_graph6(encoded)
    assert vertex_count == 14
    assert set(decoded) == edges
    assert [macaulay_raise(value, 2) for value in range(11)] == [
        0,
        1,
        2,
        4,
        5,
        7,
        10,
        11,
        13,
        16,
        20,
    ]


def test_two_product_polytope_missing_edge_graphs_survive() -> None:
    # Delta_1^3 x Delta_3^2: the three missing pairs form a matching.
    matching = encode_graph6(14, {(0, 1), (2, 3), (4, 5)})
    audit = analyze_missing_edge_graph(matching)
    assert (audit.g2, audit.g3, audit.macaulay_g3_bound) == (7, 8, 11)
    assert audit.admits_missing_four_count is True
    assert audit.missing_four_count_min <= 2 <= audit.missing_four_count_max
    assert audit.survives_graph_prefilter is True

    # Pentagon x Delta_3 x Delta_4: the missing-pair graph is a 5-cycle.
    cycle = encode_graph6(14, {(0, 1), (1, 2), (2, 3), (3, 4), (0, 4)})
    audit = analyze_missing_edge_graph(cycle)
    assert (audit.g2, audit.g3, audit.macaulay_g3_bound) == (5, 5, 7)
    assert audit.admits_missing_four_count is True
    assert audit.missing_four_count_min <= 1 <= audit.missing_four_count_max
    assert audit.survives_graph_prefilter is True


def test_g2_zero_stratum_is_rejected_after_macaulay() -> None:
    # Two disjoint K_1,5 stars and two isolated vertices have e=10, g3=0.
    edges = {(leaf, 12) for leaf in range(5)} | {
        (leaf, 13) for leaf in range(5, 10)
    }
    audit = analyze_missing_edge_graph(encode_graph6(14, edges))
    assert audit.g2 == 0
    assert audit.g3 == 0
    assert audit.passes_g3_macaulay is True
    assert audit.g4_if_no_missing_fours == 20
    assert (audit.missing_four_count_min, audit.missing_four_count_max) == (20, 20)
    assert audit.passes_g_vector_prefilter is True
    assert audit.stacked_obstruction is True
    assert audit.survives_graph_prefilter is False


def test_archived_manifest_verifies_every_graph_without_nauty() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    assert manifest["graph6_stdout_sha256"] == EXPECTED_GRAPH6_STDOUT_SHA256
    summary = verify_manifest(manifest)
    assert summary["generated_total"] == sum(EXPECTED_UNLABELED_COUNTS.values())
    assert summary["g3_macaulay_survivors_total"] == 592
    assert summary["g4_feasible_survivors_total"] == 575
    assert summary["stacked_rejections_after_g4_total"] == 131
    assert summary["graph_prefilter_survivors_total"] == 444


def test_manifest_tampering_is_detected() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    manifest["records"][0]["g3"] += 1
    with pytest.raises(ValueError, match="invariants do not verify"):
        verify_manifest(manifest)


def test_self_consistent_transcript_replacement_is_detected() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    records = manifest["records"]
    records[0], records[1] = records[1], records[0]
    replacement = "".join(f"{record['graph6']}\n" for record in records).encode(
        "ascii"
    )
    digest = sha256(replacement).hexdigest()
    manifest["graph6_stdout_sha256"] = digest
    manifest["generator"]["stdout_sha256"] = digest
    with pytest.raises(ValueError, match="frozen archive"):
        verify_manifest(manifest)

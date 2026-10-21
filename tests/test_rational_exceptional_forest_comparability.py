"""Exact regressions for exceptional-forest parity comparability."""

from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_rational_exceptional_forest_comparability import (
    OUTPUT,
    build_payload,
    digest,
    incomparable,
    proximity_children,
    weights_equalizing_incomparable_pair,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results" / "rational_exceptional_forest_comparability_theorem.md"


def test_archive_and_exhaustive_forest_regression_are_exact():
    archive = json.loads(OUTPUT.read_text())
    payload = build_payload()
    assert archive["payload"] == payload
    assert archive["payload_sha256"] == digest(payload)
    records = payload["exhaustive_chronological_forest_checks"]
    assert [record["chronological_forests"] for record in records] == [
        2,
        6,
        24,
        120,
        720,
        5040,
    ]
    assert all(record["incomparable_pairs_checked"] > 0 for record in records)


def test_branched_and_satellite_proximity_weights_equalize():
    proximity = (
        (),
        (0,),
        (0,),
        (1, 0),
        (2, 0),
        (3, 1),
        (4, 2),
        (5, 3),
        (6, 4),
    )
    for left, right in ((1, 2), (3, 4), (5, 6), (7, 8)):
        assert incomparable(left, right, proximity)
        weights = weights_equalizing_incomparable_pair(
            proximity, left, right
        )
        assert weights[left] == weights[right]
        children = proximity_children(proximity)
        assert all(
            weights[vertex]
            > sum(weights[child] for child in children[vertex])
            for vertex in range(len(proximity))
        )


def test_reflection_and_parity_width_boundaries():
    payload = build_payload()
    for record in payload["reflection_samples"]:
        assert record["surface_root_square"] == -2
        assert record["doubled_host_root_square"] == -4
        assert record["reflection_integral_on_NS"] is True
        assert record["reflection_swaps_total_exceptional_coordinates"] is True

    parity = payload["adapted_host_parity"]
    assert parity["even_pairings_sufficient"] is True
    assert parity["even_pairings_necessary"] is True
    assert parity["odd_absorption_impossible_in_integral_M"] is True

    widths = payload["parity_width_corollary"]
    assert [row["maximal_antichain_bound"] for row in widths] == [
        2 ** count for count in range(17)
    ]
    assert widths[0]["direct_lane_single_chain"] is True
    assert all(
        not row["direct_lane_single_chain"] for row in widths[1:]
    )


def test_report_states_exact_theorem_and_parity_remainder():
    text = REPORT.read_text()
    required = [
        "Lemma 2.1",
        "Theorem 3.1 (exceptional-forest comparability)",
        r"\(\lambda_i=\lambda_j\)",
        "A_t.(e_i-e_j)=0",
        "mirror crosses the interior",
        "properly contains",
        "evenness in (3.2) is not merely sufficient",
        r"\operatorname {width}(V)\le2^b",
        "the exceptional forest is one chain",
        "width theorem, not yet a length or dimension bound",
        "Without a common invariant lattice",
    ]
    for phrase in required:
        assert phrase in text

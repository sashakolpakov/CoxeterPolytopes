from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_h10_two_k3_redoubling_rays_hostile_audit import (
    OUTPUT,
    _expected_mirrors,
    _expected_relations,
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]


def test_frozen_archive_rebuilds_exactly() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert build_certificate() == archived


def test_closed_form_small_descendants_are_frozen() -> None:
    assert _expected_relations("F", 2, 0) == {
        "0": 3,
        "1": 2,
        "2": 75,
        "3": 10,
        "4": 1,
    }
    assert _expected_relations("F", 2, 1) == {
        "0": 2,
        "1": 2,
        "2": 74,
        "3": 10,
        "4": 3,
    }
    assert _expected_relations("S", 4, 0) == {
        "0": 3,
        "1": 9,
        "2": 93,
        "3": 11,
        "4": 4,
    }
    assert _expected_mirrors("F", 4, 0) == {
        "square_-4_div_2": 11,
        "square_-2_div_1": 3,
        "square_-2_div_2": 2,
    }


def test_lineage_and_scope_boundaries_are_explicit() -> None:
    archive = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert archive["all_checks_pass"] is True
    candidates = archive["bounded_lineage"]["candidates"]
    assert len(candidates["F_N2_B"]["histogram_matches"]) == 1
    assert not candidates["F_N2_B"][
        "exact_0_1_distinguished_labelled_graph_matches"
    ]
    assert candidates["S_N2_A"][
        "exact_0_1_distinguished_labelled_graph_matches"
    ]
    report = (
        ROOT / "results/h10_two_k3_redoubling_rays_hostile_audit.md"
    ).read_text(encoding="utf-8")
    for phrase in (
        "bounded project-level result",
        "period-level",
        "not\n  equation-level",
        "full nef cone",
        "surface automorphisms",
    ):
        assert phrase in report

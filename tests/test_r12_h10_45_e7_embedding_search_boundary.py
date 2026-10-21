from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_r12_h10_45_e7_embedding_search_boundary import (
    OUTPUT,
    build_certificate,
)


def archive() -> dict:
    return json.loads(Path(OUTPUT).read_text(encoding="utf-8"))


def test_archive_rebuilds_exactly() -> None:
    assert build_certificate() == archive()


def test_explicit_e7_embedding_is_primitive_but_does_not_glue() -> None:
    row = archive()["primitive_near_miss"]
    assert row["anti_lattice"] == "E7(-2)"
    assert row["e7_embedding_smith"] == [1] * 7
    assert row["eigensublattice_index"] == 8
    assert row["fixed_lattice_signature_determinant"] == [[2, 1], -16]
    assert row["fixed_lattice_norms_divisible_by_four"] is True
    assert row["fixed_q2_histogram"] == {"0": 16, "1": 0, "2": 16, "3": 0}
    assert row["equivariant_glue_exists"] is False


def test_target_picard_action_has_the_incompatible_fixed_q2_type() -> None:
    target = archive()["archived_picard_action"]
    assert target["coinvariant_lattice"] == "<-4>"
    assert target["coinvariant_divisibility"] == 2
    assert target["fixed_q2_histogram"] == {"0": 8, "1": 8, "2": 8, "3": 8}


def test_index_four_false_positive_is_rejected_by_primitivity() -> None:
    false = archive()["index_four_false_positive"]
    assert false["coarse_overlattice_invariants"]["determinant"] == 64
    assert false["coarse_overlattice_invariants"]["integral_reversal"] is True
    assert false["K_embedding_smith"] == [1] * 6 + [2]
    assert false["K_saturation_index"] == 2
    assert false["rejected_by_primitive_embedding_gate"] is True


def test_vinberg_handoff_remains_literal_and_positive() -> None:
    volume = archive()["vinberg_handoff"]
    assert (volume["dimension"], volume["facet_count"]) == (10, 45)
    assert volume["finite_volume"] is True
    assert volume["compact"] is False
    assert (volume["finite_vertex_count"], volume["ideal_vertex_count"]) == (
        1048,
        93,
    )
    assert volume["elliptic_edge_count"] == 7332
    assert volume["bad_edge_count"] == 0
    assert volume["bad_parabolic_neighborhood_count"] == 0


def test_bounded_search_is_not_misstated_as_exhaustive() -> None:
    search = archive()["bounded_discovery_search"]
    assert search["theorem_status"] == "evidence_only_not_used_for_exhaustion"
    assert search["fixed_central_run"]["fixed_q2_histograms"] == {
        "16,0,16,0": 9046
    }
    assert search["random_central_run"]["fixed_q2_histograms"] == {
        "16,0,16,0": 6598
    }


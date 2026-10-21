from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_roulleau_non2elementary_involution_folds import (
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT / "data/enumerations/roulleau_non2elementary_involution_folds.json"
)


def _cases() -> dict[str, dict]:
    artifact = json.loads(ARTIFACT.read_text())
    return {case["id"]: case for case in artifact["cases"]}


def test_frozen_certificate_rebuilds_exactly() -> None:
    assert build_certificate() == json.loads(ARTIFACT.read_text())


def test_both_divisor_folds_are_genuine_finite_volume_polytopes() -> None:
    cases = _cases()
    assert set(cases) == {"R12_UA2E8_C2_FIXED", "R13_UE8A3_C2_FIXED"}
    expected = {
        "R12_UA2E8_C2_FIXED": (10, 12, 99, 18, 1),
        "R13_UE8A3_C2_FIXED": (11, 13, 162, 26, 2),
    }
    for identifier, (dimension, facets, edges, finite, cusps) in expected.items():
        case = cases[identifier]
        geometry = case["geometry"]
        vinberg = geometry["vinberg_volume_criterion"]
        assert geometry["hyperbolic_dimension"] == dimension
        assert geometry["facet_count"] == facets
        assert vinberg["elliptic_edge_count"] == edges
        assert geometry["finite_vertex_count"] == finite
        assert geometry["ideal_vertex_count"] == cusps
        assert vinberg["decision"] == "finite_volume_noncompact"
        assert vinberg["finite_volume_bad_edges"] == []
        assert geometry["all_checks_pass"] is True


def test_fixed_lattices_and_folded_reflections_are_exact() -> None:
    cases = _cases()
    first = cases["R12_UA2E8_C2_FIXED"]
    second = cases["R13_UE8A3_C2_FIXED"]
    assert first["fixed_lattice"]["determinant"] == 2
    assert first["fixed_lattice"]["smith_invariants"] == [1] * 10 + [2]
    assert first["folded_walls"]["root_type_histogram"] == {
        "-2,1": 11,
        "-2,2": 1,
    }
    assert second["fixed_lattice"]["determinant"] == -4
    assert second["fixed_lattice"]["smith_invariants"] == [1] * 10 + [2, 2]
    assert second["folded_walls"]["root_type_histogram"] == {
        "-2,1": 11,
        "-2,2": 1,
        "-4,2": 1,
    }
    assert all(
        case["folded_walls"]["reflection_audit"][
            "all_folded_reflections_are_integral"
        ]
        for case in cases.values()
    )


def test_lineage_separates_old_abstract_chambers_from_new_realizations() -> None:
    artifact = json.loads(ARTIFACT.read_text())
    assert artifact["summary"]["new_abstract_chamber_count"] == 0
    for case in artifact["cases"]:
        lineage = case["lineage"]
        assert lineage["abstract_chamber_status"] == (
            "known_exact_CoxIter_relation_graph_match"
        )
        assert len(lineage["coxiter"]["ordinary_exact_matches"]) == 2
        assert lineage["coxiter"]["doubling_exact_matches"] == []
        assert lineage["maximal_reflection_index_status"] == "exact_index_one"
        assert lineage["full_reflection_chamber"]["index"] == 1
        assert lineage["full_reflection_chamber"][
            "wall_by_wall_signed_Gram_and_vector_conjugacy"
        ] is True
        assert lineage["surface_automorphism_required_for_polytope_proof"] is False
        assert lineage["surface_automorphism_realization"].startswith(
            "proved_by_Roulleau_C2_nontrivial_NS_action"
        )

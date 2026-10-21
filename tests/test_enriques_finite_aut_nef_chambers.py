from __future__ import annotations

from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path

import sympy as sp

from scripts.certify_enriques_finite_aut_nef_chambers import (
    OUTPUT,
    RAW_DATA,
    RAW_DATA_SHA256,
    RAY_OUTPUT,
    certify,
)


ROOT = Path(__file__).resolve().parents[1]
SUMMARY_SHA256 = "245715b7352219a7ed1bfc34a43b5b37921e289d39829df51015c347b97072bc"
RAYS_SHA256 = "b72a2a5479a773c8151ed92cd8b698e71fe4b8cb0b4bf527245847df39cb46c5"
EXPECTED = {
    "I": (12, 14, 9, 114),
    "II": (12, 38, 17, 270),
    "III": (20, 200, 44, 1301),
    "IV": (20, 720, 87, 4405),
    "V": (20, 192, 29, 1182),
    "VI": (20, 600, 57, 3420),
    "VII": (20, 450, 51, 2490),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_archived_seven_image_certificate() -> None:
    assert _sha256(RAW_DATA) == RAW_DATA_SHA256
    assert _sha256(OUTPUT) == SUMMARY_SHA256
    assert _sha256(RAY_OUTPUT) == RAYS_SHA256

    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert document["schema"] == "enriques_finite_automorphism_nef_chambers_v1"
    assert document["surface_to_polytope_map"]["image_cardinality"] == 7
    assert document["pairwise_distinction"]["all_seven_distinct"] is True
    assert (
        document["pairwise_distinction"][
            "20C_and_20D_configuration_matrices_literally_equal"
        ]
        is True
    )
    assert document["global_checks"] == {
        "all_seven_complete_nef_chambers": True,
        "all_seven_exact_cdd_cones_have_no_hyperideal_ray": True,
        "all_seven_exact_wall_lists_rebuilt_from_raw_vectors": True,
        "all_seven_finite_volume": True,
        "all_seven_literal_Vinberg_audits_pass": True,
        "all_seven_noncompact": True,
        "facet_count_multiset": [12, 12, 20, 20, 20, 20, 20],
        "no_full_nef_finite_Aut_Enriques_H9_14_member": True,
    }

    lattice_gram = sp.Matrix(document["lattice"]["gram"])
    seen_pair_signatures = set()
    for record in document["types"]:
        type_name = record["nikulin_kondo_type"]
        facets, finite, cusps, edges = EXPECTED[type_name]
        configuration = sp.Matrix(record["intersection_matrix"])
        walls = sp.Matrix(record["wall_vectors_in_L10_basis"])
        assert configuration == walls * lattice_gram * walls.T
        assert configuration.rows == configuration.cols == facets
        assert walls.rank() == 10
        assert set(configuration) <= {-2, 0, 1, 2}
        assert record["facet_count"] == facets
        assert record["normal_gram_inertia"] == {
            "positive": 9,
            "negative": 1,
            "zero": facets - 10,
        }
        cone = record["exact_cone"]
        vinberg = record["literal_vinberg"]
        assert (
            cone["finite_vertex_count"],
            cone["ideal_vertex_count"],
            vinberg["elliptic_rank_8_edge_count"],
        ) == (finite, cusps, edges)
        assert cone["hyperideal_vertex_count"] == 0
        assert cone["finite_volume"] is True
        assert cone["compact"] is False
        assert set(cone["all_facet_face_ranks"]) == {9}
        assert vinberg["decision"] == "finite_volume_noncompact"
        assert vinberg["bad_finite_volume_edge_count"] == 0
        assert vinberg["invalid_full_zero_neighborhood_count"] == 0
        assert (
            vinberg[
                "every_elliptic_edge_has_exactly_two_finite_or_ideal_endpoints"
            ]
            is True
        )
        pair_counts = record["coxeter_pair_counts"]
        signature = (
            facets,
            pair_counts["orthogonal_m2"],
            pair_counts["angle_pi_over_3_m3"],
            pair_counts["parallel_at_infinity"],
        )
        assert signature not in seen_pair_signatures
        seen_pair_signatures.add(signature)
    assert len(seen_pair_signatures) == 7


def test_archived_exact_ray_census() -> None:
    with gzip.open(RAY_OUTPUT, "rt", encoding="utf-8") as handle:
        lines = [json.loads(line) for line in handle]
    assert lines[0]["schema"] == "enriques_finite_aut_nef_chamber_rays_v1"
    records = lines[1:]
    histogram = Counter((record["type"], record["kind"]) for record in records)
    for type_name, (_, finite, cusps, _) in EXPECTED.items():
        assert histogram[(type_name, "finite")] == finite
        assert histogram[(type_name, "ideal")] == cusps
        assert histogram[(type_name, "hyperideal")] == 0
    assert len(records) == sum(finite + cusps for _, finite, cusps, _ in EXPECTED.values())


def test_full_raw_to_vinberg_and_cdd_replay() -> None:
    rebuilt, ray_bytes = certify()
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert rebuilt == archived
    assert ray_bytes == RAY_OUTPUT.read_bytes()

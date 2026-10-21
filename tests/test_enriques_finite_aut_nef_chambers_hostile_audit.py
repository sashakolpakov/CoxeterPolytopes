from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

from scripts.audit_enriques_finite_aut_nef_chambers_hostile import (
    OUTPUT,
    RAW_DATA,
    RAW_DATA_SHA256,
    certify,
)


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_enriques_finite_aut_nef_chambers_hostile.py"
ARCHIVE_SHA256 = "224ded8c22c223cf379461f54650d2f7fd560292e1de79bbefe62e5d51d55b1b"

EXPECTED = {
    # record, facets, wall index, (m2,m3,infinity),
    # (finite, ideal, hyperideal, elliptic rank-8 edges, affine audits)
    "I": ("12A", 12, 2, (53, 10, 3), (14, 9, 0, 114, 12)),
    "II": ("12B", 12, 2, (51, 15, 0), (38, 17, 0, 270, 20)),
    "III": ("20B", 20, 4, (138, 16, 36), (200, 44, 0, 1301, 84)),
    "IV": ("20F", 20, 2, (140, 40, 10), (720, 87, 0, 4405, 184)),
    "V": ("20A", 20, 2, (135, 24, 31), (192, 29, 0, 1182, 68)),
    "VI": ("20E", 20, 1, (135, 45, 10), (600, 57, 0, 3420, 124)),
    "VII": ("20C", 20, 2, (135, 30, 25), (450, 51, 0, 2490, 92)),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _document() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_hostile_archive_and_independence_boundary() -> None:
    assert _sha256(RAW_DATA) == RAW_DATA_SHA256
    assert _sha256(OUTPUT) == ARCHIVE_SHA256
    document = _document()
    assert document["schema"] == "enriques_finite_aut_nef_chambers_hostile_audit_v1"
    boundary = document["independence_boundary"]
    assert boundary["production_certifier_imported"] is False
    assert boundary["production_summary_used_as_input"] is False
    assert boundary["production_ray_archive_used_as_input"] is False

    parsed = ast.parse(SCRIPT.read_text(encoding="utf-8"))
    imported_modules = {
        node.module
        for node in ast.walk(parsed)
        if isinstance(node, ast.ImportFrom) and node.module is not None
    }
    assert "scripts.certify_enriques_finite_aut_nef_chambers" not in imported_modules


def test_exact_seven_type_census_and_all_rays() -> None:
    document = _document()
    assert len(document["types"]) == 7
    seen_signatures = set()
    for record in document["types"]:
        type_name = record["nikulin_kondo_type"]
        raw_name, facets, wall_index, pair_expected, volume_expected = EXPECTED[type_name]
        assert record["raw_shimada_record"] == raw_name
        assert record["facet_count"] == facets
        assert record["raw_reconstruction"] == {
            **record["raw_reconstruction"],
            "configuration_equals_W_GramL10_Wt": True,
            "wall_rank": 10,
            "all_walls_primitive": True,
            "all_minus_two_reflections_integral": True,
            "wall_span_index_in_L10": wall_index,
        }
        assert record["gram_inertia"] == {
            "positive": 9,
            "negative": 1,
            "zero": facets - 10,
        }

        pair_counts = record["coxeter_pair_histogram"]
        pair_actual = (
            pair_counts["orthogonal_m2"],
            pair_counts["angle_pi_over_3_m3"],
            pair_counts["parallel_infinity"],
        )
        assert pair_actual == pair_expected
        assert sum(pair_actual) == facets * (facets - 1) // 2
        signature = (facets, *pair_actual)
        assert signature not in seen_signatures
        seen_signatures.add(signature)

        finite, ideal, hyperideal, edges, affine = volume_expected
        cone = record["exact_cone"]
        vinberg = record["literal_vinberg"]
        assert (
            cone["finite_vertex_count"],
            cone["ideal_vertex_count"],
            cone["hyperideal_ray_count"],
            vinberg["elliptic_rank_8_edge_count"],
            vinberg["affine_zero_neighborhood_audit_count"],
        ) == (finite, ideal, hyperideal, edges, affine)
        assert cone["extreme_ray_count"] == finite + ideal
        assert len(cone["rays"]) == finite + ideal
        assert {ray["kind"] for ray in cone["rays"]} <= {"finite", "ideal"}
        assert all(len(ray["active_facets"]) >= 9 for ray in cone["rays"])
        assert set(cone["facet_face_ranks"]) == {9}
        assert cone["all_facets_irredundant"] is True
        assert cone["all_rays_in_one_positive_time_sheet"] is True
        assert cone["finite_volume_from_exact_rays"] is True
        assert cone["compact_from_exact_rays"] is False
        assert vinberg["decision"] == "finite_volume_noncompact"
        assert vinberg["bad_finite_volume_edge_count"] == 0
        assert vinberg["invalid_affine_zero_neighborhood_count"] == 0
        assert sum(vinberg["edge_endpoint_histogram"].values()) == edges

    assert len(seen_signatures) == 7
    cross = document["cross_type_checks"]
    assert cross["20C_configuration_equals_20D_literally"] is True
    assert cross["20C_wall_vectors_equal_20D_literally"] is False
    assert cross["20D_duplicate_raw_replay"] == {
        "configuration_equals_W_GramL10_Wt": True,
        "wall_rank": 10,
        "all_walls_primitive": True,
        "normal_gram_has_H9_inertia": True,
        "wall_span_smith_diagonal": [1, 1, 1, 1, 1, 1, 1, 1, 1, 2],
        "wall_span_index_in_L10": 2,
    }
    assert cross["all_seven_pair_histograms_distinct"] is True
    assert cross["facet_count_multiset"] == [12, 12, 20, 20, 20, 20, 20]
    assert cross["no_fourteen_facet_image"] is True


def test_full_raw_hostile_replay() -> None:
    assert certify() == _document()

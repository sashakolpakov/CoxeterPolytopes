from __future__ import annotations

from collections import Counter
import gzip
import json

from scripts.certify_alexeev_20_2_1_15_9_1_exact_cones import (
    OUTPUT as CONE_OUTPUT,
    ROOT,
)
from scripts.certify_alexeev_20_2_1_15_9_1_literal_vinberg import (
    OUTPUT as VINBERG_OUTPUT,
)
from scripts.certify_alexeev_20_2_1_15_9_1_lineage_bounds import (
    OUTPUT as LINEAGE_OUTPUT,
)
from scripts.inventory_alexeev_20_2_1_15_9_1_wall_symmetries import (
    OUTPUT as INVENTORY_OUTPUT,
)
from scripts.screen_alexeev_20_2_1_15_9_1_fixed_sections import (
    OUTPUT as SECTION_OUTPUT,
)


def _load(path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_exact_weighted_groups_and_complete_high_rank_subgroup_queue() -> None:
    payload = _load(INVENTORY_OUTPUT)
    assert payload["source"]["sha256"] == (
        "4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457"
    )
    by_key = {tuple(row["key"]): row for row in payload["examples"]}
    assert set(by_key) == {(20, 2, 1), (15, 9, 1)}

    expected = {
        (20, 2, 1): {
            "rank": 20,
            "dimension": 19,
            "determinant": -4,
            "walls": 50,
            "name": "S5",
            "group_order": 120,
            "subgroup_classes": 19,
            "high": [(1, 2, 10, 13, 12)],
        },
        (15, 9, 1): {
            "rank": 15,
            "dimension": 14,
            "determinant": 512,
            "walls": 98,
            "name": "PGL(2,7)",
            "group_order": 336,
            "subgroup_classes": 23,
            "high": [(1, 2, 21, 11, 10)],
        },
    }
    # The (20,2,1) stock has two inequivalent involution classes; spell out
    # the second after the shared first tuple above.
    expected[(20, 2, 1)]["high"] = [
        (1, 2, 10, 13, 12),
        (2, 2, 15, 12, 11),
    ]
    for key, row in by_key.items():
        want = expected[key]
        assert row["ambient"]["rank"] == want["rank"]
        assert row["ambient"]["hyperbolic_dimension"] == want["dimension"]
        assert row["ambient"]["determinant"] == want["determinant"]
        assert row["walls"]["count"] == want["walls"]
        group = row["weighted_diagram_group"]
        assert group["name_from_Alexeev"] == want["name"]
        assert group["order"] == want["group_order"]
        assert group["all_actions_integral_unimodular_Lorentz"] is True
        assert row["summary"]["subgroup_conjugacy_class_count"] == want[
            "subgroup_classes"
        ]
        high = [
            (
                subgroup["class_index"],
                subgroup["order"],
                subgroup["conjugacy_class_size"],
                subgroup["fixed_rank"],
                subgroup["hyperbolic_dimension"],
            )
            for subgroup in row["subgroup_classes"]
            if subgroup["proper_fixed_section_above_h9"]
        ]
        assert high == want["high"]


def test_all_three_sections_are_strict_integral_coxeter_chambers() -> None:
    payload = _load(SECTION_OUTPUT)
    by_id = {row["id"]: row for row in payload["cases"]}
    expected = {
        "A20_2_1_C01": {
            "rank": 13,
            "dimension": 12,
            "determinant": 256,
            "smith": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 4],
            "orbits": 30,
            "facets": 19,
            "redundant": 11,
            "edges": {"2": 144, "3": 15, "4": 6, "infinity": 6},
        },
        "A20_2_1_C02": {
            "rank": 12,
            "dimension": 11,
            "determinant": -1024,
            "smith": [1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            "orbits": 30,
            "facets": 20,
            "redundant": 10,
            "edges": {
                "2": 149,
                "3": 12,
                "4": 11,
                "infinity": 8,
                "ultraparallel": 10,
            },
        },
        "A15_9_1_C01": {
            "rank": 11,
            "dimension": 10,
            "determinant": 2048,
            "smith": [1, 1, 2, 2, 2, 2, 2, 2, 2, 4, 4],
            "orbits": 62,
            "facets": 34,
            "redundant": 28,
            "edges": {
                "2": 331,
                "3": 17,
                "4": 30,
                "infinity": 75,
                "ultraparallel": 108,
            },
        },
    }
    assert set(by_id) == set(expected)
    for case_id, row in by_id.items():
        want = expected[case_id]
        lattice = row["saturated_fixed_lattice"]
        reduction = row["exact_reduction"]
        assert row["case"]["fixed_rank"] == want["rank"]
        assert row["case"]["hyperbolic_dimension"] == want["dimension"]
        assert lattice["inertia"] == [1, want["dimension"], 0]
        assert lattice["determinant"] == want["determinant"]
        assert lattice["smith_diagonal_absolute"] == want["smith"]
        assert row["case"]["wall_orbit_count"] == want["orbits"]
        assert reduction["active_hyperbolic_facet_count"] == want["facets"]
        assert reduction["conic_redundant_count"] == want["redundant"]
        assert reduction["unresolved_orbit_indices"] == []
        assert len(reduction["strict_facet_witnesses"]) == want["facets"]
        assert row["coxeter_gate"]["coxeter_edge_counts"] == want["edges"]
        assert row["coxeter_gate"]["all_pairs_coxeter"] is True
        assert row["restricted_mirrors"][
            "all_active_reflections_integral_on_saturated_fixed_lattice"
        ] is True
        assert row["geometry_scope"][
            "geometric_automorphism_membership_certified"
        ] is False


def test_gmp_and_literal_vinberg_volume_certificates_agree() -> None:
    cones = {row["id"]: row for row in _load(CONE_OUTPUT)["cases"]}
    vinberg = {row["id"]: row for row in _load(VINBERG_OUTPUT)["cases"]}
    expected = {
        "A20_2_1_C01": (12, 19, 552, 527, 25, 3846, 59),
        "A20_2_1_C02": (11, 20, 340, 311, 29, 2229, 64),
        "A15_9_1_C01": (10, 34, 1920, 1774, 146, 12792, 426),
    }
    assert set(cones) == set(vinberg) == set(expected)
    for case_id, want in expected.items():
        dimension, facets, rays, finite, ideal, edges, affine = want
        cone = cones[case_id]
        literal = vinberg[case_id]
        exact = cone["exact_double_description"]
        endpoint = literal["literal_endpoint_audit"]
        diagram = literal["diagram_enumeration"]
        assert (cone["hyperbolic_dimension"], cone["facet_count"]) == (
            dimension,
            facets,
        )
        assert exact["extreme_ray_count"] == rays
        assert exact["finite_vertex_count"] == finite
        assert exact["ideal_vertex_count"] == ideal
        assert exact["hyperideal_vertex_count"] == 0
        assert endpoint["finite_vertex_count"] == finite
        assert endpoint["ideal_vertex_count"] == ideal
        assert endpoint["edge_count"] == edges
        assert endpoint["edge_endpoint_count_histogram"] == {"2": edges}
        assert endpoint["bad_edge_count"] == 0
        assert endpoint["bad_parabolic_neighborhood_count"] == 0
        assert diagram["elliptic_rank_n_minus_1_edge_diagram_count"] == edges
        assert diagram["elliptic_rank_n_finite_endpoint_diagram_count"] == finite
        assert diagram["connected_affine_component_count"] == affine
        assert cone["decision"]["classification"] == "finite_volume_noncompact"
        assert literal["decision"]["classification"] == "finite_volume_noncompact"

        ray_path = ROOT / cone["artifacts"]["ray_archive_path"]
        with gzip.open(ray_path, "rt", encoding="utf-8") as stream:
            ray_header = json.loads(next(stream))
            ray_records = tuple(json.loads(line) for line in stream)
        assert ray_header["ray_count"] == len(ray_records) == rays
        assert Counter(record[3] for record in ray_records) == {
            "finite": finite,
            "ideal": ideal,
        }

        edge_path = ROOT / literal["edge_archive_path"]
        with gzip.open(edge_path, "rt", encoding="utf-8") as stream:
            edge_header = json.loads(next(stream))
            edge_records = tuple(json.loads(line) for line in stream)
        assert edge_header["edge_count"] == len(edge_records) == edges
        assert all(len(record[1]) == 2 for record in edge_records)


def test_bounded_lineage_keeps_literal_and_arithmetic_novelty_separate() -> None:
    payload = _load(LINEAGE_OUTPUT)
    assert payload["all_checks_pass"] is True
    coxiter = payload["coxiter_bounded_comparison"]
    assert coxiter["ordinary_coxiter_file_count"] == 1214
    assert coxiter["doubling_coxiter_file_count"] == 2959
    assert set(coxiter["targets"]) == {
        "H10_facets_34",
        "H11_facets_20",
        "H12_facets_19",
    }
    assert all(
        row["absent_from_both_snapshots_by_dimension_and_facet_count"]
        for row in coxiter["targets"].values()
    )

    by_id = {row["id"]: row for row in payload["cases"]}
    assert all(
        row["rational_similarity_comparisons"]["f1_odd_unimodular"][
            "rationally_similar"
        ]
        for row in by_id.values()
    )
    h11 = by_id["A20_2_1_C02"]
    assert h11["two_elementary_delta"] == 1
    assert h11["integral_classification"] == {
        "identification_method": (
            "uniqueness of indefinite even 2-elementary lattices from "
            "signature and (r,a,delta)"
        ),
        "index_value": None,
        "isometry_type": "U + E8(-2) + <-2>^2",
        "known_complete_reflection_chamber": {
            "facet_count": 14,
            "source": "Alexeev arXiv:2209.09110v4, case (12,10,1)",
        },
        "nikulin_invariants": [12, 10, 1],
        "section_group_has_finite_index_in_known_full_reflection_group": True,
    }

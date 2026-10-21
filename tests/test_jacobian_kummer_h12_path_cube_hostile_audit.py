from __future__ import annotations

import ast
from functools import lru_cache
import json

from scripts.audit_jacobian_kummer_h12_path_cube_hostile import (
    CONFIGURATIONS,
    OUTPUT,
    ROOT,
    build,
)


@lru_cache(maxsize=1)
def _payload() -> dict:
    return build()


def test_hostile_rebuild_matches_frozen_archive() -> None:
    payload = _payload()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    independence = payload["independence"]
    assert independence["project_module_imported"] is False
    assert independence["exploratory_path_cube_scout_imported_or_invoked"] is False
    assert independence["production_path_or_cube_builder_imported_or_invoked"] is False
    assert independence["production_path_or_cube_JSON_used_as_input"] is False
    assert independence["finite_N_interpolation_used_for_universal_formula"] is False


def test_exact_doubling_stock_and_all_eight_path_cubes() -> None:
    payload = _payload()
    stock = payload["doubling_stock"]
    assert stock["doubling_wall_count"] == 124
    assert stock["non_doubling_wall_count"] == 48
    assert len(stock["doubling_walls_one_based"]) == 124
    rows = payload["listed_configurations"]
    assert len(rows) == len(CONFIGURATIONS) == 8
    expected_gram = [
        [-4, 4, 0, 0, 0],
        [4, -4, 0, 0, 0],
        [0, 0, -4, 0, 0],
        [0, 0, 0, -4, 0],
        [0, 0, 0, 0, -4],
    ]
    for row, (pair, cube) in zip(rows, CONFIGURATIONS, strict=True):
        assert row["path_pair_one_based"] == [index + 1 for index in pair]
        assert row["cube_walls_one_based"] == [index + 1 for index in cube]
        assert row["selected_wall_Gram"] == expected_gram
        assert row["cube_order"] == 8
        assert row["cube_first_facet_count"] == 941
        assert row["side_facet_category_counts"] == {
            "orthogonal_to_both": 17,
            "orthogonal_to_left_only": 36,
            "orthogonal_to_neither": 850,
            "orthogonal_to_right_only": 36,
        }
        assert row["all_exact_checks_pass"]
        assert all(row["checks"].values())


def test_symbolic_all_z_formula_conjugacy_and_direct_replay() -> None:
    payload = _payload()
    proof = payload["representative_exact_proof"]
    witness = proof["strict_cube_facet_witnesses"]
    collision = proof["all_integer_collision_audit"]
    assert witness["facet_count"] == 941
    assert witness["minimum_off_target_slack"] == 1
    assert witness["minimum_witness_square"] == "17/2"
    assert collision["bounded_exponent_scan_used"] is False
    assert collision["ordered_signed_side_equation_count"] == 3_526_884
    assert collision["finite_solution_count"] == 994
    assert collision["finite_solution_histogram"] == {
        "T^z*S_right;z=-1;sign=1": 36,
        "T^z*S_right;z=0;sign=1": 36,
        "T^z;z=0;sign=1": 922,
    }
    assert collision["all_z_identity_count"] == 34
    assert collision["all_z_identity_histogram"] == {
        "T^z": 17,
        "T^z*S_right": 17,
    }
    assert collision["negative_root_collision_count"] == 0
    assert collision["cross_facet_collision_count"] == 0
    assert collision["endpoint_to_side_orbit_collision_count"] == 0

    conjugacy = payload["centralizer_conjugacy"]
    assert conjugacy["published_D16_stabilizer_order"] == 23_040
    assert conjugacy["trace_nine_involution_centralizer_order"] == 768
    assert conjugacy["configuration_orbit_size"] == 8
    assert conjugacy["configuration_stabilizer_order"] == 96
    assert conjugacy["phase_swap_within_configuration_1"][
        "centralizer_element_count"
    ] == 48
    assert all(conjugacy["checks"].values())

    family = payload["universal_family"]
    assert family["facet_formula"] == "f(N)=886*N+55 for every integer N>=1"
    assert family["reflection_subgroup_index"] == "8*N"
    assert family["volume"] == "Vol(P_N)=8*N*Vol(P_seed)"
    assert family["finite_volume"] and not family["compact"]
    assert family["Smith_diagonal_cube_boundary_without_left_endpoint"] == [1] * 13
    assert family["Smith_diagonal_cube_boundary_without_right_endpoint"] == [1] * 13
    assert family["full_fixed_Picard_lattice_for_every_N"] is True
    assert family[
        "one_isometry_class_among_the_eight_listed_configurations_for_each_N"
    ] is True

    direct = payload["direct_reconstruction_N_1_through_8"]
    assert len(direct["configurations"]) == 8
    for configuration in direct["configurations"]:
        assert len(configuration["phase_reconstructions"]) == 2
        for phase in configuration["phase_reconstructions"]:
            assert [row["facet_count"] for row in phase["checks_N_1_through_8"]] == [
                886 * n + 55 for n in range(1, 9)
            ]


def test_volume_carrier_lineage_and_claim_boundaries() -> None:
    payload = _payload()
    carrier = payload["universal_family"]["carrier"]
    assert carrier["seed_classification"] == "finite_volume_noncompact"
    assert carrier["carrier_lattice"] == "U(2)+D5(-1)^2+<-4>"
    assert carrier["carrier_rank"] == 13
    assert carrier["carrier_signature"] == [1, 12, 0]
    assert carrier["carrier_moduli_dimension"] == 7
    assert all(carrier["checks"].values())

    lineage = payload["bounded_coxiter_corpus_replay"]
    assert lineage["ordinary_file_count"] == 1214
    assert lineage["ordinary_H12_facet_histogram"] == {
        "13": 6,
        "14": 3,
        "16": 1,
        "18": 1,
    }
    assert lineage["index_two_file_count"] == 2959
    assert lineage["index_two_H12_facet_histogram"] == {
        "18": 1,
        "20": 1,
        "25": 1,
        "35": 1,
    }
    assert lineage["literal_family_headers_absent"] is True
    assert lineage["worldwide_priority"] == "open"
    assert "not a publication-complete denominator" in lineage["scope_warning"]

    boundary = payload["claim_boundary"]
    assert boundary["abstract_hyperbolic_Coxeter_polytopes"] is True
    assert boundary["finite_volume_noncompact"] is True
    assert boundary["weak_full_Picard_K3_divisor_reflection_carrier"] is True
    assert boundary["full_nef_realization_claimed"] is False
    assert boundary["surface_automorphism_realization_claimed"] is False
    assert boundary["eight_listed_configurations_exhaust_all_possible_path_cubes"] is False
    assert boundary["worldwide_priority"] == "open"


def test_hostile_script_has_no_project_import() -> None:
    script = ROOT / "scripts/audit_jacobian_kummer_h12_path_cube_hostile.py"
    tree = ast.parse(script.read_text(encoding="utf-8"))
    imported: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.append(node.module)
    assert not any(name == "scripts" or name.startswith("scripts.") for name in imported)
    assert not any(
        name == "coxeter_polytopes" or name.startswith("coxeter_polytopes.")
        for name in imported
    )

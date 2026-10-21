from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
import subprocess
import sys

from scripts.certify_coble_type_i_two_boundary_plane_inverse import (
    OUTPUT,
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]


@lru_cache(maxsize=1)
def fresh_certificate() -> dict:
    return build_certificate()


def test_frozen_certificate_is_exact_replay() -> None:
    frozen = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert fresh_certificate() == frozen
    assert frozen["all_checks_pass"] is True


def test_raw_endpoint_denominator_and_intrinsic_branch_count() -> None:
    certificate = fresh_certificate()
    raw = certificate["raw_extraction"]
    assert raw["all_sequential_P2_endpoint_count"] == 52
    assert raw["endpoint_degree_signature_count"] == 9
    assert raw["selected_endpoint_count"] == 4
    assert all(branch["proximity_chain_lengths"] == [2, 3, 3, 3] for branch in raw["branches"])
    assert all(branch["degree_zero_proximity_roots"] == 7 for branch in raw["branches"])
    assert all(branch["degree_one_three_center_root_lines"] == 3 for branch in raw["branches"])
    assert raw["invariant_canonical_program_profile_count"] == 1
    programs = [branch["raw_derived_canonical_program"] for branch in raw["branches"]]
    assert all(program["coefficientwise_boundary_match"] for program in programs)
    assert all(program["coefficientwise_all_ten_actual_roots_match"] for program in programs)
    assert all(program["coefficientwise_two_crossings_match"] for program in programs)
    assert all(program["coefficientwise_two_visible_cusp_lifts_match"] for program in programs)
    assert len(
        {
            json.dumps(
                {
                    "boundaries": program["canonical_boundary_classes"],
                    "roots": program["canonical_actual_root_classes"],
                    "deep": program["canonical_deep_exceptional_classes"],
                },
                sort_keys=True,
            )
            for program in programs
        }
    ) == 1


def test_determinantal_locus_and_rejected_component() -> None:
    geometry = fresh_certificate()["canonical_determinantal_geometry"]
    matrix = geometry["quintic_matrix"]
    assert (matrix["rows"], matrix["columns"]) == (23, 21)
    assert matrix["maximal_minor"] == "-64*(a - b)*(a + b)"
    assert matrix["exact_rank_drop_locus"] == "a^2-b^2=0"
    assert matrix["rank_off_determinantal_locus"] == 21
    assert matrix["rank_on_a_equals_b"] == 20
    assert matrix["rank_on_a_equals_minus_b"] == 20
    assert "**2" in geometry["rejected_component"]["kernel_factorization"]
    assert "nonreduced" in geometry["rejected_component"]["reason"]


def test_good_quintic_and_exact_return_scope() -> None:
    certificate = fresh_certificate()
    good = certificate["canonical_determinantal_geometry"]["surviving_horizontal_component"]
    assert good["singular_audit"]["absolutely_irreducible"] is True
    assert good["singular_audit"]["singular_radical_degree"] == 2
    assert [record["analytic_type_in_characteristic_zero"] for record in good["local_resolutions"]] == ["A5", "A5"]
    assert good["geometric_genus"] == 0
    returned = certificate["strong_Coble_and_exact_return"]
    assert returned["boundary_classes"]["intersection_matrix"] == [[-4, 0], [0, -4]]
    replay = returned["exact_return_transfer_replay"]
    assert replay["all_transfer_gates_pass"] is True
    assert replay["finite_index_CM_root_gate"]["passes"] is True
    assert replay["complete_cusp_lift_gate"]["passes"] is True
    assert replay["four_raw_plane_branches_pass"] is True
    explicit = replay["explicit_ten_missing_curve_gate"]
    assert explicit["passes"] is True
    assert explicit["ten_positive_degree_curve_count"] == 10
    assert explicit["degree_multiset"] == [2, 2, 2, 2, 3, 3, 4, 4, 5, 5]
    assert explicit["all_ten_absolutely_irreducible"] is True
    assert explicit["all_ten_exact_strict_transforms_smooth_rational"] is True
    assert replay["uniform_hidden_wall_theorem_gate"]["passes"] is True
    assert replay["uniform_hidden_wall_theorem_gate"]["not_a_phrase_gate"] is True
    assert replay["literal_full_NS_cone_gate"]["passes"] is True
    assert all(
        record["explicit_minus_one_wall_count"] == 4
        and record["remaining_count"] == 10
        and record["remaining_is_two_opposite_partners_plus_four_complete_pairs"]
        for record in replay["four_raw_plane_branch_wall_accounting"]
    )
    assert returned["full_nef_polytope"] == {
        "compactness": "noncompact",
        "dimension": 11,
        "facets": 26,
        "finite_vertices": 508,
        "hyperideal_rays": 0,
        "ideal_vertices": 42,
        "volume": "finite",
    }


def test_cli_check_output() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/certify_coble_type_i_two_boundary_plane_inverse.py"),
            "--check-output",
            str(OUTPUT),
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr

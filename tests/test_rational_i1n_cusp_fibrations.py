from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_rational_i1n_cusp_fibrations import (
    A1_WORKER_OUTPUT,
    OUTPUT,
    _a0_cluster_jet_chart,
    _build_a1_worker_certificate,
    build_certificate,
    validate,
)
from scripts.certify_rational_i1n_inverse_design import (
    enumerate_simple_roots,
    lorentz_pair,
)


def _archive() -> dict[str, object]:
    return json.loads(Path(OUTPUT).read_text(encoding="utf-8"))


def test_frozen_cusp_fibration_certificate() -> None:
    archive = _archive()
    validate(archive)
    assert build_certificate() == archive


def test_frozen_a1_worker_certificate_replays_in_fresh_processes() -> None:
    archive = json.loads(Path(A1_WORKER_OUTPUT).read_text(encoding="utf-8"))
    assert _build_a1_worker_certificate() == archive


def test_every_affine_cycle_reconstructs_its_primitive_cusp() -> None:
    archive = _archive()
    for dimension, record in archive["dimensions"].items():
        roots = enumerate_simple_roots(int(dimension))
        for cusp in record["cusps"]:
            vector = tuple(cusp["primitive_nef_vector"])
            assert lorentz_pair(vector, vector) == 0
            assert cusp["active_facet_count"] == (
                int(dimension) + cusp["affine_component_count"] - 1
            )
            for component in cusp["components"]:
                labels = [label - 1 for label in component["facet_labels_one_based"]]
                weights = component["primitive_null_weights"]
                reconstructed = tuple(
                    sum(weights[position] * roots[label][coordinate]
                        for position, label in enumerate(labels))
                    for coordinate in range(int(dimension) + 1)
                )
                assert reconstructed == vector


def test_semiample_claim_partition_is_exact() -> None:
    archive = _archive()
    semiample = []
    unresolved = []
    for record in archive["dimensions"].values():
        for cusp in record["cusps"]:
            (semiample if cusp["semiample_certified"] else unresolved).append(cusp)
    assert len(semiample) == 28
    assert unresolved == []
    assert all(cusp["semiample_reason"] for cusp in semiample)
    index_four = next(
        cusp for cusp in semiample
        if cusp.get("semiample_multiple_system") is not None
    )
    assert index_four["arithmetic_genus"] == 3
    rigid = index_four["primitive_linear_system"]
    assert rigid["complete_linear_system_h0"] == 1
    assert rigid["unique_section_is_affine_wall_cycle"] is True
    assert "audited separately" in rigid["semiampleness_boundary"]
    assert index_four["semiample_multiple_system"]["semiample_index"] == 4
    assert "exact first moving multiple is |4F|" in index_four["semiample_reason"]


def test_exact_a0_cluster_chart_is_derived_by_strict_transforms() -> None:
    chart = _a0_cluster_jet_chart()
    assert chart["Q_selected_branch_coefficients_through_x15"] == {
        "3": "1",
        "9": "1",
        "11": "-1",
        "13": "1",
        "15": "2",
    }
    assert chart["S_selected_branch_coefficients_through_x15"] == {
        "3": "1",
        "9": "1",
        "11": "-1",
        "13": "1",
        "15": "1",
    }
    assert chart["post_p15_exceptional_restrictions"] == {
        "Q": "z",
        "S": "z+1",
    }
    replay = chart["blowup_replay"]
    assert [row["Q_multiplicity"] for row in replay] == [2] + [1] * 14
    assert [row["S_multiplicity"] for row in replay] == [2] * 7 + [1] * 8


def test_five_formerly_open_cusps_have_exact_coprime_second_fibres() -> None:
    archive = _archive()
    records = {}
    for dimension, row in archive["dimensions"].items():
        for cusp in row["cusps"]:
            jet = cusp.get("primitive_linear_system")
            if jet is not None:
                records[(int(dimension), cusp["arithmetic_genus"])] = jet

    expected = {
        (14, 5): (78, 627, 76, 31,
                  "03752569fc28f9fa673f4d44f36baaa616de4f2ded403c9a577a56c832c928b3"),
        (15, 5): (78, 708, 76, 31,
                  "03752569fc28f9fa673f4d44f36baaa616de4f2ded403c9a577a56c832c928b3"),
        (16, 5): (78, 786, 76, 31,
                  "03752569fc28f9fa673f4d44f36baaa616de4f2ded403c9a577a56c832c928b3"),
        (16, 10): (253, 2655, 251, 91,
                   "8e8e7208cb2b46b9d929d528a6a1b55f75c9cc2566fadef91f1e2310576df2e6"),
        (16, 11): (300, 3139, 298, 114,
                   "20110990490c2d0836421fadf56bba8e3f25032679f05bb751374aaec19c1b28"),
    }
    for key, wanted in expected.items():
        jet = records[key]
        actual = (
            jet["homogeneous_monomial_count"],
            jet["forbidden_jet_condition_count"],
            jet["exact_condition_matrix_rank"],
            jet["second_plane_section_term_count"],
            jet["second_plane_section_sparse_sha256"],
        )
        assert actual == wanted
        assert jet["complete_linear_system_h0"] == 2
        assert jet["second_section_actual_multiplicities"] == (
            jet["assigned_multiplicities"]
        )
        assert jet["second_section_has_no_exceptional_component"] is True
        assert jet["second_section_coprime_to_wall_plane_product"] is True
        assert jet["base_point_free_pencil"] is True
        assert jet["semiample_certified"] is True


def test_the_sixth_primitive_jet_kernel_is_exactly_one_dimensional() -> None:
    archive = _archive()
    jet = next(
        cusp["primitive_linear_system"]
        for cusp in archive["dimensions"]["16"]["cusps"]
        if cusp["arithmetic_genus"] == 3
        and cusp.get("primitive_linear_system") is not None
    )
    assert (
        jet["homogeneous_monomial_count"],
        jet["forbidden_jet_condition_count"],
        jet["exact_condition_matrix_rank"],
        jet["complete_linear_system_h0"],
    ) == (15, 136, 14, 1)
    assert jet["wall_cycle_plane_factorization"] == "Q"
    assert jet["primitive_linear_system_rigid"] is True
    assert jet["semiample_certified"] is False


def test_genus_three_cusp_has_exact_semiample_index_four() -> None:
    archive = _archive()
    cusp = next(
        cusp
        for cusp in archive["dimensions"]["16"]["cusps"]
        if cusp.get("semiample_multiple_system") is not None
    )
    audit = cusp["semiample_multiple_system"]
    assert audit["semiample_index"] == 4
    assert audit["fourth_multiple_base_point_free"] is True
    assert [
        (
            row["homogeneous_monomial_count"],
            row["forbidden_jet_condition_count"],
            row["exact_condition_matrix_rank"],
            row["complete_linear_system_h0"],
        )
        for row in audit["multiples"]
    ] == [
        (15, 136, 14, 1),
        (45, 408, 44, 1),
        (91, 816, 90, 1),
        (153, 1360, 151, 2),
    ]
    fourth = audit["multiples"][3]
    assert fourth["second_plane_section_term_count"] == 64
    assert fourth["second_plane_section_sparse_sha256"] == (
        "1393b70fc23154a593c95781a8990caab364716c3cbcc308296f91883657ec12"
    )
    assert fourth["second_section_actual_multiplicities"] == [4] * 16
    assert fourth["second_section_has_no_exceptional_component"] is True
    assert fourth["second_section_coprime_to_Q4"] is True
    assert fourth["base_point_free_pencil"] is True
    assert fourth["fibre_arithmetic_genus"] == 9


def test_nonisomorphic_a1_carrier_has_the_same_28_semiample_cusps() -> None:
    archive = _archive()
    audit = archive["a1_nonisomorphic_carrier_primitive_jet_audit"]
    assert audit["nonisomorphic_to_a0_j0_carrier"] is True
    assert (
        audit["semiample_cusp_count"],
        audit["primitive_base_point_free_pencil_count"],
        audit["index_four_base_point_free_pencil_count"],
    ) == (28, 27, 1)
    assert audit["all_six_nullities_equal_a0"] is True
    assert audit["all_five_second_plane_sections_differ_from_a0"] is True
    assert audit["index_four_second_plane_section_differs_from_a0"] is True
    moving = [row for row in audit["records"] if row["semiample_certified"]]
    assert len(moving) == 5
    assert all(row["second_plane_section_differs_from_a0"] for row in moving)
    multiple = audit["genus_three_multiple_audit"]
    assert multiple["semiample_index"] == 4
    fourth = multiple["multiples"][3]
    assert fourth["second_plane_section_term_count"] == 63
    assert fourth["second_plane_section_sparse_sha256"] == (
        "f8b43b63e6e0c560303c398d983884f737623cc732b61ead6c2d604b18866b0b"
    )
    assert fourth["second_section_actual_multiplicities"] == [4] * 16
    assert fourth["second_section_coprime_to_Q4"] is True

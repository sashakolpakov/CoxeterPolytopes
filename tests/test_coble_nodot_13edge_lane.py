from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
import shutil

from scripts.certify_coble_nodot_13edge_lane import (
    OUTPUT,
    _analyze_cone,
    _even_unimodular_extensions,
    _full_cusp_lift,
    _literal_vinberg,
    _matrix,
    _quotient_lattice,
)


ARCHIVE_SHA256 = "103f4867ef89bca47f96d0ee806ae4c6fa157a62058294ed5852b112bfad6ffa"


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_lane_scope_and_exact_denominator_counts() -> None:
    assert sha256(OUTPUT.read_bytes()).hexdigest() == ARCHIVE_SHA256
    payload = _archive()
    assert payload["all_checks_pass"] is True
    assert payload["scope"] == {
        "connected_ordinary_support": True,
        "not_the_full_nodot_denominator": True,
        "ordinary_support_edge_count": 13,
        "ordinary_support_maximum_degree": 3,
        "pairing_alphabet": [0, 1, 2],
        "parallel_pair_counts": [2, 3],
        "type_II_excluded_by_scope_because_it_has_15_simple_edges": True,
        "saturated_radical_quotient_determinant": -4,
        "root_span_index_in_an_E10_extension": 2,
        "vertices": 12,
    }
    denominator = payload["scout_denominator"]
    assert denominator == {
        "underlying_unlabelled_graph_count": 2888,
        "exact_lorentzian_corank_two_count": 100,
        "determinant_four_labelled_count": 38,
        "determinant_four_weighted_isomorphism_class_count": 36,
    }
    classes = payload["class_records"]
    assert sum(row["labelled_member_count"] for row in classes) == 38
    assert Counter(row["parallel_pair_count"] for row in classes) == {2: 8, 3: 28}
    assert Counter(
        row["radical_quotient"]["even_unimodular_extension_count"]
        for row in classes
    ) == {1: 26, 2: 10}
    assert all(
        row["radical_quotient"]["determinant"] == -4
        and row["radical_quotient"]["smith_diagonal_absolute"]
        == [1] * 8 + [2, 2]
        for row in classes
    )


def test_only_type_i_survives_true_volume_and_full_cusp_gates() -> None:
    payload = _archive()
    branches = [
        (row, branch)
        for row in payload["class_records"]
        for branch in row["extension_branches"]
    ]
    assert len(branches) == 46
    assert sum(branch["cm_gate_pass"] for _, branch in branches) == 1
    assert sum(branch["full_D_1_split_gate_pass"] for _, branch in branches) == 1
    assert sum(
        branch["cm_cone"][
            "complete_cdd_V_representation_has_no_hyperideal_ray"
        ]
        for _, branch in branches
    ) == 1
    row, survivor = next(
        (row, branch) for row, branch in branches if branch["cm_gate_pass"]
    )
    assert row["class_index_zero_based"] == 11
    assert row["known_root_diagram_matches"] == ["Type I"]
    assert survivor["cm_cone"]["finite_vertex_count"] == 14
    assert survivor["cm_cone"]["ideal_vertex_count"] == 9
    assert survivor["cm_cone"]["hyperideal_vertex_count"] == 0
    assert survivor["cm_literal_vinberg"]["elliptic_edge_count"] == 114
    assert survivor["cm_literal_vinberg"]["connected_affine_component_count"] == 12
    assert survivor["cm_literal_vinberg"]["literal_vinberg_passes"] is True
    full = survivor["full_cusp_lift"]
    assert full["wall_count"] == 22
    assert full["characteristic_wall_positions_zero_based"] == [0]
    assert full["wall_coordinate_span_index_in_I_1_10"] == 1
    assert full["cone"]["finite_vertex_count"] == 315
    assert full["cone"]["ideal_vertex_count"] == 38
    assert full["cone"]["hyperideal_vertex_count"] == 0
    assert full["literal_vinberg"]["elliptic_edge_count"] == 2121
    assert full["literal_vinberg"]["connected_affine_component_count"] == 83
    assert full["literal_vinberg"]["literal_vinberg_passes"] is True
    assert payload["full_survivor_ledger"] == [
        {
            "compiler_scope": (
                "inherited exact raw_22 replay after literal full-wall-Gram "
                "isomorphism"
            ),
            "full_isomorphism_class_zero_based": 0,
            "full_split_K3_compiler_pass": True,
            "known_full_chamber_matches": ["Type I"],
            "source_branches": [[11, 0]],
        }
    ]


def test_unique_survivor_replays_from_its_raw_gram() -> None:
    payload = _archive()
    row = payload["class_records"][11]
    quotient = _quotient_lattice(-_matrix(row["negative_cartan"]))
    extensions = _even_unimodular_extensions(
        quotient["gram"], quotient["root_coordinates"]
    )
    assert len(extensions) == 1
    extension = extensions[0]
    cdd = shutil.which("cddexec_gmp") or "cddexec_gmp"
    cm = _analyze_cone(extension["gram"], extension["root_coordinates"], cdd)
    vinberg = _literal_vinberg(
        extension["gram"], extension["root_coordinates"], cm
    )
    assert (cm["finite_vertex_count"], cm["ideal_vertex_count"], cm["hyperideal_vertex_count"]) == (14, 9, 0)
    assert vinberg["literal_vinberg_passes"] is True
    full = _full_cusp_lift(
        extension["gram"], extension["root_coordinates"], cm, cdd
    )
    assert full["wall_coordinate_span_index_in_I_1_10"] == 1
    assert (
        full["cone"]["finite_vertex_count"],
        full["cone"]["ideal_vertex_count"],
        full["cone"]["hyperideal_vertex_count"],
    ) == (315, 38, 0)
    assert full["literal_vinberg"]["literal_vinberg_passes"] is True


def test_no_L4_claim_is_made() -> None:
    payload = _archive()
    assert payload["actual_full_D_1_split_branch_count"] == 1
    assert payload["actual_full_D_1_split_isomorphism_class_count"] == 1
    assert payload["claims"]["strict_L4_for_D_1_split_nodot"] is False
    assert payload["claims"]["full_split_compiler_only_inferred_after_literal_known_gram_match"] is True
    assert "proper slice" in payload["claims"]["reason"]

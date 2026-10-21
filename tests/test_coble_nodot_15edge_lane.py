from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
import shutil

from scripts.certify_coble_nodot_13edge_lane import (
    _all_even_unimodular_extensions,
    _analyze_cone,
    _full_cusp_lift,
    _literal_vinberg,
    _matrix,
    _quotient_lattice,
)
from scripts.certify_coble_nodot_15edge_lane import OUTPUT


ARCHIVE_SHA256 = "c9e2da6bcc08736fbc9848d64ca3fd273c23d535394ef727b671f369e7b20dae"


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_exact_surface_free_denominator() -> None:
    assert sha256(OUTPUT.read_bytes()).hexdigest() == ARCHIVE_SHA256
    payload = _archive()
    assert payload["all_checks_pass"] is True
    assert payload["scope"] == {
        "all_even_unimodular_E10_extensions_retained": True,
        "all_exact_Lorentzian_corank_two_determinants_retained": True,
        "connected_ordinary_support": True,
        "not_the_full_D_1_split_nodot_denominator": True,
        "ordinary_support_edge_count": 15,
        "ordinary_support_maximum_degree": 3,
        "pairing_alphabet": [0, 1],
        "parallel_pair_count": 0,
        "vertices": 12,
    }
    assert payload["denominator"] == {
        "exact_Lorentzian_corank_two_labelled_count": 2,
        "integral_E10_extension_branch_count": 3,
        "radical_quotient_determinant_histogram": {
            "15:0:-1": 1,
            "15:0:-4": 1,
        },
        "underlying_unlabelled_graph_count": 5741,
        "weighted_diagram_isometry_class_count": 2,
        "weighted_trial_count": 5741,
    }
    assert sorted(
        (
            row["radical_quotient"]["determinant"],
            row["radical_quotient"]["even_unimodular_E10_extension_count"],
        )
        for row in payload["class_records"]
    ) == [(-4, 2), (-1, 1)]


def test_exact_volume_and_type_ii_return_ledger() -> None:
    payload = _archive()
    branches = [
        (row, branch)
        for row in payload["class_records"]
        for branch in row["extension_branches"]
    ]
    assert Counter(branch["classification"] for _, branch in branches) == {
        "finite_volume_noncompact": 2,
        "infinite_volume_hyperideal": 1,
    }
    finite = [(row, branch) for row, branch in branches if branch["classification"] == "finite_volume_noncompact"]
    assert len(finite) == 2
    assert all(row["known_root_diagram_matches"] == ["Type II"] for row, _ in finite)
    assert all(
        (
            branch["cm_cone"]["finite_vertex_count"],
            branch["cm_cone"]["ideal_vertex_count"],
            branch["cm_cone"]["hyperideal_vertex_count"],
        )
        == (38, 17, 0)
        and branch["cm_literal_vinberg"]["elliptic_edge_count"] == 270
        and branch["cm_literal_vinberg"]["connected_affine_component_count"] == 20
        and branch["cm_literal_vinberg"]["literal_vinberg_passes"] is True
        for _, branch in finite
    )
    for _, branch in finite:
        full = branch["full_cusp_lift"]
        assert full["gate_pass"] is True
        assert full["known_full_chamber_matches"] == ["Type II"]
        assert full["wall_count"] == 30
        assert full["wall_coordinate_span_index_in_I_1_10"] == 1
        assert full["characteristic_wall_positions_zero_based"] == [0]
        assert (
            full["cone"]["finite_vertex_count"],
            full["cone"]["ideal_vertex_count"],
            full["cone"]["hyperideal_vertex_count"],
        ) == (1333, 207, 0)
        assert full["literal_vinberg"]["elliptic_edge_count"] == 10030
        assert full["literal_vinberg"]["connected_affine_component_count"] == 436
        assert full["literal_vinberg"]["literal_vinberg_passes"] is True
    rejected = next(branch for _, branch in branches if branch["classification"] == "infinite_volume_hyperideal")
    assert (
        rejected["cm_cone"]["finite_vertex_count"],
        rejected["cm_cone"]["ideal_vertex_count"],
        rejected["cm_cone"]["hyperideal_vertex_count"],
    ) == (18, 8, 24)
    assert rejected["cm_literal_vinberg"]["stopped_at_hyperideal_endpoint"] is True
    assert rejected["full_cusp_lift"] is None
    assert payload["full_D_1_split_branch_count"] == 2
    assert payload["full_D_1_split_isomorphism_class_count"] == 1
    assert all(
        row["known_full_chamber_matches"] == ["Type II"]
        and row["full_split_K3_Coble_return_pass"] is True
        for row in payload["full_survivor_ledger"]
    )


def test_three_integral_branches_replay_at_the_root_cone() -> None:
    payload = _archive()
    cdd = shutil.which("cddexec_gmp") or "cddexec_gmp"
    replayed = []
    for row in payload["class_records"]:
        quotient = _quotient_lattice(-_matrix(row["negative_cartan"]))
        extensions = _all_even_unimodular_extensions(
            quotient["gram"], quotient["root_coordinates"]
        )
        assert len(extensions) == row["radical_quotient"][
            "even_unimodular_E10_extension_count"
        ]
        for extension in extensions:
            cone = _analyze_cone(
                extension["gram"], extension["root_coordinates"], cdd
            )
            vinberg = _literal_vinberg(
                extension["gram"], extension["root_coordinates"], cone
            )
            replayed.append(
                (
                    cone["finite_vertex_count"],
                    cone["ideal_vertex_count"],
                    cone["hyperideal_vertex_count"],
                    vinberg["literal_vinberg_passes"],
                )
            )
    assert sorted(replayed) == [
        (18, 8, 24, False),
        (38, 17, 0, True),
        (38, 17, 0, True),
    ]


def test_scope_is_not_promoted_to_the_flagship_theorem() -> None:
    claims = _archive()["claims"]
    assert claims["strict_L4_for_full_D_1_split_nodot"] is False
    assert claims["full_split_compiler_only_inferred_after_literal_known_gram_match"] is True
    assert "only the 15-edge" in claims["reason"]

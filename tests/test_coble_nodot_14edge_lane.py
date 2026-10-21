from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
import math
import shutil

from scripts.certify_coble_nodot_13edge_lane import (
    _all_even_unimodular_extensions,
    _analyze_cone,
    _literal_vinberg,
    _matrix,
    _quotient_lattice,
)
from scripts.certify_coble_nodot_14edge_lane import OUTPUT


ARCHIVE_SHA256 = "e1a14bfed652de62979baf61afaabf262bd96482e0ec7e3376994267884a03aa"


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_frozen_denominator_and_scope() -> None:
    assert sha256(OUTPUT.read_bytes()).hexdigest() == ARCHIVE_SHA256
    payload = _archive()
    assert payload["all_checks_pass"] is True
    assert payload["scope"]["ordinary_support_edge_count"] == 14
    assert payload["scope"]["ordinary_support_maximum_degree"] == 3
    assert payload["scope"]["parallel_pair_counts"] == [0, 1, 2, 3]
    assert payload["scope"]["all_exact_Lorentzian_corank_two_determinants_retained"] is True
    assert payload["scope"]["all_even_unimodular_E10_extensions_retained"] is True
    assert payload["scope"]["not_the_full_D_1_split_nodot_denominator"] is True
    denominator = payload["denominator"]
    assert denominator["underlying_unlabelled_graph_count"] == 5291
    assert denominator["weighted_trial_count"] == 2_486_770
    assert denominator["exact_Lorentzian_corank_two_labelled_count"] == 236
    assert denominator["weighted_diagram_isometry_class_count"] == 194
    assert denominator["E10_extendable_weighted_class_count"] == 148
    assert denominator["integral_E10_extension_branch_count"] == 169
    assert sum(row["labelled_member_count"] for row in payload["class_records"]) == 236


def test_every_integral_branch_is_infinite_volume() -> None:
    payload = _archive()
    branches = [
        branch
        for row in payload["class_records"]
        for branch in row["extension_branches"]
    ]
    assert len(branches) == 169
    assert Counter(branch["classification"] for branch in branches) == {
        "infinite_volume_hyperideal": 169
    }
    assert min(branch["cone"]["hyperideal_vertex_count"] for branch in branches) == 5
    assert max(branch["cone"]["hyperideal_vertex_count"] for branch in branches) == 56
    assert all(
        branch["cone"]["hyperideal_vertex_count"] > 0
        and branch["hyperideal_witness"]["type"] == "hyperideal"
        and branch["hyperideal_witness"]["square"] < 0
        and math.gcd(
            *(abs(value) for value in branch["hyperideal_witness"]["coordinates"])
        )
        == 1
        and branch["finite_volume_gate_pass"] is False
        and branch["literal_vinberg"]["stopped_at_hyperideal_endpoint"] is True
        and branch["literal_vinberg"]["literal_vinberg_passes"] is False
        and branch["full_cusp_lift_attempted"] is False
        and branch["split_K3_compiler_attempted"] is False
        for branch in branches
    )
    ledger = payload["volume_ledger"]
    assert ledger["compact_branch_count"] == 0
    assert ledger["finite_volume_noncompact_branch_count"] == 0
    assert ledger["infinite_volume_hyperideal_branch_count"] == 169
    assert ledger["finite_volume_survivors"] == []


def test_every_square_determinant_E10_lane_replays() -> None:
    payload = _archive()
    expected_branches = {-1: 75, -4: 63, -9: 6, -16: 10, -36: 12, -64: 3}
    assert {
        int(key): value
        for key, value in payload["denominator"][
            "determinant_histogram_on_E10_branches"
        ].items()
    } == expected_branches
    representatives = {}
    for row in payload["class_records"]:
        determinant = row["radical_quotient"]["determinant"]
        if row["extension_branches"] and determinant not in representatives:
            representatives[determinant] = row
    assert set(representatives) == set(expected_branches)
    for determinant, row in representatives.items():
        quotient = _quotient_lattice(-_matrix(row["negative_cartan"]))
        extensions = _all_even_unimodular_extensions(
            quotient["gram"], quotient["root_coordinates"]
        )
        assert len(extensions) == row["radical_quotient"][
            "even_unimodular_E10_extension_count"
        ]
        for extension in extensions:
            gram = extension["gram"]
            assert int(gram.det()) == -1
            assert all(int(gram[index, index]) % 2 == 0 for index in range(10))
            assert extension["index_over_root_quotient"] ** 2 == abs(determinant)


def test_minimum_and_maximum_hyperideal_witness_branches_replay() -> None:
    payload = _archive()
    selected = {}
    for row in payload["class_records"]:
        for branch in row["extension_branches"]:
            count = branch["cone"]["hyperideal_vertex_count"]
            if count in {5, 56} and count not in selected:
                selected[count] = (row, branch)
    assert set(selected) == {5, 56}
    cdd = shutil.which("cddexec_gmp") or "cddexec_gmp"
    for expected, (row, archived_branch) in selected.items():
        quotient = _quotient_lattice(-_matrix(row["negative_cartan"]))
        extensions = _all_even_unimodular_extensions(
            quotient["gram"], quotient["root_coordinates"]
        )
        extension = extensions[archived_branch["extension_index_zero_based"]]
        cone = _analyze_cone(
            extension["gram"], extension["root_coordinates"], cdd
        )
        assert cone["hyperideal_vertex_count"] == expected
        vinberg = _literal_vinberg(
            extension["gram"], extension["root_coordinates"], cone
        )
        assert vinberg["stopped_at_hyperideal_endpoint"] is True
        assert vinberg["literal_vinberg_passes"] is False


def test_claims_are_strictly_scoped() -> None:
    claims = _archive()["claims"]
    assert claims == {
        "all_nodot": False,
        "new_hyperbolic_example": False,
        "scoped_negative_classification": True,
        "strict_L4": False,
    }

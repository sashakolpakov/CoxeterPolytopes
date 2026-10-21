import json
from collections import Counter

import pytest

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import SILVER_ZERO, SQRT2, Silver
from coxeter_polytopes.search.silver_k3_vinberg import (
    run_silver_k3_vinberg_four_plus_two_sqrt2_extension,
)
from scripts.certify_silver_k3_vinberg_4_plus_2sqrt2 import (
    ARCHIVE,
    EXPECTED_CERTIFICATE_SHA256,
    verify_manifest,
)
from scripts.independent_silver_k3_vinberg_4_plus_2sqrt2 import (
    EXPECTED_REPLAY_SHA256,
    audit_block_silver_four_plus_two_shell,
)


def test_four_plus_two_sqrt2_shell_is_exact_and_adds_121_walls() -> None:
    extension = run_silver_k3_vinberg_four_plus_two_sqrt2_extension()
    level = extension.new_level

    assert extension.passes
    assert extension.completed_through_control_coordinate == 4 + 2 * SQRT2
    assert extension.next_admissible_control_coordinate == 3 + 3 * SQRT2
    assert level.target == 58 + 40 * SQRT2
    assert level.alphabet.values == (
        SILVER_ZERO,
        Silver(1),
        SQRT2,
        1 + SQRT2,
        2 + SQRT2,
        3 + SQRT2,
        2 + 2 * SQRT2,
        3 + 2 * SQRT2,
        4 + 2 * SQRT2,
        3 + 3 * SQRT2,
        4 + 3 * SQRT2,
        5 + 3 * SQRT2,
        4 + 4 * SQRT2,
        5 + 4 * SQRT2,
        6 + 4 * SQRT2,
        7 + 4 * SQRT2,
        6 + 5 * SQRT2,
        7 + 5 * SQRT2,
        8 + 5 * SQRT2,
    )
    assert level.visited_by_depth == (
        1,
        19,
        361,
        1178,
        2261,
        3211,
        22363,
        73568,
        56012,
        21299,
    )
    assert level.schur_pruned_by_depth == (
        0,
        0,
        299,
        1059,
        2092,
        2034,
        18491,
        70620,
        54891,
        15783,
    )
    assert level.exact_weight_solutions == 653
    assert level.nonintegral_reconstructions == 330
    assert level.previous_wall_rejections == 202
    assert len(level.roots) == 121
    assert len(extension.roots) == 274
    assert extension.wall_identity_inertia == Inertia(9, 1, 264)
    assert extension.wall_conjugate_inertia == Inertia(10, 0, 264)


def test_all_121_same_height_roots_survive_simultaneously() -> None:
    extension = run_silver_k3_vinberg_four_plus_two_sqrt2_extension()
    cohort = tuple(range(extension.prior_wall_count, len(extension.roots)))
    pairings = tuple(
        extension.wall_gram[i][j]
        for i in cohort
        for j in cohort
        if j < i
    )

    assert cohort == tuple(range(153, 274))
    assert len(pairings) == 7260
    assert Counter(pairing.sign(0) for pairing in pairings) == {
        -1: 7180,
        0: 80,
    }
    assert extension.same_height_negative_pair_count == 7180
    assert extension.same_height_zero_pair_count == 80
    assert all(
        extension.wall_gram[i][j].sign(0) <= 0
        for i in cohort
        for j in range(i)
    )


def test_targeted_edge_extension_audit_is_exact() -> None:
    extension = run_silver_k3_vinberg_four_plus_two_sqrt2_extension()

    assert extension.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
    assert extension.finite_endpoint_root_indices == (33,)
    assert extension.new_finite_endpoint_root_indices == ()
    assert len(extension.new_root_edge_pairings) == 121
    assert len(extension.new_root_edge_extension_determinants) == 121
    assert len(extension.new_root_edge_extension_identity_inertias) == 121
    assert len(extension.new_root_edge_extension_conjugate_inertias) == 121
    assert all(
        determinant.sign(0) < 0
        for determinant in extension.new_root_edge_extension_determinants
    )
    assert extension.new_root_edge_extension_identity_inertias == (
        Inertia(8, 1, 0),
    ) * 121
    assert extension.new_root_edge_extension_conjugate_inertias == (
        Inertia(9, 0, 0),
    ) * 121
    assert extension.ideal_endpoint_audit_status == (
        "not_run_until_complete_simple_roots"
    )
    assert extension.stopping_status == "undecidable_incomplete_simple_roots"
    assert extension.second_endpoint_status == "not_found_through_this_shell"


def test_block_factorized_replay_matches_production_root_sets() -> None:
    audit = audit_block_silver_four_plus_two_shell()

    assert audit.passes
    assert audit.component_box_counts == (130321, 130321, 19)
    assert audit.component_bounded_counts == (169, 46, 19)
    assert audit.norm_join_probe_count == 3211
    assert audit.exact_weight_solutions == 653
    assert audit.nonintegral_reconstructions == 330
    assert audit.rejected_root_count == 202
    assert len(audit.accepted_roots) == 121
    assert audit.same_height_negative_pair_count == 7180
    assert audit.same_height_zero_pair_count == 80
    assert audit.all_prefix_pairings_nonpositive
    assert audit.prefix_wall_count == 274
    assert audit.prefix_identity_inertia == Inertia(9, 1, 264)
    assert audit.prefix_conjugate_inertia == Inertia(10, 0, 264)
    assert audit.finite_endpoint_root_indices == (33,)
    assert audit.new_finite_endpoint_root_indices == ()
    assert audit.next_admissible_control_coordinate == 3 + 3 * SQRT2
    assert audit.production_accepted_root_set_matches
    assert audit.production_rejection_witness_set_matches
    assert audit.replay_sha256 == EXPECTED_REPLAY_SHA256


def test_four_plus_two_sqrt2_certificate_is_frozen() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))

    summary = verify_manifest(manifest)
    assert summary["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    assert summary["exact_weight_solutions"] == 653
    assert summary["accepted_root_count"] == 121
    assert summary["wall_count"] == 274
    assert summary["known_finite_endpoint_root_indices"] == [33]
    assert summary["independent_replay_sha256"] == EXPECTED_REPLAY_SHA256
    assert summary["next_admissible_control_coordinate"] == {"a": "3", "b": "3"}
    assert summary["stopping_status"] == "undecidable_incomplete_simple_roots"


def test_four_plus_two_sqrt2_certificate_rejects_mutation() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    manifest["new_level"]["accepted_root_count"] = 120

    with pytest.raises(ValueError, match="does not reproduce"):
        verify_manifest(manifest)

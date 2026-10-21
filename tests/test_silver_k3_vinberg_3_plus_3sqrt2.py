import json
from collections import Counter

import pytest

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import SILVER_ZERO, SQRT2, Silver
from coxeter_polytopes.search.silver_k3_vinberg import (
    run_silver_k3_vinberg_three_plus_three_sqrt2_extension,
)
from scripts.certify_silver_k3_vinberg_3_plus_3sqrt2 import (
    ARCHIVE,
    EXPECTED_CERTIFICATE_SHA256,
    verify_manifest,
)
from scripts.independent_silver_k3_vinberg_3_plus_3sqrt2 import (
    EXPECTED_REPLAY_SHA256,
    audit_block_silver_three_plus_three_shell,
)


def test_three_plus_three_sqrt2_shell_is_exact_and_adds_59_walls() -> None:
    extension = run_silver_k3_vinberg_three_plus_three_sqrt2_extension()
    level = extension.new_level

    assert extension.passes
    assert extension.completed_through_control_coordinate == 3 + 3 * SQRT2
    assert extension.next_admissible_control_coordinate == 5 + 2 * SQRT2
    assert level.target == 65 + 45 * SQRT2
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
        5 + 4 * SQRT2,
        6 + 4 * SQRT2,
        7 + 4 * SQRT2,
        6 + 5 * SQRT2,
        7 + 5 * SQRT2,
        8 + 5 * SQRT2,
        7 + 6 * SQRT2,
    )
    assert level.visited_by_depth == (
        1,
        19,
        361,
        1197,
        2318,
        3154,
        23218,
        77463,
        56506,
        18943,
    )
    assert level.schur_pruned_by_depth == (
        0,
        0,
        298,
        1075,
        2152,
        1932,
        19141,
        74489,
        55509,
        13412,
    )
    assert level.exact_weight_solutions == 154
    assert level.nonintegral_reconstructions == 0
    assert level.previous_wall_rejections == 95
    assert len(level.roots) == 59
    assert len(extension.roots) == 333
    assert extension.wall_identity_inertia == Inertia(9, 1, 323)
    assert extension.wall_conjugate_inertia == Inertia(10, 0, 323)


def test_all_59_same_height_roots_survive_simultaneously() -> None:
    extension = run_silver_k3_vinberg_three_plus_three_sqrt2_extension()
    cohort = tuple(range(extension.prior_wall_count, len(extension.roots)))
    pairings = tuple(
        extension.wall_gram[i][j]
        for i in cohort
        for j in cohort
        if j < i
    )

    assert cohort == tuple(range(274, 333))
    assert len(pairings) == 1711
    assert Counter(pairing.sign(0) for pairing in pairings) == {-1: 1689, 0: 22}
    assert extension.same_height_negative_pair_count == 1689
    assert extension.same_height_zero_pair_count == 22
    assert all(
        extension.wall_gram[i][j].sign(0) <= 0
        for i in cohort
        for j in range(i)
    )


def test_targeted_edge_remains_open_at_three_plus_three_sqrt2() -> None:
    extension = run_silver_k3_vinberg_three_plus_three_sqrt2_extension()

    assert extension.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
    assert extension.finite_endpoint_root_indices == (33,)
    assert extension.new_finite_endpoint_root_indices == ()
    assert len(extension.new_root_edge_pairings) == 59
    assert all(
        determinant.sign(0) < 0
        for determinant in extension.new_root_edge_extension_determinants
    )
    assert extension.new_root_edge_extension_identity_inertias == (
        Inertia(8, 1, 0),
    ) * 59
    assert extension.new_root_edge_extension_conjugate_inertias == (
        Inertia(9, 0, 0),
    ) * 59
    assert extension.ideal_endpoint_audit_status == (
        "not_run_until_complete_simple_roots"
    )
    assert extension.second_endpoint_status == "not_found_through_this_shell"
    assert extension.stopping_status == "undecidable_incomplete_simple_roots"


def test_block_factorized_replay_matches_three_plus_three_production() -> None:
    audit = audit_block_silver_three_plus_three_shell()

    assert audit.passes
    assert audit.component_box_counts == (130321, 130321, 19)
    assert audit.component_bounded_counts == (166, 38, 19)
    assert audit.norm_join_probe_count == 3154
    assert audit.exact_weight_solutions == 154
    assert audit.nonintegral_reconstructions == 0
    assert audit.rejected_root_count == 95
    assert len(audit.accepted_roots) == 59
    assert audit.same_height_negative_pair_count == 1689
    assert audit.same_height_zero_pair_count == 22
    assert audit.all_prefix_pairings_nonpositive
    assert audit.prefix_wall_count == 333
    assert audit.prefix_identity_inertia == Inertia(9, 1, 323)
    assert audit.prefix_conjugate_inertia == Inertia(10, 0, 323)
    assert audit.finite_endpoint_root_indices == (33,)
    assert audit.new_finite_endpoint_root_indices == ()
    assert audit.next_admissible_control_coordinate == 5 + 2 * SQRT2
    assert audit.production_accepted_root_set_matches
    assert audit.production_rejection_witness_set_matches
    assert audit.replay_sha256 == EXPECTED_REPLAY_SHA256


def test_three_plus_three_sqrt2_certificate_is_frozen() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))

    summary = verify_manifest(manifest)
    assert summary["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    assert summary["exact_weight_solutions"] == 154
    assert summary["accepted_root_count"] == 59
    assert summary["wall_count"] == 333
    assert summary["known_finite_endpoint_root_indices"] == [33]
    assert summary["independent_replay_sha256"] == EXPECTED_REPLAY_SHA256
    assert summary["next_admissible_control_coordinate"] == {"a": "5", "b": "2"}
    assert summary["stopping_status"] == "undecidable_incomplete_simple_roots"


def test_three_plus_three_sqrt2_certificate_rejects_mutation() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    manifest["new_level"]["accepted_root_count"] = 58

    with pytest.raises(ValueError, match="does not reproduce"):
        verify_manifest(manifest)

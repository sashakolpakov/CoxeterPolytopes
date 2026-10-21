import json

import pytest

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import SILVER_ZERO, SQRT2
from coxeter_polytopes.search.silver_k3_vinberg import (
    _silver_pair,
    run_silver_k3_vinberg_five_plus_two_sqrt2_extension,
    silver_k3_alexeev_field_gram,
)
from scripts.independent_silver_k3_vinberg_5_plus_2sqrt2 import (
    EXPECTED_REPLAY_SHA256,
    audit_block_silver_five_plus_two_shell,
)
from scripts.certify_silver_k3_vinberg_5_plus_2sqrt2 import (
    ARCHIVE,
    EXPECTED_CERTIFICATE_SHA256,
    verify_manifest,
)


def test_five_plus_two_sqrt2_shell_is_exact_and_empty() -> None:
    extension = run_silver_k3_vinberg_five_plus_two_sqrt2_extension()
    level = extension.new_level

    assert extension.passes
    assert extension.completed_through_control_coordinate == 5 + 2 * SQRT2
    assert extension.next_admissible_control_coordinate == 4 + 3 * SQRT2
    assert level.target == 75 + 53 * SQRT2
    assert level.alphabet.values == (
        SILVER_ZERO,
        3 + 2 * SQRT2,
        4 + 3 * SQRT2,
        7 + 5 * SQRT2,
    )
    assert (
        level.alphabet.identity_upper,
        level.alphabet.conjugate_upper,
        level.alphabet.coefficient_bound,
    ) == (18, 1, 19)
    assert level.visited_by_depth == (1, 4, 16, 16, 16, 12, 32, 44, 24, 12)
    assert level.schur_pruned_by_depth == (0, 0, 12, 12, 13, 4, 21, 38, 21, 4)
    assert level.exact_weight_solutions == 1
    assert level.nonintegral_reconstructions == 0
    assert level.previous_wall_rejections == 1
    assert level.roots == ()


def test_unique_candidate_has_its_complete_lower_wall_blocker_set() -> None:
    extension = run_silver_k3_vinberg_five_plus_two_sqrt2_extension()
    ((root, blockers),) = extension.new_level.previous_wall_rejection_witnesses
    gram = silver_k3_alexeev_field_gram()
    recomputed = tuple(
        index
        for index, wall in enumerate(extension.roots)
        if _silver_pair(root, gram, wall).sign(0) > 0
    )

    assert blockers == recomputed
    assert blockers == (9,)
    assert all(
        _silver_pair(root, gram, extension.roots[index]).sign(0) > 0
        for index in blockers
    )


def test_empty_shell_preserves_prefix_and_targeted_edge_scope() -> None:
    extension = run_silver_k3_vinberg_five_plus_two_sqrt2_extension()

    assert extension.prior_wall_count == 333
    assert len(extension.roots) == 333
    assert extension.wall_identity_inertia == Inertia(9, 1, 323)
    assert extension.wall_conjugate_inertia == Inertia(10, 0, 323)
    assert extension.same_height_negative_pair_count == 0
    assert extension.same_height_zero_pair_count == 0
    assert extension.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
    assert extension.finite_endpoint_root_indices == (33,)
    assert extension.new_finite_endpoint_root_indices == ()
    assert extension.new_root_edge_pairings == ()
    assert extension.ideal_endpoint_audit_status == (
        "not_run_until_complete_simple_roots"
    )
    assert extension.second_endpoint_status == "not_found_through_this_shell"
    assert extension.stopping_status == "undecidable_incomplete_simple_roots"


def test_block_factorized_replay_matches_empty_shell_production() -> None:
    audit = audit_block_silver_five_plus_two_shell()

    assert audit.passes
    assert audit.component_box_counts == (256, 256, 4)
    assert audit.component_bounded_counts == (3, 1, 4)
    assert audit.norm_join_probe_count == 12
    assert audit.exact_weight_solutions == 1
    assert audit.nonintegral_reconstructions == 0
    assert audit.rejected_root_count == 1
    assert audit.accepted_roots == ()
    assert audit.same_height_negative_pair_count == 0
    assert audit.same_height_zero_pair_count == 0
    assert audit.all_prefix_pairings_nonpositive
    assert audit.prefix_wall_count == 333
    assert audit.prefix_identity_inertia == Inertia(9, 1, 323)
    assert audit.prefix_conjugate_inertia == Inertia(10, 0, 323)
    assert audit.finite_endpoint_root_indices == (33,)
    assert audit.new_finite_endpoint_root_indices == ()
    assert audit.next_admissible_control_coordinate == 4 + 3 * SQRT2
    assert audit.production_accepted_root_set_matches
    assert audit.production_rejection_witness_set_matches
    assert audit.replay_sha256 == EXPECTED_REPLAY_SHA256


def test_five_plus_two_sqrt2_certificate_is_frozen() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))

    summary = verify_manifest(manifest)
    assert summary["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    assert summary["exact_weight_solutions"] == 1
    assert summary["accepted_root_count"] == 0
    assert summary["wall_count"] == 333
    assert summary["known_finite_endpoint_root_indices"] == [33]
    assert summary["independent_replay_sha256"] == EXPECTED_REPLAY_SHA256
    assert summary["next_admissible_control_coordinate"] == {"a": "4", "b": "3"}
    assert summary["stopping_status"] == "undecidable_incomplete_simple_roots"


def test_five_plus_two_sqrt2_certificate_rejects_mutation() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    manifest["new_level"]["accepted_root_count"] = 1

    with pytest.raises(ValueError, match="does not reproduce"):
        verify_manifest(manifest)

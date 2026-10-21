import json
from collections import Counter
from hashlib import sha256

import pytest

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import SQRT2
from coxeter_polytopes.search.silver_k3_vinberg import (
    run_silver_k3_vinberg_four_plus_three_sqrt2_extension,
)
from scripts.certify_silver_k3_vinberg_4_plus_3sqrt2 import (
    ARCHIVE,
    EXPECTED_CERTIFICATE_SHA256,
    verify_manifest,
)
from scripts.independent_silver_k3_vinberg_4_plus_3sqrt2 import (
    EXPECTED_REPLAY_SHA256,
    audit_block_silver_four_plus_three_shell,
)


def test_four_plus_three_sqrt2_shell_is_exact_and_adds_393_walls() -> None:
    extension = run_silver_k3_vinberg_four_plus_three_sqrt2_extension()
    level = extension.new_level

    assert extension.passes
    assert extension.completed_through_control_coordinate == 4 + 3 * SQRT2
    assert extension.next_admissible_control_coordinate == 5 + 3 * SQRT2
    assert level.target == 84 + 58 * SQRT2
    assert len(level.alphabet.values) == 26
    assert (
        level.alphabet.identity_upper,
        level.alphabet.conjugate_upper,
        level.alphabet.coefficient_bound,
    ) == (19, 2, 21)
    assert level.visited_by_depth == (
        1,
        26,
        676,
        3016,
        7462,
        12688,
        117078,
        493064,
        425750,
        156000,
    )
    assert level.schur_pruned_by_depth == (
        0,
        0,
        560,
        2729,
        6974,
        8185,
        98114,
        476689,
        419750,
        117374,
    )
    assert level.exact_weight_solutions == 2507
    assert level.nonintegral_reconstructions == 1544
    assert level.previous_wall_rejections == 570
    assert len(level.roots) == 393
    assert len(extension.roots) == 726


def test_all_blocker_sets_and_same_height_pairings_are_complete() -> None:
    extension = run_silver_k3_vinberg_four_plus_three_sqrt2_extension()
    level = extension.new_level
    cohort = tuple(range(extension.prior_wall_count, len(extension.roots)))
    pairings = tuple(
        extension.wall_gram[i][j]
        for i in cohort
        for j in cohort
        if j < i
    )

    assert Counter(
        len(blockers) for _root, blockers in level.previous_wall_rejection_witnesses
    ) == {1: 400, 2: 170}
    assert all(blockers for _root, blockers in level.previous_wall_rejection_witnesses)
    assert len(pairings) == 77028
    assert Counter(pairing.sign(0) for pairing in pairings) == {-1: 76453, 0: 575}
    assert all(
        extension.wall_gram[i][j].sign(0) <= 0
        for i in cohort
        for j in range(i)
    )


def test_targeted_edge_closes_locally_but_prefix_remains_incomplete() -> None:
    extension = run_silver_k3_vinberg_four_plus_three_sqrt2_extension()
    local_index = 622 - extension.prior_wall_count

    assert extension.wall_identity_inertia == Inertia(9, 1, 716)
    assert extension.wall_conjugate_inertia == Inertia(10, 0, 716)
    assert extension.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
    assert extension.finite_endpoint_root_indices == (33, 622)
    assert extension.new_finite_endpoint_root_indices == (622,)
    assert extension.new_root_edge_extension_determinants[local_index].sign(0) > 0
    assert extension.new_root_edge_extension_identity_inertias[local_index] == Inertia(
        9, 0, 0
    )
    assert Counter(
        determinant.sign(0)
        for determinant in extension.new_root_edge_extension_determinants
    ) == {-1: 392, 1: 1}
    assert Counter(extension.new_root_edge_extension_identity_inertias) == {
        Inertia(8, 1, 0): 392,
        Inertia(9, 0, 0): 1,
    }
    assert extension.new_root_edge_extension_conjugate_inertias == (
        Inertia(9, 0, 0),
    ) * 393
    assert extension.second_endpoint_status == "finite_endpoint_found_at_this_shell"
    assert extension.ideal_endpoint_audit_status == (
        "not_run_until_complete_simple_roots"
    )
    assert extension.stopping_status == "undecidable_incomplete_simple_roots"


def test_block_factorized_replay_matches_four_plus_three_production() -> None:
    audit = audit_block_silver_four_plus_three_shell()

    assert audit.passes
    assert audit.component_box_counts == (456976, 456976, 26)
    assert audit.component_bounded_counts == (488, 115, 26)
    assert audit.norm_join_probe_count == 12688
    assert audit.exact_weight_solutions == 2507
    assert audit.nonintegral_reconstructions == 1544
    assert audit.rejected_root_count == 570
    assert len(audit.accepted_roots) == 393
    assert audit.same_height_negative_pair_count == 76453
    assert audit.same_height_zero_pair_count == 575
    assert audit.prefix_wall_count == 726
    assert audit.prefix_identity_inertia == Inertia(9, 1, 716)
    assert audit.prefix_conjugate_inertia == Inertia(10, 0, 716)
    assert audit.finite_endpoint_root_indices == (33, 622)
    assert audit.new_finite_endpoint_root_indices == (622,)
    assert audit.native_finite_endpoint_root_indices == (33, 443)
    assert audit.native_new_finite_endpoint_root_indices == (443,)
    assert audit.next_admissible_control_coordinate == 5 + 3 * SQRT2
    assert audit.production_accepted_root_set_matches
    assert audit.production_rejection_witness_set_matches
    assert audit.replay_sha256 == EXPECTED_REPLAY_SHA256


def test_four_plus_three_sqrt2_certificate_is_frozen() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))

    summary = verify_manifest(manifest)
    assert summary["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    assert summary["exact_weight_solutions"] == 2507
    assert summary["accepted_root_count"] == 393
    assert summary["wall_count"] == 726
    assert summary["known_finite_endpoint_root_indices"] == [33, 622]
    assert summary["independent_replay_sha256"] == EXPECTED_REPLAY_SHA256
    assert summary["next_admissible_control_coordinate"] == {"a": "5", "b": "3"}
    assert summary["stopping_status"] == "undecidable_incomplete_simple_roots"
    assert len(manifest["new_level"]["accepted_roots"]) == 393
    assert manifest["new_level"]["accepted_root_order"] == (
        "production_weight_recursion_order"
    )
    assert manifest["ordered_prior_prefix"]["wall_count"] == 333
    assert manifest["ordered_prior_prefix"]["certificate_sha256"] == (
        "e2283e6330c93632035febe1b84ef1d499a7b9507931008940b956f1470e10d5"
    )
    ordered_roots = manifest["ordered_prefix_roots"]
    accepted_roots = manifest["new_level"]["accepted_roots"]
    assert manifest["ordered_prefix_root_count"] == len(ordered_roots) == 726
    assert ordered_roots[333:] == accepted_roots
    assert all(
        set(coordinate) == {"a", "b"}
        for root in ordered_roots
        for coordinate in root
    )
    assert manifest["ordered_prefix_roots_sha256"] == sha256(
        json.dumps(ordered_roots, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def test_four_plus_three_sqrt2_certificate_rejects_mutation() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    manifest["new_level"]["accepted_root_count"] = 392

    with pytest.raises(ValueError, match="does not reproduce"):
        verify_manifest(manifest)

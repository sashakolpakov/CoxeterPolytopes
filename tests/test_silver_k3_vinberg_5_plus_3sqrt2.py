import json
from collections import Counter
from hashlib import sha256

import pytest

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import SQRT2
from coxeter_polytopes.search.silver_k3_vinberg import (
    run_silver_k3_vinberg_five_plus_three_sqrt2_extension,
)
from scripts.certify_silver_k3_vinberg_5_plus_3sqrt2 import (
    ARCHIVE,
    EXPECTED_CERTIFICATE_SHA256,
    verify_manifest,
)
from scripts.independent_silver_k3_vinberg_5_plus_3sqrt2 import (
    EXPECTED_REPLAY_SHA256,
    build_census,
)


def test_five_plus_three_sqrt2_shell_is_exact_and_adds_350_walls() -> None:
    extension = run_silver_k3_vinberg_five_plus_three_sqrt2_extension()
    level = extension.new_level

    assert extension.passes
    assert extension.completed_through_control_coordinate == 5 + 3 * SQRT2
    assert extension.next_admissible_control_coordinate == 4 + 4 * SQRT2
    assert level.target == 105 + 73 * SQRT2
    assert len(level.alphabet.values) == 27
    assert (
        level.alphabet.identity_upper,
        level.alphabet.conjugate_upper,
        level.alphabet.coefficient_bound,
    ) == (21, 2, 23)
    assert level.visited_by_depth == (
        1,
        27,
        729,
        3456,
        8775,
        15120,
        146556,
        639441,
        553014,
        191916,
    )
    assert level.schur_pruned_by_depth == (
        0,
        0,
        601,
        3131,
        8215,
        9692,
        122873,
        618959,
        545906,
        142182,
    )
    assert level.exact_weight_solutions == 861
    assert level.nonintegral_reconstructions == 0
    assert level.previous_wall_rejections == 511
    assert len(level.roots) == 350
    assert len(extension.roots) == 1076


def test_five_plus_three_blockers_and_same_height_pairs_are_complete() -> None:
    extension = run_silver_k3_vinberg_five_plus_three_sqrt2_extension()
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
    ) == {1: 323, 2: 188}
    assert len(pairings) == 61075
    assert Counter(pairing.sign(0) for pairing in pairings) == {-1: 60759, 0: 316}
    assert all(
        extension.wall_gram[i][j].sign(0) <= 0
        for i in cohort
        for j in range(i)
    )


def test_targeted_edge_stays_locally_closed_but_prefix_is_incomplete() -> None:
    extension = run_silver_k3_vinberg_five_plus_three_sqrt2_extension()

    assert extension.wall_identity_inertia == Inertia(9, 1, 1066)
    assert extension.wall_conjugate_inertia == Inertia(10, 0, 1066)
    assert extension.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
    assert extension.finite_endpoint_root_indices == (33, 622)
    assert extension.new_finite_endpoint_root_indices == ()
    assert Counter(
        determinant.sign(0)
        for determinant in extension.new_root_edge_extension_determinants
    ) == {-1: 350}
    assert extension.new_root_edge_extension_identity_inertias == (
        Inertia(8, 1, 0),
    ) * 350
    assert extension.new_root_edge_extension_conjugate_inertias == (
        Inertia(9, 0, 0),
    ) * 350
    assert extension.ideal_endpoint_audit_status == (
        "not_run_until_complete_simple_roots"
    )
    assert extension.stopping_status == "undecidable_incomplete_simple_roots"


def test_block_factorized_replay_matches_five_plus_three_production() -> None:
    replay = build_census()

    assert replay["component_box_counts"] == [531441, 531441, 27]
    assert replay["component_bounded_counts"] == [560, 126, 27]
    assert replay["norm_join_probe_count"] == 15120
    assert replay["exact_weight_solutions"] == 861
    assert replay["nonintegral_reconstructions"] == 0
    assert replay["rejected_root_count"] == 511
    assert replay["accepted_root_count"] == 350
    assert replay["same_height_negative_pair_count"] == 60759
    assert replay["same_height_zero_pair_count"] == 316
    assert replay["prefix_wall_count"] == 1076
    assert replay["census_sha256"] == EXPECTED_REPLAY_SHA256
    assert replay["stopping_status"] == "undecidable_incomplete_simple_roots"


def test_five_plus_three_sqrt2_certificate_is_frozen_and_literal() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))

    summary = verify_manifest(manifest)
    assert summary["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    assert summary["exact_weight_solutions"] == 861
    assert summary["accepted_root_count"] == 350
    assert summary["wall_count"] == 1076
    assert summary["known_finite_endpoint_root_indices"] == [33, 622]
    assert summary["independent_replay_sha256"] == EXPECTED_REPLAY_SHA256
    assert summary["next_admissible_control_coordinate"] == {"a": "4", "b": "4"}
    assert summary["stopping_status"] == "undecidable_incomplete_simple_roots"
    ordered_roots = manifest["ordered_prefix_roots"]
    accepted_roots = manifest["new_level"]["accepted_roots"]
    assert manifest["ordered_prefix_root_count"] == len(ordered_roots) == 1076
    assert len(accepted_roots) == 350
    assert ordered_roots[726:] == accepted_roots
    assert all(
        set(coordinate) == {"a", "b"}
        for root in ordered_roots
        for coordinate in root
    )
    assert manifest["ordered_prefix_roots_sha256"] == sha256(
        json.dumps(ordered_roots, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def test_five_plus_three_sqrt2_certificate_rejects_mutation() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    manifest["new_level"]["accepted_root_count"] = 349

    with pytest.raises(ValueError, match="does not reproduce"):
        verify_manifest(manifest)

import json
from collections import Counter
from hashlib import sha256

import pytest

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import SQRT2
from coxeter_polytopes.search.silver_k3_vinberg import (
    run_silver_k3_vinberg_four_plus_four_sqrt2_extension,
)
from scripts.certify_silver_k3_vinberg_4_plus_4sqrt2 import (
    ARCHIVE,
    EXPECTED_CERTIFICATE_SHA256,
    verify_manifest,
)
from scripts.independent_silver_k3_vinberg_4_plus_4sqrt2 import (
    EXPECTED_REPLAY_SHA256,
    build_census,
)


def test_four_plus_four_sqrt2_shell_is_exact_and_adds_142_walls() -> None:
    extension = run_silver_k3_vinberg_four_plus_four_sqrt2_extension()
    level = extension.new_level

    assert extension.passes
    assert extension.completed_through_control_coordinate == 4 + 4 * SQRT2
    assert extension.next_admissible_control_coordinate == 6 + 3 * SQRT2
    assert level.target == 114 + 80 * SQRT2
    assert len(level.alphabet.values) == 21
    assert (
        level.alphabet.identity_upper,
        level.alphabet.conjugate_upper,
        level.alphabet.coefficient_bound,
    ) == (22, 2, 24)
    assert level.visited_by_depth == (
        1,
        21,
        441,
        1533,
        3192,
        4746,
        35049,
        122157,
        95718,
        36141,
    )
    assert level.schur_pruned_by_depth == (
        0,
        0,
        368,
        1381,
        2966,
        3077,
        29232,
        117599,
        93997,
        27184,
    )
    assert level.exact_weight_solutions == 957
    assert level.nonintegral_reconstructions == 491
    assert level.previous_wall_rejections == 324
    assert len(level.roots) == 142
    assert len(extension.roots) == 1218


def test_four_plus_four_blockers_and_pairings_are_complete() -> None:
    extension = run_silver_k3_vinberg_four_plus_four_sqrt2_extension()
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
    ) == {1: 255, 2: 69}
    assert len(pairings) == 10011
    assert Counter(pairing.sign(0) for pairing in pairings) == {-1: 10011}
    assert all(
        extension.wall_gram[i][j].sign(0) <= 0
        for i in cohort
        for j in range(i)
    )


def test_four_plus_four_targeted_edge_does_not_imply_global_stopping() -> None:
    extension = run_silver_k3_vinberg_four_plus_four_sqrt2_extension()

    assert extension.wall_identity_inertia == Inertia(9, 1, 1208)
    assert extension.wall_conjugate_inertia == Inertia(10, 0, 1208)
    assert extension.finite_endpoint_root_indices == (33, 622)
    assert extension.new_finite_endpoint_root_indices == ()
    assert Counter(
        determinant.sign(0)
        for determinant in extension.new_root_edge_extension_determinants
    ) == {-1: 142}
    assert extension.new_root_edge_extension_identity_inertias == (
        Inertia(8, 1, 0),
    ) * 142
    assert extension.new_root_edge_extension_conjugate_inertias == (
        Inertia(9, 0, 0),
    ) * 142
    assert extension.stopping_status == "undecidable_incomplete_simple_roots"


def test_four_plus_four_independent_replay_is_frozen() -> None:
    replay = build_census()

    assert replay["component_box_counts"] == [194481, 194481, 21]
    assert replay["component_bounded_counts"] == [226, 61, 21]
    assert replay["norm_join_probe_count"] == 4746
    assert replay["exact_weight_solutions"] == 957
    assert replay["nonintegral_reconstructions"] == 491
    assert replay["rejected_root_count"] == 324
    assert replay["accepted_root_count"] == 142
    assert replay["same_height_negative_pair_count"] == 10011
    assert replay["same_height_zero_pair_count"] == 0
    assert replay["census_sha256"] == EXPECTED_REPLAY_SHA256


def test_four_plus_four_certificate_is_frozen_and_literal() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))

    summary = verify_manifest(manifest)
    assert summary["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    assert summary["exact_weight_solutions"] == 957
    assert summary["accepted_root_count"] == 142
    assert summary["wall_count"] == 1218
    assert summary["known_finite_endpoint_root_indices"] == [33, 622]
    assert summary["independent_replay_sha256"] == EXPECTED_REPLAY_SHA256
    assert summary["next_admissible_control_coordinate"] == {"a": "6", "b": "3"}
    assert summary["stopping_status"] == "undecidable_incomplete_simple_roots"
    ordered_roots = manifest["ordered_prefix_roots"]
    accepted_roots = manifest["new_level"]["accepted_roots"]
    assert manifest["ordered_prefix_root_count"] == len(ordered_roots) == 1218
    assert len(accepted_roots) == 142
    assert ordered_roots[1076:] == accepted_roots
    assert manifest["ordered_prefix_roots_sha256"] == sha256(
        json.dumps(ordered_roots, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def test_four_plus_four_certificate_rejects_mutation() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    manifest["new_level"]["accepted_root_count"] = 141

    with pytest.raises(ValueError, match="does not reproduce"):
        verify_manifest(manifest)

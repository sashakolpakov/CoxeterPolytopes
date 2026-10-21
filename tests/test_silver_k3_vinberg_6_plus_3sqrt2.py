from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json

import pytest

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import SQRT2
from coxeter_polytopes.search.silver_k3_vinberg import (
    run_silver_k3_vinberg_six_plus_three_sqrt2_extension,
)
from scripts.certify_silver_k3_vinberg_6_plus_3sqrt2 import (
    ARCHIVE,
    EXPECTED_CERTIFICATE_SHA256,
    verify_manifest,
)
from scripts.independent_silver_k3_vinberg_6_plus_3sqrt2 import (
    EXPECTED_REPLAY_SHA256,
    build_census,
)


def test_six_plus_three_sqrt2_shell_is_exact_and_adds_73_walls() -> None:
    extension = run_silver_k3_vinberg_six_plus_three_sqrt2_extension()
    level = extension.new_level

    assert extension.passes
    assert extension.completed_through_control_coordinate == 6 + 3 * SQRT2
    assert extension.next_admissible_control_coordinate == 5 + 4 * SQRT2
    assert level.target == 128 + 90 * SQRT2
    assert len(level.alphabet.values) == 20
    assert (
        level.alphabet.identity_upper,
        level.alphabet.conjugate_upper,
        level.alphabet.coefficient_bound,
    ) == (23, 2, 25)
    assert level.visited_by_depth == (
        1,
        20,
        400,
        1380,
        2780,
        3920,
        28820,
        98760,
        76020,
        27780,
    )
    assert level.schur_pruned_by_depth == (
        0,
        0,
        331,
        1241,
        2584,
        2479,
        23882,
        94959,
        74631,
        20465,
    )
    assert level.exact_weight_solutions == 697
    assert level.nonintegral_reconstructions == 425
    assert level.previous_wall_rejections == 199
    assert len(level.roots) == 73
    assert len(extension.roots) == 1291


def test_six_plus_three_blockers_and_pairings_are_complete() -> None:
    extension = run_silver_k3_vinberg_six_plus_three_sqrt2_extension()
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
    ) == {1: 127, 2: 72}
    assert len(pairings) == 2628
    assert Counter(pairing.sign(0) for pairing in pairings) == {-1: 2628}
    assert all(
        extension.wall_gram[i][j].sign(0) <= 0
        for i in cohort
        for j in range(i)
    )


def test_six_plus_three_prefix_inertia_does_not_imply_stopping() -> None:
    extension = run_silver_k3_vinberg_six_plus_three_sqrt2_extension()

    assert extension.wall_identity_inertia == Inertia(9, 1, 1281)
    assert extension.wall_conjugate_inertia == Inertia(10, 0, 1281)
    assert extension.stopping_status == "undecidable_incomplete_simple_roots"


def test_six_plus_three_independent_replay_is_frozen() -> None:
    replay = build_census()

    assert replay["component_box_counts"] == [160000, 160000, 20]
    assert replay["component_bounded_counts"] == [196, 51, 20]
    assert replay["norm_join_probe_count"] == 3920
    assert replay["exact_weight_solutions"] == 697
    assert replay["nonintegral_reconstructions"] == 425
    assert replay["rejected_root_count"] == 199
    assert replay["accepted_root_count"] == 73
    assert replay["same_height_negative_pair_count"] == 2628
    assert replay["same_height_zero_pair_count"] == 0
    assert replay["census_sha256"] == EXPECTED_REPLAY_SHA256


def test_six_plus_three_certificate_is_frozen_and_literal() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))

    summary = verify_manifest(manifest)
    assert summary["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    assert summary["exact_weight_solutions"] == 697
    assert summary["accepted_root_count"] == 73
    assert summary["wall_count"] == 1291
    assert summary["independent_replay_sha256"] == EXPECTED_REPLAY_SHA256
    assert summary["wall_identity_inertia"] == {
        "positive": 9,
        "negative": 1,
        "zero": 1281,
    }
    assert summary["wall_conjugate_inertia"] == {
        "positive": 10,
        "negative": 0,
        "zero": 1281,
    }
    assert summary["next_admissible_control_coordinate"] == {"a": "5", "b": "4"}
    assert summary["stopping_status"] == "undecidable_incomplete_simple_roots"
    ordered_roots = manifest["ordered_prefix_roots"]
    accepted_roots = manifest["new_level"]["accepted_roots"]
    rejection_witnesses = manifest["new_level"][
        "production_rejection_witnesses"
    ]
    assert manifest["ordered_prefix_root_count"] == len(ordered_roots) == 1291
    assert len(accepted_roots) == 73
    assert len(rejection_witnesses) == 199
    assert ordered_roots[1218:] == accepted_roots
    assert manifest["ordered_prefix_roots_sha256"] == sha256(
        json.dumps(ordered_roots, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    assert manifest["new_level"][
        "production_rejection_witness_set_sha256"
    ] == sha256(
        json.dumps(
            rejection_witnesses, sort_keys=True, separators=(",", ":")
        ).encode()
    ).hexdigest()


def test_six_plus_three_certificate_rejects_mutation() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    manifest["new_level"]["accepted_root_count"] = 72

    with pytest.raises(ValueError, match="does not reproduce"):
        verify_manifest(manifest)

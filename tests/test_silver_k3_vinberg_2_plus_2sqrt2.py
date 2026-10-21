import json
from collections import Counter

import pytest

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import SILVER_ZERO, SQRT2, Silver
from coxeter_polytopes.search.silver_k3_vinberg import (
    run_silver_k3_vinberg_two_plus_two_sqrt2_extension,
)
from scripts.certify_silver_k3_vinberg_2_plus_2sqrt2 import (
    ARCHIVE,
    EXPECTED_CERTIFICATE_SHA256,
    verify_manifest,
)
from scripts.independent_silver_k3_vinberg_2_plus_2sqrt2 import (
    EXPECTED_REPLAY_SHA256,
    audit_reverse_silver_two_plus_two_shell,
)


def test_two_plus_two_sqrt2_shell_is_exact_and_adds_56_walls() -> None:
    extension = run_silver_k3_vinberg_two_plus_two_sqrt2_extension()
    level = extension.new_level

    assert extension.passes
    assert extension.completed_through_control_coordinate == 2 + 2 * SQRT2
    assert extension.next_admissible_control_coordinate == 3 + 2 * SQRT2
    assert level.target == 30 + 20 * SQRT2
    assert level.alphabet.values == (
        SILVER_ZERO,
        Silver(1),
        SQRT2,
        1 + SQRT2,
        2 + SQRT2,
        1 + 2 * SQRT2,
        3 + SQRT2,
        2 + 2 * SQRT2,
        3 + 2 * SQRT2,
        4 + 2 * SQRT2,
        3 + 3 * SQRT2,
        4 + 3 * SQRT2,
        5 + 3 * SQRT2,
        4 + 4 * SQRT2,
        6 + 3 * SQRT2,
        5 + 4 * SQRT2,
    )
    assert level.exact_weight_solutions == 277
    assert level.nonintegral_reconstructions == 129
    assert level.previous_wall_rejections == 92
    assert len(level.roots) == 56
    assert len(extension.roots) == 93
    assert extension.wall_identity_inertia == Inertia(9, 1, 83)
    assert extension.wall_conjugate_inertia == Inertia(10, 0, 83)


def test_all_56_same_height_roots_survive_simultaneously() -> None:
    extension = run_silver_k3_vinberg_two_plus_two_sqrt2_extension()
    cohort = tuple(range(extension.prior_wall_count, len(extension.roots)))
    pairings = tuple(
        extension.wall_gram[i][j]
        for i in cohort
        for j in cohort
        if j < i
    )

    assert cohort == tuple(range(37, 93))
    assert len(pairings) == 1540
    assert Counter(pairing.sign(0) for pairing in pairings) == {-1: 1482, 0: 58}
    assert all(
        extension.wall_gram[i][j].sign(0) <= 0
        for i in cohort
        for j in range(i)
    )


def test_targeted_edge_does_not_gain_its_second_finite_endpoint() -> None:
    extension = run_silver_k3_vinberg_two_plus_two_sqrt2_extension()

    assert extension.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
    assert extension.finite_endpoint_root_indices == (33,)
    assert extension.new_finite_endpoint_root_indices == ()
    assert len(extension.new_root_edge_pairings) == 56
    assert len(extension.new_root_edge_extension_determinants) == 56
    assert all(
        determinant.sign(0) < 0
        for determinant in extension.new_root_edge_extension_determinants
    )
    assert extension.new_root_edge_extension_identity_inertias == (
        Inertia(8, 1, 0),
    ) * 56
    assert extension.new_root_edge_extension_conjugate_inertias == (
        Inertia(9, 0, 0),
    ) * 56
    assert extension.second_endpoint_status == "not_found_through_this_shell"
    assert extension.stopping_status == "undecidable_incomplete_simple_roots"


def test_reverse_coordinate_replay_matches_production_root_sets() -> None:
    audit = audit_reverse_silver_two_plus_two_shell()

    assert audit.passes
    assert audit.exact_weight_solutions == 277
    assert audit.nonintegral_reconstructions == 129
    assert audit.rejected_root_count == 92
    assert len(audit.accepted_roots) == 56
    assert audit.same_height_negative_pair_count == 1482
    assert audit.same_height_zero_pair_count == 58
    assert audit.production_accepted_root_set_matches
    assert audit.production_rejection_witness_set_matches
    assert audit.replay_sha256 == EXPECTED_REPLAY_SHA256


def test_two_plus_two_sqrt2_certificate_is_frozen() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))

    assert verify_manifest(manifest) == {
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
        "exact_weight_solutions": 277,
        "accepted_root_count": 56,
        "wall_count": 93,
        "known_finite_endpoint_root_indices": [33],
        "independent_replay_sha256": EXPECTED_REPLAY_SHA256,
        "next_admissible_control_coordinate": {"a": "3", "b": "2"},
        "stopping_status": "undecidable_incomplete_simple_roots",
    }


def test_two_plus_two_sqrt2_certificate_rejects_mutation() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    manifest["new_level"]["accepted_root_count"] = 55

    with pytest.raises(ValueError, match="does not reproduce"):
        verify_manifest(manifest)

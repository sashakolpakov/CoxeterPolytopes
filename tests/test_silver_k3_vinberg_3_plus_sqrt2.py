import json

import pytest

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import SILVER_ZERO, SQRT2, Silver
from coxeter_polytopes.search.silver_k3_vinberg import (
    run_silver_k3_vinberg_three_plus_sqrt2_extension,
)
from scripts.certify_silver_k3_vinberg_3_plus_sqrt2 import (
    ARCHIVE,
    EXPECTED_CERTIFICATE_SHA256,
    verify_manifest,
)
from scripts.independent_silver_k3_vinberg_3_plus_sqrt2 import (
    EXPECTED_REPLAY_SHA256,
    audit_reverse_silver_shell,
)


def test_three_plus_sqrt2_shell_is_exact_and_adds_three_walls() -> None:
    extension = run_silver_k3_vinberg_three_plus_sqrt2_extension()
    level = extension.new_level

    assert extension.passes
    assert extension.completed_through_control_coordinate == 3 + SQRT2
    assert extension.next_admissible_control_coordinate == 2 + 2 * SQRT2
    assert level.target == 25 + 17 * SQRT2
    assert level.alphabet.values == (
        Silver(0),
        Silver(1),
        1 + SQRT2,
        2 + SQRT2,
        2 + 2 * SQRT2,
        3 + 2 * SQRT2,
        4 + 2 * SQRT2,
        3 + 3 * SQRT2,
        4 + 3 * SQRT2,
        5 + 3 * SQRT2,
    )
    assert level.exact_weight_solutions == 12
    assert level.nonintegral_reconstructions == 0
    assert level.previous_wall_rejections == 9
    assert len(level.roots) == 3
    assert len(extension.roots) == 37
    assert extension.wall_identity_inertia == Inertia(9, 1, 27)
    assert extension.wall_conjugate_inertia == Inertia(10, 0, 27)


def test_all_three_same_height_roots_survive_simultaneously() -> None:
    extension = run_silver_k3_vinberg_three_plus_sqrt2_extension()
    first_new = extension.prior_wall_count
    cohort = tuple(range(first_new, len(extension.roots)))
    same_height = tuple(
        tuple(extension.wall_gram[i][j] for j in cohort) for i in cohort
    )

    assert cohort == (34, 35, 36)
    assert same_height == (
        (Silver(2), -3 - 3 * SQRT2, -SQRT2),
        (-3 - 3 * SQRT2, Silver(2), -3 - 3 * SQRT2),
        (-SQRT2, -3 - 3 * SQRT2, Silver(2)),
    )
    assert all(
        extension.wall_gram[i][j].sign(0) <= 0
        for i in cohort
        for j in range(i)
    )


def test_targeted_edge_does_not_gain_its_second_finite_endpoint() -> None:
    extension = run_silver_k3_vinberg_three_plus_sqrt2_extension()

    assert extension.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
    assert extension.finite_endpoint_root_indices == (33,)
    assert extension.new_finite_endpoint_root_indices == ()
    assert extension.new_root_edge_pairings == (
        (
            SILVER_ZERO,
            SILVER_ZERO,
            SILVER_ZERO,
            SILVER_ZERO,
            SILVER_ZERO,
            -1 - SQRT2,
            SILVER_ZERO,
            -2 - SQRT2,
        ),
        (
            SILVER_ZERO,
            SILVER_ZERO,
            SILVER_ZERO,
            Silver(-1),
            Silver(-1),
            SILVER_ZERO,
            SILVER_ZERO,
            -8 - 5 * SQRT2,
        ),
        (
            SILVER_ZERO,
            SILVER_ZERO,
            -1 - SQRT2,
            SILVER_ZERO,
            SILVER_ZERO,
            SILVER_ZERO,
            SILVER_ZERO,
            -2 - SQRT2,
        ),
    )
    assert extension.new_root_edge_extension_identity_inertias == (
        Inertia(8, 1, 0),
    ) * 3
    assert extension.second_endpoint_status == "not_found_through_this_shell"
    assert extension.stopping_status == "undecidable_incomplete_simple_roots"


def test_reverse_coordinate_replay_matches_production_root_sets() -> None:
    audit = audit_reverse_silver_shell()

    assert audit.passes
    assert audit.exact_weight_solutions == 12
    assert audit.nonintegral_reconstructions == 0
    assert audit.rejected_root_count == 9
    assert len(audit.accepted_roots) == 3
    assert audit.production_accepted_root_set_matches
    assert audit.production_rejection_witness_set_matches
    assert audit.replay_sha256 == EXPECTED_REPLAY_SHA256


def test_three_plus_sqrt2_certificate_is_frozen() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))

    assert verify_manifest(manifest) == {
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
        "exact_weight_solutions": 12,
        "accepted_root_count": 3,
        "wall_count": 37,
        "known_finite_endpoint_root_indices": [33],
        "independent_replay_sha256": EXPECTED_REPLAY_SHA256,
        "next_admissible_control_coordinate": {"a": "2", "b": "2"},
        "stopping_status": "undecidable_incomplete_simple_roots",
    }


def test_three_plus_sqrt2_certificate_rejects_mutation() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    manifest["new_level"]["accepted_root_indices"] = [34, 35]

    with pytest.raises(ValueError, match="does not reproduce"):
        verify_manifest(manifest)

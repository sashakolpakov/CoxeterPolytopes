import json

import pytest

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_k3_vinberg import (
    run_silver_k3_vinberg_one_plus_two_sqrt2_extension,
)
from scripts.certify_silver_k3_vinberg_1_plus_2sqrt2 import (
    ARCHIVE,
    EXPECTED_CERTIFICATE_SHA256,
    verify_manifest,
)


def test_one_plus_two_sqrt2_shell_is_exact_and_empty() -> None:
    extension = run_silver_k3_vinberg_one_plus_two_sqrt2_extension()
    level = extension.new_level

    assert extension.passes
    assert extension.completed_through_control_coordinate == 1 + 2 * SQRT2
    assert extension.next_admissible_control_coordinate == 3 + SQRT2
    assert level.target == 19 + 13 * SQRT2
    assert level.alphabet.values == (
        Silver(0),
        Silver(1),
        1 + SQRT2,
        2 + SQRT2,
        2 + 2 * SQRT2,
        3 + 2 * SQRT2,
        4 + 3 * SQRT2,
    )
    assert level.exact_weight_solutions == 4
    assert level.nonintegral_reconstructions == 0
    assert level.previous_wall_rejections == 4
    assert level.roots == ()
    assert tuple(
        blockers for _root, blockers in level.previous_wall_rejection_witnesses
    ) == ((9,), (9,), (9,), (9,))


def test_targeted_edge_still_has_only_one_certified_finite_endpoint() -> None:
    extension = run_silver_k3_vinberg_one_plus_two_sqrt2_extension()

    assert extension.prior_wall_count == 34
    assert len(extension.roots) == 34
    assert extension.wall_identity_inertia == Inertia(9, 1, 24)
    assert extension.wall_conjugate_inertia == Inertia(10, 0, 24)
    assert extension.finite_endpoint_root_indices == (33,)
    assert extension.new_finite_endpoint_root_indices == ()
    assert extension.ideal_endpoint_audit_status == "not_run_until_complete_simple_roots"
    assert extension.second_endpoint_status == "not_found_through_this_shell"
    assert extension.stopping_status == "undecidable_incomplete_simple_roots"


def test_one_plus_two_sqrt2_certificate_is_frozen() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))

    assert verify_manifest(manifest) == {
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
        "exact_weight_solutions": 4,
        "accepted_root_count": 0,
        "wall_count": 34,
        "known_finite_endpoint_root_indices": [33],
        "next_admissible_control_coordinate": {"a": "3", "b": "1"},
        "stopping_status": "undecidable_incomplete_simple_roots",
    }


def test_one_plus_two_sqrt2_certificate_rejects_mutation() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    manifest["new_level"]["accepted_root_count"] = 1

    with pytest.raises(ValueError, match="does not reproduce"):
        verify_manifest(manifest)

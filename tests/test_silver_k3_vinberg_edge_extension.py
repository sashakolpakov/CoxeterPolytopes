import json

import pytest

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import SQRT2
from coxeter_polytopes.search.silver_k3_vinberg import (
    run_silver_k3_vinberg_edge_extension,
)
from scripts.certify_silver_k3_vinberg_edge_extension import (
    ARCHIVE,
    EXPECTED_CERTIFICATE_SHA256,
    verify_manifest,
)


def test_next_silver_shell_is_exact() -> None:
    extension = run_silver_k3_vinberg_edge_extension()

    assert extension.passes
    assert extension.completed_through_control_coordinate == 2 + SQRT2
    assert extension.next_admissible_control_coordinate == 1 + 2 * SQRT2
    assert extension.new_level.exact_weight_solutions == 72
    assert extension.new_level.nonintegral_reconstructions == 43
    assert extension.new_level.previous_wall_rejections == 11
    assert len(extension.new_level.roots) == 18
    assert len(extension.roots) == 34
    assert extension.wall_identity_inertia == Inertia(9, 1, 24)
    assert extension.wall_conjugate_inertia == Inertia(10, 0, 24)


def test_zero_endpoint_edge_gets_one_finite_endpoint() -> None:
    audit = run_silver_k3_vinberg_edge_extension().edge_endpoint

    assert audit.passes
    assert audit.finite_type == "A2 + F4 + B2"
    assert len(audit.allowable_pairing_patterns) == 11
    assert audit.new_endpoint_root_indices == (33,)
    assert audit.endpoint_gram_determinant.a == 4
    assert audit.endpoint_gram_determinant.b == 0
    assert audit.endpoint_identity_inertia == Inertia(9, 0, 0)


def test_silver_edge_extension_certificate_is_frozen() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))

    assert verify_manifest(manifest) == {
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
        "new_wall_count": 18,
        "wall_count": 34,
        "allowable_edge_pairing_patterns": 11,
        "new_endpoint_root_indices": [33],
        "stopping_status": "undecidable_incomplete_simple_roots",
    }


def test_silver_edge_extension_certificate_rejects_mutation() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    manifest["wall_count"] = 35

    with pytest.raises(ValueError, match="does not reproduce"):
        verify_manifest(manifest)

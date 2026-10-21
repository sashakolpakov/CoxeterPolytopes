import json

import pytest

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_k3_vinberg import (
    run_silver_k3_vinberg_extension,
)
from scripts.certify_silver_k3_vinberg_extension import (
    ARCHIVE,
    EXPECTED_CERTIFICATE_SHA256,
    verify_manifest,
)


def test_first_silver_k3_extension_is_exact() -> None:
    extension = run_silver_k3_vinberg_extension()

    assert extension.passes
    assert extension.candidate_coordinates == (
        Silver(1),
        SQRT2,
        Silver(2),
        1 + SQRT2,
    )
    assert extension.new_level.control_coordinate == 1 + SQRT2
    assert extension.new_level.exact_weight_solutions == 6
    assert extension.new_level.nonintegral_reconstructions == 0
    assert extension.new_level.previous_wall_rejections == 2
    assert len(extension.new_level.roots) == 4
    assert len(extension.roots) == 16
    assert extension.wall_identity_inertia == Inertia(9, 1, 6)
    assert extension.wall_conjugate_inertia == Inertia(10, 0, 6)


def test_first_silver_k3_extension_does_not_stop_vinberg() -> None:
    diagram = run_silver_k3_vinberg_extension().diagram

    assert diagram.elliptic_rank_eight_count == 267
    assert diagram.elliptic_rank_nine_count == 41
    assert diagram.finite_endpoint_histogram == ((0, 1), (1, 163), (2, 103))
    assert diagram.open_edge_count == 164
    assert diagram.connected_affine_subdiagrams == ()
    assert diagram.compact_stopping_available is False
    assert diagram.finite_volume_stopping_available is False


def test_first_silver_k3_extension_certificate_is_frozen() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))

    assert verify_manifest(manifest) == {
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
        "new_wall_count": 4,
        "wall_count": 16,
        "open_edge_count": 164,
        "stopping_status": "undecidable_incomplete_simple_roots",
    }


def test_first_silver_k3_extension_certificate_rejects_mutation() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    manifest["wall_count"] = 17

    with pytest.raises(ValueError, match="does not reproduce"):
        verify_manifest(manifest)

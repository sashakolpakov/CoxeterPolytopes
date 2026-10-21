import json

import pytest

from scripts.certify_silver_k3_vinberg_prefix import (
    ARCHIVE,
    EXPECTED_CERTIFICATE_SHA256,
    verify_manifest,
)


def test_silver_k3_vinberg_prefix_certificate_is_frozen() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))

    assert verify_manifest(manifest) == {
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
        "height_zero_root_count": 82,
        "wall_count": 12,
        "stopping_status": "undecidable_incomplete_simple_roots",
    }


def test_silver_k3_vinberg_prefix_certificate_rejects_mutation() -> None:
    manifest = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    manifest["wall_count"] = 13

    with pytest.raises(ValueError, match="does not reproduce"):
        verify_manifest(manifest)

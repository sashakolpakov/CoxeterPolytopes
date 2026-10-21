import json
from pathlib import Path

from scripts.certify_silver_k3_carrier import verify


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "data/candidates/silver_k3_carrier.json"


def test_frozen_silver_k3_carrier_certificate_replays_exactly() -> None:
    manifest = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    result = verify(manifest)

    assert result == {
        "seed_saturation": "E8(-1) + U",
        "carrier_types": [
            "U + E8(-1) + D10(-1)",
            "U + E8(-1)^2 + A1(-1)^2",
        ],
        "k3_complement": [[2, 0], [0, 2]],
        "certificate_sha256": (
            "f26310dc2c9475653f8d64fdabe7ba15abe312953621f5f757fc7fb344dec0b9"
        ),
    }
    assert manifest["c5_saturation"]["unique_stable_nonzero_even_support"] == [
        3,
        5,
        9,
    ]
    assert manifest["c5_saturation"]["extra_saturation_root"]["field_norm"] == {
        "a": "3",
        "b": "-1",
    }
    assert manifest["k3_embedding"]["primitive_embedding_verified"] is True
    assert manifest["scope"]["wall_closure_required"] is True

import json
from pathlib import Path

from scripts.certify_silver_k3_forced_cores import verify


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "data/candidates/silver_k3_forced_cores.json"


def test_frozen_silver_k3_certificate_replays_exactly() -> None:
    manifest = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    result = verify(manifest)

    assert result == {
        "finite_trace_types": ["D4(-1)", "D6(-1)", "E8(-1)"],
        "lanner_orders": [4, 5],
        "certificate_sha256": (
            "8f814d5db8a8073191c708879fe701cb489e69999122b66910e300995b2ee1e0"
        ),
    }
    assert manifest["lanner_enumeration_stages"] == [
        {
            "order": 4,
            "graph_masks": 64,
            "graph_prefilter_survivors": 19,
            "decorated_candidates": 176,
            "labelled_exact_survivors": 12,
            "isomorphism_types": 1,
        },
        {
            "order": 5,
            "graph_masks": 1024,
            "graph_prefilter_survivors": 137,
            "decorated_candidates": 2384,
            "labelled_exact_survivors": 60,
            "isomorphism_types": 1,
        },
    ]
    assert [row["trace_type"] for row in manifest["lanner_blocks"]] == [
        "U + A6(-1)",
        "E8(-1) + U(2)",
    ]
    assert all(row["trace_isometry_verified"] for row in manifest["lanner_blocks"])


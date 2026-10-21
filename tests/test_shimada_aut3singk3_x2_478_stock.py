from __future__ import annotations

import hashlib
import json

from scripts.ingest_shimada_aut3singk3_x2_478_stock import (
    OUTPUT,
    SCHEMA,
    validate_stock,
)


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_complete_x2_wall_and_surface_automorphism_stock_is_frozen() -> None:
    assert hashlib.sha256(OUTPUT.read_bytes()).hexdigest() == (
        "10d1c767d21d7194c6ccb603c5dff29520201954b4a07fd3c5f214d452ea2f74"
    )
    payload = _archive()
    validate_stock(payload)
    assert payload["schema"] == SCHEMA
    assert payload["complete_induced_chamber"]["wall_orbit_sizes"] == [
        36, 12, 40, 90, 30, 30, 120, 120
    ]
    assert payload["complete_induced_chamber"]["wall_count"] == 478
    assert payload["ample_stabilizer"]["abstract_group"] == "PGL_2(F_9)"
    assert payload["ample_stabilizer"]["order"] == 720
    assert payload["ample_stabilizer"]["identity_zero_based_index"] == 599
    assert len(
        payload["ample_stabilizer"][
            "all_integral_20_by_20_matrices_in_supplement_order"
        ]
    ) == 720
    assert payload["surface"] == {
        **payload["surface"],
        "name": "X_2",
        "kind": "complex singular K3 surface",
        "transcendental_gram": [[2, 1], [1, 8]],
        "neron_severi_rank": 20,
        "neron_severi_signature": [1, 19],
        "neron_severi_determinant": -15,
    }


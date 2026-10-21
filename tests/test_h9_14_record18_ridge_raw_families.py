from __future__ import annotations

import json

from scripts.certify_h9_14_record18_a1_ridge_raw_families import (
    DEFAULT_OUTPUT as A1_OUTPUT,
    verify_record as verify_a1_record,
)
from scripts.certify_h9_14_record18_a2_add4_ridge_raw_families import (
    DEFAULT_OUTPUT as A2_OUTPUT,
    REPRESENTATIVES as A2_REPRESENTATIVES,
    verify_record as verify_a2_record,
)


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def test_record18_a1_literal_raw_family_crosscheck() -> None:
    document = _load(A1_OUTPUT)
    verify_a1_record(document)
    assert document["production_independent_literal_set_equality"] is True
    assert document["a1_unsat_representatives"] == [23, 15]
    assert document["a1_target_representative"] == 53
    assert document["a1_target_orbit_count"] == 8
    assert document["a1_target_literal_count"] == 1080
    assert document["normalized_raw_target_count"] == 1080
    assert document["compact_only"] is True
    assert document["finite_volume_noncompact_checked"] is False


def test_record18_all_a2_add4_literal_raw_family_crosschecks() -> None:
    document = _load(A2_OUTPUT)
    verify_a2_record(document)
    assert document["production_independent_literal_set_equality"] is True
    assert document["representative_count"] == 12
    assert document["zero_target_representative_count"] == 10
    assert document["literal_target_count"] == 8
    assert document["target_representatives"] == [[53, 198], [53, 202]]
    assert document["compact_only"] is True
    assert document["finite_volume_noncompact_checked"] is False

    rows = {
        tuple(entry["representative_A4_masks"]): entry
        for entry in document["normalized_raw_target_sets"]
    }
    assert tuple(rows) == A2_REPRESENTATIVES
    assert rows[(53, 198)]["target_count"] == 6
    assert rows[(53, 202)]["target_count"] == 2
    for representative, entry in rows.items():
        if representative not in {(53, 198), (53, 202)}:
            assert entry["target_count"] == 0
            assert entry["target_families"] == []

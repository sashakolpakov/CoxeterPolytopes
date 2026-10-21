from __future__ import annotations

import json

from scripts.certify_h9_14_record30_a0_candidate import ARCHIVE, fresh_record, verify


def test_record30_a0_candidate_archive_and_scope() -> None:
    document = json.loads(ARCHIVE.read_text())
    assert verify(document) == {
        "facet_count": 132,
        "ridge_count": 594,
        "base_multiwedge_facet_count": 36,
        "compact_label_terminal_states": 0,
        "convex_polytopality_proved": False,
        "finite_volume_noncompact_chamber_excluded": False,
    }


def test_record30_a0_candidate_fresh_independent_replay() -> None:
    assert fresh_record() == json.loads(ARCHIVE.read_text())

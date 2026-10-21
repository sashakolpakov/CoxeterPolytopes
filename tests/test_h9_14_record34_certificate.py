import json

from scripts.certify_h9_14_record34 import (
    A0_ARCHIVE,
    EXPECTED_CERTIFICATE_SHA256,
    OUTPUT,
    verify,
)


def test_frozen_record34_solver_independent_exclusion() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(manifest) == {
        "record_index": 34,
        "graph6": "M?????????O?O?_@?",
        "theorem_status": "certified_unsat",
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }
    assert manifest["updated_theorem_partition"] == {
        "independently_certified_impossible": 121,
        "certified_explicit_ordinary_polytopal": 6,
        "theorem_undecided": 15,
        "total_graph_shadows": 142,
    }


def test_record34_a0_archive_has_no_ridge_valid_target() -> None:
    document = json.loads(A0_ARCHIVE.read_text(encoding="utf-8"))
    assert document["totals"] == {
        "genbg_representatives": 1_288_153,
        "active_colourings": 26_598_204,
        "distinct_row_states": 25_501_701,
        "no_small_transversal_states": 3_101_794,
        "target_114_states": 18,
        "pure_target_states": 18,
        "ridge_valid_states": 0,
    }
    assert document["canonical_orbit_count"] == 9
    assert all(
        set(map(int, state["ridge_degree_histogram"])) != {2}
        for state in document["target_states"]
    )

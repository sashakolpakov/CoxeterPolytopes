import json

from scripts.certify_h9_14_record37 import OUTPUT, verify
from scripts.independent_h9_14_record37_audit import replay


EXPECTED_CERTIFICATE_SHA256 = (
    "4da8617bd4dfc77889c3c3f987e7988efd87d6d366d7efd29c7403370536cc5c"
)


def test_frozen_record37_solver_independent_exclusion() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(manifest) == {
        "record_index": 37,
        "graph6": "M?????????O?O?o@?",
        "theorem_status": "certified_unsat",
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }
    assert manifest["updated_theorem_partition"] == {
        "independently_certified_impossible": 120,
        "certified_explicit_ordinary_polytopal": 5,
        "theorem_undecided": 17,
        "total_graph_shadows": 142,
    }


def test_record37_independent_five_blocker_replay() -> None:
    result = replay()
    assert result["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    assert result["raw_target_state_count"] == 12
    assert result["unique_family_count"] == 5
    assert result["passing_family_count"] == 0
    assert {
        tuple(sorted(family["ridge_degree_histogram"].items()))
        for family in result["families"]
    } == {
        ((1, 28), (2, 403), (3, 28)),
        ((1, 24), (2, 411), (3, 24)),
    }
    assert {family["dual_component_count"] for family in result["families"]} == {
        2,
        3,
    }

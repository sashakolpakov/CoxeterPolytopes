import json
from pathlib import Path

from coxeter_polytopes.combinatorics.h9_14_record16 import TARGET_GRAPH6
from coxeter_polytopes.combinatorics.h9_14_record16_positive import (
    audit_record16_positive_branches,
)
from scripts.certify_h9_14_record16 import (
    EXPECTED_CERTIFICATE_SHA256,
    OUTPUT,
    X_ARCHIVE,
    verify,
)
from scripts.independent_h9_14_record16_audit import (
    EXPECTED_SUMMARY_SHA256,
    audit as independent_audit,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"


def _schedules() -> list[dict[str, object]]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    survivors = [
        record for record in source["records"] if record["survives_second_layer"]
    ]
    record = survivors[16]
    assert record["graph6"] == TARGET_GRAPH6
    return record["dehn_sommerville_schedule"]


def test_record16_positive_branches_are_exhaustively_excluded() -> None:
    audit = audit_record16_positive_branches(_schedules(), X_ARCHIVE)
    assert audit.verified is True
    assert audit.support.verified is True
    residual = [
        (
            case.missing_four_count,
            case.smaller_cut_side,
            case.larger_cut_side,
            case.required_y_count,
        )
        for case in audit.cut_cases
        if case.survives_coarse_bound
    ]
    assert residual == [
        (3, 2, 8, 72),
        (3, 3, 7, 62),
        (4, 2, 8, 70),
        (4, 3, 7, 60),
        (4, 4, 6, 54),
    ]
    archive = audit.x_archive
    assert archive.verified is True
    assert dict(archive.stream_counts) == {6: 38, 7: 107, 8: 410}
    assert dict(archive.admissible_counts) == {6: 22, 7: 41, 8: 82}
    assert dict(archive.pairwise_nonisomorphic) == {6: True, 7: True, 8: True}
    target_counts = {
        key: len(value) for key, value in archive.target_graphs.items()
    }
    assert target_counts == {
        (8, 72): 0,
        (7, 62): 1,
        (8, 70): 1,
        (7, 60): 2,
        (6, 54): 1,
    }
    assert archive.target_graphs[(6, 54)][0].graph6 == "EFz_"
    assert archive.target_graphs[(6, 54)][0].complete_bipartite_3_3 is True
    assert len(audit.forced_a4_masks) == 4
    assert audit.forced_five_shadow_size == 30
    assert audit.forced_additional_minimal_a5_count == 2


def test_record16_final_case_has_two_independent_contradictions() -> None:
    audit = audit_record16_positive_branches(_schedules(), X_ARCHIVE)
    a5 = audit.final_a5_audit
    assert a5.verified is True
    assert a5.candidate_minimal_a5_count == 762
    assert a5.unordered_two_a5_families == 289_941
    assert a5.tau_z_at_least_two_families == 126
    assert a5.z_edge_count_histogram == ((25, 126),)
    assert a5.z_part_size_histogram == (((5, 5), 126),)
    blocker = audit.full_blocker_audit
    assert blocker.verified is True
    assert blocker.facet_complement_count == 102
    assert blocker.histogram == ((2, 5), (4, 4), (5, 2), (6, 2), (7, 2))


def test_independent_and_frozen_record16_certificates_reproduce() -> None:
    hostile = independent_audit()
    assert hostile["summary_sha256"] == EXPECTED_SUMMARY_SHA256
    assert hostile["final_case"]["contradiction"] == "K5,5 versus K4,6"
    assert hostile["final_case"]["direct_facet_count"] == 102
    assert hostile["final_case"]["direct_minimal_transversal_histogram"] == {
        "2": 5,
        "4": 4,
        "5": 2,
        "6": 2,
        "7": 2,
    }

    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert {source["path"] for source in manifest["proof_sources"]} == {
        "scripts/certify_h9_14_record16_a0.py",
        "scripts/enumerate_h9_14_record16_x_graphs.py",
        "scripts/independent_h9_14_record16_audit.py",
        "src/coxeter_polytopes/combinatorics/h9_14_record16.py",
        "src/coxeter_polytopes/combinatorics/h9_14_record16_positive.py",
        "src/coxeter_polytopes/combinatorics/h9_14_record92.py",
        "src/coxeter_polytopes/combinatorics/local_link.py",
        "src/coxeter_polytopes/combinatorics/missing_edges.py",
        "src/coxeter_polytopes/combinatorics/nerve_support.py",
    }
    assert verify(manifest) == {
        "record_index": 16,
        "graph6": TARGET_GRAPH6,
        "theorem_status": "certified_unsat",
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
        "updated_theorem_partition": {
            "independently_certified_impossible": 122,
            "certified_explicit_ordinary_polytopal": 7,
            "theorem_undecided": 13,
            "total_graph_shadows": 142,
        },
    }
    assert manifest["backend_corroboration"]["trusted_for_proof"] is False

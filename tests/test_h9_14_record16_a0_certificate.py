import json
from pathlib import Path

from coxeter_polytopes.combinatorics.h9_14_record16 import (
    TARGET_GRAPH6,
    audit_record16_a0,
)
from scripts.certify_h9_14_record16_a0 import (
    EXPECTED_CERTIFICATE_SHA256,
    OUTPUT,
    verify,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"


def _source_record() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    survivors = [
        record for record in source["records"] if record["survives_second_layer"]
    ]
    return survivors[16]


def test_record16_a0_symbolic_ridge_reduction_is_complete() -> None:
    record = _source_record()
    assert record["graph6"] == TARGET_GRAPH6
    audit = audit_record16_a0(
        record["graph6"], record["dehn_sommerville_schedule"]
    )
    assert audit.verified is True
    assert [
        tuple(item) for item in audit.to_record()["minimal_graph_cover_vertices"]
    ] == [(0, 1, 12), (0, 1, 13), (12, 13)]
    assert audit.base_vertex_link_g3 == (
        4,
        4,
        6,
        6,
        6,
        6,
        6,
        6,
        6,
        6,
        6,
        6,
        1,
        1,
    )
    assert audit.base_vertex_link_g4 == (
        5,
        5,
        7,
        7,
        7,
        7,
        7,
        7,
        7,
        7,
        7,
        7,
        1,
        1,
    )
    assert audit.ridge_system.to_record() == {
        "total_six_sets": 3_003,
        "nonempty_ridge_groups": 735,
        "two_special_groups": 495,
        "two_special_groups_with_doubled_z": 45,
        "one_special_three_isolate_groups": 240,
        "y_facet_six_set_incidences": 1_980,
        "z_facet_six_set_incidences": 810,
        "verified": True,
    }
    assert audit.blocker_types.to_record() == {
        "graph_cover_five_sets": 310,
        "y_type_five_sets": 220,
        "forbidden_y_01i_five_sets": 10,
        "z_type_five_sets": 90,
        "paired_z_variables": 45,
        "every_graph_cover_classified_once": True,
        "verified": True,
    }
    assert audit.z_cut_replay == {
        "seed_edge_count": 9,
        "seed_count": 512,
        "distinct_solution_count": 512,
        "all_solutions_complete_bipartite": True,
        "edge_count_histogram": {
            "0": 1,
            "9": 10,
            "16": 45,
            "21": 120,
            "24": 210,
            "25": 126,
        },
    }


def test_record16_a0_every_cut_case_has_a_strict_gap() -> None:
    record = _source_record()
    audit = audit_record16_a0(
        record["graph6"], record["dehn_sommerville_schedule"]
    )
    assert [case.z_edge_count for case in audit.cases] == [0, 9, 16, 21, 24, 25]
    assert [case.required_y_count for case in audit.cases] == [
        110,
        92,
        78,
        68,
        62,
        60,
    ]
    assert [case.cauchy_y_upper_bound for case in audit.cases] == [
        90,
        81,
        72,
        63,
        54,
        45,
    ]
    assert [case.contradiction_gap for case in audit.cases] == [
        20,
        11,
        6,
        5,
        8,
        15,
    ]
    assert all(case.excluded for case in audit.cases)


def test_frozen_record16_a0_certificate_reproduces() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(manifest) == {
        "record_index": 16,
        "graph6": TARGET_GRAPH6,
        "theorem_status": "certified_branch_unsat",
        "excluded_missing_four_counts": [0, 5, 6],
        "remaining_missing_four_counts": [1, 2, 3, 4],
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }
    assert manifest["backend_corroboration"]["trusted_for_proof"] is False
    assert manifest["frontier_effect"]["graph_shadow_partition_changed"] is False
    assert manifest["frontier_effect"]["current_global_partition"] == {
        "independently_certified_impossible": 121,
        "certified_explicit_ordinary_polytopal": 7,
        "theorem_undecided": 14,
        "total_graph_shadows": 142,
    }

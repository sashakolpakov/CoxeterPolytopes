import json
from pathlib import Path

from coxeter_polytopes.combinatorics.h9_14_record46 import (
    TARGET_GRAPH6,
    audit_record46_structure,
)
from scripts.certify_h9_14_record46 import EXPECTED_HISTOGRAM, verify


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "enumerations" / "h9_14_nerve_support.json"
CERTIFICATE = (
    ROOT / "data" / "enumerations" / "h9_14_record46_certificate.json"
)


def _source_record() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    survivors = [
        record for record in source["records"] if record["survives_second_layer"]
    ]
    return survivors[46]


def test_record46_has_exact_four_a5_blocker_certificate() -> None:
    record = _source_record()
    assert record["graph6"] == TARGET_GRAPH6
    audit = audit_record46_structure(
        record["graph6"], record["dehn_sommerville_schedule"]
    )
    assert audit.verified is True
    assert audit.schedule_link_g4_sums == (
        (0, 24),
        (1, 18),
        (2, 12),
        (3, 6),
    )
    assert audit.excluded_high_a_branches == ((2, 12, 14), (3, 6, 14))
    assert [item.active_vertices for item in audit.a4_types] == [(), (1,), (11,)]
    assert [item.multiplicity for item in audit.a4_types] == [126, 84, 84]
    assert [item.shadow_size for item in audit.a4_types] == [10, 8, 9]
    assert all(item.impossible for item in audit.a4_types)
    assert audit.a0_forbidden_vertices == (0, 12, 13)
    assert audit.a0_individually_allowed_a5_count == 378
    assert audit.a0_missing_five_count == 4
    assert audit.a0_required_facet_count == 90

    manifest = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    proof = verify(manifest, recompute=True)
    assert proof == audit.to_record()
    assert manifest["certificate_sha256"] == (
        "5f28acdd31ae4c043e799c4dfb62c63dd4b22f95af7032c6b8af83a5fa04b3f3"
    )
    assert manifest["full_census_sha256"] == (
        "2cbafdc1d17c6690aa1aa92b2c22d943245f23053c203ffc43b64070097a6af0"
    )
    assert manifest["theorem_status"] == "certified_unsat"

    census = manifest["a0_blocker_census"]
    assert census["active_assignment_row_count_cells"] == 15
    assert census["raw_active_isolate_joins"] == 55_458
    assert census["full_local_states_before_row_quotient"] == 55_458
    assert census["distinct_row_states_before_row_quotient"] == 53_902
    assert census["canonical_local_states"] == 2_523
    assert census["blocker_valid_states"] == 330
    histogram = {
        int(key): value
        for key, value in census["facet_count_histogram"].items()
    }
    assert histogram == EXPECTED_HISTOGRAM
    assert sum(census["facet_count_histogram"].values()) == 330
    assert max(map(int, census["facet_count_histogram"])) == 83
    assert census["required_facet_count"] == 90
    assert census["target_state_count"] == 0
    assert census["target_states"] == []

    assert manifest["backend_corroboration"]["status"] == "unknown"
    assert manifest["backend_corroboration"]["trusted_for_proof"] is False
    assert manifest["updated_theorem_partition"] == {
        "independently_certified_impossible": 119,
        "certified_explicit_ordinary_polytopal": 2,
        "theorem_undecided": 21,
        "total_graph_shadows": 142,
    }

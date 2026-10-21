import json
from pathlib import Path

from coxeter_polytopes.combinatorics.h9_14_record34 import (
    TARGET_GRAPH6,
    TARGET_RECORD_INDEX,
    audit_record34_positive_branches,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "enumerations" / "h9_14_nerve_support.json"


def _target_record() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    survivors = [
        record for record in source["records"] if record["survives_second_layer"]
    ]
    return survivors[TARGET_RECORD_INDEX]


def test_record34_all_positive_a_branches_are_excluded_exactly() -> None:
    record = _target_record()
    assert record["graph6"] == TARGET_GRAPH6
    audit = audit_record34_positive_branches(
        record["graph6"], record["dehn_sommerville_schedule"]
    )
    assert audit.verified
    assert audit.summed_link_excluded_branches == ((4, 12, 14), (5, 6, 14))
    assert min(dict(audit.a2_forced_shadow_histogram)) > (
        audit.a2_scheduled_missing_fives
    )
    assert min(dict(audit.a3_forced_shadow_histogram)) > (
        audit.a3_scheduled_missing_fives
    )
    assert dict(audit.a1_facet_count_histogram) == {
        90: 480,
        96: 600,
        102: 490,
        108: 75,
    }
    assert audit.a1_required_facet_count not in dict(
        audit.a1_facet_count_histogram
    )

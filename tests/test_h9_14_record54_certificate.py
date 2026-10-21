import json
from pathlib import Path

from coxeter_polytopes.combinatorics.h9_14_record54 import (
    EXPECTED_A0_CANONICAL_ORBIT_HISTOGRAM,
    EXPECTED_A0_HISTOGRAM,
    EXPECTED_A1_AGGREGATE_HISTOGRAM,
    TARGET_GRAPH6,
    audit_record54,
)
from scripts.certify_h9_14_record54 import verify


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "enumerations" / "h9_14_nerve_support.json"
CERTIFICATE = (
    ROOT / "data" / "enumerations" / "h9_14_record54_certificate.json"
)


def _source_record() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    survivors = [
        record for record in source["records"] if record["survives_second_layer"]
    ]
    return survivors[54]


def test_record54_has_solver_independent_blocker_certificate() -> None:
    record = _source_record()
    assert record["graph6"] == TARGET_GRAPH6
    audit = audit_record54(
        record["graph6"], record["dehn_sommerville_schedule"]
    )
    assert audit.verified is True
    assert audit.schedule_link_g4_sums == ((0, 24), (1, 18), (2, 12), (3, 6))
    assert audit.excluded_high_a_branches == ((2, 12, 14), (3, 6, 14))

    assert audit.a0.active_mask_count == 21
    assert audit.a0.active_part_multiset_count == 235
    assert audit.a0.column_distribution_signature_count == 35
    assert audit.a0.compressed_state_count == 1_422
    assert audit.a0.canonical_orbit_count == 159
    assert dict(audit.a0.canonical_orbit_histogram) == (
        EXPECTED_A0_CANONICAL_ORBIT_HISTOGRAM
    )
    assert audit.a0.no_small_transversal_family_count == 2_229_164
    assert dict(audit.a0.five_transversal_histogram) == EXPECTED_A0_HISTOGRAM
    assert 96 not in dict(audit.a0.five_transversal_histogram)

    assert [orbit.multiplicity for orbit in audit.missing_four_orbits] == [
        70,
        112,
        224,
        32,
        224,
        56,
    ]
    assert [item.missing_four_orbit_size for item in audit.a1_types] == [224, 32]
    assert all(item.no_small_transversal_triples == 840 for item in audit.a1_types)
    assert all(
        item.active_free_no_small_transversal_triples == 840
        for item in audit.a1_types
    )
    assert dict(audit.a1_aggregate_histogram) == EXPECTED_A1_AGGREGATE_HISTOGRAM
    assert 94 not in dict(audit.a1_aggregate_histogram)

    manifest = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    proof = verify(manifest)
    assert proof == audit.to_record()
    assert manifest["theorem_status"] == "certified_unsat"
    assert manifest["backend_corroboration"]["trusted_for_proof"] is False
    assert manifest["updated_theorem_partition"] == {
        "independently_certified_impossible": 116,
        "certified_explicit_witness": 2,
        "theorem_undecided": 24,
        "total_graph_shadows": 142,
    }

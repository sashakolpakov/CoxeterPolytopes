import json
from pathlib import Path

from coxeter_polytopes.combinatorics.h9_14_record93 import (
    TARGET_GRAPH6,
    audit_record93,
)
from scripts.certify_h9_14_record93 import verify


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "enumerations" / "h9_14_nerve_support.json"
CERTIFICATE = (
    ROOT / "data" / "enumerations" / "h9_14_record93_certificate.json"
)


def _source_record() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    survivors = [
        record for record in source["records"] if record["survives_second_layer"]
    ]
    return survivors[93]


def test_record93_has_solver_independent_blocker_certificate() -> None:
    record = _source_record()
    assert record["graph6"] == TARGET_GRAPH6
    audit = audit_record93(
        record["graph6"], record["dehn_sommerville_schedule"]
    )
    assert audit.verified is True
    assert audit.schedule_link_g4_sums == ((0, 18), (1, 12), (2, 6))
    assert audit.a0_missing_five_count == 3
    assert audit.a0_required_facet_count == 78
    assert audit.a0_pairwise_pattern_compatible_triples == 12_605_586
    assert audit.a0_no_small_transversal_triples == 32_256
    assert dict(audit.a0_five_transversal_histogram) == {
        54: 7_560,
        57: 22_680,
        63: 1_260,
        75: 756,
    }
    assert {
        (pairwise, triple, count): multiplicity
        for pairwise, triple, count, multiplicity in (
            audit.a0_zero_pattern_signature_histogram
        )
    } == {
        ((0, 0, 3), 0, 63): 1_260,
        ((0, 1, 2), 0, 57): 15_120,
        ((1, 1, 1), 0, 54): 7_560,
        ((1, 1, 2), 1, 57): 7_560,
    }
    assert sum(family.multiplicity for family in audit.a0_families) == 32_256

    manifest = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    proof = verify(manifest)
    assert proof == audit.to_record()
    assert manifest["theorem_status"] == "certified_unsat"
    assert manifest["backend_corroboration"]["trusted_for_proof"] is False
    assert manifest["updated_theorem_partition"] == {
        "independently_certified_impossible": 115,
        "certified_explicit_witness": 2,
        "theorem_undecided": 25,
        "total_graph_shadows": 142,
    }


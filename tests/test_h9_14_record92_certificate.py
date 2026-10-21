import json
from pathlib import Path

from coxeter_polytopes.combinatorics.h9_14_record92 import (
    TARGET_GRAPH6,
    audit_record92,
)
from scripts.certify_h9_14_record92 import verify


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "enumerations" / "h9_14_nerve_support.json"
CERTIFICATE = (
    ROOT / "data" / "enumerations" / "h9_14_record92_certificate.json"
)


def _source_record() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    survivors = [
        record for record in source["records"] if record["survives_second_layer"]
    ]
    return survivors[92]


def test_record92_has_solver_independent_blocker_certificate() -> None:
    record = _source_record()
    assert record["graph6"] == TARGET_GRAPH6
    audit = audit_record92(
        record["graph6"], record["dehn_sommerville_schedule"]
    )
    assert audit.verified is True
    assert audit.schedule_link_g4_sums == ((0, 18), (1, 12), (2, 6))
    assert audit.a0_missing_five_count == 3
    assert audit.a0_required_facet_count == 78
    assert audit.a0_pairwise_pattern_compatible_triples == 6_773_256
    assert audit.a0_no_small_transversal_triples == 1_176
    assert dict(audit.a0_five_transversal_histogram) == {57: 840, 75: 336}
    assert {
        (family.pattern_multiset, family.five_transversal_count): (
            family.multiplicity
        )
        for family in audit.a0_families
    } == {
        ((0, 0, 0), 57): 840,
        ((0, 3, 4), 75): 336,
    }

    manifest = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    proof = verify(manifest)
    assert proof == audit.to_record()
    assert manifest["theorem_status"] == "certified_unsat"
    assert manifest["proof_backend"] == (
        "finite integer enumeration; no SAT/SMT backend"
    )
    assert manifest["backend_corroboration"]["trusted_for_proof"] is False
    assert manifest["updated_theorem_partition"] == {
        "independently_certified_impossible": 114,
        "certified_explicit_witness": 2,
        "theorem_undecided": 26,
        "total_graph_shadows": 142,
    }


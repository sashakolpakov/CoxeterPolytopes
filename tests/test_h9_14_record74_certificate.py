import json
from pathlib import Path

from coxeter_polytopes.combinatorics.h9_14_record74 import (
    ACTIVE_VERTICES,
    TARGET_GRAPH6,
    audit_record74,
)
from scripts.certify_h9_14_record74 import verify


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "enumerations" / "h9_14_nerve_support.json"
CERTIFICATE = (
    ROOT / "data" / "enumerations" / "h9_14_record74_certificate.json"
)


def _source_record() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    survivors = [
        record for record in source["records"] if record["survives_second_layer"]
    ]
    return survivors[74]


def test_record74_has_exact_local_divisibility_certificate() -> None:
    record = _source_record()
    assert record["graph6"] == TARGET_GRAPH6
    audit = audit_record74(
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
    assert len(audit.minimal_graph_covers) == len(ACTIVE_VERTICES) == 4

    assert [item.orbit_multiplicity for item in audit.a4_types] == [210, 480]
    assert [item.independent_five_shadow_size for item in audit.a4_types] == [
        10,
        7,
    ]
    assert audit.a4_types[0].strict_g4_zero_vertices_before_additional_a5
    assert audit.a4_types[1].strict_g3_zero_vertices
    assert all(item.locally_impossible for item in audit.a4_types)

    assert audit.a0_missing_five_count == 4
    assert audit.a0_locally_allowed_missing_five_count == 252
    assert audit.a0_forbidden_active_missing_five_count == 840
    assert audit.a0_every_active_a5_forces_nonpositive_g4 is True
    assert audit.a0_three_active_two_isolate_cover_count == 180
    assert audit.a0_four_active_one_isolate_cover_count == 10
    assert audit.a0_three_active_class_multiplicity == 4
    assert audit.a0_four_active_survivors_under_empty_intersection == 0
    assert audit.a0_required_facet_count == 90
    assert audit.a0_required_facet_count % audit.a0_facet_divisor == 2

    manifest = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    proof = verify(manifest)
    assert proof == audit.to_record()
    assert manifest["certificate_sha256"] == (
        "0a8812c4dbade45cabe3c544c968d265581e4d3f0fdd6b48673ff7847148ceb4"
    )
    assert manifest["theorem_status"] == "certified_unsat"
    assert manifest["backend_corroboration"]["status"] == "unknown"
    assert manifest["backend_corroboration"]["trusted_for_proof"] is False
    assert manifest["updated_theorem_partition"] == {
        "independently_certified_impossible": 118,
        "certified_explicit_witness": 2,
        "theorem_undecided": 22,
        "total_graph_shadows": 142,
    }

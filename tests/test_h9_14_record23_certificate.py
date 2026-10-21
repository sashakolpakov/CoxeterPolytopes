import json
from pathlib import Path

from coxeter_polytopes.combinatorics.h9_14_record23 import (
    TARGET_GRAPH6,
    audit_record23_structure,
)
from scripts.certify_h9_14_record23 import verify


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "enumerations" / "h9_14_nerve_support.json"
CERTIFICATE = (
    ROOT / "data" / "enumerations" / "h9_14_record23_certificate.json"
)


def _source_record() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    survivors = [
        record for record in source["records"] if record["survives_second_layer"]
    ]
    return survivors[23]


def test_record23_has_solver_independent_incidence_certificate() -> None:
    record = _source_record()
    assert record["graph6"] == TARGET_GRAPH6
    structure = audit_record23_structure(
        record["graph6"], record["dehn_sommerville_schedule"]
    )
    assert structure.verified is True
    assert structure.schedule_link_g4_sums == (
        (0, 36),
        (1, 30),
        (2, 24),
        (3, 18),
        (4, 12),
        (5, 6),
    )
    assert [layer.all_orbit_count for layer in structure.a4_orbit_layers] == [
        1,
        6,
        72,
        1_357,
    ]
    assert [
        layer.vertex_interval_feasible_orbit_count
        for layer in structure.a4_orbit_layers
    ] == [1, 3, 23, 86]

    manifest = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    proof = verify(manifest, recompute=False)
    assert proof == structure.to_record()
    assert manifest["certificate_sha256"] == (
        "7b13a6b215c4cf6c2c871361c984f150f044a7107bce5b5fb9f484096360da5c"
    )
    assert manifest["theorem_status"] == "certified_unsat"
    assert manifest["backend_corroboration"]["trusted_for_proof"] is False

    sources = {
        item["path"]: item["sha256"]
        for item in manifest["enumerator_provenance"]["sources"]
    }
    assert sources == {
        "scripts/certify_h9_14_record23_a0.cpp": (
            "b376f5d00e70d6c7e11138324e05b109f1d89ffb095fdfd92e8d1bffb7840a50"
        ),
        "scripts/certify_h9_14_record23_positive.cpp": (
            "ba3ef1a4bbaadb680f6de96e98e610d9fa3cb312804fce04f51131b9d3ab36f5"
        ),
        "scripts/independent_record23_positive_a.cpp": (
            "fb674b757d97c68d2522f30b3c93845312450e1126635d36bc8952331ad8aac1"
        ),
        "scripts/independent_record23_a1_audit.cpp": (
            "0fb98ae8fba3bf495507f5e5006ed7dcdf42af35d837ef661fd97a14e01fe55a"
        ),
    }

    censuses = manifest["censuses"]
    assert {
        name: census["raw_output_sha256"]
        for name, census in censuses.items()
    } == {
        "a0_exact_genbg": (
            "bfd872ff8d7edd07271eabfd8575b035793d393ba4ae3843c96a776ecdad0d68"
        ),
        "a1_exact_genbg": (
            "8f1e607254fe6b03ee46769ac85534d13f0580bec0602a4494161d86bae8b08d"
        ),
        "a1_independent_replay": (
            "99b7df9981c7c16cb441c1eac1fb94137925117c0c661fa1f9c6c7410d2c795e"
        ),
        "positive_a_relaxed_independent": (
            "9203373d6be07b7b10689600b4452592a3850685af974faeba43b887d2202806"
        ),
    }

    a0 = censuses["a0_exact_genbg"]
    assert a0["generated_graph_count"] == 2_807_714
    assert a0["distinct_row_graph_count"] == 2_765_889
    assert a0["locally_valid_left_choice_count"] == 11_456_748
    assert a0["maximum_facet_count"] == 94
    assert a0["target_facet_count_occurrences"] == 0
    assert "98" not in a0["facet_count_histogram"]

    a1 = censuses["a1_exact_genbg"]
    assert a1["generated_graph_count"] == 75_890
    assert a1["schedule_compatible_graph_count"] == 73_919
    assert a1["blocker_valid_left_choice_count"] == 86_241
    assert a1["maximum_facet_count"] == 90
    assert a1["target_facet_count_occurrences"] == 0
    assert a1["facet_count_histogram"] == (
        censuses["a1_independent_replay"]["facet_count_histogram"]
    )

    relaxed = censuses["positive_a_relaxed_independent"]
    assert relaxed["a1_maximum_facet_count"] == 92
    assert relaxed["a1_target_facet_count_occurrences"] == 0
    later = {
        item["missing_four_count"]: item for item in relaxed["later_branches"]
    }
    assert later[2]["a4_orbit_count"] == 17
    assert later[2]["valid_family_count"] == 16_778
    assert later[2]["target_facet_count_occurrences"] == 0
    assert "94" not in later[2]["facet_count_histogram"]
    assert later[3]["a4_orbit_count"] == 144
    assert later[3]["valid_family_count"] == 1_935
    assert later[3]["target_facet_count_occurrences"] == 0
    assert "92" not in later[3]["facet_count_histogram"]

    assert manifest["updated_theorem_partition"] == {
        "independently_certified_impossible": 117,
        "certified_explicit_witness": 2,
        "theorem_undecided": 23,
        "total_graph_shadows": 142,
    }

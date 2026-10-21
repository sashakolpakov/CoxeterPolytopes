import json

from scripts.ingest_generic_hessian_quartic_84_stock import OUTPUT, build_payload


def test_generic_hessian_84_wall_archive_replays_exactly() -> None:
    archived = json.loads(OUTPUT.read_text())
    rebuilt = build_payload()
    assert rebuilt == archived
    assert rebuilt["all_checks_pass"]
    assert rebuilt["completeness"]["wall_count"] == 84
    assert rebuilt["exact_ingestion"]["wall_gram_rank"] == 16
    assert rebuilt["symmetry_census"]["source_exact_group"] == "C2 x S5"
    assert rebuilt["symmetry_census"]["gap_subgroup_conjugacy_class_count"] == 57
    assert len(rebuilt["symmetry_census"]["proper_element_classes_fixed_rank_at_least_11"]) == 1


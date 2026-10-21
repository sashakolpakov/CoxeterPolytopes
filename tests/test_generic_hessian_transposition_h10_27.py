import json

from scripts.certify_generic_hessian_transposition_h10_27 import OUTPUT, build_payload


def test_generic_hessian_transposition_h10_27_archive_replays_exactly() -> None:
    archived = json.loads(OUTPUT.read_text())
    rebuilt = build_payload()
    assert rebuilt == archived
    assert rebuilt["all_checks_pass"]
    assert rebuilt["literal_vinberg"]["finite_volume"]
    assert rebuilt["literal_vinberg"]["compact"] is False
    assert rebuilt["restricted_chamber"]["negative_genuine_facet_count"] == 27
    assert rebuilt["lineage"]["known_id"] == "D16-D4-C7-H10-27"
    assert rebuilt["lineage"]["new_abstract_polytope"] is False


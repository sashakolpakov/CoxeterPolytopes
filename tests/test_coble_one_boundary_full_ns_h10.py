import json
from pathlib import Path

from scripts.certify_coble_one_boundary_full_ns_h10 import (
    OUTPUT,
    build_certificate,
)


def test_archived_coble_h10_theorem_certificate() -> None:
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"]
    assert payload["cases"]["I"]["exact_double_description"][
        "finite_vertex_count"
    ] == 315
    assert payload["cases"]["I"]["literal_vinberg"]["bad_edge_count"] == 0
    assert payload["cases"]["II"]["exact_double_description"][
        "ideal_vertex_count"
    ] == 207
    assert payload["lineage"]["Type_II"]["entrywise_oriented_gram_equality"]
    assert payload["polytope_first_recognition"][
        "D16_D5_H10_30_input_only_audit"
    ]["successful_reconstructions"][0]["recognized_kondo_type"] == "II"


def test_exact_coble_h10_archive_rebuild() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    # Compare the JSON data model, rather than Python-only distinctions such
    # as integer dictionary keys that JSON canonically serializes as strings.
    rebuilt = json.loads(json.dumps(build_certificate(), sort_keys=True))
    assert rebuilt == archived

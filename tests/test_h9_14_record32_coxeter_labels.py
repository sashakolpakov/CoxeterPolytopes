import json

from scripts.certify_h9_14_record32_coxeter_labels import (
    DEPENDENCIES,
    OUTPUT,
    ROOT,
    verify,
)


def test_frozen_record32_unrestricted_coxeter_no_go() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(document) == {
        "a4_compatible_states": 546,
        "a4_passing_states": 48,
        "sparse_states": 20_160,
        "cross_tests": 403_200,
        "cross_survivors": 0,
        "unrestricted_compact_coxeter_no_go": True,
        "certificate_sha256": (
            "15bf62a37b5f65fdf7118e80be347482d30b0e47b84eeab2dd831f089b5b4412"
        ),
    }
    overlap = document["theorem"]["overlapping_order_four_gate"]
    assert overlap["passing_state_count"] == 48
    assert overlap["distinct_first_support_projection_count"] == 24
    assert overlap["first_support_projection_multiplicity_histogram"] == {"2": 24}


def test_record32_coxeter_manifest_hashes_every_imported_gate() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    archived = document["source"]["imported_dependency_sha256"]
    assert set(archived) == {str(path.relative_to(ROOT)) for path in DEPENDENCIES}
    assert set(archived) == {
        "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py",
        "src/coxeter_polytopes/combinatorics/h9_14_disjoint_lanner_gate.py",
        "src/coxeter_polytopes/combinatorics/h9_14_record32.py",
    }

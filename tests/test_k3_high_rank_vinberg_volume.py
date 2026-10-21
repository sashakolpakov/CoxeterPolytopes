import hashlib
import json

from scripts.certify_k3_high_rank_vinberg_volume import (
    ARCHIVE_PATH,
    DATA_PATH,
    analyze_case,
)
from scripts.verify_k3_high_rank_nef_chambers import load_artifact


def _archive() -> dict:
    return json.loads(ARCHIVE_PATH.read_text())


def test_all_five_complete_k3_diagrams_have_frozen_vinberg_certificates() -> None:
    archived = _archive()
    assert archived["source_sha256"] == hashlib.sha256(DATA_PATH.read_bytes()).hexdigest()
    assert [
        (
            case["rho"],
            case["hyperbolic_dimension"],
            case["facet_count"],
            case["elliptic_edge_count"],
            case["finite_vertex_count"],
            case["ideal_vertex_count"],
            case["connected_affine_audit_count"],
        )
        for case in archived["cases"]
    ] == [
        (15, 14, 21, 2676, 328, 15, 34),
        (16, 15, 19, 1563, 186, 7, 13),
        (17, 16, 19, 1087, 125, 3, 6),
        (18, 17, 19, 729, 80, 2, 3),
        (19, 18, 24, 5160, 525, 10, 22),
    ]
    assert all(case["simple_roots_complete"] for case in archived["cases"])
    assert all(case["decision"] == "finite_volume_noncompact" for case in archived["cases"])
    assert all(case["finite_volume_bad_edge_count"] == 0 for case in archived["cases"])
    assert all(
        case["bad_full_parabolic_neighborhood_count"] == 0
        for case in archived["cases"]
    )
    assert [case["canonical_full_audit_sha256"] for case in archived["cases"]] == [
        "1fd946e0fcdfb09814a397de2e2c755cddf15372e8b1cc34a8591070bcba2566",
        "3a0e32d14a4856536bb3205f8b21a5752eea49b3f8934e443a0487c64b332434",
        "dff839af85d94c42ceb378a664d4576167e99aa842f47d946435800c68ff06bb",
        "6e77968ffeb645ad008cc52bcc4995c9a983d928a88d6a4a34e5dc35ad4d0139",
        "dfead85c337bc4e6e58a2925a4f905949c73f0e271c559333fd48adf6d0715be",
    ]
    assert all(case["all_checks_pass"] for case in archived["cases"])


def test_rank18_case_replays_the_archived_direct_diagram_gate() -> None:
    source_case = next(
        case for case in load_artifact()["cases"] if int(case["rho"]) == 18
    )
    archived_case = next(case for case in _archive()["cases"] if case["rho"] == 18)
    assert analyze_case(source_case) == archived_case

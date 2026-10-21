import hashlib
import json

from scripts.certify_k3_rank11_14_literal_vinberg_volume import (
    CASES,
    CPP_SOURCE,
    INVENTORY,
    OUTPUT,
    _canonical_sha256,
    _reconstruct_ud4d4a1,
    build_certificate,
)


EXPECTED_COUNTS = [
    ("R11_UE8A1", 10, 12, 99, 18, 1, 2),
    ("R11_UD8A1", 10, 14, 218, 35, 4, 9),
    ("R11_UD4D4A1", 10, 22, 1343, 185, 25, 72),
    ("R11_UD4A1x5", 10, 90, 55176, 6490, 603, 2571),
    ("R12_UE8A1x2", 11, 14, 216, 35, 2, 4),
    ("R12_UD8A1x2", 11, 19, 954, 136, 12, 31),
    ("R12_UD4x2A1x2", 11, 59, 25470, 3070, 261, 957),
    ("R12_UA2E8", 11, 13, 162, 27, 1, 2),
    ("R13_UE8A1x3", 12, 17, 647, 93, 5, 12),
    ("R13_UD8A1x3", 12, 39, 12266, 1534, 101, 316),
    ("R13_UE8A3", 12, 14, 234, 35, 2, 3),
    ("R14_UE8D4", 13, 15, 315, 42, 4, 5),
    ("R14_UD8D4", 13, 20, 3832, 436, 61, 96),
    ("R14_UE8A1x4", 13, 27, 5698, 718, 37, 97),
]


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_all_fourteen_literal_vinberg_counts_are_frozen() -> None:
    archive = _archive()
    assert archive["all_checks_pass"]
    assert [
        (
            case["key"],
            case["hyperbolic_dimension"],
            case["facet_count"],
            case["elliptic_edge_count"],
            case["finite_vertex_count"],
            case["ideal_vertex_count"],
            case["connected_affine_component_count"],
        )
        for case in archive["cases"]
    ] == EXPECTED_COUNTS
    assert archive["totals"] == {
        "bad_edge_count": 0,
        "bad_parabolic_neighborhood_count": 0,
        "case_count": 14,
        "connected_affine_component_count": 4177,
        "elliptic_edge_count": 106630,
        "facet_count": 375,
        "finite_vertex_count": 12854,
        "ideal_vertex_count": 1119,
    }
    assert all(case["simple_roots_complete"] for case in archive["cases"])
    assert all(
        case["decision"] == "finite_volume_noncompact"
        for case in archive["cases"]
    )
    assert all(case["bad_edge_count"] == 0 for case in archive["cases"])
    assert all(
        case["bad_parabolic_neighborhood_count"] == 0
        for case in archive["cases"]
    )


def test_rank11_22_standard_lattice_reconstruction_replays() -> None:
    gram, reconstruction = _reconstruct_ud4d4a1()
    archived = _archive()
    case = next(
        record for record in archived["cases"] if record["key"] == "R11_UD4D4A1"
    )
    assert reconstruction == archived["rank11_22_reconstruction"]
    assert _canonical_sha256(gram) == case["curve_gram_sha256"]
    assert [
        (
            level["control"],
            level["dominant_candidate_count"],
            level["accepted_root_count"],
        )
        for level in reconstruction["levels"]
    ] == [(1, 3, 3), (2, 12, 0), (3, 18, 0), (4, 54, 9)]
    diagnostic = archived["rank11_22_source_transcription_diagnostic"]
    assert diagnostic["literal_basis_graph_absolute_determinant"] == 128
    assert diagnostic["printed_dependent_curve_squares"]["A13"] == -14
    assert diagnostic["obstruction_is_terminal"] is False


def test_source_hashes_and_full_fourteen_case_replay() -> None:
    archive = _archive()
    assert len(CASES) == 14
    assert archive["source_artifacts"][str(INVENTORY.relative_to(INVENTORY.parents[2]))] == hashlib.sha256(
        INVENTORY.read_bytes()
    ).hexdigest()
    assert archive["source_artifacts"][str(CPP_SOURCE.relative_to(CPP_SOURCE.parents[1]))] == hashlib.sha256(
        CPP_SOURCE.read_bytes()
    ).hexdigest()
    assert build_certificate() == archive

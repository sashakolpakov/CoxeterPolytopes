import json

from scripts.certify_h9_14_record0_a2_source29_triple012_cut import (
    EXPECTED_MATHEMATICAL_DIGEST,
    EXPECTED_RADIUS_TWO,
    OUTPUT,
    _analytic_cut,
    _near_seed,
    verify_document,
)


def test_frozen_source29_triple012_localization() -> None:
    document = json.loads(OUTPUT.read_text())
    result = verify_document(document)
    assert result == {
        "status": "PASS",
        "source29_status": "OPEN",
        "exact_branch_count": 6,
        "radius_two_replayed": False,
        "peak_rss_bytes": result["peak_rss_bytes"],
    }
    assert document["mathematical_digest"] == EXPECTED_MATHEMATICAL_DIGEST


def test_triple012_cut_has_six_exact_branches() -> None:
    cut = _analytic_cut()
    assert cut["source_stabilizer"]["order"] == 161280
    assert cut["candidate_A5_orbit_count"] == 22
    assert cut["forbidden_containing_triple_orbit"]["size"] == 28
    assert sum(row["size"] for row in cut["surviving_killer_orbits"]) == 168
    assert cut["canonical_killer_masks"] == [227, 229, 230]
    assert cut["normalized_active_subset_orbit_count"] == 5
    assert cut["exact_branch_count"] == 6
    assert cut["branches"]["zero_support"]["undecided_candidates"] == 1404


def test_near_seed_and_radius_two_scope() -> None:
    seed = _near_seed()
    assert seed["scheduled_missing_face_counts"] == {"6": 147, "7": 420, "8": 629}
    assert seed["strict_failure_count"] == 1
    assert seed["strict_failure"] == {
        "face": [1, 2, 3],
        "g": [2, 1, -1],
        "selected_A5_incidence": 2,
    }
    assert seed["compact_incidence_witness"] is False
    assert EXPECTED_RADIUS_TWO["all_ground_one_replacement_shadow_exact"] == 599
    assert EXPECTED_RADIUS_TWO["triple012_forbidden_one_replacement_shadow_exact"] == 3
    assert EXPECTED_RADIUS_TWO["post_cut_one_replacement_shadow_exact"] == 596
    assert EXPECTED_RADIUS_TWO["post_cut_two_replacement_shadow_exact"] == 515354
    assert EXPECTED_RADIUS_TWO["post_cut_two_replacement_strict_pass"] == 0

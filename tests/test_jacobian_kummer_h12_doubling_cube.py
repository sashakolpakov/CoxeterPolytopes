from __future__ import annotations

from functools import lru_cache
import json

from scripts.certify_jacobian_kummer_h12_doubling_cube import OUTPUT, build


@lru_cache(maxsize=1)
def _payload():
    return build()


def test_cube_archive_rebuilds_exactly() -> None:
    payload = _payload()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True


def test_complete_c2_four_census() -> None:
    cube = _payload()["orthogonal_doubling_cube"]
    assert cube["walls_one_based"] == [31, 58, 59, 68]
    assert cube["group"] == "C2^4"
    assert cube["faithful"] is True
    assert cube["complete_nonempty_subset_count"] == 15
    assert cube["facet_count_histogram"] == {
        "296": 3,
        "317": 1,
        "523": 3,
        "559": 3,
        "941": 1,
        "1004": 3,
        "1827": 1,
    }
    assert cube["pairwise_nonisomorphic_type_lower_bound"] == 7


def test_seven_type_representatives_have_strict_facet_witnesses() -> None:
    payload = _payload()
    representatives = [
        row
        for row in payload["constructions"]
        if row["strict_facet_witnesses"].get("validated_for_every_facet")
    ]
    assert [row["facet_count"] for row in representatives] == [
        296,
        317,
        523,
        559,
        941,
        1004,
        1827,
    ]
    assert all(
        row["all_exact_checks_pass"]
        and row["finite_volume"]
        and not row["compact"]
        and row["full_linear_rank"] == 13
        and row["strict_facet_witnesses"]["count"] == row["facet_count"]
        and row["strict_facet_witnesses"]["minimum_off_target_slack"] == 1
        for row in representatives
    )
    assert payload["decision"]["new_to_project_types_beyond_single_doublings"] == 5
    assert payload["decision"]["largest_certified_descendant"] == "H12/1827"


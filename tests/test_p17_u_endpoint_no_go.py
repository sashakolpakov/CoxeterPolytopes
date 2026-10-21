import json
from pathlib import Path

from scripts.certify_p17_u_endpoint_no_go import build_certificate


def test_all_eight_u_endpoints_have_the_forced_six_curve_profile():
    certificate = build_certificate()
    frozen = json.loads(
        Path("data/enumerations/p17_u_endpoint_no_go.json").read_text()
    )
    assert certificate == frozen
    assert certificate["rank_two_state_count"] == 16
    assert certificate["u_endpoint_count"] == 8
    assert certificate["f1_type_endpoint_count"] == 8
    assert certificate["all_u_endpoints_have_one_common_profile"]

    profile = certificate["common_u_endpoint_profile"]
    assert {
        label: tuple(entry["base_class_f_s"])
        for label, entry in profile.items()
    } == {
        "s": (0, 1),
        "f": (1, 0),
        "C": (4, 2),
        "Q": (5, 2),
        "A": (6, 3),
        "S": (8, 4),
    }
    assert tuple(profile["f"]["multiplicities_p1_through_p16"]) == (
        1,
        1,
        *((0,) * 14),
    )
    assert tuple(profile["C"]["multiplicities_p1_through_p16"]) == (
        *((1,) * 10),
        *((0,) * 6),
    )
    assert tuple(profile["A"]["multiplicities_p1_through_p16"]) == (
        2,
        *((1,) * 15),
    )


def test_the_exact_local_order_is_the_forbidden_ruled_parity_gap():
    certificate = build_certificate()
    assert certificate["forced_intersections"] == {
        "C_dot_f": 2,
        "C_dot_A": 12,
        "cluster_local_Ip_C_A": 11,
    }
    obstruction = certificate["geometrically_ruled_obstruction"]
    assert obstruction["required_local_order"] == 11
    assert obstruction["required_base_orders_for_order_at_least_11"] == {
        "P": 6,
        "Q": 5,
    }
    assert obstruction["degree_bounds"] == {"deg_L3": 6, "deg_L2": 4}


from __future__ import annotations

import json

from scripts.certify_jacobian_kummer_h10_reflection_index import (
    OUTPUT,
    _connected_masks_by_nonright_components,
    _connected_masks_direct,
)


def _payload() -> dict[str, object]:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_exact_euler_ratio_and_index() -> None:
    payload = _payload()
    assert payload["all_checks_pass"] is True
    subgroup = payload["subgroup_84_wall_chamber"]
    base = payload["vinberg_f2_h10_13_wall_chamber"]
    assert subgroup["rational_euler_characteristic"] == "-527/11520"
    assert base["rational_euler_characteristic"] == "-31/11147673600"
    assert payload["index"] == {
        "formula": "chi(W_84) / chi(W_f2_H10)",
        "prime_factorization": {"2": 10, "3": 3, "5": 1, "7": 1, "17": 1},
        "value": 16_450_560,
    }


def test_finite_special_counts_recover_vinberg_incidence() -> None:
    payload = _payload()
    counts = payload["subgroup_84_wall_chamber"][
        "finite_special_subset_count_by_rank"
    ]
    assert counts == [
        1,
        84,
        1_350,
        10_650,
        49_643,
        146_599,
        281_244,
        347_634,
        263_403,
        107_233,
        15_698,
    ]
    assert payload["hyperbolic_covolume"]["value"] == "pi^5 * 527/340200"
    assert payload["independent_coxiter_1_1_8_replay"]["f_vector"] == [
        16_560,
        107_233,
        263_403,
        347_634,
        281_244,
        146_599,
        49_643,
        10_650,
        1_350,
        84,
        1,
    ]


def test_componentwise_connected_elliptic_enumerator_matches_direct() -> None:
    # A3, a B2 tail, and an isolated node, with one parallel obstruction.
    relations = [
        ["1", "3", "2", "2", "0", "2"],
        ["3", "1", "3", "2", "2", "2"],
        ["2", "3", "1", "2", "2", "2"],
        ["2", "2", "2", "1", "4", "2"],
        ["0", "2", "2", "4", "1", "2"],
        ["2", "2", "2", "2", "2", "1"],
    ]
    assert _connected_masks_by_nonright_components(relations) == _connected_masks_direct(
        relations
    )

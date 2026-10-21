from __future__ import annotations

import json

from scripts.certify_ss_fermat_char3_c2_h10_mcleod_index import OUTPUT


def _payload() -> dict[str, object]:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_exact_mcleod_lineage_and_index() -> None:
    payload = _payload()
    assert payload["all_checks_pass"] is True
    assert payload["integral_lineage"] == {
        "half_scaled_fixed_lattice": "<3> + <-1>^10",
        "explicit_unimodular_isometry": True,
        "primitive_root_type_histogram": {
            "-6,3": 80,
            "-2,1": 36,
            "-1,1": 80,
        },
        "all_196_reflections_integral": True,
        "reflection_subgroup_inclusion": "W_196 < W(f_3^10)",
    }
    assert payload["index"] == {
        "formula": "chi(W_196) / chi(W_f3_H10)",
        "value": 75_202_560,
        "prime_factorization": {"2": 15, "3": 3, "5": 1, "17": 1},
    }


def test_euler_characteristics_volume_and_face_vector() -> None:
    payload = _payload()
    assert payload["subgroup_196_wall_chamber"][
        "rational_euler_characteristic"
    ] == "-5797/3780"
    assert payload["mcleod_f3_h10_14_wall_chamber"][
        "rational_euler_characteristic"
    ] == "-341/16721510400"
    assert payload["hyperbolic_covolume"]["value"] == "pi^5 * 46376/893025"
    assert payload["full_face_vector"] == [
        223_384,
        1_427_720,
        3_435_180,
        4_355_400,
        3_290_232,
        1_543_260,
        447_360,
        76_260,
        6_680,
        196,
        1,
    ]


def test_lineage_claim_is_bounded_correctly() -> None:
    decision = _payload()["lineage_decision"]
    assert decision["maximal_reflection_group"] == "McLeod f_3^10"
    assert decision["maximal_chamber_facet_count"] == 14
    assert decision["subgroup_index_is_exact"] is True
    assert decision["literal_196_wall_subgroup_priority_remains_open"] is True
    assert decision["new_lattice_or_commensurability_class"] is False

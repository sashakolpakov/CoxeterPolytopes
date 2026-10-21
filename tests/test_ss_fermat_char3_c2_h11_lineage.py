from __future__ import annotations

import json

from scripts.certify_ss_fermat_char3_c2_h11_lineage import OUTPUT


def _payload() -> dict[str, object]:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_c2_h11_lattice_and_alexeev_maximal_chamber() -> None:
    payload = _payload()
    assert payload["all_checks_pass"] is True
    identification = payload["fixed_lattice_identification"]
    assert identification["target_lattice"] == "U + E8(-2) + A1(-2)^2"
    assert identification["nikulin_invariants"] == {"r": 12, "a": 10, "delta": 1}
    assert identification["explicit_isometry_determinant"] == 1
    assert identification["rational_comparison_with_I_1_11"]["rationally_similar"] is True

    maximal = payload["known_maximal_reflection_chamber"]
    assert len(maximal["coxeter_gram"]) == 14
    assert maximal["coxeter_gram_sha256"] == (
        "a12618b59ad911429dc4dfe1a18a1d7b6b64a7f5dae53fcf0c677f28c43eb9e5"
    )
    assert maximal["primitive_square_divisibility_histogram"] == {
        "-4,2": 11,
        "-2,1": 2,
        "-2,2": 1,
    }
    assert maximal["all_reflections_integral_involutive_isometries"] is True


def test_all_206_mirrors_are_integral_reflections_of_the_known_lattice() -> None:
    subgroup = _payload()["induced_206_wall_reflection_subgroup"]
    assert subgroup["facet_count"] == 206
    assert subgroup["primitive_square_divisibility_histogram"] == {
        "-4,2": 30,
        "-2,1": 12,
        "-2,2": 164,
    }
    assert subgroup["all_reflections_integral_involutive_isometries"] is True
    assert len(
        {tuple(row["fixed_basis_coordinates"]) for row in subgroup["mirrors"]}
    ) == 206
    assert all(
        row["reflection_integral_involutive_isometry"]
        for row in subgroup["mirrors"]
    )


def test_bounded_h11_weighted_diagram_audit() -> None:
    payload = _payload()
    bounded = payload["bounded_weighted_diagram_audit"]
    assert bounded["literal_206_facet_record_present"] is False
    snapshot = bounded["coxiter_snapshot"]
    assert snapshot["ordinary_file_count"] == 1_214
    assert snapshot["ordinary_H11_facet_histogram"] == {"12": 14, "13": 6, "15": 2}
    assert snapshot["doubling_file_count"] == 2_959
    assert snapshot["doubling_H11_facet_histogram"] == {
        "15": 2,
        "16": 2,
        "18": 2,
        "22": 2,
        "30": 2,
        "46": 2,
    }
    assert payload["lineage_decision"] == {
        "new_lattice": False,
        "new_maximal_reflection_group": False,
        "new_arithmetic_commensurability_class": False,
        "positive_algebraic_construction": True,
        "subgroup_chamber_priority": "open after bounded corpus exclusion",
    }

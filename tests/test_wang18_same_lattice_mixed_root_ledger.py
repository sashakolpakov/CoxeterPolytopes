import json

from scripts.certify_wang18_same_lattice_mixed_root_ledger import OUTPUT_PATH, analyze


def _rows():
    return {row["wang_number"]: row for row in analyze()["rows"]}


def test_frozen_artifact_rebuilds_exactly() -> None:
    assert json.loads(OUTPUT_PATH.read_text()) == analyze()


def test_all_eighteen_rows_are_decided_and_exactly_seventeen_enlarge() -> None:
    report = analyze()
    assert all(report["checks"].values())
    assert report["status"] == "closed_18_of_18"
    assert report["counts"] == {
        "rows": 18,
        "strict_same_lattice_mixed_reflection_enlargements": 17,
        "proper_integral_overclosures": 15,
        "arithmetic_reflection_maximal_rows": 1,
        "undecided_rows": 0,
    }


def test_exact_mixed_root_witness_table() -> None:
    rows = _rows()
    expected = {
        1: (None, None, None),
        2: ("D6", -4, 2),
        3: ("D4", -4, 2),
        4: ("D4", -4, 2),
        5: ("D8", -4, 2),
        6: ("2A1", -4, 2),
        7: ("2A1", -4, 2),
        8: ("D8", -4, 2),
        9: ("A3", -4, 2),
        10: ("D8", -4, 2),
        11: ("D4", -4, 2),
        12: ("A2", -6, 3),
        13: ("A2", -6, 3),
        14: ("D4", -4, 2),
        15: ("D5", -4, 2),
        16: ("D4", -4, 2),
        17: ("A2", -6, 3),
        18: ("A3", -4, 2),
    }
    observed = {}
    for number, row in rows.items():
        witness = row["same_lattice_mixed_root"]
        observed[number] = (
            (witness["component"], witness["root_norm"], witness["root_divisibility"])
            if witness["exists"]
            else (None, None, None)
        )
        if witness["exists"]:
            replay = witness["exact_coordinate_replay"]
            assert replay["root_coordinate_gcd"] == 1
            assert replay["multiple_is_integral"]
            assert replay["reflection_coefficients_integral"]
            assert witness["outside_W2"]
    assert observed == expected


def test_integral_overclosure_ledger_and_squarefree_exceptions() -> None:
    rows = _rows()
    without = {
        number for number, row in rows.items()
        if not row["proper_integral_overclosure"]["exists"]
    }
    assert without == {1, 12, 17}
    assert {rows[number]["absolute_determinant"] for number in without} == {2, 3, 15}
    for number, row in rows.items():
        closure = row["proper_integral_overclosure"]
        if closure["exists"]:
            assert closure["W2_stable"]
            assert (
                closure["index"] ** 2 * closure["target_absolute_determinant"]
                == row["absolute_determinant"]
            )
            assert closure["exact_coordinate_replay"]["glue_pairings_integral"]
            assert closure["exact_coordinate_replay"]["proper"]


def test_UE8E7_is_the_unique_negative_maximal_row() -> None:
    report = analyze()
    rows = report["rows"]
    maximal = [row for row in rows if row["arithmetic_reflection_maximal"]]
    assert [row["lattice"] for row in maximal] == ["U+E8+E7"]
    discriminant = report["UE8E7_discriminant_check"]
    assert discriminant["E7_determinant"] == 2
    assert discriminant["negative_discriminant_q_mod_2"] == "1/2"
    assert discriminant["mixed_root_required_q_mod_2"] == "1"
    assert discriminant["values_are_distinct"]


def test_global_69_corollary_is_one_way_and_has_one_exception() -> None:
    corollary = analyze()["global_69_corollary"]
    assert (
        corollary["with_Wang_automorphic_correction"]
        + corollary["without_Wang_automorphic_correction"]
        == corollary["complete_norm_two_reflective_lattices_of_rank_at_least_6"]
        == 69
    )
    assert "or S is U+E8(-1)+E7(-1)" in corollary["theorem"]
    assert corollary["converse_claimed"] is False

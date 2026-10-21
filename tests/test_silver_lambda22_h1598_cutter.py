import json
from hashlib import sha256
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/silver_lambda22_h1598_cutter.json"
H738 = ROOT / "data/enumerations/silver_lambda22_h738_cutter.json"
REPORT = ROOT / "results/silver_lambda22_h1598_cutter.md"


def payload():
    return json.loads(ARTIFACT.read_text())


def test_explicit_h1598_root_is_primitive_reflective_and_cuts_old_ray():
    root = payload()["explicit_cutter"]
    assert root["Gram_d_e"] == [[4, 0], [0, 10]]
    assert root["coefficient_gram"] == [[1610, 1130], [1130, 794]]
    assert root["control"]["display"] == "62 + 44*sqrt(2)"
    assert root["norm"]["display"] == "2"
    assert root["height"]["display"] == "1598 + 1130*sqrt(2)"
    assert root["primitive"] is True
    assert root["reflection_integral"] is True
    assert root["old_32d_over_99_ray_pairing"]["display"] == "-14/99 + 10/99*sqrt(2)"
    assert root["old_32d_over_99_ray_pairing"]["signs"][0] > 0


def test_complete_h14_packets_and_prior_cutters_are_strictly_acute():
    acute = payload()["complete_acuteness"]
    assert [
        (
            row["packet"],
            row["PSD_extension_count"],
            row["minimum_four_extension_count"],
            row["maximizer"],
            row["maximum_pairing"]["display"],
        )
        for row in acute["complete_H14_packets"]
    ] == [
        ("ordinary_q2", 201, 165, {"d_dot_w": 5, "e_dot_w": 5}, "-45 - 32*sqrt(2)"),
        ("dyadic_q4_S2", 103, 89, {"d_dot_w": 8, "e_dot_w": 4}, "-28 - 20*sqrt(2)"),
    ]
    assert acute["H164_cutter_pairing"]["display"] == "-34 - 24*sqrt(2)"
    assert acute["H738_cutter_pairing"]["display"] == "-37 - 26*sqrt(2)"
    assert acute["all_strict"] is True


def test_chronological_census_has_only_one_viable_positive_root():
    scan = payload()["chronological_acceptance"]
    assert scan["control_count"] == 284
    assert scan["gate_compatible_post_H14_type_count"] == 1134
    assert scan["nonnegative_row_count"] == len(scan["rows"]) == 32
    extensions = [extension for row in scan["rows"] for extension in row["extensions"]]
    assert scan["nonnegative_extension_count"] == len(extensions) == 34
    positives = [
        (row, extension)
        for row in scan["rows"]
        for extension in row["extensions"]
        if extension["pairing"]["signs"][0] > 0
    ]
    assert scan["positive_extension_count"] == len(positives) == 5
    assert [(row["height"], row["norm"], row["gram"]) for row, _ in positives] == [
        ("21 + 15*sqrt(2)", "6", [60, 36, 45]),
        ("21 + 15*sqrt(2)", "6", [60, 36, 45]),
        ("274 + 194*sqrt(2)", "4", [580, 260, 388]),
        ("1598 + 1130*sqrt(2)", "2", [1610, 794, 1130]),
        ("1598 + 1130*sqrt(2)", "4", [3156, 1620, 2260]),
    ]
    assert "unique realizable positive extension" in scan["conclusion"]


def test_same_height_q4_impostor_has_exact_norm_two_contradiction():
    rows = payload()["chronological_acceptance"]["rows"]
    q4 = next(row for row in rows if row["height"] == "1598 + 1130*sqrt(2)" and row["gram"] == [3156, 1620, 2260])
    extension = q4["extensions"][0]
    assert [extension[key] for key in ("d_dot_u", "d_dot_v", "e_dot_u", "e_dot_v")] == [112, 80, 14, 14]
    assert 3156 + 1620 - 2 * 2260 == 8 * 8 * 4
    assert 1620 + 20 * 20 * 4 - 2 * 20 * 80 == 20
    assert 20 + 10 - 2 * 14 == 2


def test_10d_over_31_is_strict_hyperideal_replacement_through_h1598():
    ray = payload()["replacement_hyperideal_ray"]
    assert ray["scale"] == "10/31"
    assert ray["norm"]["display"] == "1361/961 - 1*sqrt(2)"
    assert ray["norm"]["signs"][0] > 0
    assert ray["control_count"] == 284
    assert ray["retained_necessary_type_count"] == 1138
    assert ray["relaxed_nonnegative_row_count"] == len(ray["relaxed_nonnegative_rows"]) == 10
    assert all(row["resolution"] for row in ray["relaxed_nonnegative_rows"])
    assert [row["height"] for row in ray["relaxed_nonnegative_rows"]][-3:] == [
        "1393 + 985*sqrt(2)",
        "1476 + 1044*sqrt(2)",
        "1598 + 1130*sqrt(2)",
    ]
    assert {
        name: value["display"] for name, value in ray["accepted_wall_pairings"].items()
    } == {
        "H164": "-32/31 + 12/31*sqrt(2)",
        "H738": "-18/31 + 8/31*sqrt(2)",
        "H1598": "-6/31 + 2/31*sqrt(2)",
    }
    assert all(value["signs"][0] < 0 for value in ray["accepted_wall_pairings"].values())


def test_parent_pin_and_bounded_scope_are_explicit():
    data = payload()
    assert data["base_H738_artifact"]["sha256"] == sha256(H738.read_bytes()).hexdigest()
    assert data["discovery"] is None
    assert "bounded prefix only" in data["scope"]
    assert "not finite volume" in data["scope"]
    report = REPORT.read_text()
    assert "third accepted cutter" in report
    assert "still has infinite volume" in report
    assert "not a" in report and "global nonreflectivity" in report

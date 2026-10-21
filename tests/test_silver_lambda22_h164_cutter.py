import json
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/silver_lambda22_h164_cutter.json"


def payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def test_explicit_root_is_primitive_reflective_and_cuts_old_ray() -> None:
    root = payload()["explicit_cutter"]
    assert root["Gram_d_e"] == [[4, -2], [-2, 12]]
    assert root["coefficient_gram"] == [[300, 232], [232, 180]]
    assert root["d_pairings"] == [34, 26]
    assert root["norm"]["display"] == "4"
    assert root["height"]["display"] == "164 + 116*sqrt(2)"
    assert root["primitive"] is True
    assert root["primitive_unit_coordinate"] == "-1 - 1*sqrt(2)"
    assert root["reflection_integral"] is True
    assert root["fixed_ray_pairing"]["display"] == "-2/3 + 2/3*sqrt(2)"
    assert root["fixed_ray_pairing"]["signs"][0] > 0


def test_cutter_is_strictly_acute_to_both_complete_h14_packets() -> None:
    audit = payload()["complete_H14_packet_acuteness"]
    assert audit["both_strict"] is True
    ordinary, dyadic = audit["packets"]
    assert ordinary["max_spatial_pairing"]["display"] == "124 + 86*sqrt(2)"
    assert ordinary["maximizer"] == {"d_dot_w": 5, "e_dot_w": -8}
    assert ordinary["max_root_pairing"]["display"] == "-12 - 10*sqrt(2)"
    assert dyadic["max_spatial_pairing"]["display"] == "190 + 132*sqrt(2)"
    assert dyadic["maximizer"] == {"d_dot_w": 8, "e_dot_w": -10}
    assert dyadic["max_root_pairing"]["display"] == "-2 - 4*sqrt(2)"
    assert all(row["max_root_pairing"]["signs"][0] < 0 for row in audit["packets"])


def test_chronological_cross_scan_has_only_tangencies_before_cutter() -> None:
    audit = payload()["chronological_acceptance"]
    assert audit["nonnegative_row_count"] == 12
    assert audit["pre_H164_tangent_row_count"] == 9
    assert audit["H164_q2_tangent_row_count"] == 2
    assert audit["H164_positive_row_count"] == 1
    rows = audit["rows"]
    assert all(
        extension["pairing"]["signs"][0] == 0
        for row in rows[:-1]
        for extension in row["extensions"]
    )
    final = rows[-1]
    assert final["gram"] == [300, 180, 232]
    assert final["extensions"] == [
        {
            "pairing": {
                "exact": Silver(4).to_json(),
                "display": "4",
                "signs": [1, 1],
            },
            "d_dot_u": 34,
            "d_dot_v": 26,
            "e_dot_u": -28,
            "e_dot_v": -24,
            "Schur_A": 0,
            "Schur_B": 0,
            "Schur_N": 0,
        }
    ]


def test_exact_pre_h164_tangent_table() -> None:
    rows = payload()["chronological_acceptance"]["rows"][:9]
    assert [(row["height"]["display"], row["norm"]["display"], row["gram"]) for row in rows] == [
        ("21 + 15*sqrt(2)", "6", [60, 36, 45]),
        ("41 + 29*sqrt(2)", "2", [40, 22, 29]),
        ("41 + 29*sqrt(2)", "4", [88, 40, 58]),
        ("56 + 40*sqrt(2)", "4", [124, 52, 80]),
        ("82 + 58*sqrt(2)", "2", [86, 40, 58]),
        ("82 + 58*sqrt(2)", "2", [78, 44, 58]),
        ("103 + 73*sqrt(2)", "2", [104, 52, 73]),
        ("103 + 73*sqrt(2)", "4", [208, 104, 146]),
        ("137 + 97*sqrt(2)", "2", [140, 68, 97]),
    ]


def test_replacement_ray_is_hyperideal_and_strict_through_h164() -> None:
    ray = payload()["replacement_hyperideal_ray"]
    assert ray["scale"] == "13/40"
    assert ray["norm"]["display"] == "569/400 - 1*sqrt(2)"
    assert ray["norm"]["signs"][0] > 0
    assert ray["control_count"] == 91
    assert ray["retained_necessary_type_count"] == 112
    assert [row["gram"] for row in ray["relaxed_nonnegative_rows"]] == [
        [28, 16, 20],
        [60, 36, 45],
        [228, 144, 180],
        [100, 64, 80],
    ]
    assert ray["explicit_cutter_pairing"]["display"] == "-19/20 + 9/20*sqrt(2)"
    assert ray["explicit_cutter_pairing"]["signs"][0] < 0
    assert ray["alternate_H164_relaxed_margin"]["signs"][0] < 0


def test_scope_is_bounded_and_does_not_claim_discovery() -> None:
    data = payload()
    assert "bounded prefix" in data["scope"]
    assert "not Vinberg termination" in data["scope"]
    assert data["discovery"] is None

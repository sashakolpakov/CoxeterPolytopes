import json
from hashlib import sha256
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_lambda22 import (
    LAMBDA22_COEFFICIENT,
    lambda22_gram,
    lambda22_lorentz_pair_silver,
    lambda22_pair_silver,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT / "data/enumerations/silver_lambda22_h164_independent_audit.json"
)


def payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def silver(row: dict[str, object]) -> Silver:
    return Silver.from_json(row["exact"])


def test_explicit_candidate_is_primitive_reflective_norm_four() -> None:
    data = payload()
    candidate = data["candidate"]
    d = tuple(candidate["d"])
    e = tuple(candidate["e"])
    u = tuple(candidate["u"])
    v = tuple(candidate["v"])
    x = silver(candidate["control"])
    assert u == tuple(8 * d[index] - e[index] for index in range(22))
    assert v == tuple(6 * d[index] - e[index] for index in range(22))
    assert (
        lambda22_pair_silver(u, u),
        lambda22_pair_silver(v, v),
        lambda22_pair_silver(u, v),
    ) == (300, 180, 232)
    root = (x, *(Silver(u[index], v[index]) for index in range(22)))
    assert lambda22_lorentz_pair_silver(root, root) == 4
    assert LAMBDA22_COEFFICIENT * x * x / 4 == Silver(164, 116)
    assert root[12] == Silver(-1, -1)
    assert abs(root[12].norm()) == 1
    assert candidate["reflection_integral"] is True


def test_reflection_integrality_is_checked_on_the_full_catalogue_basis() -> None:
    candidate = payload()["candidate"]
    u = tuple(candidate["u"])
    v = tuple(candidate["v"])
    gram = lambda22_gram()
    assert silver(candidate["time_reflection_coefficient"]) == Silver(6, 4)
    expected = []
    for row in range(22):
        gu = sum(gram[row][column] * u[column] for column in range(22))
        gv = sum(gram[row][column] * v[column] for column in range(22))
        assert gu % 2 == 0 and gv % 2 == 0
        expected.append(Silver(gu, gv) / 2)
    actual = [silver(row) for row in candidate["spatial_reflection_coefficients"]]
    assert actual == expected
    assert all(value.is_integral for value in actual)


def test_candidate_cuts_the_old_ray_but_is_strictly_acute_to_H14() -> None:
    data = payload()
    old_pairing = silver(data["candidate"]["old_ray"]["pairing"])
    assert old_pairing == Silver(-2, 2) / 3
    assert old_pairing.sign(0) > 0
    maxima = data["H14_packet_cross_maxima"]
    assert [silver(row["spatial_maximum"]) for row in maxima] == [
        Silver(124, 86),
        Silver(190, 132),
    ]
    assert [silver(row["lorentz_threshold"]) for row in maxima] == [
        Silver(136, 96),
        Silver(192, 136),
    ]
    assert [silver(row["pairing_maximum"]) for row in maxima] == [
        Silver(-12, -10),
        Silver(-2, -4),
    ]


def test_complete_necessary_chronology_has_no_earlier_positive_blocker() -> None:
    chronology = payload()["chronology"]
    assert chronology["control_count"] == 91
    assert chronology["post_H14_necessary_binary_type_count"] == 108
    assert chronology["strictly_positive_pre_H164_count"] == 0
    assert chronology["pre_H164_tangency_count"] == 9
    assert all(
        silver(row["pairing"]) == 0
        for row in chronology["pre_H164_tangencies"]
    )
    assert [row["height"]["display"] for row in chronology["pre_H164_tangencies"]] == [
        "21 + 15*sqrt(2)",
        "41 + 29*sqrt(2)",
        "41 + 29*sqrt(2)",
        "56 + 40*sqrt(2)",
        "82 + 58*sqrt(2)",
        "82 + 58*sqrt(2)",
        "103 + 73*sqrt(2)",
        "103 + 73*sqrt(2)",
        "137 + 97*sqrt(2)",
    ]


def test_unique_positive_same_height_extension_is_the_root_itself() -> None:
    chronology = payload()["chronology"]
    assert chronology["H164_tangent_q2_count"] == 2
    assert [row["gram"] for row in chronology["H164_tangent_q2_rows"]] == [
        [162, 84, 116],
        [158, 86, 116],
    ]
    assert chronology["H164_positive_count"] == 1
    row = chronology["H164_positive_rows"][0]
    assert row["norm"]["display"] == "4"
    assert row["gram"] == [300, 180, 232]
    assert row["cross"] == [34, 26, -28, -24]
    assert row["schur_scaled22"] == [0, 0, 0]
    assert silver(row["pairing"]) == 4
    # Solve against inverse([[4,-2],[-2,12]]) exactly.
    p, q, r, s = row["cross"]
    assert ((12 * p + 2 * r) // 44, (2 * p + 4 * r) // 44) == (8, -1)
    assert ((12 * q + 2 * s) // 44, (2 * q + 4 * s) // 44) == (6, -1)


def test_replacement_ray_is_hyperideal_and_strict_on_the_full_superset() -> None:
    replacement = payload()["replacement_hyperideal_ray"]
    norm = silver(replacement["norm"])
    assert norm == Silver(569) / 400 - Silver(0, 1)
    assert norm.sign(0) > 0
    assert silver(replacement["pairing_with_H164_candidate"]) == Silver(-19, 9) / 20
    assert replacement["complete_binary_gate_type_count"] == 62
    maximum = replacement["global_necessary_superset_maximum"]
    assert maximum["height"]["display"] == "164 + 116*sqrt(2)"
    assert maximum["norm"]["display"] == "2"
    assert maximum["gram"] == [162, 84, 116]
    assert maximum["d_cross"] == [25, 18]
    assert silver(maximum["pairing_times_40"]) == Silver(5, -6)
    assert silver(replacement["global_pairing_maximum"]) == Silver(5, -6) / 40
    assert silver(replacement["global_pairing_maximum"]).sign(0) < 0
    assert replacement["strict_against_every_necessary_type"] is True


def test_independent_base_and_scope_are_pinned() -> None:
    data = payload()
    base = data["base_H14_independent_audit"]
    path = ROOT / base["path"]
    assert sha256(path.read_bytes()).hexdigest() == base["sha256"]
    assert data["status"] == "EXACT_INDEPENDENT_PASS"
    assert "bounded Vinberg prefix" in data["scope"]
    assert "no global termination claim" in data["scope"]

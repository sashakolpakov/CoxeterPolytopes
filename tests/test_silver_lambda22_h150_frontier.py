import json
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_lambda22 import (
    lambda22_gram,
    lambda22_pairing_radical_basis,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/silver_lambda22_h150_frontier.json"


def payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def test_complete_relaxed_frontier_census_and_seven_resolutions() -> None:
    data = payload()
    scan = data["necessary_superset_scan"]
    assert scan["control_count"] == 58
    assert scan["retained_level_count"] == 18
    assert scan["retained_necessary_type_count"] == 72
    assert scan["post_H14_nonnegative_type_count"] == 7
    assert len(data["exact_resolutions"]) == 7
    assert [row["gram"] for row in scan["post_H14_nonnegative_types"]] == [
        [60, 36, 45],
        [354, 198, 264],
        [88, 40, 58],
        [228, 144, 180],
        [124, 52, 80],
        [100, 64, 80],
        [230, 126, 170],
    ]


def test_h36_difference_dichotomy_is_exact() -> None:
    triples = [(406 - 4 * index, 172 + 2 * index, 264) for index in range(17)]
    for u_square, v_square, u_dot_v in triples:
        w_square = u_square + v_square - 2 * u_dot_v
        s_square = u_square + 4 * v_square - 4 * u_dot_v
        assert w_square == 222 - v_square
        assert s_square == 2 * v_square - 306
        if v_square <= 198:
            assert 0 < s_square < 96
        else:
            assert 0 < w_square < 24


def test_s2_radical_classes_have_only_zero_or_four_norm_mod8() -> None:
    data = payload()["S2_radical_norm_residues_mod8"]
    assert data["class_count"] == 4
    assert data["residue_multiset"] == [0, 4, 4, 4]
    gram = lambda22_gram()
    basis = lambda22_pairing_radical_basis(2)
    for row in data["classes"]:
        vector = row["representative"]
        assert all(
            sum(gram[i][j] * vector[j] for j in range(22)) % 2 == 0
            for i in range(22)
        )
        norm = sum(
            vector[i] * gram[i][j] * vector[j]
            for i in range(22)
            for j in range(22)
        )
        assert norm % 8 == row["norm_mod_8"]


def test_h70_all_possible_u_norms_violate_the_radical_residue_rule() -> None:
    possible = [250, 246, 242, 238, 234, 230]
    assert {value % 8 for value in possible} == {2, 6}
    assert {value % 8 for value in possible}.isdisjoint({0, 4})


def test_cross_pairing_uses_lcm_and_kills_the_whole_h72_level() -> None:
    rejected = payload()["necessary_superset_scan"]["H72_cross_gate_rejections"]
    assert rejected == [
        [780, 354, 525],
        [768, 360, 525],
        [756, 366, 525],
        [744, 372, 525],
        [732, 378, 525],
        [720, 384, 525],
        [708, 390, 525],
    ]
    assert all(u % 6 == 0 for u, _, _ in rejected)
    assert all(v % 3 == 0 for _, v, _ in rejected)
    assert all(dot % 6 == 3 for _, _, dot in rejected)


def test_complete_pari_short_shell_pairing_spectra() -> None:
    shells = payload()["pari_direction_shell_audit"]["pairing_radical_shells"]
    assert shells["2"] == {
        "signed_count": 58_344,
        "maximum": 16,
        "representatives": 29_172,
        "norm_12_representatives": 4_224,
        "norm_12_abs_pairings": [0, 2, 4],
        "norm_16_representatives": 24_948,
        "norm_16_abs_pairings": [0, 2, 4, 8],
    }
    assert shells["6"] == {
        "signed_count": 1_782,
        "maximum": 96,
        "representatives": 891,
        "norm_96_representatives": 891,
        "norm_96_abs_pairings": [0, 6, 12],
    }


def test_scope_remains_bounded_and_hyperideal() -> None:
    data = payload()
    assert data["height_bound"]["display"] == "150"
    assert data["hyperideal_witness"]["norm"]["display"] == "13/9 - 1*sqrt(2)"
    assert data["hyperideal_witness"]["norm"]["signs"][0] > 0
    assert data["discovery"] is None
    assert "not global nonreflectivity" in data["scope"]

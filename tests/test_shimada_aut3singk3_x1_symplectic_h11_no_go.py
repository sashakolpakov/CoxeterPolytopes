from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
STOCK = ROOT / "data/enumerations/shimada_aut3singk3_x1_1094_wall_stock.json"
CERTIFICATE = (
    ROOT / "data/enumerations/shimada_aut3singk3_x1_symplectic_h11_no_go.json"
)
STOCK_SHA256 = "8fab28bffb6c617b09d81fd9735581abd39a28f89cd6d760bedaefbbdf8b330e"


def _dot(left: list[Fraction], gram: list[list[Fraction]], right: list[Fraction]) -> Fraction:
    return sum(
        left[i] * gram[i][j] * right[j]
        for i in range(len(left))
        for j in range(len(right))
    )


def _fractions(values: list[int | str]) -> list[Fraction]:
    return [Fraction(value) for value in values]


def test_x1_complete_source_and_denominator_are_frozen() -> None:
    assert hashlib.sha256(STOCK.read_bytes()).hexdigest() == STOCK_SHA256
    stock = json.loads(STOCK.read_text(encoding="utf-8"))
    assert stock["all_checks_pass"] is True
    assert stock["complete_induced_chamber"]["wall_count"] == 1094
    assert stock["complete_induced_chamber"]["wall_orbit_sizes"] == [
        45, 45, 10, 30, 72, 60, 12, 40, 120, 120, 120, 180, 120, 120
    ]
    assert stock["ample_stabilizer"]["order"] == 720


def test_x1_high_source_group_queue_is_closed() -> None:
    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    census = certificate["high_subgroup_census"]
    assert census["fixed_rank_histogram"] == {
        "2": 144,
        "4": 324,
        "6": 90,
        "8": 80,
        "10": 36,
        "12": 45,
        "20": 1,
    }
    assert census["unique_high_element_class"]["centralizer_order"] == 16
    assert census["unique_high_element_class"]["conjugacy_class_size"] == 45
    assert census["high_closed_V4_count"] == 30
    assert census["high_closed_V4_common_fixed_rank_histogram"] == {"8": 30}
    assert census["only_C2_can_be_a_proper_fixed_rank_at_least_11_subgroup"]


def test_independent_strict_forbidden_pair_replay() -> None:
    """Replay one decisive pair without importing the production certifier."""

    stock = json.loads(STOCK.read_text(encoding="utf-8"))
    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    gram = [[Fraction(value) for value in row] for row in stock["surface"]["published_neron_severi_gram"]]
    fold = certificate["unique_high_fixed_section"]
    restrictions = [
        _fractions(row["primitive_fixed_vector_in_published_NS_basis"])
        for row in fold["unique_restrictions"]
    ]
    facets = [
        vector
        for vector, record in zip(restrictions, fold["unique_restrictions"], strict=True)
        if record["lorentz_role"] == "strict_hyperbolic_facet"
    ]
    assert len(restrictions) == 580
    assert len(facets) == 162

    obstruction = fold["pair_gate"]["first_forbidden_pairs"][0]
    assert obstruction["facet_pair_one_based"] == [7, 34]
    left, right = (facets[index - 1] for index in obstruction["facet_pair_one_based"])
    left_norm = _dot(left, gram, left)
    right_norm = _dot(right, gram, right)
    pairing = _dot(left, gram, right)
    assert (left_norm, right_norm, pairing) == (Fraction(-4), Fraction(-6), Fraction(4))
    assert pairing * pairing / (left_norm * right_norm) == Fraction(2, 3)

    witness_by_restriction = {
        row["unique_restriction_index_one_based"]: _fractions(
            row["witness_in_published_NS_basis"]
        )
        for row in fold["strict_facet_witnesses"]
    }
    strict_positions = [
        index + 1
        for index, record in enumerate(fold["unique_restrictions"])
        if record["lorentz_role"] == "strict_hyperbolic_facet"
    ]
    for facet_number in obstruction["facet_pair_one_based"]:
        restriction_position = strict_positions[facet_number - 1]
        witness = witness_by_restriction[restriction_position]
        pairings = [_dot(witness, gram, vector) for vector in restrictions]
        assert _dot(witness, gram, witness) > 0
        assert pairings[restriction_position - 1] == 0
        assert all(
            value > 0
            for index, value in enumerate(pairings)
            if index != restriction_position - 1
        )

    assert fold["pair_gate"]["pair_count"] == 13_041
    assert fold["pair_gate"]["coxeter_edge_counts"]["forbidden"] == 915
    assert certificate["conclusion"][
        "Coxeter_polytope_above_H9_in_source_group_denominator"
    ] is False


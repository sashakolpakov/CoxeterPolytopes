from __future__ import annotations

import json
from pathlib import Path

from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_a2_ternary_golay import (
    SQRT2_A2_GOLAY_SEED,
    TWO_SQRT2_A2_GOLAY_SEED,
    all_seed_h17_searches,
    all_seed_searches,
    control_decisions,
    first_h24_cutter,
    first_packet_geometry,
    first_packet_is_exact,
    h24_candidate_census,
    h24_uniform_witness,
    h100_uniform_witness,
    niemeier_construction_is_even_unimodular,
    ternary_golay_codewords,
    ternary_golay_is_self_dual,
    ternary_golay_pair_contingencies,
    ternary_golay_weight_enumerator,
)


ROOT = Path(__file__).resolve().parents[1]


def test_extended_ternary_golay_code_is_exact() -> None:
    assert len(ternary_golay_codewords()) == 729
    assert ternary_golay_is_self_dual()
    assert ternary_golay_weight_enumerator() == {0: 1, 6: 264, 9: 440, 12: 24}
    assert niemeier_construction_is_even_unimodular()


def test_first_packets_are_exact_affine_a2_caps() -> None:
    for seed in (SQRT2_A2_GOLAY_SEED, TWO_SQRT2_A2_GOLAY_SEED):
        assert first_packet_is_exact(seed)
        decisions = control_decisions(seed)
        assert [decision.status for decision in decisions].count("accepted") == 1
        assert [decision.status for decision in decisions][-1] == "rejected"


def test_product_triangle_geometry_is_complete_and_nonfinite_volume() -> None:
    for seed in (SQRT2_A2_GOLAY_SEED, TWO_SQRT2_A2_GOLAY_SEED):
        geometry = first_packet_geometry(seed)
        assert geometry.extreme_ray_count == 3**12 == 531_441
        assert geometry.finite_vertex_count == 25
        assert geometry.cusp_count == 0
        assert geometry.hyperideal_ray_count == 531_416
        assert not geometry.finite_volume
        assert not geometry.compact
        assert geometry.exact_norms_by_nonzero_blocks[1].sign(0) < 0
        assert geometry.exact_norms_by_nonzero_blocks[2].sign(0) > 0


def test_control_census_through_h10() -> None:
    searches = all_seed_searches()
    assert [len(search.controls) for search in searches] == [5, 3]
    assert all(search.height_zero_wall_count == 24 for search in searches)
    assert all(search.accepted_positive_wall_count == 12 for search in searches)
    assert all(search.total_wall_count == 36 for search in searches)
    assert searches[0].seed.timelike_coefficient == SQRT2
    assert searches[1].seed.timelike_coefficient == 2 * SQRT2
    assert searches[0].controls[-1].target == Silver(10, 6)
    assert searches[1].controls[-1].target == Silver(10, 6)


def test_second_packet_contingency_census_and_witness() -> None:
    contingencies = ternary_golay_pair_contingencies()
    assert len(contingencies) == 288
    assert sum(multiplicity for _, multiplicity in contingencies) == 3**12
    searches = all_seed_h17_searches()
    assert [len(search.controls) for search in searches] == [7, 5]
    for search in searches:
        assert search.accepted_positive_wall_count == 3_972
        assert search.total_wall_count == 3_996
        assert search.second_packet is not None
        assert search.second_packet.height == Silver(8, 6)
        assert search.second_packet.norm == Silver(2)
        assert search.second_packet.target == Silver(18, 12)
        assert search.second_packet.candidate_count == 3_960
        assert search.second_packet.surviving_contingency_count == 6
        assert search.second_packet.maximum_uniform_weight_sum == Silver(9, 9)
        assert search.strict_hyperideal_witness is not None
        assert search.strict_hyperideal_witness.strict
        assert search.strict_hyperideal_witness.lorentz_norm.sign(0) > 0
        assert search.strict_hyperideal_witness.first_packet_maximum_pairing.sign(0) < 0
        assert search.strict_hyperideal_witness.second_packet_maximum_pairing.sign(0) < 0


def test_second_packet_witnesses_have_the_expected_exact_norms() -> None:
    silver, twice_silver = all_seed_h17_searches()
    assert silver.strict_hyperideal_witness is not None
    assert twice_silver.strict_hyperideal_witness is not None
    assert silver.strict_hyperideal_witness.block_dynkin_weight == Silver(0, "1/5")
    assert silver.strict_hyperideal_witness.lorentz_norm == Silver("48/25", -1)
    assert twice_silver.strict_hyperideal_witness.block_dynkin_weight == Silver("2/5")
    assert twice_silver.strict_hyperideal_witness.lorentz_norm == Silver("96/25", -2)


def test_explicit_h24_cutter_is_reflective_and_acute() -> None:
    cutters = [
        first_h24_cutter(SQRT2_A2_GOLAY_SEED),
        first_h24_cutter(TWO_SQRT2_A2_GOLAY_SEED),
    ]
    assert [cutter.control for cutter in cutters] == [Silver(4, 3), Silver(3, 2)]
    for cutter in cutters:
        assert cutter.height == Silver(24, 17)
        assert cutter.norm == Silver(2)
        assert cutter.target == Silver(50, 34)
        assert cutter.maximum_second_packet_pairing == -1
        assert cutter.primitive
        assert cutter.reflection_integral
        assert cutter.acute_to_previous_packets


def test_h24_candidate_superset_and_uniform_witnesses() -> None:
    searches = (SQRT2_A2_GOLAY_SEED, TWO_SQRT2_A2_GOLAY_SEED)
    censuses = [h24_candidate_census(seed) for seed in searches]
    for census in censuses:
        assert census.height == Silver(24, 17)
        assert census.norm == Silver(2)
        assert census.target == Silver(50, 34)
        assert census.candidate_count == 4_733_256
        assert census.surviving_contingency_count == 106
        assert census.maximum_uniform_weight_sum == Silver(21, 15)
    witnesses = [h24_uniform_witness(seed) for seed in searches]
    assert witnesses[0].block_dynkin_weight == Silver("1/4")
    assert witnesses[0].lorentz_norm == Silver("3/2", -1)
    assert witnesses[0].h24_candidate_maximum_pairing == Silver("-3/4", "-1/4")
    assert witnesses[1].block_dynkin_weight == Silver("3/8")
    assert witnesses[1].lorentz_norm == Silver("27/8", -2)
    assert witnesses[1].h24_candidate_maximum_pairing == Silver("-1/8", "-3/8")
    assert all(witness.strict for witness in witnesses)


def test_frozen_json_matches_search() -> None:
    path = ROOT / "data/enumerations/silver_a2_ternary_golay_h24_h10.json"
    payload = json.loads(path.read_text())
    assert payload["schema"] == "silver-a2-ternary-golay-h24-h10-v1"
    assert payload["construction"]["weight_enumerator"] == {
        "0": 1,
        "6": 264,
        "9": 440,
        "12": 24,
    }
    assert payload["construction"]["even_unimodular_lattice"] is True
    assert len(payload["records"]) == 2
    assert all(
        record["geometry"]["hyperideal_ray_count"] == 531_416
        for record in payload["records"]
    )
    assert payload["outcome"]["compact_survivor_count"] == 0


def test_h24_cutter_json_matches_frontier() -> None:
    path = ROOT / "data/enumerations/silver_a2_ternary_golay_h24_cutter.json"
    payload = json.loads(path.read_text())
    assert payload["schema"] == "silver-a2-ternary-golay-h24-cutter-v1"
    assert payload["glue_dynamic_program"]["ordinary_contingency_type_count"] == 288
    assert payload["glue_dynamic_program"]["two_pointed_contingency_type_count"] == 7758
    assert payload["glue_dynamic_program"]["h8_same_height_pairing_census"] == {
        "maximum": {"a": "0", "b": "0"},
        "method": "independent exact 3960-by-3960 integral coefficient replay",
        "negative": 7_839_114,
        "positive": 0,
        "unordered_pair_count": 7_839_180,
        "zero": 66,
    }
    assert len(payload["records"]) == 2
    for record in payload["records"]:
        assert record["h24_candidate_superset"]["candidate_count"] == 4_733_256
        assert record["explicit_accepted_cutter"]["maximum_second_packet_pairing"] == {
            "a": "-1",
            "b": "0",
        }
        assert record["strict_uniform_witness"]["strict"] is True
    assert payload["outcome"]["compact_survivor_count"] == 0


def test_h100_uniform_witnesses_survive_all_candidate_supersets() -> None:
    silver = h100_uniform_witness(SQRT2_A2_GOLAY_SEED)
    twice_silver = h100_uniform_witness(TWO_SQRT2_A2_GOLAY_SEED)
    assert silver.block_dynkin_weight == Silver("1/4")
    assert silver.lorentz_norm == Silver("3/2", -1)
    assert twice_silver.block_dynkin_weight == Silver("7/20")
    assert twice_silver.lorentz_norm == Silver("147/50", -2)
    for witness in (silver, twice_silver):
        assert witness.strict
        assert [height for height, _ in witness.nonempty_level_maximum_pairings] == [
            Silver(8, 6),
            Silver(24, 17),
            Silver(41, 29),
            Silver(48, 34),
        ]
        assert all(
            pairing.sign(0) < 0
            for _, pairing in witness.nonempty_level_maximum_pairings
        )


def test_h200_frozen_frontier_has_eight_nonempty_levels_and_strict_witnesses() -> None:
    path = ROOT / "data/enumerations/silver_a2_ternary_golay_h200.json"
    payload = json.loads(path.read_text())
    assert payload["schema"] == "silver-a2-ternary-golay-h24-h200-v1"
    assert payload["real_height_bound"] == 200
    assert len(payload["records"]) == 2
    expected_counts = [
        3_960,
        4_733_256,
        1_724_736,
        720_156_144,
        27_720,
        27_720,
        45_875_104,
        103_188_357_648,
    ]
    for record in payload["records"]:
        assert record["nonempty_packet_count"] == 8
        assert [
            row["candidate_count"]
            for row in record["post_affine_censuses"]
            if row["candidate_count"]
        ] == expected_counts
        witness = record["strict_uniform_witness"]
        assert witness["strict"] is True
        assert len(witness["nonempty_level_maximum_pairings"]) == 8
    assert payload["outcome"]["compact_survivor_count"] == 0

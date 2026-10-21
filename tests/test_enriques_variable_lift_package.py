from __future__ import annotations

import hashlib
import json

from scripts.certify_enriques_variable_lift_package import (
    OUTPUT,
    RAW_DATA,
    RAW_DATA_SHA256,
    TYPE_TO_RECORD,
    analyze,
)


ARCHIVE_SHA256 = "27115e714ef4acd5a8bb7801a12795668f9c9b1f0309816f336454c591583a6a"


def _document() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_archive_replays_from_pinned_shimada_input() -> None:
    assert hashlib.sha256(RAW_DATA.read_bytes()).hexdigest() == RAW_DATA_SHA256
    assert analyze() == _document()
    assert hashlib.sha256(OUTPUT.read_bytes()).hexdigest() == ARCHIVE_SHA256


def test_fixed_enriques_cover_lattice_and_glue() -> None:
    document = _document()
    assert document["downstairs_lattice_M"] == {
        "description": "L10 = U + E8(-1)",
        "determinant": -1,
        "inertia": {"negative": 9, "positive": 1, "zero": 0},
        "rank": 10,
    }
    model = document["fixed_K3_cover_model"]
    assert model["ambient"]["inertia"] == {
        "negative": 19,
        "positive": 3,
        "zero": 0,
    }
    assert model["invariant_lattice_S"]["inertia"] == {
        "negative": 9,
        "positive": 1,
        "zero": 0,
    }
    assert model["anti_invariant_lattice_T"]["inertia"] == {
        "negative": 10,
        "positive": 2,
        "zero": 0,
    }
    assert model["invariant_lattice_S"]["absolute_discriminant"] == 2**10
    assert model["anti_invariant_lattice_T"]["absolute_discriminant"] == 2**10
    assert model["invariant_lattice_S"]["primitive_embedding_smith_diagonal"] == [1] * 10
    assert model["anti_invariant_lattice_T"]["primitive_embedding_smith_diagonal"] == [1] * 12
    assert model["glue"] == {
        "S_perp_T": True,
        "half_sum_formula": "(s_i+t_i)/2=(0,e_i,0)",
        "index_of_S_plus_T_in_Lambda": 2**10,
        "quotient": "(Z/2Z)^10",
        "ten_half_sum_generators_integral": True,
    }


def test_all_seven_canonical_simultaneous_lifts_are_period_obstructed() -> None:
    document = _document()
    obstruction = document["canonical_simultaneous_lift_obstruction"]
    assert obstruction["all_seven_wall_systems_span_L10"] is True
    assert obstruction["all_seven_canonical_simultaneous_lifts_obstructed"] is True
    assert len(document["types"]) == 7
    for record in document["types"]:
        assert record["raw_shimada_record"] == TYPE_TO_RECORD[
            record["nikulin_kondo_type"]
        ]
        assert record["wall_rank"] == 10
        assert record["all_downstairs_walls_primitive_minus_two_roots"] is True
        lift = record["canonical_lift"]
        assert lift["all_individual_pairs_are_disjoint_exchanged_minus_two_classes"] is True
        assert lift["b_i_dot_b_j_equals_downstairs_configuration"] is True
        assert lift["b_i_dot_tau_b_j_is_zero"] is True
        assert lift["delta_gram_equals_twice_downstairs_configuration"] is True
        assert lift["delta_span_rank"] == 10
        assert lift["delta_span_signature"] == {"positive": 1, "negative": 9}
        assert lift["delta_gram_inertia"] == {
            "positive": 1,
            "negative": 9,
            "zero": record["facet_count"] - 10,
        }
        assert lift["orthogonal_complement_in_T_inertia"] == {
            "positive": 1,
            "negative": 1,
            "zero": 0,
        }
        assert lift["positive_directions_left_for_period_plane"] == 1
        assert lift["required_positive_directions_for_K3_period_plane"] == 2
        assert lift["simultaneous_period_exists"] is False


def test_two_sheet_intersection_distribution_is_exact() -> None:
    identities = _document()["general_lift_identities"]
    assert identities["nonnegative_curve_distribution"] == {
        "m=0": [0],
        "m=1": [-2, 2],
        "m=2": [-4, 0, 4],
    }
    assert "a+c=m" in identities["two_wall_formula"]
    assert "2(2a-m)" in identities["two_wall_formula"]

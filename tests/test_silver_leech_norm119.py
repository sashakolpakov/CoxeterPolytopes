from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_leech_manifest import (
    _control_levels,
    fixed_direction_pairing_bound,
    has_short_binary_combination,
    mirror_norm_classes,
)
from coxeter_polytopes.search.silver_leech_norm119 import (
    INVERSE_UNIT_SQUARE,
    NORM119_COEFFICIENT,
    P7_DIVISOR,
    P17_DIVISOR,
    P119_DIVISOR,
    SIGN_CHANGING_UNIT,
    _fast_fixed_direction_pairing_bound,
    norm119_leech_prefix_audit,
    norm119_local_classes,
    strongest_norm119_prefix_witness,
)


def test_norm119_local_factorization_and_normalized_mirror_classes() -> None:
    assert NORM119_COEFFICIENT == P119_DIVISOR * SIGN_CHANGING_UNIT
    assert P119_DIVISOR == P7_DIVISOR * P17_DIVISOR == Silver(19, 11)
    assert P7_DIVISOR.norm() == 7
    assert P17_DIVISOR.norm() == 17
    assert NORM119_COEFFICIENT.norm() == -119
    assert INVERSE_UNIT_SQUARE == Silver(1, 1) ** -2 == Silver(3, -2)

    rows = norm119_local_classes()
    assert tuple(row.canonical_norm for row in rows) == (
        Silver(1),
        Silver(2),
        Silver(3, 1),
        Silver(13, -5),
        Silver(5, 2),
        Silver(6, 2),
        Silver(26, -10),
        Silver(10, 4),
    )
    assert tuple(row.odd_spatial_divisor for row in rows) == (
        Silver(1),
        Silver(1),
        P7_DIVISOR,
        P119_DIVISOR,
        P17_DIVISOR,
        P7_DIVISOR,
        P119_DIVISOR,
        P17_DIVISOR,
    )
    assert tuple(row.dyadic_factor for row in rows) == (1, 2, 1, 1, 1, 2, 2, 2)
    assert all(
        row.canonical_norm
        == row.dyadic_factor
        * row.odd_spatial_divisor
        * row.unit_square_multiplier
        for row in rows
    )


def test_fast_pairing_oracle_matches_the_rectangular_reference() -> None:
    audit = norm119_leech_prefix_audit()
    representatives = tuple(
        row.surviving_triples[0]
        for row in audit.controls
        if row.surviving_triples
    )
    # One triple from every nonempty component exercises all seven packets.
    for triple in representatives:
        for direction_square in (4, 6, 8, 10):
            assert _fast_fixed_direction_pairing_bound(
                triple, direction_square
            ) == fixed_direction_pairing_bound(triple, direction_square)


def test_h10000_control_enumeration_is_complete_and_normalized() -> None:
    audit = norm119_leech_prefix_audit()
    reference = _control_levels(
        NORM119_COEFFICIENT,
        mirror_norm_classes(NORM119_COEFFICIENT),
        10_000,
    )
    assert tuple(
        (
            row.height,
            row.norm,
            row.control,
            row.target_before_gate,
            row.odd_spatial_divisor,
        )
        for row in audit.controls
    ) == reference
    assert audit.control_count == 358
    assert audit.empty_control_count == 347
    assert audit.nonempty_component_count == 11
    assert audit.nonempty_height_count == 7
    assert audit.raw_candidate_triple_count == 176
    assert audit.short_combination_exclusion_count == 12
    assert audit.surviving_triple_count == 164
    assert all(
        row.target_after_gate * row.odd_spatial_divisor**2
        == row.target_before_gate
        for row in audit.controls
    )
    assert all(
        row.normalized_wall_bound * row.odd_spatial_divisor
        == NORM119_COEFFICIENT * row.control
        for row in audit.controls
    )
    assert all(
        not row.surviving_triples
        for row in audit.controls
        if row.odd_spatial_divisor != 1
    )


def test_h10000_nonempty_packet_census_and_gram_filter() -> None:
    rows = tuple(
        row for row in norm119_leech_prefix_audit().controls if row.surviving_triples
    )
    assert tuple(
        (
            row.ordinal,
            row.nonempty_packet_ordinal,
            row.height,
            row.norm,
            row.control,
            len(row.raw_candidate_triples),
            len(row.surviving_triples),
        )
        for row in rows
    ) == (
        (75, 0, Silver(243, 172), Silver(1), Silver(3, 2), 4, 3),
        (76, 0, Silver(243, 172), Silver(2), Silver(4, 3), 9, 8),
        (110, 1, Silver(486, 344), Silver(2), Silver(6, 4), 11, 9),
        (157, 2, Silver(972, 688), Silver(2), Silver(8, 6), 3, 1),
        (187, 3, Silver(1417, 1002), Silver(1), Silver(7, 5), 13, 12),
        (188, 3, Silver(1417, 1002), Silver(2), Silver(10, 7), 26, 26),
        (267, 4, Silver(2834, 2004), Silver(2), Silver(14, 10), 36, 36),
        (302, 5, Silver(3563, 2520), Silver(1), Silver(11, 8), 10, 9),
        (303, 5, Silver(3563, 2520), Silver(2), Silver(16, 11), 19, 18),
        (347, 6, Silver(4737, 3350), Silver(1), Silver(13, 9), 15, 13),
        (348, 6, Silver(4737, 3350), Silver(2), Silver(18, 13), 30, 29),
    )
    assert all(
        not has_short_binary_combination(triple)
        for row in rows
        for triple in row.surviving_triples
    )


def test_only_norm_six_fixed_direction_survives_h10000_strictly() -> None:
    audits = {
        row.direction_square: row
        for row in norm119_leech_prefix_audit().direction_audits
    }
    assert tuple(
        audits[direction].strict_hyperideal_witness_exists
        for direction in (4, 6, 8, 10)
    ) == (False, True, False, False)
    assert audits[4].limiting_hyperideal_excess == Silver(
        Fraction(823, 36), Fraction(-97, 6)
    )
    assert audits[8].limiting_hyperideal_excess == audits[4].limiting_hyperideal_excess
    assert audits[10].limiting_hyperideal_excess == Silver(
        Fraction(4_240_917, 829_921), Fraction(-2_999_448, 829_921)
    )
    assert all(
        audits[direction].limiting_hyperideal_excess.sign(0) < 0
        for direction in (4, 8, 10)
    )

    witness = strongest_norm119_prefix_witness()
    assert witness.direction_square == 6
    assert witness.limiting_control_ordinal == 267
    assert witness.limiting_height == Silver(2834, 2004)
    assert witness.limiting_norm == 2
    assert witness.limiting_control == Silver(14, 10)
    assert witness.limiting_coordinate == Silver(
        Fraction(118, 219), Fraction(52, 73)
    )
    assert witness.limiting_hyperideal_excess == Silver(
        Fraction(77_231, 15_987), Fraction(-18_088, 5_329)
    )
    assert witness.shrink_denominator == 936
    assert witness.witness_coordinate == Silver(
        Fraction(55_165, 102_492), Fraction(935, 1_314)
    )
    assert witness.minimum_wall_margin == Silver(
        Fraction(57_515, 51_246), Fraction(-19_541, 34_164)
    )
    assert witness.witness_hyperideal_excess == Silver(
        Fraction(8_428_441_993, 1_750_768_344),
        Fraction(-38_203_717, 11_222_874),
    )
    assert all(margin.sign(0) > 0 for margin in witness.component_wall_margins)
    assert witness.witness_hyperideal_excess.sign(0) > 0


def test_frozen_norm119_artifact_is_pinned_and_scoped() -> None:
    root = Path(__file__).resolve().parents[1]
    artifact = root / "data/enumerations/silver_leech_norm119_h10000.json"
    assert sha256(artifact.read_bytes()).hexdigest() == (
        "3e845cada57b93a1aeecdaaf15310ca5fdaf1a97a05f4e9948e96bd9d0c52b13"
    )
    payload = json.loads(artifact.read_text())
    assert payload["schema"] == "coxeter.silver_leech_norm119_necessary_prefix.v1"
    assert payload["audit"]["complete_height_bound"] == 10_000
    assert payload["audit"]["control_count"] == 358
    assert payload["audit"]["surviving_triple_count"] == 164
    assert payload["strongest_fixed_direction_witness"]["direction_square"] == 6
    assert "not a global nonreflectivity theorem" in payload["scope"]

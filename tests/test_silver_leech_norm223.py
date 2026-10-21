from fractions import Fraction

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_leech import candidate_leech_gram_triples
from coxeter_polytopes.search.silver_leech_manifest import (
    fixed_direction_pairing_bound,
    has_short_binary_combination,
)
from coxeter_polytopes.search.silver_leech_norm223 import (
    NORM223_COEFFICIENT,
    P223_DIVISOR,
    SIGN_CHANGING_UNIT,
    _candidate_triples,
    d4_pairing_bound,
    norm223_control_skeleton,
    norm223_leech_prefix_audit,
    norm223_local_classes,
)


def test_norm223_factorization_and_local_classes() -> None:
    assert NORM223_COEFFICIENT == P223_DIVISOR * SIGN_CHANGING_UNIT
    assert P223_DIVISOR.norm() == 223
    assert NORM223_COEFFICIENT.norm() == -223
    rows = norm223_local_classes()
    assert tuple(row.canonical_norm for row in rows) == (
        Silver(1),
        Silver(2),
        Silver(15, 1),
        Silver(30, 2),
    )
    assert tuple(row.odd_spatial_divisor for row in rows) == (
        Silver(1),
        Silver(1),
        P223_DIVISOR,
        P223_DIVISOR,
    )


def test_fast_binary_and_pairing_oracles_match_references() -> None:
    audit = norm223_leech_prefix_audit()
    representatives = audit.components[::5]
    for row in representatives:
        assert _candidate_triples(row.target_after_gate) == candidate_leech_gram_triples(
            row.target_after_gate
        )
        survivors = tuple(
            triple
            for triple in _candidate_triples(row.target_after_gate)
            if not has_short_binary_combination(triple)
        )
        assert d4_pairing_bound(survivors[0]) == fixed_direction_pairing_bound(
            survivors[0], 4
        )


def test_h250000_complete_census_and_strict_witness() -> None:
    audit = norm223_leech_prefix_audit()
    assert len(norm223_control_skeleton()) == audit.control_count == 1129
    assert audit.control_skeleton_sha256 == (
        "07d22789e741b59f956e22829fea5fdd51f743ee1e2c5245c000b94a8bd36a0d"
    )
    assert audit.empty_control_count == 1094
    assert audit.nonempty_component_count == 35
    assert audit.nonempty_height_count == 24
    assert audit.raw_candidate_triple_count == 2875
    assert audit.short_combination_exclusion_count == 20
    assert audit.surviving_triple_count == 2855
    assert audit.limiting_control_ordinal == 1084
    assert audit.limiting_height == Silver(114998, 81316)
    assert audit.limiting_norm == 2
    assert audit.limiting_control == Silver(130, 92)
    assert audit.limiting_coordinate == Silver(
        Fraction(3802, 2159), Fraction(-698, 2159)
    )
    assert audit.limiting_hyperideal_excess == Silver(
        Fraction(122315101, 4661281), Fraction(-86488302, 4661281)
    )
    assert audit.shrink_denominator == 29193
    assert audit.witness_coordinate == Silver(
        Fraction(110987984, 63027687), Fraction(-20376016, 63027687)
    )
    assert audit.minimum_wall_margin == Silver(
        Fraction(886, 29193), Fraction(208, 9731)
    )
    assert audit.witness_hyperideal_excess == Silver(
        Fraction(34745715955068223, 1324163109523323),
        Fraction(-73706794102313518, 3972489328569969),
    )
    assert audit.witness_hyperideal_excess.sign(0) > 0

from collections import Counter
from fractions import Fraction

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_odd_leech import (
    EPSILON,
    FIRST_PACKET_CONTROL,
    FIRST_PACKET_HEIGHT,
    FIRST_PACKET_NORM,
    FIRST_PACKET_TARGET,
    HEIGHT_ZERO_CANONICAL_NORMS,
    HEIGHT_ZERO_TRACE_MINIMUM,
    ODD_LEECH_BASIS_SCALED,
    ODD_LEECH_COEFFICIENT,
    ODD_LEECH_GRAM,
    ODD_LEECH_MIRROR_NORMS,
    ODD_LEECH_NORM4_WITNESS_SCALED,
    H45_ENDPOINT,
    OddLeechGramTriple,
    bareiss_determinant,
    catalog_golay_codewords,
    iter_minimal_odd_leech_vectors_scaled,
    odd_leech_controls_through_first_packet,
    odd_leech_first_packet_audit,
    odd_leech_h45_prefix_audit,
    odd_leech_leading_principal_minors,
    odd_leech_model_audit,
    odd_leech_pair_scaled,
    odd_leech_raw_member,
    psd_extension_joint_maximum,
)


def test_catalogue_basis_is_positive_odd_unimodular() -> None:
    assert len(ODD_LEECH_BASIS_SCALED) == 24
    assert all(len(row) == 24 for row in ODD_LEECH_BASIS_SCALED)
    assert all(odd_leech_raw_member(row) for row in ODD_LEECH_BASIS_SCALED)
    assert ODD_LEECH_GRAM == tuple(zip(*ODD_LEECH_GRAM))
    assert bareiss_determinant(ODD_LEECH_GRAM) == 1
    assert odd_leech_leading_principal_minors() == (
        8,
        16,
        32,
        64,
        128,
        256,
        512,
        256,
        512,
        1024,
        2048,
        1024,
        2048,
        1024,
        512,
        256,
        512,
        256,
        128,
        64,
        32,
        16,
        8,
        1,
    )
    assert ODD_LEECH_GRAM[-1][-1] == 3


def test_catalogue_code_is_the_extended_golay_code() -> None:
    words = catalog_golay_codewords()
    assert len(words) == len(set(words)) == 4096
    assert Counter(word.bit_count() for word in words) == Counter(
        {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}
    )
    assert (1 << 24) - 1 in words


def test_complete_norm_three_shell_is_proved_not_guessed() -> None:
    vectors = tuple(iter_minimal_odd_leech_vectors_scaled())
    assert len(vectors) == len(set(vectors)) == 4096
    assert all(odd_leech_raw_member(vector) for vector in vectors)
    assert {odd_leech_pair_scaled(vector, vector) for vector in vectors} == {3}
    # In the shifted construction coset all 24 raw coordinates are odd, so
    # norm three forces all of them to be +/-1.  Their residues recover one
    # and only one codeword; the iterator realizes every such codeword.
    assert all(set(vector) <= {-1, 1} for vector in vectors)


def test_model_audit_and_independent_witness_pairings() -> None:
    audit = odd_leech_model_audit()
    assert audit.determinant == 1
    assert audit.minimal_shell_count == 4096
    assert audit.minimal_shell_digest == (
        "96bd0b37ddd6c627b5082988834a6a16fbdfa98b057da60a05b434c5ba20ad22"
    )
    assert audit.witness_square == 4
    assert audit.witness_pairing_histogram == ((-1, 1024), (0, 2048), (1, 1024))
    assert odd_leech_raw_member(ODD_LEECH_NORM4_WITNESS_SCALED)


def test_complete_local_mirror_norm_candidates_and_height_zero_exclusion() -> None:
    assert ODD_LEECH_COEFFICIENT.norm() == -23
    assert ODD_LEECH_COEFFICIENT.sign(0) > 0
    assert ODD_LEECH_COEFFICIENT.sign(1) < 0
    assert ODD_LEECH_MIRROR_NORMS == (
        Silver(1),
        Silver(2, 1),
        Silver(2),
        Silver(5, 1),
        Silver(12, 7),
        Silver(10, 2),
    )
    assert tuple(int(norm.trace()) for norm in HEIGHT_ZERO_CANONICAL_NORMS) == (
        2,
        4,
        4,
    )
    assert HEIGHT_ZERO_TRACE_MINIMUM == 6


def test_exact_control_census_through_first_packet() -> None:
    controls = odd_leech_controls_through_first_packet()
    assert len(controls) == 9
    assert sum(bool(level.necessary_gram_triples) for level in controls) == 1
    assert all(
        not level.reduced_target.is_integral
        for level in controls[:-1]
    )
    first = controls[-1]
    assert first.height == FIRST_PACKET_HEIGHT == Silver(7, 6) / 2
    assert first.norm == FIRST_PACKET_NORM == 2
    assert first.control == FIRST_PACKET_CONTROL == EPSILON
    assert first.spatial_square_target == FIRST_PACKET_TARGET == Silver(9, 6)
    assert first.necessary_gram_triples == (OddLeechGramTriple(3, 3, 3),)


def test_first_packet_is_primitive_reflective_and_vinberg_acute() -> None:
    audit = odd_leech_first_packet_audit()
    # Every root is (epsilon,epsilon*b), b^2=3.
    assert EPSILON ** 2 * (-ODD_LEECH_COEFFICIENT + 3) == 2
    assert (2 * ODD_LEECH_COEFFICIENT * EPSILON / 2).is_integral
    assert audit.first_packet_root_count == 4096
    assert audit.distinct_spatial_pairings == (1, 0, -1, -3)
    assert audit.distinct_root_pairings == (
        Silver(-4, -4),
        Silver(-7, -6),
        Silver(-10, -8),
        Silver(-16, -12),
    )
    assert all(pairing.sign(0) < 0 for pairing in audit.distinct_root_pairings)


def test_first_prefix_has_a_strict_hyperideal_patch() -> None:
    audit = odd_leech_first_packet_audit()
    assert audit.witness_square == 4
    assert audit.witness_maximum_pairing == 1
    assert audit.witness_minimum_wall_margin == Silver(-4, 4)
    assert audit.witness_minimum_wall_margin.sign(0) > 0
    assert audit.witness_hyperideal_excess == Silver(7, -4)
    assert audit.witness_hyperideal_excess.sign(0) > 0


def test_complete_h45_control_and_gram_census() -> None:
    audit = odd_leech_h45_prefix_audit()
    assert len(audit.controls) == 22
    assert len(audit.nonempty_necessary_levels) == 2
    assert tuple(level.height for level in audit.nonempty_necessary_levels) == (
        FIRST_PACKET_HEIGHT,
        H45_ENDPOINT,
    )
    assert audit.terminal_height == Silver(45, 32) / 2
    assert audit.terminal_norm == 2
    assert audit.terminal_control == Silver(3, 2)
    assert audit.terminal_target == Silver(47, 32)
    assert audit.terminal_triples == (
        OddLeechGramTriple(27, 10, 16),
        OddLeechGramTriple(25, 11, 16),
        OddLeechGramTriple(23, 12, 16),
    )


def test_norm_four_witness_survives_the_h45_necessary_superset() -> None:
    audit = odd_leech_h45_prefix_audit()
    assert audit.witness_scale == Fraction(5, 6)
    assert tuple(bound.fixed_direction_maximum_pairing for bound in audit.fixed_direction_bounds) == (
        Silver(Fraction(25, 3), 5),
        Silver(Fraction(15, 2), 5),
        Silver(Fraction(15, 2), 5),
    )
    assert {bound.wall_bound for bound in audit.fixed_direction_bounds} == {Silver(7, 6)}
    assert tuple(bound.strict_margin for bound in audit.fixed_direction_bounds) == (
        Silver(Fraction(-4, 3), 1),
        Silver(Fraction(-1, 2), 1),
        Silver(Fraction(-1, 2), 1),
    )
    assert all(bound.strict_margin.sign(0) > 0 for bound in audit.fixed_direction_bounds)
    assert audit.minimum_wall_margin == Silver(Fraction(-4, 3), 1)
    assert audit.minimum_wall_margin.sign(0) > 0
    assert audit.hyperideal_excess == Silver(Fraction(52, 9), -4)
    assert audit.hyperideal_excess.sign(0) > 0


def test_h45_psd_maxima_are_joint_not_coordinatewise() -> None:
    assert psd_extension_joint_maximum(3, 10, -4, 4, Silver(2, 1)) == (
        -2,
        6,
        Silver(10, 6),
    )
    assert psd_extension_joint_maximum(4, 11, 5, 4, EPSILON) == (
        3,
        6,
        Silver(9, 6),
    )
    assert psd_extension_joint_maximum(3, 12, 4, 4, EPSILON) == (
        3,
        6,
        Silver(9, 6),
    )

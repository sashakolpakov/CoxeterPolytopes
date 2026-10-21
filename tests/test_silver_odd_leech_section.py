from fractions import Fraction

from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_odd_leech import OddLeechGramTriple
from coxeter_polytopes.search.silver_odd_leech_section import (
    CANONICAL_MIRROR_NORMS,
    norm_six_section_shell_audit,
    sqrt2_section_first_packet_audit,
    unit_section_first_packet_audit,
)


def test_complete_section_local_and_height_zero_census():
    audit = unit_section_first_packet_audit()
    assert audit.model.determinant == 3
    assert audit.model.smith_invariants == (1,) * 22 + (3,)
    assert audit.model.minimum == 3
    assert audit.model.minimal_shell_count == 2576
    assert CANONICAL_MIRROR_NORMS == (
        Silver(1), Silver(2, 1), Silver(2),
        Silver(3), Silver(6, 3), Silver(6),
    )
    assert [len(types) for _norm, types in audit.height_zero_ambient_types] == [
        0, 0, 0, 1, 0, 2
    ]
    assert all(not types for _norm, types in audit.height_zero_gated_types)


def test_complete_norm_six_shell_census():
    shell = norm_six_section_shell_audit()
    assert shell.odd_coset_count == 0
    assert shell.even_weight_zero_count == 0
    assert shell.even_weight_eight_count == 1_360_128
    assert shell.even_weight_twelve_count == 2_380_224
    assert shell.total_count == 3_740_352
    assert shell.witness_pairing_histogram == (
        (-3, 12320), (-2, 229152), (-1, 924000), (0, 1409408),
        (1, 924000), (2, 229152), (3, 12320),
    )


def test_unit_first_packet_and_strict_hyperideal_witness():
    audit = unit_section_first_packet_audit()
    nonempty = [level for level in audit.controls if level.necessary_gram_triples]
    assert len(audit.controls) == 6
    assert len(nonempty) == 1
    assert nonempty[0].necessary_gram_triples == (OddLeechGramTriple(3, 3, 3),)
    assert audit.first_height == Silver(2, Fraction(3, 2))
    assert audit.first_packet_norm_counts == ((Silver(2, 1), 2576),)
    assert audit.accepted_wall_count == 2576
    assert audit.witness_spatial_scale == 1
    assert audit.witness_square == Silver(3, -1)
    assert audit.closest_wall_pairing == Silver(-2, -1)
    assert audit.strict_hyperideal and not audit.finite_volume and not audit.compact


def test_sqrt2_simultaneous_first_packet_and_strict_witness():
    audit = sqrt2_section_first_packet_audit()
    nonempty = [level for level in audit.controls if level.necessary_gram_triples]
    assert len(audit.controls) == 10
    assert len(nonempty) == 2
    assert {level.norm for level in nonempty} == {Silver(1), Silver(2)}
    assert all(level.height == Silver(8, 6) for level in nonempty)
    assert audit.first_packet_norm_counts == (
        (Silver(1), 2576),
        (Silver(2), 3_740_352),
    )
    assert audit.imprimitive_candidate_count == 2576
    assert audit.accepted_wall_count == 3_742_928
    assert audit.witness_spatial_scale == Fraction(2, 3)
    assert audit.witness_square == Silver(Fraction(16, 9), -1)
    assert audit.closest_wall_pairing == Silver(-2)
    assert all(bound.sign(0) < 0 for _label, bound in audit.packet_pairing_upper_bounds)
    assert audit.coefficient == SQRT2
    assert audit.strict_hyperideal and not audit.finite_volume and not audit.compact

from collections import Counter
from fractions import Fraction

from coxeter_polytopes.exact.silver import SILVER_ZERO, Silver
from coxeter_polytopes.search.bugaenko_vinberg import _fraction_inverse
from coxeter_polytopes.search.niemeier_host_triage import cartan_gram, highest_root_marks
from coxeter_polytopes.search.silver_a15_d9 import (
    A15_COMPONENT,
    COMPLETE_UNIT_NORM_PROOF,
    D9_COMPONENT,
    FIRST_PACKET_HEIGHT,
    SECOND_PACKET_HEIGHT,
    a15_d9_basis_audit,
    a15_d9_glue_membership,
    a15_d9_pair,
    axis_witness,
    enumerate_control_level,
    enumerate_dominant_a15_weights,
    enumerate_dominant_d9_weights,
    first_packet_geometry,
    glue_coset_minima,
    height_zero_dominant_roots,
    run_a15_d9_prefix,
)


def test_index_eight_glue_is_explicit_even_and_unimodular():
    assert a15_d9_basis_audit()
    a = tuple(Silver(1) if index == 1 else SILVER_ZERO for index in range(15))
    d = tuple(Silver(1) if index == 8 else SILVER_ZERO for index in range(9))
    assert a15_d9_glue_membership(a, d)
    inverse_a = _fraction_inverse(cartan_gram(A15_COMPONENT))
    inverse_d = _fraction_inverse(cartan_gram(D9_COMPONENT))
    assert inverse_a[1][1] == Fraction(7, 4)
    assert inverse_d[8][8] == Fraction(9, 4)
    assert inverse_a[1][1] + inverse_d[8][8] == 4
    assert all((8 * inverse_a[index][1]).denominator == 1 for index in range(15))
    assert 8 * inverse_a[14][1] == 1
    assert all((8 * inverse_d[index][8]).denominator == 1 for index in range(9))
    glue = {(2 * multiple % 16, multiple % 4) for multiple in range(8)}
    assert len(glue) == 8
    assert glue_coset_minima() == (4, 4, 6, 4, 6, 4, 4)
    assert "complete canonical positive norms are 1 and 2" in COMPLETE_UNIT_NORM_PROOF


def test_height_zero_and_controls_are_complete():
    assert highest_root_marks(A15_COMPONENT) == (1,) * 15
    assert highest_root_marks(D9_COMPONENT) == (1, 2, 2, 2, 2, 2, 2, 1, 1)
    roots = height_zero_dominant_roots()
    assert len(roots) == 2
    assert all(root.norm == 2 and root.control == 0 for root in roots)
    assert {(root.a15.square, root.d9.square) for root in roots} == {
        (Silver(2), Silver(0)),
        (Silver(0), Silver(2)),
    }
    prefix = run_a15_d9_prefix(15)
    assert len(prefix.controls) == 8
    assert [
        (packet.height, packet.candidate_count, packet.rejected_count, len(packet.roots))
        for packet in prefix.packets
    ] == [
        (Silver(Fraction(1, 2), Fraction(1, 2)), 0, 0, 0),
        (Silver(1, 1), 0, 0, 0),
        (FIRST_PACKET_HEIGHT, 2, 0, 2),
        (Silver(Fraction(7, 2), Fraction(5, 2)), 0, 0, 0),
        (SECOND_PACKET_HEIGHT, 15, 0, 15),
    ]


def test_second_packet_glue_and_dynamic_program_are_frozen():
    prefix = run_a15_d9_prefix(15)
    second = next(packet for packet in prefix.packets if packet.height == SECOND_PACKET_HEIGHT)
    level = next(level for level in second.controls if level.norm == 2)
    d9 = enumerate_dominant_d9_weights(level.target, level.norm, level.control)
    d9_vectors = tuple(
        vector
        for vector in d9.vectors
        if all((level.target - vector.square).sign(place) >= 0 for place in (0, 1))
    )
    targets = tuple(sorted({level.target - vector.square for vector in d9_vectors}))
    a15 = enumerate_dominant_a15_weights(targets, level.norm, level.control)
    assert (len(d9.vectors), len(d9_vectors), len(targets)) == (48, 48, 23)
    assert len(a15.vectors) == 24
    assert a15.leaf_count == 222
    assert sum(a15.visited_by_depth) == 10259
    assert sum(a15.cap_pruned_by_depth) == 19849
    assert sum(a15.energy_pruned_by_depth) == 7300
    assert sum(root.glue_nontrivial for root in second.roots) == 9
    assert sum(root.cross_component_glue for root in second.roots) == 8
    assert Counter((root.a15.residue, root.d9.residue) for root in second.roots) == {
        ((0, 0), (0, 0)): 6,
        ((14, 14), (3, 3)): 2,
        ((12, 12), (2, 2)): 2,
        ((4, 4), (2, 2)): 2,
        ((2, 2), (1, 1)): 2,
        ((8, 8), (0, 0)): 1,
    }
    assert enumerate_control_level(level, True) == second.roots


def test_all_walls_are_exact_and_chronologically_acute():
    prefix = run_a15_d9_prefix(15)
    assert len(prefix.accepted_roots) == 17
    assert all(a15_d9_pair(root, root) == root.norm for root in prefix.accepted_roots)
    assert all(
        a15_d9_pair(left, right).sign(0) <= 0
        for index, left in enumerate(prefix.accepted_roots)
        for right in prefix.accepted_roots[index + 1 :]
    )


def test_geometry_and_strict_d9_axis_witness():
    geometry = first_packet_geometry()
    assert (
        geometry.extreme_ray_count,
        geometry.finite_count,
        geometry.cusp_count,
        geometry.hyperideal_count,
    ) == (160, 2, 0, 158)
    witness = axis_witness("D9", 4, Fraction(7, 10))
    assert witness.norm == Silver(Fraction(29, 20), -1)
    assert witness.height_zero_nonpositive
    assert witness.all_positive_walls_strict
    assert witness.hyperideal
    assert len(witness.positive_wall_pairings) == 17
    assert max(witness.positive_wall_pairings) == Silver(Fraction(1, 5), Fraction(-1, 5))


def test_complete_h50_third_packet_and_surviving_witness():
    prefix = run_a15_d9_prefix(50)
    assert len(prefix.controls) == 13
    assert [len(packet.roots) for packet in prefix.packets if packet.roots] == [2, 15, 191]
    third = next(packet for packet in prefix.packets if packet.height == Silver(14, 10))
    assert (third.candidate_count, third.rejected_count, len(third.roots)) == (197, 6, 191)
    assert sum(root.glue_nontrivial for root in third.roots) == 140
    assert sum(root.cross_component_glue for root in third.roots) == 112
    pairings = tuple(
        a15_d9_pair(left, right)
        for index, left in enumerate(third.roots)
        for right in third.roots[index + 1 :]
    )
    assert Counter(pairing.sign(0) for pairing in pairings) == {-1: 18121, 0: 24}
    final = prefix.packets[-1]
    assert final.height == Silver(Fraction(41, 2), Fraction(29, 2))
    assert final.candidate_count == len(final.roots) == 0
    witness = axis_witness("D9", 4, Fraction(7, 10), 50)
    assert witness.norm == Silver(Fraction(29, 20), -1)
    assert len(witness.positive_wall_pairings) == 208
    assert witness.all_positive_walls_strict and witness.hyperideal

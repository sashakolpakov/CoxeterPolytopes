from collections import Counter
from fractions import Fraction

from coxeter_polytopes.exact.silver import SILVER_ZERO, Silver
from coxeter_polytopes.search.bugaenko_vinberg import _fraction_inverse
from coxeter_polytopes.search.niemeier_host_triage import cartan_gram
from coxeter_polytopes.search.silver_a17_e7 import (
    A17_COMPONENT,
    COMPLETE_UNIT_NORM_PROOF,
    E7_COMPONENT,
    FIRST_PACKET_HEIGHT,
    SECOND_PACKET_HEIGHT,
    a17_e7_basis_audit,
    a17_e7_glue_membership,
    a17_e7_pair,
    enumerate_control_level,
    enumerate_dominant_a17_weights,
    enumerate_dominant_e7_weights,
    first_packet_geometry,
    height_zero_dominant_roots,
    run_a17_e7_prefix,
    second_packet_hyperideal_witness,
)


def test_index_six_glue_is_even_unimodular_and_explicit():
    assert a17_e7_basis_audit()
    a_weights = tuple(Silver(1) if index == 2 else SILVER_ZERO for index in range(17))
    e_weights = tuple(Silver(1) if index == 5 else SILVER_ZERO for index in range(7))
    assert a17_e7_glue_membership(a_weights, e_weights)
    inverse_a = _fraction_inverse(cartan_gram(A17_COMPONENT))
    inverse_e = _fraction_inverse(cartan_gram(E7_COMPONENT))
    assert inverse_a[2][2] == Fraction(5, 2)
    assert inverse_e[5][5] == Fraction(3, 2)
    assert inverse_a[2][2] + inverse_e[5][5] == 4
    assert all((6 * inverse_a[index][2]).denominator == 1 for index in range(17))
    assert (6 * inverse_a[16][2]) == 1
    assert all((6 * inverse_e[index][5]).denominator == 1 for index in range(7))
    glue = {(3 * multiple % 18, multiple % 2) for multiple in range(6)}
    assert len(glue) == 6
    assert all((Fraction(17 * a * a, 18) + Fraction(3 * e * e, 2)) % 2 == 0 for a, e in glue)
    assert "complete canonical positive norms are 1 and 2" in COMPLETE_UNIT_NORM_PROOF


def test_height_zero_and_controls_are_complete():
    roots = height_zero_dominant_roots()
    assert len(roots) == 2
    assert all(root.norm == 2 and root.control == 0 for root in roots)
    assert {(root.a17.square, root.e7.square) for root in roots} == {(Silver(2), Silver(0)), (Silver(0), Silver(2))}
    prefix = run_a17_e7_prefix(15)
    assert len(prefix.controls) == 8
    assert [
        (packet.height, packet.candidate_count, packet.rejected_count, len(packet.roots))
        for packet in prefix.packets
    ] == [
        (Silver(Fraction(1, 2), Fraction(1, 2)), 0, 0, 0),
        (Silver(1, 1), 0, 0, 0),
        (FIRST_PACKET_HEIGHT, 2, 0, 2),
        (Silver(Fraction(7, 2), Fraction(5, 2)), 0, 0, 0),
        (SECOND_PACKET_HEIGHT, 12, 0, 12),
    ]


def test_rank17_dynamic_program_and_second_packet_glue_are_frozen():
    prefix = run_a17_e7_prefix(15)
    second = next(packet for packet in prefix.packets if packet.height == SECOND_PACKET_HEIGHT)
    level = next(level for level in second.controls if level.norm == 2)
    e7 = enumerate_dominant_e7_weights(level.target, level.norm, level.control)
    e7_vectors = tuple(
        vector
        for vector in e7.vectors
        if all((level.target - vector.square).sign(place) >= 0 for place in (0, 1))
    )
    targets = tuple(sorted({level.target - vector.square for vector in e7_vectors}))
    a17 = enumerate_dominant_a17_weights(targets, level.norm, level.control)
    assert len(e7.vectors) == 18
    assert len(a17.vectors) == 12
    assert a17.leaf_count == 205
    assert sum(a17.visited_by_depth) == 11992
    assert sum(a17.cap_pruned_by_depth) == 22824
    assert sum(a17.energy_pruned_by_depth) == 8622
    assert sum(root.glue_nontrivial for root in second.roots) == 6
    assert sum(root.cross_component_glue for root in second.roots) == 4
    assert Counter((root.a17.residue, root.e7.residue) for root in second.roots) == {
        ((0, 0), (0, 0)): 6,
        ((15, 15), (1, 1)): 2,
        ((3, 3), (1, 1)): 2,
        ((12, 12), (0, 0)): 1,
        ((6, 6), (0, 0)): 1,
    }
    assert enumerate_control_level(level, True) == second.roots


def test_all_walls_are_exact_roots_and_chronologically_acute():
    prefix = run_a17_e7_prefix(15)
    assert len(prefix.accepted_roots) == 14
    assert all(a17_e7_pair(root, root) == root.norm for root in prefix.accepted_roots)
    assert all(
        a17_e7_pair(left, right).sign(0) <= 0
        for index, left in enumerate(prefix.accepted_roots)
        for right in prefix.accepted_roots[index + 1 :]
    )


def test_geometry_and_strict_hyperideal_witness():
    geometry = first_packet_geometry()
    assert (
        geometry.extreme_ray_count,
        geometry.finite_count,
        geometry.cusp_count,
        geometry.hyperideal_count,
    ) == (144, 2, 0, 142)
    witness = second_packet_hyperideal_witness()
    assert witness.norm == Silver(Fraction(529, 200), -1)
    assert witness.height_zero_nonpositive
    assert witness.all_positive_walls_strict
    assert witness.hyperideal
    assert len(witness.positive_wall_pairings) == 14
    assert Counter(str(pairing) for pairing in witness.positive_wall_pairings) == {
        "-2 - 2*sqrt(2)": 1,
        "-11/10 - 11/10*sqrt(2)": 1,
        "-4 - 3*sqrt(2)": 1,
        "-53/20 - 33/20*sqrt(2)": 2,
        "-13/10 - 6/5*sqrt(2)": 2,
        "-7/4 - 33/20*sqrt(2)": 2,
        "-2/5 - 3/10*sqrt(2)": 2,
        "-11/5 - 21/10*sqrt(2)": 2,
        "-31/10 - 21/10*sqrt(2)": 1,
    }

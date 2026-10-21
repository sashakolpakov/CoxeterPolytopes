from collections import Counter
from fractions import Fraction

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_a15_d9 import a15_d9_basis_audit
from coxeter_polytopes.search.silver_p7_a15_d9 import (
    COMPLETE_P7_A15_D9_NORM_PROOF,
    FIRST_PACKET_HEIGHT,
    P7_NORM_CLASSES,
    SECOND_PACKET_HEIGHT,
    _same_height_acute,
    p7_a15_d9_pair,
    p7_d9_axis_witness,
    run_p7_a15_d9_prefix,
)
from coxeter_polytopes.search.silver_p7_d24_plus import P7_DELTA, P7_TIMELIKE


def test_prime_seven_local_data_and_controls_are_complete():
    assert a15_d9_basis_audit()
    assert P7_TIMELIKE == Silver(-1, 2)
    assert P7_DELTA == Silver(3, 1)
    assert P7_NORM_CLASSES == (Silver(1), Silver(2), P7_DELTA, 2 * P7_DELTA)
    assert "complete canonical norms are 1, 2, delta, and 2*delta" in COMPLETE_P7_A15_D9_NORM_PROOF
    prefix = run_p7_a15_d9_prefix(63)
    assert len(prefix.controls) == 21
    assert len(prefix.packets) == 14
    assert prefix.packets[-1].height == SECOND_PACKET_HEIGHT


def test_first_and_dense_second_packets_are_frozen():
    prefix = run_p7_a15_d9_prefix(63)
    nonempty = [packet for packet in prefix.packets if packet.roots]
    assert [packet.height for packet in nonempty] == [FIRST_PACKET_HEIGHT, SECOND_PACKET_HEIGHT]
    assert [
        (packet.candidate_count, packet.rejected_count, len(packet.roots))
        for packet in nonempty
    ] == [(11, 0, 11), (1670, 15, 1655)]
    assert all(
        packet.candidate_count == 0
        for packet in prefix.packets
        if FIRST_PACKET_HEIGHT < packet.height < SECOND_PACKET_HEIGHT
    )
    second = nonempty[1]
    assert sum(root.glue_nontrivial for root in second.roots) == 1412
    assert sum(root.cross_component_glue for root in second.roots) == 1252


def test_all_1666_walls_are_exact_and_simultaneously_acute():
    prefix = run_p7_a15_d9_prefix(63)
    assert len(prefix.accepted_roots) == 1666
    assert all(p7_a15_d9_pair(root, root) == root.norm for root in prefix.accepted_roots)
    first = next(packet for packet in prefix.packets if packet.height == FIRST_PACKET_HEIGHT)
    second = next(packet for packet in prefix.packets if packet.height == SECOND_PACKET_HEIGHT)
    old_signs = Counter(
        p7_a15_d9_pair(root, wall).sign(0)
        for root in second.roots
        for wall in first.roots
    )
    assert old_signs == {-1: 18187, 0: 18}
    assert _same_height_acute(second.roots)


def test_strict_hyperideal_witness_survives_dense_packet():
    witness = p7_d9_axis_witness(height_bound=63)
    assert witness.norm == Silver(Fraction(1808, 625), -2)
    assert witness.height_zero_nonpositive
    assert witness.all_positive_walls_strict and witness.hyperideal
    assert len(witness.positive_wall_pairings) == 1666
    assert Counter(pairing.sign(0) for pairing in witness.positive_wall_pairings) == {-1: 1666}
    assert max(witness.positive_wall_pairings) == Silver(Fraction(-18, 25), Fraction(1, 5))


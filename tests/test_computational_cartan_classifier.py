"""Exact regression tests for the finite Cartan-cycle classifier."""

from __future__ import annotations

from fractions import Fraction

import pytest

from coxeter_polytopes.classification.cartan_packets import (
    CyclePacket,
    PacketError,
    analyze_simplex_gram,
    cartan_cycle_packet,
    cartan_gram,
    diagonal_rescale,
    enumerate_q_simplex_candidates,
    integral_cartan_gauge,
    reconstruct_candidate,
    reconstruct_cartan,
    symmetrizer,
    verify_q_simplex_candidate,
)
from coxeter_polytopes.exact.matrix import Inertia


def _matrix(rows):
    return tuple(tuple(Fraction(value) for value in row) for row in rows)


def test_complete_cycle_packet_recovers_cartan_up_to_diagonal_gauge():
    cartan = _matrix(
        (
            (2, -1, -1, 0),
            (-1, 2, -1, -1),
            (-1, -1, 2, -1),
            (0, -1, -1, 2),
        )
    )
    changed = diagonal_rescale(cartan, (2, 3, 5, 7))
    packet = cartan_cycle_packet(changed)
    recovered = reconstruct_cartan(packet)
    assert cartan_cycle_packet(recovered) == packet
    assert symmetrizer(recovered)
    assert analyze_simplex_gram(cartan_gram(recovered)).full_inertia == Inertia(3, 1, 0)


def test_integral_cycle_products_construct_an_integral_cartan_gauge():
    integral = _matrix(
        (
            (2, -1, -1, 0),
            (-2, 2, -2, -1),
            (-2, -1, 2, -1),
            (0, -2, -2, 2),
        )
    )
    rational = diagonal_rescale(integral, (2, 3, 5, 7))
    gauged, scales = integral_cartan_gauge(rational)
    assert all(value.denominator == 1 for row in gauged for value in row)
    assert cartan_cycle_packet(gauged) == cartan_cycle_packet(rational)
    assert all(scale > 0 for scale in scales)


def test_nonintegral_closed_cycle_is_an_exact_integral_gauge_obstruction():
    cartan = _matrix(((2, -1, -1), (-1, 2, -1), (-1, Fraction(-1, 2), 2)))
    with pytest.raises(PacketError, match="nonintegral closed cycle"):
        integral_cartan_gauge(cartan)


def test_finite_area_rational_triangle_emits_and_replays_a_chamber_certificate():
    # Coxeter triangle (2,3,infinity): edge products 0,1,4.
    packet = CyclePacket.from_edge_products(3, {(0, 2): 1, (1, 2): 4})
    candidate = reconstruct_candidate(packet)
    assert candidate.analysis.finite_volume is True
    assert candidate.analysis.compact is False
    assert sorted(candidate.analysis.vertex_kinds) == ["finite", "finite", "ideal"]
    replay = verify_q_simplex_candidate(candidate)
    assert replay.dimension == 2
    assert len(replay.cusps) == 1
    assert not replay.hyperideal_rays


def test_compact_and_euclidean_triangle_are_separated_by_exact_inertia():
    compact = reconstruct_candidate(
        CyclePacket.from_edge_products(3, {(0, 2): 2, (1, 2): 3})
    )  # (2,4,6)
    euclidean = reconstruct_candidate(
        CyclePacket.from_edge_products(3, {(0, 2): 1, (1, 2): 3})
    )  # (2,3,6)
    assert compact.analysis.compact is True
    assert compact.analysis.full_inertia == Inertia(2, 1, 0)
    assert euclidean.analysis.lorentzian is False
    assert euclidean.analysis.full_inertia == Inertia(2, 0, 1)


def test_rank_three_box_is_exhaustive_modulo_gauge_and_facet_permutation():
    census = enumerate_q_simplex_candidates(3, (1, 2, 3, 4))
    assert census.assignment_count == 125
    assert census.connected_assignment_count == 112
    assert census.integral_q_packet_assignment_count == 68
    assert census.unlabelled_packet_count == 18
    assert census.lorentzian_packet_count == 13
    assert len(census.finite_volume_candidates) == 13
    assert census.independently_verified_count == 13
    assert census.compact_count == 4
    assert census.noncompact_count == 9

    edge_to_label = {Fraction(0): 2, Fraction(1): 3, Fraction(2): 4,
                     Fraction(3): 6, Fraction(4): 10**9}
    label_triples = {
        tuple(
            sorted(
                edge_to_label[candidate.packet.edge_products.get((i, j), Fraction(0))]
                for i in range(3)
                for j in range(i + 1, 3)
            )
        )
        for candidate in census.finite_volume_candidates
    }
    infinity = 10**9
    assert label_triples == {
        (2, 3, infinity),
        (2, 4, 6),
        (2, 4, infinity),
        (2, 6, 6),
        (2, 6, infinity),
        (2, infinity, infinity),
        (3, 3, infinity),
        (3, 4, 4),
        (3, 6, 6),
        (3, infinity, infinity),
        (4, 4, infinity),
        (6, 6, infinity),
        (infinity, infinity, infinity),
    }


def test_rank_four_rational_label_box_has_eighteen_verified_tetrahedra():
    census = enumerate_q_simplex_candidates(4, (1, 2, 3, 4))
    assert census.assignment_count == 15625
    assert census.connected_assignment_count == 15104
    assert census.integral_q_packet_assignment_count == 2976
    assert census.unlabelled_packet_count == 207
    assert census.lorentzian_packet_count == 200
    assert len(census.finite_volume_candidates) == 18
    assert census.compact_count == 1
    assert census.noncompact_count == 17
    assert census.independently_verified_count == 18


def test_packet_canonicalization_identifies_all_facet_permutations():
    packet = CyclePacket.from_edge_products(
        4,
        {(0, 1): 1, (1, 2): 2, (2, 3): 3},
    )
    relabelled = packet.relabel((2, 0, 3, 1))
    assert packet.canonicalized() == relabelled.canonicalized()


def test_nonsquare_cycle_product_is_routed_out_of_the_Q_field_lane():
    with pytest.raises(PacketError, match="field is not Q"):
        CyclePacket.from_edge_products(3, {(0, 1): 1, (0, 2): 1, (1, 2): 2})


def test_search_cap_requires_explicit_partition_for_large_boxes():
    with pytest.raises(ValueError, match="partition the edge box"):
        enumerate_q_simplex_candidates(5, (1, 2, 3, 4), max_assignments=1000)

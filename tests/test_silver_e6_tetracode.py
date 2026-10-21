from fractions import Fraction

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.niemeier_host_triage import highest_root_marks
from coxeter_polytopes.search.silver_e6_tetracode import (
    COMPLETE_UNIT_NORM_PROOF,
    E6_HIGHEST_ROOT,
    SECOND_PACKET_HEIGHT,
    e6_cartan_gram,
    e6_inverse_cartan,
    first_packet_geometry,
    height_zero_dominant_roots,
    run_e6_tetracode_prefix,
    second_packet_hyperideal_witness,
    silver_e6_pair,
    ternary_tetracode,
    tetracode_membership,
)


def _determinant(matrix):
    work = [list(row) for row in matrix]
    result = Fraction(1)
    for column in range(len(work)):
        pivot = next(row for row in range(column, len(work)) if work[row][column])
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result *= -1
        value = work[column][column]
        result *= value
        for row in range(column + 1, len(work)):
            factor = work[row][column] / value
            for index in range(column + 1, len(work)):
                work[row][index] -= factor * work[column][index]
    return result


def test_e6_discriminant_and_tetracode_data_are_exact():
    gram = e6_cartan_gram()
    inverse = e6_inverse_cartan()
    assert _determinant(gram) == 3
    assert tuple(
        tuple(sum(gram[i][k] * inverse[k][j] for k in range(6)) for j in range(6))
        for i in range(6)
    ) == tuple(
        tuple(Fraction(i == j) for j in range(6)) for i in range(6)
    )
    code = ternary_tetracode()
    assert len(code) == 9
    assert sorted(sum(value != 0 for value in word) for word in code) == [0] + [3] * 8
    assert all(sum(x * y for x, y in zip(left, right)) % 3 == 0 for left in code for right in code)
    assert E6_HIGHEST_ROOT == highest_root_marks(("E", 6))
    assert "complete canonical positive norms are 1 and 2" in COMPLETE_UNIT_NORM_PROOF


def test_height_zero_is_exactly_e6_four():
    roots = height_zero_dominant_roots()
    assert len(roots) == 4
    assert all(root.norm == 2 and root.control == 0 for root in roots)
    assert all(tetracode_membership(tuple(block.weights for block in root.blocks)) for root in roots)


def test_two_nonempty_packets_and_coupled_second_packet():
    prefix = run_e6_tetracode_prefix(15)
    assert len(prefix.controls) == 8
    assert [
        (packet.height, packet.candidate_count, packet.rejected_count, len(packet.roots))
        for packet in prefix.packets
    ] == [
        (Silver(Fraction(1, 2), Fraction(1, 2)), 0, 0, 0),
        (Silver(1, 1), 0, 0, 0),
        (Silver(2, 2), 4, 0, 4),
        (Silver(Fraction(7, 2), Fraction(5, 2)), 0, 0, 0),
        (SECOND_PACKET_HEIGHT, 44, 0, 44),
    ]
    assert len(prefix.accepted_roots) == 48
    second = next(packet for packet in prefix.packets if packet.height == SECOND_PACKET_HEIGHT)
    assert sum(root.glue_coupled for root in second.roots) == 24
    assert sum(not root.glue_coupled for root in second.roots) == 20
    assert all(silver_e6_pair(root.vector, root.vector) == root.norm for root in prefix.accepted_roots)
    assert all(
        silver_e6_pair(left.vector, right.vector).sign(0) <= 0
        for index, left in enumerate(second.roots)
        for right in second.roots[index + 1 :]
    )


def test_first_geometry_and_second_packet_witness_are_hyperideal():
    geometry = first_packet_geometry()
    assert (
        geometry.extreme_ray_count,
        geometry.finite_vertex_count,
        geometry.cusp_count,
        geometry.hyperideal_ray_count,
        geometry.compact,
    ) == (2401, 5, 0, 2396, False)
    witness = second_packet_hyperideal_witness()
    assert witness.norm == Silver(Fraction(47, 25), -1)
    assert witness.height_zero_nonpositive
    assert witness.all_positive_walls_strict
    assert witness.hyperideal
    assert len(witness.positive_wall_pairings) == 48


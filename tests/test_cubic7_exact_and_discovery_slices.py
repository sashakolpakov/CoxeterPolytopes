from fractions import Fraction

from coxeter_polytopes.classification.cubic7_packets import (
    PRISM_VERTICES,
    enumerate_equal_norm_compact_simplices,
    enumerate_symmetric_compact_prisms,
)
from coxeter_polytopes.exact.cubic7 import (
    CUBIC7_FIELD_DISCRIMINANT,
    CUBIC7_ONE,
    LABEL7_PRODUCT,
    LABEL14_PRODUCT,
    RHO7,
    THETA7,
    Cubic7,
    cubic7_determinant,
    cubic7_inertia,
)
from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.search.cubic7_diagonal import (
    bounded_cubic7_integers,
    diagonal_pair,
    enumerate_diagonal_root_pool,
    run_diagonal_vinberg_search,
)


def _pair(vector, gram, other):
    return sum(
        (
            vector[i] * gram[i][j] * other[j]
            for i in range(len(vector))
            for j in range(len(other))
        ),
        start=Cubic7(0),
    )


def _submatrix(matrix, indices):
    return tuple(tuple(matrix[i][j] for j in indices) for i in indices)


def test_cubic7_field_model_and_coxeter_invariants_are_exact():
    assert THETA7**3 + THETA7**2 - 2 * THETA7 - 1 == 0
    assert THETA7.conjugate() == THETA7**2 - 2
    assert THETA7.conjugate(2) == 1 - THETA7 - THETA7**2
    assert THETA7.trace() == -1
    assert THETA7.norm() == 1
    assert [THETA7.sign(i) for i in range(3)] == [1, -1, -1]

    assert RHO7 == THETA7 + THETA7**2 - 1
    assert RHO7**2 == LABEL7_PRODUCT == 2 + THETA7
    assert LABEL14_PRODUCT == 2 + RHO7
    assert RHO7.norm() == -1
    assert CUBIC7_FIELD_DISCRIMINANT == 49


def test_cubic7_inverse_json_and_exact_signs():
    values = (
        THETA7,
        RHO7,
        1 + THETA7,
        Cubic7(Fraction(2, 3), Fraction(-5, 7), Fraction(11, 13)),
    )
    for value in values:
        assert value * value.inverse() == CUBIC7_ONE
        assert Cubic7.from_json(value.to_json()) == value
        assert all(value.sign(embedding) in (-1, 0, 1) for embedding in range(3))


def test_equal_norm_cubic_simplex_box_is_exhausted_in_ranks_three_and_four():
    triangles = enumerate_equal_norm_compact_simplices(3)
    assert [packet.labels for packet in triangles] == [
        (2, 3, 7),
        (2, 7, 7),
        (3, 3, 7),
        (7, 7, 7),
    ]
    assert enumerate_equal_norm_compact_simplices(4) == ()
    for packet in triangles:
        assert cubic7_inertia(packet.gram, embedding=0) == Inertia(2, 1, 0)
        assert all(
            cubic7_inertia(packet.gram, embedding=embedding) == Inertia(3, 0, 0)
            for embedding in (1, 2)
        )


def test_symmetric_prism_slice_has_two_exact_arithmetic_doubles():
    packets = enumerate_symmetric_compact_prisms()
    assert [
        (packet.base_labels, packet.lateral_labels, packet.dotted_entry)
        for packet in packets
    ] == [
        ((2, 2, 3), (7, 2, 3), (1 + THETA7) ** 2),
        ((2, 2, 3), (7, 3, 3), 1 + THETA7),
    ]

    for packet in packets:
        gram = packet.gram
        assert cubic7_determinant(gram) == 0
        assert cubic7_inertia(gram, embedding=0) == Inertia(3, 1, 1)
        assert all(
            cubic7_inertia(gram, embedding=embedding) == Inertia(4, 0, 1)
            for embedding in (1, 2)
        )
        assert all(
            cubic7_inertia(_submatrix(gram, vertex), embedding=0)
            == Inertia(3, 0, 0)
            for vertex in PRISM_VERTICES
        )
        assert all(value.is_integral for value in packet.radical_vector)
        assert all(
            sum(
                (gram[i][j] * packet.radical_vector[j] for j in range(5)),
                start=Cubic7(0),
            )
            == 0
            for i in range(5)
        )

        # The extra integral reflection is explicit: it interchanges the two
        # bases and fixes every lateral root.  Hence these are doubles, not
        # candidates for new maximal reflection groups.
        swap = packet.base_swap_root
        swap_norm = _pair(swap, gram, swap)
        assert swap_norm == 2 * (2 + packet.dotted_entry)
        reflection_coefficients = tuple(
            2
            * _pair(
                tuple(Cubic7(int(i == j)) for i in range(5)),
                gram,
                swap,
            )
            / swap_norm
            for j in range(5)
        )
        assert reflection_coefficients == (
            CUBIC7_ONE,
            -CUBIC7_ONE,
            Cubic7(0),
            Cubic7(0),
            Cubic7(0),
        )


def test_diagonal_cubic_lattice_root_pool_is_exact_through_height_ten():
    assert bounded_cubic7_integers((1, 1, 1)) == (
        Cubic7(-1),
        Cubic7(0),
        Cubic7(1),
    )

    pool = enumerate_diagonal_root_pool(dimension=3, height_bound=10)
    assert pool.control_counts_by_norm == ((1, 1), (2, 4))
    assert pool.spatial_alphabet_counts_by_norm == ((1, 5), (2, 9))
    assert len(pool.cone_roots) == 3
    assert len(pool.roots) == 5
    assert {root.control for root in pool.roots} == {
        RHO7,
        Cubic7(1),
        1 + THETA7,
        2 + THETA7,
    }
    assert {root.norm for root in pool.roots} == {1, 2}
    for root in pool.roots:
        assert diagonal_pair(root.vector, root.vector) == root.norm
        assert root.height == THETA7 * root.control**2 / root.norm
        assert (Cubic7(10) - root.height).sign(0) >= 0
        assert root.control.sign(0) > 0
        assert all(coordinate.sign(0) >= 0 for coordinate in root.vector[1:])
        assert all(
            (root.vector[i] - root.vector[i + 1]).sign(0) >= 0
            for i in range(1, 3)
        )


def test_diagonal_rank_four_lattice_has_a_complete_compact_vinberg_prism():
    search = run_diagonal_vinberg_search(dimension=3, height_bound=7)
    assert len(search.pool.roots) == 5
    assert len(search.accepted_roots) == 2
    assert len(search.rejected_roots) == 3
    assert len(search.walls) == 5
    assert search.compact is True
    assert len(search.compact_prism_vertices) == 6

    norms = tuple(diagonal_pair(root, root) for root in search.walls)
    assert norms == (2, 2, 1, 2, 2)

    def edge_product(i, j):
        return (
            4
            * search.wall_gram[i][j] ** 2
            / (search.wall_gram[i][i] * search.wall_gram[j][j])
        )

    assert edge_product(0, 1) == 1  # label 3
    assert edge_product(1, 2) == 2  # label 4
    assert edge_product(0, 3) == LABEL7_PRODUCT  # label 7
    assert edge_product(2, 4) == 2 * (1 + THETA7) ** 2
    assert (edge_product(2, 4) - 4).sign(0) > 0  # divergent bases
    assert all(
        search.wall_gram[i][j] == 0
        for i, j in ((0, 2), (0, 4), (1, 3), (1, 4), (2, 3), (3, 4))
    )

    assert search.accepted_roots[-1].height == (
        Fraction(1, 2) + 3 * THETA7 + Fraction(3, 2) * THETA7**2
    )
    assert (Cubic7(7) - search.accepted_roots[-1].height).sign(0) > 0

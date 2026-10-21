import sympy as sp

from scripts.certify_project_native_minimal_surface_bound import (
    WEHLER_GRAM,
    WEHLER_ROOTS,
    certify,
    complex_dimension_ceiling,
    pairing,
    reflection_matrix,
    trace_signature,
)


def test_full_exact_certificate():
    certify()


def test_restriction_of_scalars_rank_and_signature():
    for degree in range(1, 9):
        for dimension in range(2, 12):
            positive, negative = trace_signature(degree, dimension)
            assert positive == 1
            assert positive + negative == degree * (dimension + 1)
            assert complex_dimension_ceiling(degree) == 20 // degree - 1


def test_wehler_deck_actions_are_integral_reflections():
    identity = sp.eye(3)
    for root in WEHLER_ROOTS:
        action = reflection_matrix(root)
        assert action**2 == identity
        assert action.T * WEHLER_GRAM * action == WEHLER_GRAM
        assert action * root == -root


def test_wehler_mirrors_have_three_distinct_ideal_vertices():
    basis = [sp.eye(3).col(i) for i in range(3)]
    for index, vertex in enumerate(basis):
        assert pairing(vertex, vertex) == 0
        vanishing_mirrors = {
            root_index
            for root_index, root in enumerate(WEHLER_ROOTS)
            if pairing(vertex, root) == 0
        }
        assert vanishing_mirrors == {0, 1, 2} - {index}


def test_wehler_lattice_has_no_minus_two_vector():
    # The exact symbolic identity makes an unbounded integer search unnecessary.
    a, b, c = sp.symbols("a b c", integer=True)
    vector = sp.Matrix([a, b, c])
    assert sp.expand(pairing(vector, vector)) == 4 * (a * b + a * c + b * c)
    assert -2 % 4 != 0

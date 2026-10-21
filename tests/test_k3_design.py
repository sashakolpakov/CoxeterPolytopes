from fractions import Fraction

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.geometric.k3_design import (
    audit_k3_coxeter_design,
    correspondence_gram,
    paired_trace_curve_gram,
    recover_field_gram,
    trace_curve_gram,
)
from coxeter_polytopes.known_examples.bugaenko import (
    BUGAENKO_VERTICES,
    bugaenko_wall_gram,
)


def test_three_and_five_edges_have_the_same_curve_intersection() -> None:
    field_gram = [
        [2, -1, -PHI],
        [-1, 2, 0],
        [-PHI, 0, 2],
    ]
    assert trace_curve_gram(field_gram) == (
        (Fraction(-2), Fraction(1), Fraction(1)),
        (Fraction(1), Fraction(-2), Fraction(0)),
        (Fraction(1), Fraction(0), Fraction(-2)),
    )
    assert correspondence_gram(field_gram) == (
        (Fraction(-2), Fraction(1), Fraction(2)),
        (Fraction(1), Fraction(-2), Fraction(0)),
        (Fraction(2), Fraction(0), Fraction(-2)),
    )
    assert recover_field_gram(
        trace_curve_gram(field_gram), correspondence_gram(field_gram)
    ) == tuple(tuple(Golden.coerce(entry) for entry in row) for row in field_gram)


def test_full_paired_curve_lift_distinguishes_three_and_five_edges() -> None:
    field_gram = [
        [2, -1, -PHI],
        [-1, 2, 0],
        [-PHI, 0, 2],
    ]
    paired = paired_trace_curve_gram(field_gram)

    # Each restricted wall lifts to the orthogonal pair r_i, (phi-1)r_i.
    assert all(paired[i][i] == -2 for i in range(6))
    assert all(paired[i][i + 3] == 0 for i in range(3))

    # A 3-edge lifts to A2 + A2: r_0 meets r_1 and u*r_0 meets u*r_1.
    assert paired[0][1] == paired[3][4] == 1
    assert paired[0][4] == paired[3][1] == 0

    # A 5-edge lifts to A4: three of the four cross-pairings are one.
    assert paired[0][2] == paired[0][5] == paired[3][2] == 1
    assert paired[3][5] == 0


def test_bugaenko_passes_inverse_k3_design_constraints() -> None:
    audit = audit_k3_coxeter_design(bugaenko_wall_gram(), BUGAENKO_VERTICES)

    assert audit.passes is True
    assert audit.wall_count == 11
    assert audit.field_rank == 9
    assert audit.finite_labels == (2, 3, 5)
    assert audit.trace_off_diagonal_values == (
        Fraction(0),
        Fraction(1),
        Fraction(2),
    )
    assert audit.paired_effective_curve_pairings is True
    assert audit.paired_trace_off_diagonal_values == (
        Fraction(0),
        Fraction(1),
        Fraction(2),
    )
    assert audit.identity_inertia == Inertia(8, 1, 2)
    assert audit.conjugate_inertia == Inertia(9, 0, 2)
    assert recover_field_gram(
        trace_curve_gram(bugaenko_wall_gram()),
        correspondence_gram(bugaenko_wall_gram()),
    ) == tuple(tuple(row) for row in bugaenko_wall_gram())


def test_non_effective_trace_orientation_is_rejected() -> None:
    # This is a valid ultraparallel field pair at the identity embedding, but
    # its trace intersection is -1, impossible for distinct effective curves.
    gram = [[Golden(2), Golden(6, -5)], [Golden(6, -5), Golden(2)]]
    audit = audit_k3_coxeter_design(gram, [])
    assert audit.effective_curve_pairings is False
    assert audit.passes is False

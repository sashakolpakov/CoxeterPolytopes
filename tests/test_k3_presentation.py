from sympy import Matrix as SympyMatrix

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.geometric.k3_design import (
    correspondence_gram,
    trace_curve_gram,
)
from coxeter_polytopes.geometric.k3_presentation import (
    audit_all_vertex_anchors,
    audit_k3_root_presentation,
    build_k3_root_presentation,
    doubled_trace_gram,
    golden_action_matrix,
)
from coxeter_polytopes.known_examples.bugaenko import (
    BUGAENKO_VERTICES,
    bugaenko_wall_gram,
)


def _buga_two_grams() -> tuple[tuple[tuple[object, ...], ...], ...]:
    field_gram = bugaenko_wall_gram()
    return trace_curve_gram(field_gram), correspondence_gram(field_gram)


def test_bugaenko_two_gram_presents_the_saturated_trace_lattice() -> None:
    q, b = _buga_two_grams()
    presentation = build_k3_root_presentation(q, b)

    assert presentation.wall_count == 11
    assert presentation.field_rank == 9
    assert presentation.trace_rank == 18
    assert presentation.radical_rank == 4
    assert presentation.trace_inertia == Inertia(1, 17, 0)
    assert presentation.doubled_inertia == Inertia(1, 17, 4)
    assert presentation.is_even is True
    assert presentation.discriminant == -4
    assert presentation.discriminant_group_invariants == (2, 2)

    doubled = SympyMatrix(doubled_trace_gram(q, b))
    change = SympyMatrix(presentation.generator_change)
    quotient_lifts = SympyMatrix(presentation.quotient_lifts)
    radical = SympyMatrix(presentation.saturated_radical_basis)
    quotient_gram = SympyMatrix(presentation.quotient_gram)
    quotient_action = SympyMatrix(presentation.quotient_action)

    # A unimodular completion proves that these four integral relations are
    # saturated; a denominator-cleared rational nullspace would not do so.
    assert abs(change.det()) == 1
    assert doubled * radical == SympyMatrix.zeros(22, 4)
    assert quotient_lifts.T * doubled * quotient_lifts == quotient_gram
    assert quotient_action**2 == quotient_action + SympyMatrix.eye(18)
    assert quotient_action.T * quotient_gram == quotient_gram * quotient_action


def test_bugaenko_presentation_passes_both_place_signatures() -> None:
    q, b = _buga_two_grams()
    audit = audit_k3_root_presentation(q, b, hyperbolic_dimension=8)

    assert audit.passes is True
    assert audit.identity_inertia == Inertia(8, 1, 2)
    assert audit.conjugate_inertia == Inertia(9, 0, 2)
    assert audit.expected_trace_inertia == Inertia(1, 17, 0)
    assert audit.expected_doubled_inertia == Inertia(1, 17, 4)


def test_every_bugaenko_vertex_anchor_has_zero_schur_complement() -> None:
    q, b = _buga_two_grams()
    audits = audit_all_vertex_anchors(q, b, BUGAENKO_VERTICES)

    assert len(audits) == 41 * 3
    assert all(audit.passes for audit in audits)
    assert all(len(audit.remainder) == 2 for audit in audits)
    assert all(
        all(value.sign(0) > 0 for _, value in audit.outer_anchor_coordinates)
        for audit in audits
    )


def test_one_wrong_correspondence_entry_breaks_the_buga_rank_gate() -> None:
    q, b_immutable = _buga_two_grams()
    b = [list(row) for row in b_immutable]
    b[0][1] += 1
    b[1][0] += 1

    audit = audit_k3_root_presentation(q, b, hyperbolic_dimension=8)
    assert audit.passes is False


def test_standard_golden_action_is_self_adjoint_for_the_doubled_gram() -> None:
    q, b = _buga_two_grams()
    doubled = SympyMatrix(doubled_trace_gram(q, b))
    action = SympyMatrix(golden_action_matrix(11))

    assert action**2 == action + SympyMatrix.eye(22)
    assert action.T * doubled == doubled * action

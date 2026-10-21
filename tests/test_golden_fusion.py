from coxeter_polytopes.exact.golden import Golden, PHI, ZERO
from coxeter_polytopes.exact.matrix import Inertia, determinant
from coxeter_polytopes.geometric.golden_fusion import (
    GoldenPairDomain,
    audit_lanner_overlaps,
    fusion_program_from_field_gram,
    golden_lanner_placement_options,
)
from coxeter_polytopes.geometric.golden_lanner import (
    golden_block_product_counterexample,
)
from coxeter_polytopes.geometric.k3_design import (
    correspondence_gram,
    trace_curve_gram,
)
from coxeter_polytopes.geometric.k3_presentation import audit_vertex_anchor
from coxeter_polytopes.known_examples.bugaenko import (
    BUGAENKO_VERTICES,
    bugaenko_wall_gram,
)


def test_labelled_lanner_placements_form_finite_exact_tables() -> None:
    order_four = golden_lanner_placement_options(range(4))
    order_five = golden_lanner_placement_options(range(5))

    assert len(order_four) == 54
    assert len(order_five) == 240
    assert all(len(option.requirements) == 6 for option in order_four)
    assert all(len(option.requirements) == 10 for option in order_five)
    assert {option.block_name for option in order_four} == {
        "A4(3-5-3)",
        "A4(5-3-5)",
        "D4(5,3,3)",
        "C4(one 5)",
        "C4(opposite 5s)",
    }
    assert {option.block_name for option in order_five} == {
        "A5(5-3-3-3)",
        "A5(5-3-3-5)",
        "D5(terminal 5)",
    }


def test_bugaenko_passes_every_compiled_vertex_fusion_constraint() -> None:
    program = fusion_program_from_field_gram(
        bugaenko_wall_gram(), BUGAENKO_VERTICES
    )
    result = program.propagate()

    assert len(program.constraints()) == 41
    assert result.compiled_vertex_count == 41
    assert result.ready_vertex_count == 41
    assert result.pruned_value_count == 0
    assert result.contradiction is None
    assert result.unresolved_vertices == ()
    assert result.passes is True

    first_system = program.symbolic_systems()[0]
    assert len(first_system.complement) == 3
    assert len(first_system.equations) == 3
    assert first_system.rational_equation_count == 6
    assert len(first_system.signs) == 9
    assert first_system.evaluate({}).passes is True


def test_unsquared_fusion_sign_is_the_vertex_anchor_same_side_test() -> None:
    field_gram = bugaenko_wall_gram()
    vertex = BUGAENKO_VERTICES[0]
    program = fusion_program_from_field_gram(field_gram, BUGAENKO_VERTICES)
    constraint = next(
        constraint
        for constraint in program.constraints()
        if constraint.vertex == vertex
    )
    propagated = constraint.propagate(program.domain_map)
    square_classes = constraint.determinant_square_class_audit(program.domain_map)
    diagonal = dict(propagated.diagonal_numerators)
    residual = {
        pair: values[0] for pair, values in propagated.residual_numerators
    }
    q = trace_curve_gram(field_gram)
    b = correspondence_gram(field_gram)

    assert propagated.contradiction is None
    assert square_classes.passes is True
    assert len(square_classes.relations) == 3
    assert all(
        relation.square_root * relation.square_root
        == relation.determinant_ratio
        for relation in square_classes.relations
        if relation.square_root is not None
    )
    for outside, determinant_value in square_classes.resolved_determinants:
        support = (*vertex, outside)
        ten_wall_gram = [
            [field_gram[left][right] for right in support]
            for left in support
        ]
        assert determinant(ten_wall_gram) == determinant_value
    assert all(value.sign(0) < 0 for value in diagonal.values())
    assert all(value.sign(0) < 0 for value in residual.values())
    for anchor in constraint.complement:
        anchor_audit = audit_vertex_anchor(q, b, vertex, anchor)
        assert anchor_audit.same_vertex_side is True
        for wall, coordinate in anchor_audit.outer_anchor_coordinates:
            pair = tuple(sorted((anchor, wall)))
            assert coordinate == residual[pair] / diagonal[anchor]
            assert coordinate.sign(0) > 0


def test_vertex_fusion_prunes_a_false_bugaenko_dotted_product() -> None:
    program = fusion_program_from_field_gram(
        bugaenko_wall_gram(),
        BUGAENKO_VERTICES,
        pair_domain_overrides={(9, 10): (-2 * PHI, -(PHI * PHI))},
    )
    symbolic = program.symbolic_systems()[0]

    assert symbolic.evaluate({(9, 10): -2 * PHI}).passes is True
    assert symbolic.evaluate({(9, 10): -(PHI * PHI)}).passes is False

    result = program.propagate()
    assert result.passes is True
    assert result.pruned_value_count == 1
    assert result.domain_map[(9, 10)].values == (-2 * PHI,)


def test_countercertificate_is_rejected_by_vertex_fusion() -> None:
    counterexample = golden_block_product_counterexample()
    program = fusion_program_from_field_gram(
        counterexample.field_gram, counterexample.vertices
    )
    result = program.propagate()

    assert set(program.missing_faces) == set(counterexample.missing_faces)
    assert result.contradiction is not None
    assert "nonsquare ratio" in result.contradiction
    assert "10/11 + 3/11*phi" in result.contradiction
    assert result.passes is False


def test_square_class_gate_can_reject_before_all_rows_resolve() -> None:
    counterexample = golden_block_product_counterexample()
    program = fusion_program_from_field_gram(
        counterexample.field_gram, counterexample.vertices
    )
    constraint = program.constraints()[0]
    domains = program.domain_map
    unresolved_pair = (0, 13)
    old = domains[unresolved_pair]
    domains[unresolved_pair] = GoldenPairDomain(
        unresolved_pair, old.kind, (ZERO, Golden(-1))
    )

    audit = constraint.determinant_square_class_audit(domains)
    propagated = constraint.propagate(domains)

    assert audit.unresolved_facets == (13,)
    assert audit.passes is False
    assert audit.contradiction is not None
    assert "nonsquare ratio" in audit.contradiction
    assert propagated.ready is False
    assert propagated.contradiction == audit.contradiction


def test_two_resolved_lanner_seeds_with_two_time_directions_are_rejected() -> None:
    counterexample = golden_block_product_counterexample()
    placements = []
    for support in (range(4), range(4, 8)):
        matches = [
            option
            for option in golden_lanner_placement_options(support)
            if all(
                counterexample.field_gram[left][right] == value
                for (left, right), value in option.requirements
            )
        ]
        assert len(matches) == 1
        placements.append(matches[0])

    program = fusion_program_from_field_gram(
        counterexample.field_gram,
        counterexample.vertices,
        lanner_placements=placements,
    )
    audits = audit_lanner_overlaps(placements, program.domain_map)

    assert len(audits) == 1
    assert audits[0].resolved is True
    assert audits[0].identity_inertia == Inertia(6, 2, 0)
    assert audits[0].conjugate_inertia == Inertia(8, 0, 0)
    assert audits[0].passes is False

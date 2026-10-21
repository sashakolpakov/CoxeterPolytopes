from dataclasses import replace
from fractions import Fraction

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.geometric.golden_halfspaces import (
    audit_golden_halfspace_certificate,
    decide_lorentzian_closed_cone,
    decide_lorentzian_strict_face,
    find_nonnegative_span_phase_one,
    find_positive_span,
    golden_dot,
    golden_quadratic_value,
    inverse_golden_matrix,
    positive_proportional_classes,
    positive_proportional_ratio,
    positive_span_on_support,
    project_to_halfspace_wall,
)


ROWS = {
    "x": (1, 0, 0),
    "phi*x": (PHI, 0, 0),
    "y": (0, 1, 0),
    "redundant": (1 + PHI, PHI - 1, 0),
}
GRAM = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
INTERIOR = (1, 1, 1)


def test_exact_positive_proportional_partition_and_ratio() -> None:
    classes = positive_proportional_classes(ROWS)

    assert tuple(item.members for item in classes) == (
        ("x", "phi*x"),
        ("y",),
        ("redundant",),
    )
    assert classes[0].ratios_to_first == (Golden(1), PHI)
    assert positive_proportional_ratio((1, 0), (-1, 0)) is None


def test_exact_projection_matrix_inverse_and_quadratic_value() -> None:
    inverse = inverse_golden_matrix(((1, PHI), (PHI, 1 + PHI * PHI)))
    assert inverse[0][0] == 1 + PHI * PHI
    assert inverse[0][1] == -PHI
    assert inverse[1][0] == -PHI
    assert inverse[1][1] == 1

    face = project_to_halfspace_wall(INTERIOR, ROWS["x"], GRAM)
    assert face == (Golden(0), Golden(1), Golden(1))
    assert golden_dot(ROWS["x"], face) == 0
    assert golden_quadratic_value(GRAM, face) == 2


def test_positive_span_generator_finds_exact_golden_dependence() -> None:
    terms = find_positive_span(
        ROWS["redundant"],
        {"x": ROWS["x"], "y": ROWS["y"]},
    )

    assert terms == ((1 + PHI, "x"), (PHI - 1, "y"))
    assert positive_span_on_support(
        ROWS["redundant"],
        {"x": ROWS["x"], "y": ROWS["y"]},
    ) == terms


def test_phase_one_cone_membership_is_exact_and_complete() -> None:
    terms, feasible = find_nonnegative_span_phase_one(
        ROWS["redundant"],
        {"x": ROWS["x"], "y": ROWS["y"]},
    )
    impossible, infeasible = find_nonnegative_span_phase_one(
        (-1, -1, 0),
        {"x": ROWS["x"], "y": ROWS["y"]},
    )

    assert terms == ((1 + PHI, "x"), (PHI - 1, "y"))
    assert feasible.feasible and feasible.objective == 0 and feasible.passes
    assert impossible is None
    assert not infeasible.feasible
    assert infeasible.objective.sign() > 0
    assert infeasible.passes


def test_lorentzian_active_set_distinguishes_causal_and_spacelike_faces() -> None:
    gram = ((1, 0, 0), (0, -1, 0), (0, 0, -1))
    interior = (2, 0, 1)
    target = (0, 0, 1)

    # On target=0 the remaining inequality says x>2t, so the entire strict
    # future linear face is spacelike.  The exact normalized maximum is < 0.
    inactive = decide_lorentzian_strict_face(
        target,
        {"forcing": (-2, 1, 5)},
        quadratic_gram=gram,
        interior=interior,
    )
    assert not inactive.is_facet
    assert inactive.reason == (
        "normalized_face_has_nonpositive_global_quadratic_maximum"
    )
    assert inactive.maximum_quadratic_value == Golden(Fraction(-27, 4))
    assert inactive.active_labels == ("forcing",)
    assert inactive.active_multipliers == (Golden(6),)

    # Here phase I happens to return a spacelike basic point, but the concave
    # Lorentzian QP has positive maximum and mixing with Slater yields a strict
    # timelike Q(phi)-point.  This pins the nontrivial fallback rather than only
    # its easy projection path.
    causal = decide_lorentzian_strict_face(
        target,
        {"forcing": (Fraction(-1, 2), 1, 2)},
        quadratic_gram=gram,
        interior=interior,
    )
    assert causal.is_facet
    assert causal.reason == "positive_normalized_quadratic_maximum"
    assert causal.maximum_quadratic_value == Golden(Fraction(27, 16))
    assert causal.witness is not None
    assert golden_dot(target, causal.witness) == 0
    assert golden_dot((Fraction(-1, 2), 1, 2), causal.witness).sign() > 0
    assert golden_quadratic_value(gram, causal.witness) == Golden(Fraction(63, 64))


def test_strict_future_face_phase_one_can_prove_infeasibility() -> None:
    decision = decide_lorentzian_strict_face(
        (1, 0),
        {},
        quadratic_gram=((1, 0), (0, -1)),
        interior=(1, 0),
    )

    assert not decision.is_facet
    assert decision.reason == "strict_future_linear_face_infeasible"
    assert not decision.strict_linear_certificate.feasible
    assert decision.strict_linear_certificate.objective.sign() > 0
    assert decision.strict_linear_certificate.passes


def test_closed_lorentzian_cone_accepts_an_open_timelike_section() -> None:
    decision = decide_lorentzian_closed_cone(
        {"future": (1, 0), "right": (0, 1)},
        quadratic_gram=((1, 0), (0, -1)),
        future_reference=(1, 0),
    )

    assert decision.has_future_timelike_point
    assert decision.witness is not None
    assert decision.closed_linear_certificate.passes
    assert golden_quadratic_value(
        ((1, 0), (0, -1)), decision.witness
    ).sign() > 0


def test_closed_lorentzian_cone_keeps_a_lower_dimensional_carrier_cell() -> None:
    decision = decide_lorentzian_closed_cone(
        {"right": (0, 1), "left": (0, -1)},
        quadratic_gram=((1, 0), (0, -1)),
        future_reference=(1, 0),
    )

    assert decision.has_future_timelike_point
    assert decision.witness == (Golden(1), Golden(0))
    assert golden_dot((0, 1), decision.witness) == 0
    assert golden_dot((0, -1), decision.witness) == 0


def test_closed_lorentzian_cone_rejects_a_spacelike_closed_slice_by_kkt() -> None:
    decision = decide_lorentzian_closed_cone(
        {
            "lower": (-2, 1),
            "lower-copy": (-4, 2),
            "upper": (2, -1),
        },
        quadratic_gram=((1, 0), (0, -1)),
        future_reference=(1, 0),
    )

    assert not decision.has_future_timelike_point
    assert decision.reason == (
        "normalized_closed_slice_has_nonpositive_global_quadratic_maximum"
    )
    assert decision.maximum_point == (Golden(1), Golden(2))
    assert decision.maximum_quadratic_value == Golden(-3)
    assert decision.equality_multiplier is not None
    assert decision.retained_linear_labels
    assert any(label == "lower-copy" for label, _ in decision.linear_redundancies)


def test_closed_lorentzian_cone_detects_an_infeasible_time_slice() -> None:
    decision = decide_lorentzian_closed_cone(
        {"past-only": (-1, 0)},
        quadratic_gram=((1, 0), (0, -1)),
        future_reference=(1, 0),
    )

    assert not decision.has_future_timelike_point
    assert decision.reason == "normalized_closed_slice_infeasible"
    assert not decision.closed_linear_certificate.feasible
    assert decision.closed_linear_certificate.objective.sign() > 0
    assert decision.closed_linear_certificate.passes


def test_halfspace_audit_accepts_exact_causal_inactivity_separately() -> None:
    gram = ((1, 0, 0), (0, -1, 0), (0, 0, -1))
    interior = (2, 0, 1)
    rows = {"inactive": (0, 0, 1), "facet": (-2, 1, 5)}
    inactivity = decide_lorentzian_strict_face(
        rows["inactive"],
        {"facet": rows["facet"]},
        quadratic_gram=gram,
        interior=interior,
    )
    facet = decide_lorentzian_strict_face(
        rows["facet"],
        {"inactive": rows["inactive"]},
        quadratic_gram=gram,
        interior=interior,
    )
    assert facet.witness is not None

    audit = audit_golden_halfspace_certificate(
        rows,
        facet_classes=(("facet",),),
        facet_witnesses={"facet": facet.witness},
        redundancies={},
        causal_inactivities={"inactive": inactivity},
        interior=interior,
        quadratic_gram=gram,
    )

    assert audit.passes
    assert audit.redundant_class_count == 1
    assert len(audit.causal_inactivities) == 1
    assert audit.causal_inactivities[0].passes

    tampered = audit_golden_halfspace_certificate(
        rows,
        facet_classes=(("facet",),),
        facet_witnesses={"facet": facet.witness},
        redundancies={},
        causal_inactivities={
            "inactive": replace(
                inactivity,
                maximum_quadratic_value=Golden(-1),
            )
        },
        interior=interior,
        quadratic_gram=gram,
    )
    assert not tampered.passes
    assert not tampered.causal_inactivities[0].passes


def test_complete_halfspace_certificate_passes_exactly() -> None:
    witnesses = {
        "x": project_to_halfspace_wall(INTERIOR, ROWS["x"], GRAM),
        "y": project_to_halfspace_wall(INTERIOR, ROWS["y"], GRAM),
    }
    dependence = find_positive_span(
        ROWS["redundant"],
        {"x": ROWS["x"], "y": ROWS["y"]},
    )
    assert dependence is not None

    audit = audit_golden_halfspace_certificate(
        ROWS,
        facet_classes=(("x", "phi*x"), ("y",)),
        facet_witnesses=witnesses,
        redundancies={"redundant": dependence},
        interior=INTERIOR,
        quadratic_gram=GRAM,
    )

    assert audit.passes
    assert audit.primitive_wall_count == 4
    assert audit.facet_class_count == 2
    assert audit.redundant_class_count == 1
    assert all(item.passes for item in audit.facet_witnesses)
    assert all(item.passes for item in audit.positive_dependences)


def test_certificate_rejects_circular_or_incomplete_redundancy_proofs() -> None:
    witnesses = {
        "x": (0, 1, 1),
        "y": (1, 0, 1),
    }
    circular = audit_golden_halfspace_certificate(
        ROWS,
        facet_classes=(("x", "phi*x"), ("y",)),
        facet_witnesses=witnesses,
        redundancies={"redundant": ((Golden(1), "redundant"),)},
        interior=INTERIOR,
        quadratic_gram=GRAM,
    )
    missing = audit_golden_halfspace_certificate(
        ROWS,
        facet_classes=(("x", "phi*x"), ("y",)),
        facet_witnesses=witnesses,
        redundancies={},
        interior=INTERIOR,
        quadratic_gram=GRAM,
    )

    assert not circular.passes
    assert not circular.positive_dependences[0].sources_are_facets
    assert not missing.passes
    assert not missing.redundancy_partition_exact


def test_certificate_requires_the_geometric_timelike_check() -> None:
    audit = audit_golden_halfspace_certificate(
        ROWS,
        facet_classes=(("x", "phi*x"), ("y",)),
        facet_witnesses={"x": (0, 1, 1), "y": (1, 0, 1)},
        redundancies={
            "redundant": ((1 + PHI, "x"), (PHI - 1, "y")),
        },
        interior=INTERIOR,
    )

    assert not audit.passes
    assert not audit.interior_timelike
    assert all(not witness.timelike for witness in audit.facet_witnesses)


def test_certificate_rejects_a_timelike_witness_on_the_opposite_sheet() -> None:
    audit = audit_golden_halfspace_certificate(
        {"x": (0, 1)},
        facet_classes=(("x",),),
        facet_witnesses={"x": (-1, 0)},
        redundancies={},
        interior=(2, 1),
        quadratic_gram=((1, 0), (0, -1)),
    )

    assert not audit.passes
    assert audit.facet_witnesses[0].timelike
    assert not audit.facet_witnesses[0].same_time_component

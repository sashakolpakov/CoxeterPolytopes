from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.geometric.golden_lanner import (
    audit_golden_lanner,
    audit_vertex_complement_schur,
    golden_block_product_counterexample,
    golden_lanner_blocks,
)
from coxeter_polytopes.known_examples.bugaenko import (
    BUGAENKO_VERTICES,
    bugaenko_wall_gram,
)


def test_exact_golden_lanner_census_has_eight_blocks() -> None:
    blocks = golden_lanner_blocks()
    expected = {
        "A4(3-5-3)": (Golden(5, -4), -11, (11,), "A4(-1)"),
        "A4(5-3-5)": (Golden(6, -5), -19, (19,), "A4(-1)"),
        "D4(5,3,3)": (Golden(4, -4), -16, (2, 2, 2, 2), "D4(-1)"),
        "C4(one 5)": (Golden(2, -5), -31, (31,), "affine A3(-1)"),
        "C4(opposite 5s)": (
            Golden(1, -7),
            -55,
            (55,),
            "affine A3(-1)",
        ),
        "A5(5-3-3-3)": (Golden(6, -4), -4, (2, 2), "A5(-1)"),
        "A5(5-3-3-5)": (Golden(8, -6), -20, (2, 10), "A5(-1)"),
        "D5(terminal 5)": (Golden(4, -4), -16, (4, 4), "D5(-1)"),
    }

    assert len(blocks) == 8
    assert {block.name for block in blocks} == set(expected)
    for block in blocks:
        assert block.name is not None
        assert block.is_arithmetic_lanner is True
        assert (
            block.field_determinant,
            block.restriction_trace_determinant,
            block.discriminant_group_invariants,
            block.curve_type,
        ) == expected[block.name]


def test_only_order_four_cycles_have_affine_curve_trace() -> None:
    blocks = golden_lanner_blocks()
    affine = [block for block in blocks if block.curve_inertia.zero]
    assert [block.name for block in affine] == [
        "C4(one 5)",
        "C4(opposite 5s)",
    ]
    assert all(block.curve_inertia == Inertia(0, 3, 1) for block in affine)
    assert all(
        block.curve_inertia == Inertia(0, block.order, 0)
        for block in blocks
        if not block.curve_inertia.zero
    )


def test_a_tampered_lanner_block_fails_the_exact_audit() -> None:
    gram = [
        [2, -1, 0, 0],
        [-1, 2, -1, 0],
        [0, -1, 2, -1],
        [0, 0, -1, 2],
    ]
    audit = audit_golden_lanner(gram)
    assert audit.name is None
    assert audit.is_arithmetic_lanner is False


def test_bugaenko_vertex_has_the_forced_rank_one_schur_complement() -> None:
    audit = audit_vertex_complement_schur(
        bugaenko_wall_gram(), BUGAENKO_VERTICES[0]
    )
    assert audit.passes is True
    assert len(audit.vertex) == 8
    assert len(audit.complement) == 3
    assert audit.schur_rank == 1
    assert audit.schur_identity_inertia == Inertia(0, 1, 2)
    assert audit.schur_conjugate_inertia == Inertia(1, 0, 2)
    assert audit.all_two_by_two_minors_zero is True


def test_block_product_passes_every_local_gate_but_fails_global_fusion() -> None:
    counterexample = golden_block_product_counterexample()
    assert tuple(map(len, counterexample.missing_faces)) == (4, 4, 2, 2, 2)
    assert counterexample.passes_local_constraints is True
    assert counterexample.global_identity_inertia == Inertia(9, 5, 0)
    assert counterexample.global_conjugate_inertia == Inertia(14, 0, 0)
    assert counterexample.global_trace_inertia == Inertia(5, 23, 0)
    assert counterexample.every_vertex_schur_rank_five is True
    assert counterexample.sample_schur_audit.schur_rank == 5
    assert counterexample.sample_schur_audit.passes is False
    assert counterexample.passes_global_h9_gram is False

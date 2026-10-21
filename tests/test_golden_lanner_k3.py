from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.geometric.golden_lanner import golden_lanner_blocks
from coxeter_polytopes.geometric.golden_lanner_k3 import (
    golden_lanner_k3_seed_audits,
    h11_matching_seed_completion,
    h19_golden_seed_completions,
)
from coxeter_polytopes.geometric.k3_design import correspondence_gram


def test_all_eight_lanner_trace_seeds_embed_primitively_in_k3() -> None:
    audits = golden_lanner_k3_seed_audits()

    assert len(audits) == 8
    assert all(audit.primitive_k3_embedding_exists for audit in audits)
    assert all(audit.primitive_k3_embedding_unique for audit in audits)
    assert {
        audit.name: (
            audit.k3_complement_signature,
            audit.k3_complement_rank,
            audit.discriminant_length,
            audit.nikulin_length_bound,
        )
        for audit in audits
    } == {
        "A4(3-5-3)": (Inertia(2, 12, 0), 14, 1, 12),
        "A4(5-3-5)": (Inertia(2, 12, 0), 14, 1, 12),
        "D4(5,3,3)": (Inertia(2, 12, 0), 14, 4, 12),
        "C4(one 5)": (Inertia(2, 12, 0), 14, 1, 12),
        "C4(opposite 5s)": (Inertia(2, 12, 0), 14, 1, 12),
        "A5(5-3-3-3)": (Inertia(2, 10, 0), 12, 2, 10),
        "A5(5-3-3-5)": (Inertia(2, 10, 0), 12, 2, 10),
        "D5(terminal 5)": (Inertia(2, 10, 0), 12, 2, 10),
    }


def test_no_lanner_seed_has_a_proper_tau_stable_even_saturation() -> None:
    audits = golden_lanner_k3_seed_audits()

    expected = {
        "A4(3-5-3)": [(1, True, (11,))],
        "A4(5-3-5)": [(1, True, (19,))],
        "D4(5,3,3)": [
            (1, True, (2, 2, 2, 2)),
            (2, False, (2, 2)),
            (2, False, (2, 2)),
            (2, False, (2, 2)),
        ],
        "C4(one 5)": [(1, True, (31,))],
        "C4(opposite 5s)": [(1, True, (55,))],
        "A5(5-3-3-3)": [(1, True, (2, 2)), (2, False, ())],
        "A5(5-3-3-5)": [(1, True, (2, 10)), (2, False, (5,))],
        "D5(terminal 5)": [
            (1, True, (4, 4)),
            (2, False, (2, 2)),
            (4, False, ()),
            (4, False, ()),
        ],
    }

    for audit in audits:
        assert [
            (
                overlattice.index,
                overlattice.tau_stable,
                overlattice.discriminant_group_invariants,
            )
            for overlattice in audit.even_overlattices
        ] == expected[audit.name]
        assert audit.nontrivial_tau_stable_overlattices == ()
        assert audit.primitive_in_every_tau_stable_completion is True


def test_discriminant_actions_are_computed_in_exact_smith_coordinates() -> None:
    assert {
        audit.name: (audit.discriminant_moduli, audit.tau_basis_images)
        for audit in golden_lanner_k3_seed_audits()
    } == {
        "A4(3-5-3)": ((11,), ((4,),)),
        "A4(5-3-5)": ((19,), ((5,),)),
        "D4(5,3,3)": (
            (2, 2, 2, 2),
            ((0, 1, 0, 0), (1, 1, 0, 0), (0, 0, 1, 1), (0, 0, 1, 0)),
        ),
        "C4(one 5)": ((31,), ((19,),)),
        "C4(opposite 5s)": ((55,), ((8,),)),
        "A5(5-3-3-3)": ((2, 2), ((0, 1), (1, 1))),
        "A5(5-3-3-5)": ((2, 10), ((1, 5), (1, 8))),
        "D5(terminal 5)": ((4, 4), ((2, 1), (3, 3))),
    }


def test_h19_contains_primitive_tau_stable_completions_of_all_eight_seeds() -> None:
    completions = h19_golden_seed_completions()

    assert len(completions) == 8
    assert all(completion.passes for completion in completions)
    assert {
        completion.seed_name: (
            completion.complement_rank,
            completion.complement_determinant,
            completion.complement_discriminant_group_invariants,
            completion.direct_sum_index,
        )
        for completion in completions
    } == {
        "A4(3-5-3)": (12, 209, (209,), 11),
        "A4(5-3-5)": (12, 361, (19, 19), 19),
        "D4(5,3,3)": (12, 304, (2, 2, 2, 38), 16),
        "C4(one 5)": (12, 589, (589,), 31),
        "C4(opposite 5s)": (12, 1045, (1045,), 55),
        "A5(5-3-3-3)": (10, 76, (2, 38), 4),
        "A5(5-3-3-5)": (10, 380, (2, 190), 20),
        "D5(terminal 5)": (10, 304, (4, 76), 16),
    }
    assert all(
        completion.direct_sum_index == abs(completion.seed_trace_determinant)
        for completion in completions
    )
    assert all(
        len(completion.gluing_subgroup) == completion.direct_sum_index
        for completion in completions
    )


def test_matching_prime_seeds_have_explicit_h11_and_h19_gluings() -> None:
    h11 = h11_matching_seed_completion()
    h19 = next(
        completion
        for completion in h19_golden_seed_completions()
        if completion.seed_name == "A4(5-3-5)"
    )

    assert (
        h11.seed_trace_determinant,
        h11.complement_determinant,
        h11.complement_discriminant_group_invariants,
        h11.direct_sum_index,
        h11.ambient_determinant,
        h11.ambient_tau_basis_images,
    ) == (-11, 121, (11, 11), 11, -11, ((8,),))
    assert h11.seed_tau_basis_images == ((4,),)
    assert h11.complement_tau_basis_images == ((8, 0), (1, 4))
    assert (
        h19.seed_trace_determinant,
        h19.complement_determinant,
        h19.complement_discriminant_group_invariants,
        h19.direct_sum_index,
        h19.ambient_determinant,
        h19.ambient_tau_basis_images,
    ) == (-19, 361, (19, 19), 19, -19, ((15,),))
    assert h19.seed_tau_basis_images == ((5,),)
    assert h19.complement_tau_basis_images == ((15, 0), (18, 5))
    assert h11.passes is True
    assert h19.passes is True


def test_affine_cycle_fibers_are_isotropic_but_not_tau_eigenvectors() -> None:
    affine = {
        block.name: block
        for block in golden_lanner_blocks()
        if block.name in {"C4(one 5)", "C4(opposite 5s)"}
    }

    for name, expected in (("C4(one 5)", 2), ("C4(opposite 5s)", 4)):
        block = affine[name]
        curve_square = sum(sum(row) for row in block.curve_gram)
        correspondence = correspondence_gram(block.field_gram)
        fiber_tau_intersection = sum(sum(row) for row in correspondence)

        assert curve_square == 0
        assert fiber_tau_intersection == expected
        # Self-adjointness and tau^2=tau+I give
        # q(tau*F,tau*F)=q(F,tau*F)+q(F,F).
        assert fiber_tau_intersection + curve_square == expected

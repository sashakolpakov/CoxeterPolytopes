from fractions import Fraction

from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_leech import (
    EPSILON,
    MOG_NORM10_SPAN_WITNESS,
    P7_COEFFICIENT,
    P7_EXCEPTIONAL_MOG_RAW,
    P7_GENERATOR,
    P7_LEECH_SEED,
    SQRT2_LEECH_SEED,
    TWO_SQRT2_LEECH_SEED,
    UNIT_LEECH_SEED,
    candidate_leech_gram_triples,
    integer_determinant,
    nonunit_leech_second_packet_audit,
    p7_leech_first_packet_audit,
    p7_leech_post_first_prefix_h82,
    p7_leech_minimal_direction_prefix_h82,
    p7_leech_second_packet_audit,
    reflective_candidate_leech_gram_triples,
    sqrt2_leech_first_packet_audit,
    two_sqrt2_leech_first_packet_audit,
    unit_leech_first_packet_audit,
    unit_leech_minimal_direction_prefix_h82,
    unit_leech_second_packet_audit,
)


def test_unit_norm_classes_first_packet_and_hyperideal_open_patch() -> None:
    assert UNIT_LEECH_SEED.admissible
    assert UNIT_LEECH_SEED.canonical_mirror_norms == (Silver(1), Silver(2))
    audit = unit_leech_first_packet_audit()
    assert audit.earlier_control_count == 10
    assert audit.first_height == 14 + 10 * SQRT2
    assert audit.first_norm == 2
    assert audit.first_control == 2 + 2 * SQRT2
    assert audit.first_target == 30 + 20 * SQRT2
    assert audit.shell_norm == 10
    assert audit.wall_count == 4_629_381_120
    assert audit.same_height_pairing_bound == -4 - 4 * SQRT2
    assert audit.span_witness_determinant == 2
    assert integer_determinant(MOG_NORM10_SPAN_WITNESS) == 2
    assert audit.affine_section_bounded
    assert audit.strict_witness_spatial_square == Silver(64) / 25
    assert audit.strict_witness_spatial_square > EPSILON
    assert audit.hyperideal_open_patch
    assert not audit.finite_volume and not audit.compact


def test_unit_second_packet_and_h82_persistent_minimal_direction() -> None:
    second = unit_leech_second_packet_audit()
    assert second.height == 28 + 20 * SQRT2
    assert second.norm == 2
    assert second.control == 4 + 2 * SQRT2
    assert second.target == 58 + 40 * SQRT2
    assert second.decomposition == "y=w+epsilon*v, w^2=4, v^2=14, w.v=6"
    assert second.old_boundary_witness_is_cut
    assert second.replacement_witness_coordinate == Silver(4, 10) / 23
    assert second.replacement_witness_spatial_square - EPSILON == Silver(335, -209) / 529
    assert second.replacement_first_packet_margin.sign(0) > 0
    assert second.replacement_hyperideal

    bounds = unit_leech_minimal_direction_prefix_h82()
    assert [(row.height, row.maximum_pairing) for row in bounds] == [
        (28 + 20 * SQRT2, 10 + 6 * SQRT2),
        (41 + 29 * SQRT2, 11 + 8 * SQRT2),
        (56 + 40 * SQRT2, 14 + 10 * SQRT2),
        (63 + 45 * SQRT2, 15 + 11 * SQRT2),
        (82 + 58 * SQRT2, 17 + 12 * SQRT2),
    ]
    assert all(row.wall_margin.sign(0) >= 0 for row in bounds)


def test_nonunit_first_shells_and_pi3_reflectivity_gate() -> None:
    sqrt2 = sqrt2_leech_first_packet_audit()
    assert SQRT2_LEECH_SEED.admissible
    assert sqrt2.first_height == 8 + 6 * SQRT2
    assert sqrt2.first_norm == 2
    assert sqrt2.first_control == 2 + 2 * SQRT2
    assert sqrt2.scalar_shell_norm == 6
    assert sqrt2.wall_count == 16_773_120
    assert sqrt2.strict_hyperideal_witness

    two_sqrt2 = two_sqrt2_leech_first_packet_audit()
    assert TWO_SQRT2_LEECH_SEED.admissible
    assert two_sqrt2.first_height == 8 + 6 * SQRT2
    assert two_sqrt2.first_norm == 2
    assert two_sqrt2.first_control == 2 + SQRT2
    assert two_sqrt2.scalar_shell_norm == 6
    assert two_sqrt2.strict_hyperideal_witness

    # Without the local y in pi*Lambda condition this is the spurious
    # epsilon times minimal shell.  Dividing its target by pi^2 kills it.
    target = 12 + 8 * SQRT2
    assert candidate_leech_gram_triples(target)
    assert not reflective_candidate_leech_gram_triples(
        TWO_SQRT2_LEECH_SEED, 4 + 2 * SQRT2, target
    )

    for seed, expected_square_difference in (
        (SQRT2_LEECH_SEED, Silver(3971, 0) / 4900 - Silver(0, 1367) / 2450),
        (TWO_SQRT2_LEECH_SEED, Silver(3971, 0) / 2450 - Silver(0, 1367) / 1225),
    ):
        second = nonunit_leech_second_packet_audit(seed)
        assert second.complete_height == 24 + 17 * SQRT2
        assert second.surviving_triple.u_square == 26
        assert second.first_packet_margin.sign(0) > 0
        assert second.second_packet_margin.sign(0) > 0
        assert second.witness_spatial_square - second.timelike_threshold == expected_square_difference
        assert second.hyperideal_prefix_witness


def test_p7_first_packet_is_exact_contact_polar_and_hyperideal() -> None:
    assert P7_COEFFICIENT.norm() == -7
    assert P7_GENERATOR.norm() == 7
    assert P7_LEECH_SEED.admissible
    assert P7_LEECH_SEED.canonical_mirror_norms == (
        Silver(1),
        Silver(2),
        3 + SQRT2,
        6 + 2 * SQRT2,
    )
    audit = p7_leech_first_packet_audit()
    assert audit.earlier_control_count == 7
    assert audit.earlier_surviving_triple_count == 0
    assert audit.first_height == 5 + 4 * SQRT2
    assert audit.first_norm == 2
    assert audit.first_control == 2 + SQRT2
    assert audit.first_target == 12 + 8 * SQRT2
    assert audit.scalar_shell_norm == 4
    assert audit.wall_count == 196_560
    assert audit.normalized_slab_bound == 4 - SQRT2
    assert audit.contact_polar_vertex_orbit_count == 232
    assert audit.exceptional_mog_raw == P7_EXCEPTIONAL_MOG_RAW
    assert len(audit.exceptional_mog_raw) == 24
    assert audit.exceptional_mog_raw_square == 192
    assert audit.exceptional_mog_divisor == 3
    assert audit.exceptional_mog_coordinate_scale_square == 8
    assert audit.exceptional_contact_polar_square == Fraction(8, 3)
    assert audit.scaled_exceptional_square == 12 - 16 * SQRT2 / 3
    assert audit.exceptional_lorentzian_norm == 13 - 22 * SQRT2 / 3
    assert audit.exceptional_lorentzian_norm.sign(0) > 0
    assert audit.hyperideal_open_patch
    assert not audit.finite_volume and not audit.compact


def test_p7_exceptional_ray_cut_and_corrected_h82_replacement() -> None:
    audit = p7_leech_post_first_prefix_h82()
    assert audit.complete_height == 82 + 58 * SQRT2
    assert audit.control_level_count == 27
    assert audit.empty_necessary_control_count == 23
    assert audit.earliest_later_height == 31 + 22 * SQRT2
    assert audit.earliest_later_norm == 2
    assert audit.earliest_later_control == 4 + 3 * SQRT2
    assert audit.earliest_later_target == 64 + 44 * SQRT2
    assert audit.earliest_later_triples == (
        # The first type has the explicit archived cutter below.
        type(audit.cutter_u_v_triple)(36, 14, 22),
        type(audit.cutter_u_v_triple)(32, 16, 22),
    )
    assert audit.cutter_direction_square == 6
    assert audit.cutter_direction_pairing == 14 + 8 * SQRT2
    assert audit.cutter_first_packet_pairing_bound == 10 + 6 * SQRT2
    assert audit.cutter_first_packet_margin == 2 * SQRT2
    assert audit.exceptional_ray_coordinate == (4 - SQRT2) / 3
    assert audit.exceptional_ray_cut_excess == Silver(16, 3) / 3
    assert audit.replacement_coordinate == Fraction(37, 67)
    assert audit.replacement_first_packet_margin == Silver(157) / 67 - SQRT2
    assert audit.replacement_spatial_square == Fraction(8214, 4489)
    assert audit.replacement_lorentzian_norm == Silver(12703) / 4489 - 2 * SQRT2
    assert audit.replacement_lorentzian_norm.sign(0) > 0

    rows = audit.nonempty_necessary_levels
    assert [(row.height, row.maximum_norm_six_pairing) for row in rows] == [
        (31 + 22 * SQRT2, 14 + 8 * SQRT2),
        (56 + 40 * SQRT2, 21 + 15 * SQRT2),
        (62 + 44 * SQRT2, 17 + 13 * SQRT2),
        (82 + 58 * SQRT2, 21 + 15 * SQRT2),
    ]
    # The p-divisible rows must use x, not a*x, after y/rho normalization.
    assert rows[1].spatial_normalizer.startswith("rho=")
    assert rows[1].normalized_wall_bound == 12 + 8 * SQRT2
    assert rows[1].replacement_margin == Silver(27, -19) / 67
    assert rows[1].replacement_margin.sign(0) > 0
    assert rows[3].spatial_normalizer.startswith("rho=")
    assert all(row.replacement_margin.sign(0) > 0 for row in rows)
    assert audit.strict_hyperideal_prefix_witness
    assert not audit.finite_volume_through_prefix


def test_p7_second_packet_and_h82_norm_six_direction_certificate() -> None:
    second = p7_leech_second_packet_audit()
    assert second.complete_through_height == 31 + 22 * SQRT2
    assert second.norm == 2
    assert second.control == 4 + 3 * SQRT2
    assert second.target == 64 + 44 * SQRT2
    assert second.surviving_triples == (
        type(second.explicit_realized_triple)(36, 14, 22),
        type(second.explicit_realized_triple)(32, 16, 22),
    )
    assert second.explicit_realized_triple == second.surviving_triples[0]
    assert second.old_witness_is_cut
    assert second.replacement_direction_square == 6
    assert second.replacement_witness_coordinate == Silver(37) / 67
    assert second.first_packet_margin == Silver(157) / 67 - SQRT2
    assert second.first_packet_margin.sign(0) > 0
    assert second.replacement_spatial_square == Silver(8214) / 4489
    assert second.replacement_spatial_square - P7_COEFFICIENT == Silver(12703) / 4489 - 2 * SQRT2
    assert second.replacement_hyperideal

    bounds = p7_leech_minimal_direction_prefix_h82()
    assert [(row.height, row.norm, row.maximum_pairing) for row in bounds] == [
        (31 + 22 * SQRT2, Silver(2), 14 + 8 * SQRT2),
        (56 + 40 * SQRT2, 6 + 2 * SQRT2, 21 + 15 * SQRT2),
        (62 + 44 * SQRT2, Silver(2), 17 + 13 * SQRT2),
        (82 + 58 * SQRT2, 6 + 2 * SQRT2, 21 + 15 * SQRT2),
    ]
    assert all(row.wall_margin.sign(0) > 0 for row in bounds)

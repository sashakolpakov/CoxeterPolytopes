from coxeter_polytopes.exact.golden import PHI
from coxeter_polytopes.search.golden_a1_golay import (
    a1_golay_pair,
    audit_uniform_inequality_level,
    explicit_first_uniform_cutter,
    first_uniform_ray_cutter,
    replacement_uniform_ray_audit,
    octad_incidence_rank,
    uniform_ray_audit,
    x1_packet,
    xphi_packet,
)


def test_first_two_packets_and_uniform_hyperideal_extreme_ray() -> None:
    assert len(x1_packet()) == 24
    assert len(xphi_packet()) == 276 + 759 == 1035
    assert all(a1_golay_pair(root.vector, root.vector) == 2 for root in (*x1_packet(), *xphi_packet()))
    assert octad_incidence_rank() == 24
    audit = uniform_ray_audit()
    assert audit.norm == 3 + PHI
    assert audit.all_octad_pairings_zero
    assert audit.x1_pairings_nonpositive and audit.xphi_pairings_nonpositive
    assert audit.hyperideal and audit.extreme
    assert first_uniform_ray_cutter(tuple(root.vector for root in (*x1_packet(), *xphi_packet()))) is None


def test_uniform_ray_first_exact_cutter_is_an_acute_golay_hexadecad_family() -> None:
    before = audit_uniform_inequality_level(PHI * PHI)
    assert before.violating_pattern_count == 0
    first = audit_uniform_inequality_level(1 + 2 * PHI)
    assert first.golay_admissible_pattern_count == 33
    assert first.violating_pattern_count == 16
    assert first.maximum_excess == 4 + PHI

    cutter = explicit_first_uniform_cutter()
    assert cutter.control == 1 + 2 * PHI
    assert cutter.height == 8 + 13 * PHI
    assert cutter.norm == 2
    assert cutter.spatial_sum == 8 + 9 * PHI
    assert cutter.uniform_pairing == (1 + 5 * PHI) / 2
    assert cutter.golay_member and cutter.accepted_through_earlier_levels
    assert cutter.labeled_family_count == 759 * 8 == 6072
    assert cutter.octad_intersection_sizes == (0, 2, 4)
    assert cutter.labeled_family_mutually_acute
    assert all(margin.sign(0) < 0 for margin in cutter.earlier_relaxed_pairing_margins)
    assert all(margin.sign(0) < 0 for _case, margin in cutter.same_height_pairing_margins)

    assert len(first.maximum_patterns) == 1
    replacement = replacement_uniform_ray_audit()
    assert replacement.spatial_coordinate == (7 + 6 * PHI) / 55
    assert replacement.norm == (816 - 58 * PHI) / 605
    assert replacement.active_orbit_count == 6072
    assert replacement.active_spatial_rank == 24
    assert replacement.all_active_pairings_zero
    assert replacement.all_earlier_shells_feasible
    assert replacement.full_active_shell_feasible
    assert replacement.hyperideal and replacement.extreme
    assert replacement.later_root_slope == -6 + 7 * PHI

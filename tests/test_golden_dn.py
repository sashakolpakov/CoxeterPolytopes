from coxeter_polytopes.exact.golden import ONE, PHI, ZERO, Golden
from coxeter_polytopes.search.golden_dn import (
    _in_d24_plus_spinor_coset,
    analyze_golden_dn_cone,
    golden_dn_first_cap_geometry,
    golden_dn_pair,
    run_golden_dn_search,
)


def test_d24_plus_accepts_all_three_nonzero_spinor_residues():
    half = Golden(1) / 2
    for residue in (ONE, PHI, ONE + PHI):
        vector = tuple(residue * half for _ in range(24))
        assert _in_d24_plus_spinor_coset(vector)
    assert not _in_d24_plus_spinor_coset(tuple(half if index else ZERO for index in range(24)))


def test_golden_d24_first_two_packets_match_exact_regressions():
    search = run_golden_dn_search(24, (ONE, PHI), d24_plus=True)
    assert [(level.candidate_count, len(level.roots)) for level in search.levels] == [
        (1, 1),
        (3, 2),
    ]
    assert len(search.walls) == 27
    supports = {
        tuple(value for value in root.vector[1:] if value)
        for root in search.levels[1].roots
    }
    assert supports == {(2 * PHI,), (PHI, PHI, PHI, PHI)}
    assert all(
        golden_dn_pair(left, right).sign(0) <= 0
        for index, left in enumerate(search.walls)
        for right in search.walls[index + 1 :]
    )
    geometry = golden_dn_first_cap_geometry(24)
    assert (
        geometry.extreme_ray_count,
        geometry.finite_vertex_count,
        geometry.cusp_count,
        geometry.hyperideal_ray_count,
    ) == (25, 3, 0, 22)


def test_golden_dn_h22_h29_stabilize_through_phi_squared():
    controls = (ONE, PHI, PHI * PHI)
    low = run_golden_dn_search(22, controls)
    high = run_golden_dn_search(29, controls)
    assert [(level.candidate_count, len(level.roots)) for level in low.levels] == [
        (1, 1),
        (3, 2),
        (4, 1),
    ]
    assert [(level.candidate_count, len(level.roots)) for level in high.levels] == [
        (1, 1),
        (3, 2),
        (4, 1),
    ]
    padded = {
        (
            root.control,
            (*root.vector[1:], *(Golden(0) for _ in range(29 - 22))),
        )
        for root in low.accepted_roots
    }
    assert padded == {
        (root.control, root.vector[1:]) for root in high.accepted_roots
    }
    low_geometry = analyze_golden_dn_cone(low)
    high_geometry = analyze_golden_dn_cone(high)
    assert (
        len(low_geometry.extreme_rays),
        low_geometry.finite_vertex_count,
        low_geometry.cusp_count,
        low_geometry.hyperideal_ray_count,
    ) == (254, 28, 0, 226)
    assert (
        len(high_geometry.extreme_rays),
        high_geometry.finite_vertex_count,
        high_geometry.cusp_count,
        high_geometry.hyperideal_ray_count,
    ) == (352, 28, 0, 324)
    assert low_geometry.irredundant and high_geometry.irredundant

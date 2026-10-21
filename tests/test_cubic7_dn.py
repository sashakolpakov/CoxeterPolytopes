from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_dn import (
    FIRST_HEIGHT,
    d24_plus_basis_audit,
    d24_plus_integral_basis_coordinates,
    d24_plus_membership,
    first_level_dominant_patterns,
    run_dn_first_packet,
)


def test_d24_plus_first_packet_and_complete_norm_audit() -> None:
    patterns = first_level_dominant_patterns(24, include_spinor_coset=True)
    assert len(patterns) == 1
    assert patterns[0].coset == "D"
    assert tuple(value for value in patterns[0].coordinates if value) == (
        Cubic7(1, 1), Cubic7(1, 1)
    )
    result = run_dn_first_packet(24, include_spinor_coset=True)
    assert result.canonical_mirror_norms == (Cubic7(1), Cubic7(2))
    assert result.complete_norm_class_proof is not None
    invalid_mixed = (Cubic7(1, 1) / 2,) * 8 + (Cubic7(0),) * 16
    assert d24_plus_membership(invalid_mixed) is None
    assert d24_plus_membership(patterns[0].coordinates) == "D"
    assert d24_plus_basis_audit()
    assert all(
        value.is_integral
        for value in d24_plus_integral_basis_coordinates(patterns[0].coordinates)
    )
    assert result.wall_rank == 25
    assert len(result.walls) == 25
    assert len({ray.vector for ray in result.rays}) == len(result.rays)
    assert all(ray.active_rank_verified for ray in result.rays)
    assert len(result.rays) == 25
    assert (result.finite_vertex_count, result.cusp_count, result.hyperideal_ray_count) == (2, 0, 23)
    assert not result.finite_volume


def test_dn_dimension_sweep_is_scoped_first_shell() -> None:
    for dimension in (22, 23, 25, 26, 27, 28, 29):
        result = run_dn_first_packet(dimension)
        assert result.complete_norm_class_proof is None
        assert len(result.patterns) == 1
        assert len(result.walls) == dimension + 1
        assert len(result.rays) == dimension + 1
        assert all(ray.active_rank_verified for ray in result.rays)
        assert result.cusp_count == 0
        assert result.hyperideal_ray_count == dimension - 1

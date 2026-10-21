import pytest

from coxeter_polytopes.known_examples.bugaenko import bugaenko_wall_gram
from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.verify.vinberg_volume import (
    _connected_affine_subdiagrams,
    verify_vinberg_volume_criterion,
)


def test_bugaenko_compact_h8_passes_exact_endpoint_criterion() -> None:
    result = verify_vinberg_volume_criterion(
        bugaenko_wall_gram(),
        8,
        simple_roots_complete=True,
    )
    assert result.decision == "compact"
    assert result.compact is True
    assert result.finite_volume is True
    assert result.elliptic_edge_count == 164
    assert result.finite_vertex_count == 41
    assert result.ideal_vertex_count == 0
    assert result.parabolic_neighborhood_audits == ()
    assert result.compact_bad_edges == ()
    assert result.finite_volume_bad_edges == ()


def test_noncompact_finite_area_triangle_is_calibration() -> None:
    # The (3,3,infinity) Coxeter triangle has two finite vertices and the
    # affine-A1 endpoint {0,2}.  Every one of its three edges has two ends.
    gram = [[2, -1, -2], [-1, 2, -1], [-2, -1, 2]]
    result = verify_vinberg_volume_criterion(
        gram,
        2,
        simple_roots_complete=True,
    )
    assert result.decision == "finite_volume_noncompact"
    assert result.compact is False
    assert result.finite_volume is True
    assert result.elliptic_edge_count == 3
    assert result.finite_vertices == ((0, 1), (1, 2))
    assert tuple(endpoint.facets for endpoint in result.ideal_vertices) == ((0, 2),)
    assert result.finite_volume_bad_edges == ()


def test_full_zero_neighborhood_rejects_affine_cusp_with_finite_contaminant() -> None:
    # J={0,1} is affine A1, but node 2 is orthogonal to all of J.  Therefore
    # N(J)={0,1,2}=affine-A1 + finite-A1 is not a parabolic cusp link.
    gram = [
        [2, -2, 0, -1],
        [-2, 2, 0, -1],
        [0, 0, 2, -1],
        [-1, -1, -1, 2],
    ]
    result = verify_vinberg_volume_criterion(
        gram,
        3,
        simple_roots_complete=True,
    )
    assert result.decision == "infinite_volume"
    assert result.finite_volume is False
    assert result.ideal_vertices == ()
    assert len(result.parabolic_neighborhood_audits) == 1
    audit = result.parabolic_neighborhood_audits[0]
    assert audit.affine_component == (0, 1)
    assert audit.zero_neighborhood == (0, 1, 2)
    assert audit.is_rank_dimension_minus_one_parabolic is False


def test_incomplete_vinberg_prefix_is_explicitly_undecidable() -> None:
    gram = [[2, -1, -2], [-1, 2, -1], [-2, -1, 2]]
    result = verify_vinberg_volume_criterion(
        gram,
        2,
        simple_roots_complete=False,
    )
    assert result.decision == "undecidable_incomplete_simple_roots"
    assert result.compact is None
    assert result.finite_volume is None
    assert result.edge_audits == ()
    assert result.finite_vertices == ()
    assert result.ideal_vertices == ()


def test_affine_audit_does_not_assume_the_desired_cusp_rank_bound() -> None:
    # Affine A2 has order three, larger than hyperbolic dimension two.  It is
    # invalid as a rank-one cusp and must be found before it can be rejected.
    affine_a2 = [
        [Golden(2), Golden(-1), Golden(-1)],
        [Golden(-1), Golden(2), Golden(-1)],
        [Golden(-1), Golden(-1), Golden(2)],
    ]
    assert _connected_affine_subdiagrams(affine_a2, 2) == ((0, 1, 2),)


def test_decomposable_global_diagram_is_not_a_hyperbolic_chamber() -> None:
    with pytest.raises(ValueError, match="indecomposable"):
        verify_vinberg_volume_criterion(
            [[2, 0], [0, 2]],
            2,
            simple_roots_complete=True,
        )

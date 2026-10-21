from __future__ import annotations

from coxeter_polytopes.search.silver_k3_vinberg import (
    run_silver_k3_vinberg_extension,
    run_silver_k3_vinberg_prefix,
)
from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.verify.silver_vinberg_volume import (
    audit_silver_vinberg_wall_gram,
)
from coxeter_polytopes.verify.vinberg_volume import verify_vinberg_volume_criterion


def test_sparse_audit_reproduces_twelve_wall_prefix() -> None:
    prefix = run_silver_k3_vinberg_prefix()
    audit = audit_silver_vinberg_wall_gram(prefix.wall_gram, 9)
    assert len(audit.elliptic_edges) == prefix.diagram.elliptic_rank_eight_count
    assert len(audit.finite_vertices) == prefix.diagram.elliptic_rank_nine_count
    assert audit.finite_endpoint_histogram == prefix.diagram.finite_endpoint_histogram
    assert audit.connected_affine_subdiagrams == ()
    assert audit.ideal_vertices == ()
    assert audit.parabolic_neighborhood_audits == ()
    assert audit.ideal_endpoint_histogram == ((0, len(audit.elliptic_edges)),)
    assert (
        audit.finite_plus_ideal_endpoint_histogram
        == audit.finite_endpoint_histogram
    )
    assert audit.invalid_finite_angle_pairs == ()
    assert audit.compact_stopping_available is False
    assert audit.finite_volume_stopping_available is False


def test_sparse_audit_reproduces_sixteen_wall_extension() -> None:
    extension = run_silver_k3_vinberg_extension()
    audit = audit_silver_vinberg_wall_gram(extension.wall_gram, 9)
    assert len(audit.elliptic_edges) == extension.diagram.elliptic_rank_eight_count
    assert len(audit.finite_vertices) == extension.diagram.elliptic_rank_nine_count
    assert audit.finite_endpoint_histogram == extension.diagram.finite_endpoint_histogram
    assert audit.connected_affine_subdiagrams == ()
    assert audit.ideal_vertices == ()
    assert audit.parabolic_neighborhood_audits == ()
    assert audit.ideal_endpoint_histogram == ((0, len(audit.elliptic_edges)),)
    assert (
        audit.finite_plus_ideal_endpoint_histogram
        == audit.finite_endpoint_histogram
    )
    assert audit.invalid_finite_angle_pairs == ()
    assert audit.compact_stopping_available is False
    assert audit.finite_volume_stopping_available is False


def test_sparse_audit_detects_connected_affine_a2() -> None:
    affine_a2 = tuple(
        tuple(Silver(2 if left == right else -1) for right in range(3))
        for left in range(3)
    )
    audit = audit_silver_vinberg_wall_gram(affine_a2, 2)
    assert audit.connected_affine_subdiagrams == ((0, 1, 2),)
    assert audit.ideal_vertices == ()
    assert len(audit.bad_parabolic_neighborhoods) == 1
    assert audit.compact_stopping_available is False
    assert audit.finite_volume_stopping_available is False


def test_sparse_finite_volume_matches_generic_ideal_triangle_semantics() -> None:
    gram = (
        (Silver(2), Silver(-1), Silver(-2)),
        (Silver(-1), Silver(2), Silver(-1)),
        (Silver(-2), Silver(-1), Silver(2)),
    )
    sparse = audit_silver_vinberg_wall_gram(gram, 2)
    generic = verify_vinberg_volume_criterion(
        [[2, -1, -2], [-1, 2, -1], [-2, -1, 2]],
        2,
        simple_roots_complete=True,
    )

    assert sparse.connected_affine_subdiagrams == ((0, 2),)
    assert tuple(endpoint.facets for endpoint in sparse.ideal_vertices) == ((0, 2),)
    assert sparse.parabolic_neighborhood_audits[0].zero_neighborhood == (0, 2)
    assert (
        sparse.parabolic_neighborhood_audits[0].is_rank_dimension_minus_one_parabolic
        is True
    )
    assert sparse.finite_endpoint_histogram == ((1, 2), (2, 1))
    assert sparse.ideal_endpoint_histogram == ((0, 1), (1, 2))
    assert sparse.finite_plus_ideal_endpoint_histogram == ((2, 3),)
    assert tuple(
        (edge.facets, edge.finite_endpoints, edge.ideal_endpoints)
        for edge in sparse.edge_audits
    ) == tuple(
        (edge.facets, edge.finite_endpoints, edge.ideal_endpoints)
        for edge in generic.edge_audits
    )
    assert sparse.compact_stopping_available is False
    assert sparse.finite_volume_stopping_available is True


def test_sparse_full_zero_neighborhood_rejects_finite_contaminant() -> None:
    gram = (
        (Silver(2), Silver(-2), Silver(0), Silver(-1)),
        (Silver(-2), Silver(2), Silver(0), Silver(-1)),
        (Silver(0), Silver(0), Silver(2), Silver(-1)),
        (Silver(-1), Silver(-1), Silver(-1), Silver(2)),
    )
    sparse = audit_silver_vinberg_wall_gram(gram, 3)
    assert sparse.connected_affine_subdiagrams == ((0, 1),)
    assert sparse.ideal_vertices == ()
    assert len(sparse.parabolic_neighborhood_audits) == 1
    neighborhood = sparse.parabolic_neighborhood_audits[0]
    assert neighborhood.zero_neighborhood == (0, 1, 2)
    assert neighborhood.neighborhood_components == ((0, 1), (2,))
    assert neighborhood.is_rank_dimension_minus_one_parabolic is False
    assert sparse.bad_parabolic_neighborhoods == (neighborhood,)
    assert sparse.finite_volume_stopping_available is False


def test_sparse_reducible_cusp_is_one_endpoint_not_one_per_component() -> None:
    gram = (
        (Silver(2), Silver(-2), Silver(0), Silver(0)),
        (Silver(-2), Silver(2), Silver(0), Silver(0)),
        (Silver(0), Silver(0), Silver(2), Silver(-2)),
        (Silver(0), Silver(0), Silver(-2), Silver(2)),
    )
    sparse = audit_silver_vinberg_wall_gram(gram, 3)
    assert sparse.connected_affine_subdiagrams == ((0, 1), (2, 3))
    assert len(sparse.parabolic_neighborhood_audits) == 2
    assert all(
        row.is_rank_dimension_minus_one_parabolic
        for row in sparse.parabolic_neighborhood_audits
    )
    assert len(sparse.ideal_vertices) == 1
    endpoint = sparse.ideal_vertices[0]
    assert endpoint.facets == (0, 1, 2, 3)
    assert endpoint.affine_components == ((0, 1), (2, 3))
    assert sparse.finite_endpoint_histogram == ((0, 4),)
    assert sparse.ideal_endpoint_histogram == ((1, 4),)
    assert sparse.finite_plus_ideal_endpoint_histogram == ((1, 4),)
    assert all(edge.ideal_endpoints == ((0, 1, 2, 3),) for edge in sparse.edge_audits)
    assert sparse.bad_parabolic_neighborhoods == ()
    assert sparse.finite_volume_stopping_available is False

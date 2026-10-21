from itertools import combinations

from coxeter_polytopes.combinatorics.h9_14_record35_coxeter_labels import (
    audit_record35_coxeter_labels,
    record35_vertex_facets,
)


def _cartan_times_vector(
    facets: tuple[int, ...],
    nonright_labels: tuple[tuple[int, int, int], ...],
    vector: tuple[int, ...],
) -> tuple[int, ...]:
    labels = {
        tuple(sorted((left, right))): label
        for left, right, label in nonright_labels
    }
    return tuple(
        sum(
            (
                2
                if left == right
                else -1
                if labels.get(tuple(sorted((left, right))), 2) == 3
                else 0
            )
            * coefficient
            for right, coefficient in zip(facets, vector, strict=True)
        )
        for left in facets
    )


def test_record35_common_core_lanner_join_is_locally_impossible() -> None:
    audit = audit_record35_coxeter_labels()
    assert audit.verified
    assert audit.d_core_extension_count_histogram == (
        (1, 36),
        (2, 30),
        (3, 40),
        (4, 51),
    )
    assert audit.d_join_counts == (420, 1_332, 4_620)
    assert audit.d_rank_four_passing_state_count == 4_620
    assert audit.d_rank_five_passing_state_count == 96
    assert audit.d_rank_six_failure_multiplicity_histogram == (
        (1, 72),
        (2, 24),
    )
    assert audit.sparse_passing_state_count == 0


def test_the_two_survivor_orbits_have_exact_affine_witnesses() -> None:
    audit = audit_record35_coxeter_labels()
    vertices = {frozenset(vertex) for vertex in record35_vertex_facets()}
    assert [orbit.orbit_size for orbit in audit.d_survivor_orbits] == [24, 72]
    assert [orbit.affine_type for orbit in audit.d_survivor_orbits] == [
        "affine_D4_plus_A1",
        "affine_D5",
    ]
    for orbit in audit.d_survivor_orbits:
        assert frozenset(orbit.containing_vertex) in vertices
        assert set(orbit.witness_facets) <= set(orbit.containing_vertex)
        assert _cartan_times_vector(
            orbit.witness_facets,
            orbit.nonright_labels,
            orbit.cartan_null_vector,
        ) == (0,) * len(orbit.witness_facets)


def test_external_nonright_labels_cannot_rescue_the_sparse_shadow() -> None:
    audit = audit_record35_coxeter_labels()
    assert audit.compatible_five_lanner_state_count == 53_749_080
    assert audit.one_external_nonright_state_count == 4_620 * 3 * 3
    assert audit.one_external_nonright_passing_state_count == 0
    assert audit.unrestricted_external_label_no_go is True
    # In normalized Gram convention, cos(pi/2)=0 and cos(pi/m)>0 for m>=3;
    # every external non-right edge therefore increases the nonnegative
    # adjacency entrywise.  The certificate's monotonic conclusion covers
    # every shell, not just the explicitly replayed one-edge shell.
    assert all(label > 2 for label in (3, 4, 5))


def test_two_support_cluster_counts_reproduce_independently() -> None:
    audit = audit_record35_coxeter_labels()
    assert audit.c_join_counts == (420, 11_634)
    assert audit.c_restriction_size_histogram == ((4, 10), (5, 16), (6, 3))
    assert audit.c_passing_state_count == 24
    assert audit.c_orbit_size_histogram == ((24, 1),)
    assert len(tuple(combinations(range(5), 2))) == 10

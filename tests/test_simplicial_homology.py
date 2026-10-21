from itertools import combinations

from coxeter_polytopes.combinatorics.simplicial_homology import (
    first_mod2_homology_sphere_link_failure,
    gf2_rank,
    link_face_masks,
    reduced_betti_numbers_mod2,
)


def _face_table(vertex_count: int, facets: tuple[int, ...]) -> tuple[bool, ...]:
    table = [False] * (1 << vertex_count)
    for facet in facets:
        subset = facet
        while True:
            table[subset] = True
            if subset == 0:
                break
            subset = (subset - 1) & facet
    return tuple(table)


def test_gf2_rank_and_circle_homology() -> None:
    assert gf2_rank((0b011, 0b101, 0b110)) == 2
    cycle_facets = tuple(
        (1 << first) | (1 << ((first + 1) % 5)) for first in range(5)
    )
    table = _face_table(5, cycle_facets)
    assert reduced_betti_numbers_mod2(link_face_masks(table, 0)) == (0, 1)
    assert first_mod2_homology_sphere_link_failure(table, 1) is None


def test_simplex_boundary_and_disconnected_union() -> None:
    # Boundary of a tetrahedron.
    sphere_facets = tuple(
        sum(1 << vertex for vertex in face)
        for face in combinations(range(4), 3)
    )
    sphere_table = _face_table(4, sphere_facets)
    assert first_mod2_homology_sphere_link_failure(sphere_table, 2) is None

    # Two disjoint triangles have Euler characteristic two, just like S^2,
    # but mod-two reduced H_0 detects the disconnected complex immediately.
    disconnected = _face_table(
        6,
        (
            sum(1 << vertex for vertex in (0, 1, 2)),
            sum(1 << vertex for vertex in (3, 4, 5)),
        ),
    )
    failure = first_mod2_homology_sphere_link_failure(disconnected, 2)
    assert failure is not None
    assert failure.face_mask == 0
    assert failure.observed_reduced_betti[0] == 1


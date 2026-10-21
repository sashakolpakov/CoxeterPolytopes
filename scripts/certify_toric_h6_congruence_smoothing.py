#!/usr/bin/env python3
"""Independent exact audit of the toric H6 congruence refinement.

This deliberately does not use cddlib or the production toric certificate.
All matrix calculations are over the integers in ``I_(1,6)``.  The only
external theorem used in the group identification is the theorem of
Everitt--Ratcliffe--Tschantz that the 27 exceptional-curve reflections
generate the principal level-two reflection group with chamber ``P^6``.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import factorial
from typing import Iterable


J = (1, -1, -1, -1, -1, -1, -1)
NAMES = (
    "E1-E4",
    "E4",
    "H-E1-E2-E4",
    "E2-E5",
    "E5",
    "H-E2-E3-E5",
    "E3-E6",
    "E6",
    "H-E3-E1-E6",
)
WALLS = (
    (0, 1, 0, 0, -1, 0, 0),
    (0, 0, 0, 0, 1, 0, 0),
    (1, -1, -1, 0, -1, 0, 0),
    (0, 0, 1, 0, 0, -1, 0),
    (0, 0, 0, 0, 0, 1, 0),
    (1, 0, -1, -1, 0, -1, 0),
    (0, 0, 0, 1, 0, 0, -1),
    (0, 0, 0, 0, 0, 0, 1),
    (1, -1, 0, -1, 0, 0, -1),
)
IDENTITY = tuple(
    1 if row == column else 0 for row in range(7) for column in range(7)
)


def pair(left: Iterable[int], right: Iterable[int]) -> int:
    return sum(sign * x * y for sign, x, y in zip(J, left, right, strict=True))


def matrix_product(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(left[7 * row + k] * right[7 * k + column] for k in range(7))
        for row in range(7)
        for column in range(7)
    )


def matrix_vector(matrix: tuple[int, ...], vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(matrix[7 * row + column] * vector[column] for column in range(7))
        for row in range(7)
    )


def reflection(root: tuple[int, ...]) -> tuple[int, ...]:
    norm = pair(root, root)
    assert norm in (-1, -2)
    multiplier = -2 // norm
    covector = tuple(root[index] * J[index] for index in range(7))
    return tuple(
        (1 if row == column else 0) + multiplier * root[row] * covector[column]
        for row in range(7)
        for column in range(7)
    )


def finite_a2_cubed_group() -> set[tuple[int, ...]]:
    generators = [reflection(root) for root in WALLS if pair(root, root) == -2]
    group = {IDENTITY}
    queue = [IDENTITY]
    while queue:
        current = queue.pop()
        for generator in generators:
            product = matrix_product(current, generator)
            if product not in group:
                group.add(product)
                queue.append(product)
    return group


def exceptional_classes() -> tuple[tuple[int, ...], ...]:
    result: list[tuple[int, ...]] = []
    for index in range(6):
        vector = [0] * 7
        vector[index + 1] = 1
        result.append(tuple(vector))
    for first, second in combinations(range(6), 2):
        vector = [1, 0, 0, 0, 0, 0, 0]
        vector[first + 1] = vector[second + 1] = -1
        result.append(tuple(vector))
    for missing in range(6):
        vector = [2, -1, -1, -1, -1, -1, -1]
        vector[missing + 1] = 0
        result.append(tuple(vector))
    assert len(result) == len(set(result)) == 27
    return tuple(result)


def _components(subset: tuple[int, ...]) -> list[list[int]]:
    unseen = set(subset)
    result: list[list[int]] = []
    while unseen:
        seed = min(unseen)
        unseen.remove(seed)
        queue = [seed]
        component: list[int] = []
        while queue:
            current = queue.pop()
            component.append(current)
            for neighbour in ((current - 1) % 9, (current + 1) % 9):
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    queue.append(neighbour)
        result.append(component)
    return result


EDGE_LABELS = (4, 4, 3, 4, 4, 3, 4, 4, 3)


def _path_edge_labels(component: list[int]) -> tuple[int, ...] | None:
    if len(component) == 1:
        return ()
    vertices = set(component)
    endpoints = [
        vertex
        for vertex in component
        if sum(
            neighbour in vertices
            for neighbour in ((vertex - 1) % 9, (vertex + 1) % 9)
        )
        == 1
    ]
    if not endpoints:
        return None
    previous: int | None = None
    current = endpoints[0]
    labels: list[int] = []
    while True:
        following = [
            neighbour
            for neighbour in ((current - 1) % 9, (current + 1) % 9)
            if neighbour in vertices and neighbour != previous
        ]
        if not following:
            break
        neighbour = following[0]
        labels.append(
            EDGE_LABELS[current]
            if neighbour == (current + 1) % 9
            else EDGE_LABELS[neighbour]
        )
        previous, current = current, neighbour
    return tuple(labels)


def _finite_path_order(labels: tuple[int, ...]) -> int | None:
    rank = len(labels) + 1
    if all(label == 3 for label in labels):
        return factorial(rank + 1)  # A_rank
    if (
        labels.count(4) == 1
        and 4 in (labels[0], labels[-1])
        and all(label in (3, 4) for label in labels)
    ):
        return 2**rank * factorial(rank)  # B_rank=C_rank
    if labels == (3, 4, 3):
        return 1152  # F4
    return None


def spherical_parabolic_order(subset: tuple[int, ...]) -> int | None:
    order = 1
    for component in _components(subset):
        labels = _path_edge_labels(component)
        if labels is None:
            return None
        component_order = _finite_path_order(labels)
        if component_order is None:
            return None
        order *= component_order
    return order


def toric_solomon_characteristic() -> tuple[Fraction, int]:
    characteristic = Fraction(0)
    spherical_count = 0
    for mask in range(1 << 9):
        subset = tuple(index for index in range(9) if mask & (1 << index))
        order = spherical_parabolic_order(subset)
        if order is not None:
            spherical_count += 1
            characteristic += Fraction((-1) ** len(subset), order)
    return characteristic, spherical_count


def p6_clique_vector_and_characteristic() -> tuple[list[int], Fraction]:
    roots = exceptional_classes()
    counts: list[int] = []
    for rank in range(7):
        counts.append(
            sum(
                all(pair(roots[i], roots[j]) == 0 for i, j in combinations(subset, 2))
                for subset in combinations(range(27), rank)
            )
        )
    characteristic = sum(
        (Fraction((-1) ** rank * count, 2**rank) for rank, count in enumerate(counts)),
        Fraction(0),
    )
    return counts, characteristic


def boundary_subset_representations() -> dict[tuple[int, ...], tuple[int, ...]]:
    """Write every cubic line class as a 0/1 sum of toric boundary classes."""

    subset_sums: dict[tuple[int, ...], tuple[int, ...]] = {}
    for mask in range(1 << 9):
        vector = tuple(
            sum(WALLS[index][coordinate] for index in range(9) if mask & (1 << index))
            for coordinate in range(7)
        )
        subset_sums.setdefault(
            vector, tuple(index for index in range(9) if mask & (1 << index))
        )
    return {root: subset_sums[root] for root in exceptional_classes()}


def build_certificate() -> dict[str, object]:
    norms = tuple(pair(root, root) for root in WALLS)
    assert norms == (-2, -1, -2) * 3
    assert tuple(
        {2: 4, 4: 3}[abs(norms[index] * norms[(index + 1) % 9])]
        for index in range(9)
    ) == EDGE_LABELS

    finite_group = finite_a2_cubed_group()
    mod_two_image = {tuple(entry % 2 for entry in matrix) for matrix in finite_group}
    assert len(finite_group) == len(mod_two_image) == 216

    norm_one_roots = tuple(root for root in WALLS if pair(root, root) == -1)
    orbits = [
        {matrix_vector(matrix, root) for matrix in finite_group}
        for root in norm_one_roots
    ]
    assert [len(orbit) for orbit in orbits] == [9, 9, 9]
    assert all(not (orbits[i] & orbits[j]) for i, j in combinations(range(3), 2))
    assert set().union(*orbits) == set(exceptional_classes())

    norm_one_reflections = [reflection(root) for root in norm_one_roots]
    assert all(
        tuple(entry % 2 for entry in matrix)
        == tuple(entry % 2 for entry in IDENTITY)
        for matrix in norm_one_reflections
    )

    representations = boundary_subset_representations()
    assert len(representations) == 27
    assert all(
        tuple(
            sum(WALLS[index][coordinate] for index in subset)
            for coordinate in range(7)
        )
        == root
        for root, subset in representations.items()
    )

    toric_chi, spherical_count = toric_solomon_characteristic()
    p6_cliques, p6_chi = p6_clique_vector_and_characteristic()
    assert toric_chi == Fraction(-1, 1728)
    assert p6_chi == Fraction(-1, 8)
    assert p6_chi / toric_chi == 216
    assert p6_cliques == [1, 27, 216, 720, 1080, 648, 72]

    # In H^6, Gauss--Bonnet gives Vol/pi^3 = -(8/15) chi.
    toric_volume = -Fraction(8, 15) * toric_chi
    p6_volume = -Fraction(8, 15) * p6_chi
    assert toric_volume == Fraction(1, 3240)
    assert p6_volume == Fraction(1, 15)

    return {
        "endpoint_diagram": "cyclic (4,4,3)^3",
        "finite_parabolic": {
            "type": "A2^3",
            "order": len(finite_group),
            "mod_2_image_order": len(mod_two_image),
            "mod_2_reduction_is_injective": True,
        },
        "exceptional_curve_orbits": [len(orbit) for orbit in orbits],
        "orbit_union_is_exactly_all_27_cubic_line_classes": True,
        "all_27_are_zero_one_boundary_sums": True,
        "congruence_argument": {
            "external_input": (
                "Everitt--Ratcliffe--Tschantz: the 27 line reflections generate "
                "the principal level-two group Gamma_2^6 with chamber P^6"
            ),
            "endpoint_group": "Gamma_2^6 semidirect W(A2)^3",
            "endpoint_is_full_mod_2_preimage": True,
            "index_over_Gamma_2_6": 216,
            "P6_endpoint_chamber_count": 216,
        },
        "solomon": {
            "endpoint_spherical_subset_count": spherical_count,
            "endpoint_characteristic": str(toric_chi),
            "P6_right_angled_clique_vector": p6_cliques,
            "P6_characteristic": str(p6_chi),
        },
        "volume": {
            "gauss_bonnet_volume_over_pi_cubed": "-(8/15)*chi",
            "endpoint_volume_over_pi_cubed": str(toric_volume),
            "P6_volume_over_pi_cubed": str(p6_volume),
        },
        "three_A2_smoothing": {
            "central_cubic": "XYZ=W^3",
            "singular_points_char_not_3": [
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
            ],
            "local_equation": "uv=w^3 (A2)",
            "central_lines": ["W=X=0", "W=Y=0", "W=Z=0"],
            "explicit_smoothing": "XYZ-W^3-t(X^3+Y^3+Z^3)=0",
            "generic_smoothness_condition": "t!=0 and 27*t^3!=1 in characteristic zero",
            "marked_line_class_specialization": "three W(A2)^3-orbits of size 9",
            "caveat": (
                "W(A2)^3 is the deck group of the product local simultaneous-"
                "resolution marking cover; a single one-parameter pencil need not "
                "have the full product as its monodromy image"
            ),
        },
    }


if __name__ == "__main__":
    import json

    print(json.dumps(build_certificate(), indent=2, sort_keys=True))

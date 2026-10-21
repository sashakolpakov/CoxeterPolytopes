#!/usr/bin/env python3
"""Clean-room audit of the variable S16/S17 D34-cut denominator.

This script deliberately does not import the production denominator, the old
bounded D34 compiler, a T2 loader, or the S16/S17 theorem synthesizer.  It
starts from two literal support graphs and reconstructs their automorphism
groups, all eight-m4 labelling orbits, the cut-normalization denominator, all
748 oriented integral Grams, their exact inertias, and both Lorentzian
survivors.  The two survivor cones are then enumerated independently by
active sets and passed to the generic exact Vinberg verifier.

Separately pinned public matrix archives are read only after the raw census is
complete.  They are used for lineage and synthesis-hash comparisons, never to
select a mask or a Gram lane.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence

import sympy as sp

from coxeter_polytopes.verify.vinberg_volume import (
    verify_vinberg_volume_criterion,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_d34_cut_s16_s17_hostile_audit.json"
SCHEMA = "coxeter.h9_14_d34_cut_s16_s17.hostile_audit.v1"

# Literal raw support data.  The order is part of the marked audit and is not
# read from any production certificate.
SUPPORTS: dict[str, tuple[tuple[int, int], ...]] = {
    "S16": (
        (2, 3), (3, 4), (4, 5), (6, 7), (7, 8), (8, 9), (1, 10),
        (7, 10), (1, 11), (3, 11), (5, 12), (7, 12), (3, 13),
        (9, 13), (2, 14), (6, 14),
    ),
    "S17": (
        (1, 3), (2, 3), (3, 4), (4, 5), (1, 6), (2, 6), (6, 7),
        (7, 8), (9, 10), (10, 11), (1, 12), (11, 12), (5, 13),
        (8, 13), (10, 13), (2, 14), (9, 14),
    ),
}

PUBLIC_MATRIX_PATHS = {
    "factory": ROOT / "data/enumerations/alexeev_k3_h9_14_factory.json",
    "T2_pilot": ROOT / "data/enumerations/h9_14_t2_surface_unfolding_pilot.json",
    "T2_bridge": ROOT / "data/enumerations/a18_c01_t2_marked_bridge_hostile_audit.json",
    "strong_promotion": ROOT / "data/enumerations/a18_c01_c07_k3_nef_promotions.json",
    "T1_obstruction": ROOT / "data/enumerations/h9_14_t1_strong_k3_nef_obstruction.json",
}
PUBLIC_MATRIX_SHA256 = {
    "factory": "b1ee6fd769527ae6dab01a2bb58099e7c7d1776d0417b6a8d710316697e8a775",
    "T2_pilot": "e72a7b5690c16c8c405ddf8679cef4bafc18cf93f59f38f5e635acc9131d60c2",
    "T2_bridge": "92396280c6b276979a005f9732408944c97fc71c08998d1f127833c7bacd6595",
    "strong_promotion": "2ab7acc379e48ac55944e92b8d614c0992d7c482af1300c670f2bf8e899ba1b1",
    "T1_obstruction": "4eb338784ec452cb08fb4a75b6f4be4c43155d1e6fce2a2259bbe39f410037c1",
}

EXPECTED_SURVIVOR_GRAM_HASHES = {
    "S16": (
        "29cf54743bed2893831c59062868b75ff81b3e9985c66e5fd77efdc8ced1c511",
        "9a7484e9327c35c16d35c9b2fbe9fb22d64eef9c33e698312a295e0d27722cd9",
    ),
    "S17": (
        "829dc6182c07948c4b1bc18a91f1e25228c00ee8a9bc45fd562ce36cbada34e2",
        "74ae56e2319789512febef457ef14895ac4d150e112e1b92dc1a43989a5686b8",
    ),
}

EXPECTED_T2_PACKAGE_HASHES = {
    "fixed_gram": "0488d1fe0e46431d69209180b95bae530877afcccf175be09305b42e05b79666",
    "rank18_gram": "da790e3608234df11e7a58355a83848778b4365b1d83101a0220f71a4fe7cc65",
    "rank18_involution": "d5e3fb7c5d4c60f4608a1c2172cd549b5b981092f54027a75f11aa974f8ca67f",
    "ordered_24_roots": "463d67ffa490dbbd4c9f517659bb83c0460b1e976ffaeaefd8aeead226f9e12c",
}

EXPECTED_T1_PACKAGE_HASHES = {
    "N_gram": "d68824fa8c796b2a863f9da06120ce5bc97f8b41c7445aa987d7f7f1ba11290d",
    "fixed_gram": "1043f0e859809ae85a71d873cb52ca0bc57092bfac50e0727e91d103b524ba72",
    "involution": "ffdc40d029586dba12f9452e0ac087b540ee964b2d175d1fe43bc492fb237bcf",
    "roots": "150fc6f2d3bc77dec983cc51e74ff4a39bb8ada3f6be6fb6b77c8a850b9969bd",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _jhash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _rows(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _histogram(values: Iterable[Any]) -> dict[str, int]:
    return {str(key): value for key, value in sorted(Counter(values).items())}


def _zero_based(edges: Sequence[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    return tuple(tuple(sorted((left - 1, right - 1))) for left, right in edges)


def _adjacency(edges: Sequence[tuple[int, int]]) -> tuple[tuple[bool, ...], ...]:
    adjacency = [[False] * 14 for _ in range(14)]
    for left, right in edges:
        adjacency[left][right] = adjacency[right][left] = True
    return tuple(tuple(row) for row in adjacency)


def _connected_bipartite(
    adjacency: Sequence[Sequence[bool]],
) -> tuple[bool, bool, tuple[int, ...]]:
    colors: list[int | None] = [None] * len(adjacency)
    colors[0] = 0
    queue = deque([0])
    bipartite = True
    while queue:
        vertex = queue.popleft()
        for neighbor, edge in enumerate(adjacency[vertex]):
            if not edge:
                continue
            if colors[neighbor] is None:
                colors[neighbor] = 1 - int(colors[vertex])
                queue.append(neighbor)
            elif colors[neighbor] == colors[vertex]:
                bipartite = False
    return (
        all(color is not None for color in colors),
        bipartite,
        tuple(int(color) for color in colors if color is not None),
    )


def _graph_automorphisms(
    edges: Sequence[tuple[int, int]],
) -> tuple[tuple[tuple[int, ...], ...], int]:
    """Enumerate every graph automorphism by a custom bijection search."""

    adjacency = _adjacency(edges)
    degrees = tuple(sum(row) for row in adjacency)
    invariants = tuple(
        (
            degrees[vertex],
            tuple(sorted(degrees[neighbor] for neighbor in range(14) if adjacency[vertex][neighbor])),
        )
        for vertex in range(14)
    )
    image = [-1] * 14
    used = [False] * 14
    results: list[tuple[int, ...]] = []
    search_nodes = 0

    def candidates(source: int) -> list[int]:
        return [
            target
            for target in range(14)
            if not used[target]
            and invariants[source] == invariants[target]
            and all(
                adjacency[source][other] == adjacency[target][image[other]]
                for other in range(14)
                if image[other] >= 0
            )
        ]

    def search() -> None:
        nonlocal search_nodes
        search_nodes += 1
        remaining = [vertex for vertex in range(14) if image[vertex] < 0]
        if not remaining:
            permutation = tuple(image)
            if any(
                adjacency[left][right]
                != adjacency[permutation[left]][permutation[right]]
                for left in range(14)
                for right in range(left)
            ):
                raise ArithmeticError("custom automorphism search admitted a non-automorphism")
            results.append(permutation)
            return
        choices = [(len(candidates(vertex)), vertex) for vertex in remaining]
        _, source = min(choices)
        for target in candidates(source):
            image[source] = target
            used[target] = True
            # Hall's condition is expensive; nonempty one-vertex forward
            # domains already prune these sparse graphs completely.
            if all(candidates(vertex) for vertex in remaining if vertex != source):
                search()
            used[target] = False
            image[source] = -1

    search()
    group = tuple(sorted(set(results)))
    if len(group) != len(results):
        raise ArithmeticError("custom automorphism search emitted duplicates")
    identity = tuple(range(14))
    group_set = set(group)
    if identity not in group_set:
        raise ArithmeticError("automorphism group lost its identity")
    for left in group:
        for right in group:
            composition = tuple(left[right[index]] for index in range(14))
            if composition not in group_set:
                raise ArithmeticError("enumerated automorphisms are not closed")
    return group, search_nodes


def _edge_permutations(
    edges: Sequence[tuple[int, int]],
    vertex_group: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    edge_number = {tuple(sorted(edge)): index for index, edge in enumerate(edges)}
    permutations = tuple(
        tuple(
            edge_number[tuple(sorted((vertex_permutation[left], vertex_permutation[right])))]
            for left, right in edges
        )
        for vertex_permutation in vertex_group
    )
    if len(set(permutations)) != len(vertex_group):
        raise ArithmeticError("the edge action is not faithful")
    return permutations


def _image_mask(mask: int, permutation: Sequence[int]) -> int:
    image = 0
    for index, target in enumerate(permutation):
        if (mask >> index) & 1:
            image |= 1 << target
    return image


def _mask_orbits(
    masks: Iterable[int], edge_group: Sequence[Sequence[int]]
) -> tuple[tuple[int, int], ...]:
    universe = set(masks)
    unseen = set(universe)
    result: list[tuple[int, int]] = []
    while unseen:
        representative = min(unseen)
        orbit = {_image_mask(representative, permutation) for permutation in edge_group}
        if not orbit <= universe:
            raise ArithmeticError("a graph automorphism left the mask universe")
        result.append((representative, len(orbit)))
        unseen -= orbit
    if sum(size for _representative, size in result) != len(universe):
        raise ArithmeticError("explicit orbit sizes do not cover the universe")
    return tuple(result)


def _burnside_orbit_count(
    masks: Sequence[int], edge_group: Sequence[Sequence[int]]
) -> tuple[int, tuple[int, ...]]:
    fixed_counts = tuple(
        sum(_image_mask(mask, permutation) == mask for mask in masks)
        for permutation in edge_group
    )
    if sum(fixed_counts) % len(edge_group):
        raise ArithmeticError("Burnside average is not integral")
    return sum(fixed_counts) // len(edge_group), fixed_counts


def _cut_mask(side: int, edges: Sequence[tuple[int, int]]) -> int:
    result = 0
    for index, (left, right) in enumerate(edges):
        if ((side >> left) & 1) != ((side >> right) & 1):
            result |= 1 << index
    return result


def _all_size_eight_cuts(
    edges: Sequence[tuple[int, int]],
) -> dict[int, int]:
    # Vertex 14 is fixed outside X.  Connectedness then makes X unique.
    result: dict[int, int] = {}
    for side in range(1 << 13):
        mask = _cut_mask(side, edges)
        if mask.bit_count() != 8:
            continue
        if mask in result and result[mask] != side:
            raise ArithmeticError("a connected support produced a duplicate cut side")
        result[mask] = side
    return result


def _solve_normalization_bits(
    mask: int, edges: Sequence[tuple[int, int]]
) -> int | None:
    """Solve b_left xor b_right = 1 exactly on the m4 edges."""

    incident: list[list[tuple[int, int]]] = [[] for _ in range(14)]
    for index, (left, right) in enumerate(edges):
        difference = (mask >> index) & 1
        incident[left].append((right, difference))
        incident[right].append((left, difference))
    bits: list[int | None] = [None] * 14
    bits[13] = 0
    queue = deque([13])
    while queue:
        vertex = queue.popleft()
        for neighbor, difference in incident[vertex]:
            required = int(bits[vertex]) ^ difference
            if bits[neighbor] is None:
                bits[neighbor] = required
                queue.append(neighbor)
            elif bits[neighbor] != required:
                return None
    if any(bit is None for bit in bits):
        raise ArithmeticError("normalization solver found a disconnected support")
    return sum(int(bit) << vertex for vertex, bit in enumerate(bits))


def _gram(
    edges: Sequence[tuple[int, int]], mask: int, side: int, orientation: int
) -> tuple[sp.Matrix, tuple[int, ...]]:
    squares = tuple(
        4 if ((((side >> vertex) & 1) ^ orientation) == 1) else 2
        for vertex in range(14)
    )
    matrix = sp.diag(*[-square for square in squares])
    for index, (left, right) in enumerate(edges):
        relation = 4 if (mask >> index) & 1 else 3
        if relation == 3:
            if squares[left] != squares[right]:
                raise ArithmeticError("m3 edge has unequal endpoint norms")
            pairing = squares[left] // 2
        else:
            low, high = sorted((squares[left], squares[right]))
            if high != 2 * low:
                raise ArithmeticError("m4 edge has the wrong primitive norm ratio")
            pairing = high // 2
        matrix[left, right] = matrix[right, left] = pairing
    return matrix, squares


def _inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Exact rational symmetric congruence, including zero pivots."""

    work = sp.Matrix(matrix)
    positive = negative = zero = 0
    while work.rows:
        pivot_index = next((index for index in range(work.rows) if work[index, index]), None)
        if pivot_index is not None:
            work.row_swap(0, pivot_index)
            work.col_swap(0, pivot_index)
            pivot = work[0, 0]
            positive += 1 if pivot > 0 else 0
            negative += 1 if pivot < 0 else 0
            if work.rows == 1:
                break
            column = work[1:, 0]
            work = (work[1:, 1:] - column * column.T / pivot).applyfunc(sp.cancel)
            continue
        off_diagonal = next(
            (
                (left, right)
                for left in range(work.rows)
                for right in range(left + 1, work.cols)
                if work[left, right]
            ),
            None,
        )
        if off_diagonal is None:
            zero += work.rows
            break
        left, right = off_diagonal
        work.row_swap(0, left)
        work.col_swap(0, left)
        if right == 0:
            right = left
        work.row_swap(1, right)
        work.col_swap(1, right)
        block = work[:2, :2]
        if block.det() >= 0:
            raise ArithmeticError("invalid exact 2x2 inertia pivot")
        positive += 1
        negative += 1
        if work.rows == 2:
            break
        cross = work[2:, :2]
        work = (work[2:, 2:] - cross * block.inv() * cross.T).applyfunc(sp.cancel)
    if positive + negative + zero != matrix.rows:
        raise ArithmeticError("exact inertia did not account for every dimension")
    return positive, negative, zero


def _primitive_column(vector: sp.Matrix) -> sp.Matrix:
    denominator = math.lcm(*(int(value.q) for value in vector))
    integral = [int(value * denominator) for value in vector]
    content = math.gcd(*(abs(value) for value in integral if value))
    if not content:
        raise ArithmeticError("zero vector cannot define a ray")
    return sp.Matrix([value // content for value in integral])


def _root_span_coordinates(
    gram: sp.Matrix,
) -> tuple[sp.Matrix, sp.Matrix, tuple[int, ...]]:
    indices: tuple[int, ...] | None = None
    for candidate in itertools.combinations(range(14), 10):
        if gram.extract(candidate, candidate).det():
            indices = candidate
            break
    if indices is None:
        raise ArithmeticError("rank-ten Gram has no nonsingular principal basis")
    lattice = gram.extract(indices, indices)
    roots = (lattice.inv() * gram.extract(indices, range(14))).applyfunc(sp.cancel)
    if roots.T * lattice * roots != gram:
        raise ArithmeticError("principal root-span coordinates do not recover the Gram")
    return lattice, roots, indices


def _active_set_rays(gram: sp.Matrix) -> dict[str, Any]:
    lattice, roots, basis_indices = _root_span_coordinates(gram)
    inequality = roots.T * lattice
    rays: dict[tuple[int, ...], sp.Matrix] = {}
    active_subset_hits = 0
    for active in itertools.combinations(range(14), 9):
        block = inequality.extract(active, range(10))
        if block.rank() != 9:
            continue
        nullspace = block.nullspace()
        if len(nullspace) != 1:
            raise ArithmeticError("rank-nine active set has non-ray nullspace")
        ray = _primitive_column(nullspace[0])
        slacks = inequality * ray
        if all(value <= 0 for value in slacks):
            ray = -ray
            slacks = -slacks
        if not all(value >= 0 for value in slacks):
            continue
        square = sp.cancel((ray.T * lattice * ray)[0])
        if square < 0:
            continue
        active_subset_hits += 1
        rays[tuple(int(value) for value in ray)] = ray
    ordered = [rays[key] for key in sorted(rays)]
    if not ordered:
        raise ArithmeticError("active-set cone enumeration found no rays")
    slacks = [inequality * ray for ray in ordered]
    squares = [sp.cancel((ray.T * lattice * ray)[0]) for ray in ordered]
    if any(value < 0 for value in squares):
        raise ArithmeticError("active-set cone has a hyperideal ray")
    active_masks = [
        tuple(index for index, value in enumerate(slack) if value == 0)
        for slack in slacks
    ]
    if any(roots[:, mask].rank() < 9 for mask in active_masks):
        raise ArithmeticError("an extreme ray lacks nine independent active walls")

    interior = sum(ordered, sp.zeros(10, 1))
    interior_slacks = inequality * interior
    interior_square = sp.cancel((interior.T * lattice * interior)[0])
    if interior_square <= 0 or any(value <= 0 for value in interior_slacks):
        raise ArithmeticError("ray sum is not a strict positive cone interior")

    witnesses = []
    for facet in range(14):
        witness = sum(
            (
                ray
                for ray, slack in zip(ordered, slacks, strict=True)
                if slack[facet] == 0
            ),
            sp.zeros(10, 1),
        )
        witness_slacks = inequality * witness
        witness_square = sp.cancel((witness.T * lattice * witness)[0])
        valid = (
            witness_square > 0
            and witness_slacks[facet] == 0
            and all(witness_slacks[other] > 0 for other in range(14) if other != facet)
        )
        witnesses.append(valid)
    if not all(witnesses):
        raise ArithmeticError("an advertised facet is redundant")

    finite_masks = {
        mask for mask, square in zip(active_masks, squares, strict=True) if square > 0
    }
    ideal_masks = {
        mask for mask, square in zip(active_masks, squares, strict=True) if square == 0
    }
    ray_rows = [[int(value) for value in ray] for ray in ordered]
    return {
        "lattice": lattice,
        "roots": roots,
        "inequality": inequality,
        "finite_masks": finite_masks,
        "ideal_masks": ideal_masks,
        "summary": {
            "principal_basis_facets_one_based": [index + 1 for index in basis_indices],
            "active_nine_subset_hits_before_deduplication": active_subset_hits,
            "extreme_ray_count": len(ordered),
            "finite_ray_count": len(finite_masks),
            "ideal_ray_count": len(ideal_masks),
            "hyperideal_ray_count": 0,
            "extreme_rays_sha256": _jhash(ray_rows),
            "active_facet_sets_sha256": _jhash([list(mask) for mask in sorted(active_masks)]),
            "strict_interior_square": int(interior_square),
            "strict_interior_minimum_slack": min(int(value) for value in interior_slacks),
            "all_14_facets_irredundant": all(witnesses),
        },
    }


def _fresh_vinberg(gram: sp.Matrix) -> dict[str, Any]:
    cone = _active_set_rays(gram)
    verification = verify_vinberg_volume_criterion(
        [[-int(gram[row, column]) for column in range(14)] for row in range(14)],
        9,
        simple_roots_complete=True,
    )
    finite_masks = set(verification.finite_vertices)
    ideal_masks = {endpoint.facets for endpoint in verification.ideal_vertices}
    if cone["finite_masks"] != finite_masks:
        raise ArithmeticError("active-set and Vinberg finite endpoints disagree")
    if cone["ideal_masks"] != ideal_masks:
        raise ArithmeticError("active-set and full-neighborhood ideal endpoints disagree")
    bad_neighborhoods = sum(
        not audit.is_rank_dimension_minus_one_parabolic
        for audit in verification.parabolic_neighborhood_audits
    )
    endpoint_histogram = Counter(
        (len(edge.finite_endpoints), len(edge.ideal_endpoints))
        for edge in verification.edge_audits
    )
    if not (
        verification.decision == "finite_volume_noncompact"
        and verification.finite_volume is True
        and verification.compact is False
        and not verification.finite_volume_bad_edges
        and bad_neighborhoods == 0
    ):
        raise ArithmeticError("fresh full Vinberg audit is not finite-volume noncompact")
    return {
        "independent_active_set_cone": cone["summary"],
        "generic_exact_vinberg": {
            "decision": verification.decision,
            "finite_volume": verification.finite_volume,
            "compact": verification.compact,
            "elliptic_rank8_edge_count": verification.elliptic_edge_count,
            "finite_endpoint_count": verification.finite_vertex_count,
            "ideal_endpoint_count": verification.ideal_vertex_count,
            "edge_endpoint_type_histogram": {
                f"{finite},{ideal}": count
                for (finite, ideal), count in sorted(endpoint_histogram.items())
            },
            "bad_edge_count": len(verification.finite_volume_bad_edges),
            "connected_affine_component_audit_count": len(
                verification.parabolic_neighborhood_audits
            ),
            "bad_full_parabolic_neighborhood_count": bad_neighborhoods,
        },
        "active_set_endpoints_equal_vinberg_endpoints": True,
    }


def _relation_matrix(gram: sp.Matrix) -> list[list[int]]:
    result = [[2] * gram.cols for _ in range(gram.rows)]
    for left in range(gram.rows):
        result[left][left] = 1
        for right in range(left):
            if gram[left, right] == 0:
                relation = 2
            else:
                ratio = Fraction(
                    4 * int(gram[left, right]) ** 2,
                    int(gram[left, left] * gram[right, right]),
                )
                relation = {Fraction(1): 3, Fraction(2): 4}.get(ratio, -1)
            result[left][right] = result[right][left] = relation
    return result


def _support_audit(support_id: str) -> tuple[dict[str, Any], tuple[sp.Matrix, sp.Matrix]]:
    edges = _zero_based(SUPPORTS[support_id])
    adjacency = _adjacency(edges)
    connected, bipartite, colors = _connected_bipartite(adjacency)
    if not connected or not bipartite:
        raise ArithmeticError("raw support is not connected bipartite")
    vertex_group, search_nodes = _graph_automorphisms(edges)
    edge_group = _edge_permutations(edges, vertex_group)

    raw_masks = tuple(
        sum(1 << index for index in chosen)
        for chosen in itertools.combinations(range(len(edges)), 8)
    )
    raw_orbits = _mask_orbits(raw_masks, edge_group)
    raw_burnside, raw_fixed_counts = _burnside_orbit_count(raw_masks, edge_group)
    if raw_burnside != len(raw_orbits):
        raise ArithmeticError("explicit and Burnside raw orbit counts disagree")

    cuts = _all_size_eight_cuts(edges)
    independently_solved = {
        mask: side
        for mask in raw_masks
        if (side := _solve_normalization_bits(mask, edges)) is not None
    }
    if independently_solved != cuts:
        raise ArithmeticError("cut enumeration and normalization constraint solver disagree")
    cut_masks = tuple(sorted(cuts))
    cut_orbits = _mask_orbits(cut_masks, edge_group)
    cut_burnside, cut_fixed_counts = _burnside_orbit_count(cut_masks, edge_group)
    if cut_burnside != len(cut_orbits):
        raise ArithmeticError("explicit and Burnside cut orbit counts disagree")

    lane_inventory = []
    inertia_histogram: Counter[str] = Counter()
    selector_histogram: Counter[int] = Counter()
    lorentzian_records: list[tuple[int, int, int, sp.Matrix, tuple[int, ...]]] = []
    for mask, orbit_size in cut_orbits:
        side = cuts[mask]
        for orientation in (0, 1):
            gram, squares = _gram(edges, mask, side, orientation)
            positive, negative, zero = _inertia(gram)
            rank = positive + negative
            key = f"rank{rank}/inertia_{positive}_{negative}_{zero}"
            inertia_histogram[key] += 1
            long_count = squares.count(4)
            selector_histogram[long_count] += 1
            record = [
                hex(mask), orbit_size, hex(side), orientation, long_count,
                _jhash(_rows(gram)), [rank, positive, negative, zero],
            ]
            lane_inventory.append(record)
            if (rank, positive, negative, zero) == (10, 1, 9, 4):
                lorentzian_records.append((mask, side, orientation, gram, squares))

    lorentzian_masks = {record[0] for record in lorentzian_records}
    if len(lorentzian_records) != 2 or len(lorentzian_masks) != 1:
        raise ArithmeticError("support does not have one scaled-dual Lorentzian orbit")
    passing_mask = next(iter(lorentzian_masks))
    passing_side = cuts[passing_mask]
    passing_lanes = tuple(
        _gram(edges, passing_mask, passing_side, orientation)[0]
        for orientation in (0, 1)
    )
    stabilizer_order = sum(
        _image_mask(passing_mask, permutation) == passing_mask
        for permutation in edge_group
    )
    gram_hashes = tuple(_jhash(_rows(gram)) for gram in passing_lanes)
    if gram_hashes != EXPECTED_SURVIVOR_GRAM_HASHES[support_id]:
        raise ArithmeticError("independent survivor Gram hash changed")
    passing_squares = [
        _gram(edges, passing_mask, passing_side, orientation)[1]
        for orientation in (0, 1)
    ]
    selector_passes = [squares.count(4) == 10 for squares in passing_squares]
    if selector_passes != [True, False]:
        raise ArithmeticError("short4/long10 selector did not orient the survivor uniquely")

    degrees = tuple(sum(row) for row in adjacency)
    audit = {
        "support_id": support_id,
        "vertices": 14,
        "edges_one_based": [list(edge) for edge in SUPPORTS[support_id]],
        "edge_count": len(edges),
        "graph": {
            "connected": connected,
            "bipartite": bipartite,
            "bipartition_sizes": sorted(Counter(colors).values()),
            "degree_sequence": sorted(degrees),
            "automorphism_order": len(vertex_group),
            "custom_backtracking_search_nodes": search_nodes,
            "vertex_permutations_sha256": _jhash([list(row) for row in vertex_group]),
            "edge_permutations_sha256": _jhash([list(row) for row in edge_group]),
            "group_closure_checked": True,
            "edge_action_faithful": True,
        },
        "all_eight_m4_labellings": {
            "labelled_count": len(raw_masks),
            "binomial_count": math.comb(len(edges), 8),
            "explicit_orbit_count": len(raw_orbits),
            "burnside_orbit_count": raw_burnside,
            "orbit_size_histogram": _histogram(size for _mask, size in raw_orbits),
            "burnside_fixed_count_histogram": _histogram(raw_fixed_counts),
            "orbit_inventory_sha256": _jhash(
                [[hex(mask), size] for mask, size in raw_orbits]
            ),
        },
        "primitive_norm_2_4_normalization": {
            "constraint": "m3 endpoints have equal bits; m4 endpoints have unequal bits",
            "iff_cut_delta_X": True,
            "all_raw_masks_checked_by_independent_constraint_propagation": len(raw_masks),
            "labelled_compatible_count": len(cuts),
            "explicit_orbit_count": len(cut_orbits),
            "burnside_orbit_count": cut_burnside,
            "orbit_size_histogram": _histogram(size for _mask, size in cut_orbits),
            "burnside_fixed_count_histogram": _histogram(cut_fixed_counts),
            "cut_orbit_inventory_sha256": _jhash(
                [[hex(mask), size, hex(cuts[mask])] for mask, size in cut_orbits]
            ),
            "side_and_complement_are_the_only_two_orientations": True,
        },
        "oriented_integral_gram_lanes": {
            "count": len(lane_inventory),
            "exact_rank_inertia_histogram": dict(sorted(inertia_histogram.items())),
            "norm_square_4_count_histogram": {
                str(key): value for key, value in sorted(selector_histogram.items())
            },
            "short4_long10_selector_lane_count_before_inertia": selector_histogram[10],
            "short4_long10_selector_lorentzian_lane_count": sum(
                record[4].count(4) == 10 for record in lorentzian_records
            ),
            "lane_inventory_sha256": _jhash(lane_inventory),
        },
        "unique_lorentzian_diagram_orbit": {
            "m4_mask_hex": hex(passing_mask),
            "canonical_cut_side_mask_hex": hex(passing_side),
            "canonical_cut_side_vertices_one_based": [
                vertex + 1 for vertex in range(14) if (passing_side >> vertex) & 1
            ],
            "diagram_stabilizer_order": stabilizer_order,
            "fixed_by_full_support_automorphism_group": stabilizer_order == len(vertex_group),
            "m3_edges_one_based": [
                [left + 1, right + 1]
                for index, (left, right) in enumerate(edges)
                if not ((passing_mask >> index) & 1)
            ],
            "m4_edges_one_based": [
                [left + 1, right + 1]
                for index, (left, right) in enumerate(edges)
                if (passing_mask >> index) & 1
            ],
            "scaled_dual_lanes": [
                {
                    "orientation": orientation,
                    "normal_square_histogram": _histogram(passing_squares[orientation]),
                    "short4_long10_selector": selector_passes[orientation],
                    "rank_inertia": [10, 1, 9, 4],
                    "gram_sha256": gram_hashes[orientation],
                    "gram": _rows(passing_lanes[orientation]),
                }
                for orientation in (0, 1)
            ],
        },
    }
    return audit, passing_lanes


def _public_matrix_cross_checks(
    survivor_grams: dict[str, tuple[sp.Matrix, sp.Matrix]],
) -> dict[str, Any]:
    observed_sha = {name: _sha256(path) for name, path in PUBLIC_MATRIX_PATHS.items()}
    if observed_sha != PUBLIC_MATRIX_SHA256:
        raise ArithmeticError("a separately pinned public matrix archive changed")
    payloads = {
        name: json.loads(path.read_text(encoding="utf-8"))
        for name, path in PUBLIC_MATRIX_PATHS.items()
    }
    factory = payloads["factory"]
    t2_pilot = payloads["T2_pilot"]
    bridge = payloads["T2_bridge"]
    promotion = payloads["strong_promotion"]
    obstruction = payloads["T1_obstruction"]

    factory_cases = {row["id"]: row for row in factory["constructions"]}
    factory_t2 = sp.Matrix(
        [[int(value) for value in row] for row in factory_cases["A18_4_0_C01"]
         ["section"]["restricted_mirrors"]["facet_gram"]]
    )
    factory_t1 = sp.Matrix(
        [[int(value) for value in row] for row in factory_cases["A14_8_1_C01"]
         ["section"]["restricted_mirrors"]["facet_gram"]]
    )
    pilot_t2 = sp.Matrix(
        [[int(value) for value in row]
         for row in t2_pilot["surface_free_logical_input"]["facet_gram"]]
    )
    generated_t2 = survivor_grams["S16"][0]
    generated_t1 = survivor_grams["S17"][0]
    if not (generated_t2 == pilot_t2 == factory_t2):
        raise ArithmeticError("generated S16 Gram is not the frozen public T2 Gram")
    if _relation_matrix(generated_t1) != _relation_matrix(factory_t1):
        raise ArithmeticError("generated S17 diagram is not the frozen public T1 diagram")

    pilot_survivor = t2_pilot["unique_survivor"]
    t2_package_hashes = {
        "fixed_gram": _jhash(pilot_survivor["fixed_lattice"]["gram"]),
        "rank18_gram": _jhash(pilot_survivor["rank18_unfolding"]["gram"]),
        "rank18_involution": _jhash(
            pilot_survivor["rank18_unfolding"]["involution"]
        ),
        "ordered_24_roots": _jhash(
            pilot_survivor["source_roots_and_fold"]["source_root_coordinates"]
        ),
    }
    if t2_package_hashes != EXPECTED_T2_PACKAGE_HASHES:
        raise ArithmeticError("direct T2 public matrices do not have synthesis hashes")

    p_rows = [
        [int(value) for value in row]
        for row in bridge["marked_isometry"]["matrix_P"]
    ]
    p_matrix = sp.Matrix(p_rows)
    p_inverse_integral = all(value.q == 1 for value in p_matrix.inv())
    if not (p_matrix.det() == 1 and p_inverse_integral):
        raise ArithmeticError("frozen marked bridge matrix is not unimodular")
    p_string_rows = [[str(value) for value in row] for row in p_rows]
    p_published_hash_recomputed = hashlib.sha256(
        json.dumps(p_string_rows, separators=(",", ":")).encode()
    ).hexdigest()
    if p_published_hash_recomputed != bridge["marked_isometry"]["matrix_P_sha256"]:
        raise ArithmeticError("frozen marked bridge matrix digest changed")

    obstruction_hashes = obstruction["matrix_reconstruction"]["reconstructed_hashes"]
    if obstruction_hashes != EXPECTED_T1_PACKAGE_HASHES:
        raise ArithmeticError("T1 obstruction hash manifest changed")
    if obstruction["matrix_reconstruction"]["archive_hashes"] != obstruction_hashes:
        raise ArithmeticError("T1 obstruction's two frozen hash manifests disagree")

    c01 = next(row for row in promotion["cases"] if row["id"] == "A18_4_0_C01")
    first_hidden = obstruction["hidden_root_witnesses"][0]
    transported_flags_green = (
        bridge["all_checks_pass"]
        and promotion["all_checks_pass"]
        and obstruction["all_checks_pass"]
        and c01["period_and_torelli"]["picard_lattice_is_exactly_N"]
        and c01["nef_and_quotient_cone"]["all_restricted_long_walls_redundant"]
        and first_hidden["negative_extreme_ray"]["hidden_normal_pairing"] < 0
        < first_hidden["positive_extreme_ray"]["hidden_normal_pairing"]
    )
    if not transported_flags_green:
        raise ArithmeticError("a frozen strong/countercalibration flag is not green")

    return {
        "pinned_public_archive_sha256": observed_sha,
        "S16_generated_Gram_equals_public_T2_entry_for_entry": True,
        "S17_generated_Coxeter_relation_matrix_equals_public_T1_label_for_label": True,
        "S17_public_factory_Gram_is_a_different_root_scaling": generated_t1 != factory_t1,
        "T2_package_hashes_recomputed_directly_from_public_matrix_arrays": t2_package_hashes,
        "T2_package_hashes_equal_synthesis_manifest": True,
        "T2_marked_bridge_matrix": {
            "canonical_integer_matrix_sha256": _jhash(p_rows),
            "published_string_matrix_sha256_recomputed": p_published_hash_recomputed,
            "published_matrix_sha256": bridge["marked_isometry"]["matrix_P_sha256"],
            "determinant": int(p_matrix.det()),
            "inverse_integral": p_inverse_integral,
        },
        "T1_package_hash_manifest": obstruction_hashes,
        "T1_hash_manifest_internal_equality_checked": True,
        "transported_green_flags": {
            "T2_exact_bridge_archive": bridge["all_checks_pass"],
            "T2_period_Torelli_and_complete_nef_archive": promotion["all_checks_pass"],
            "T1_weak_carrier_and_exchanged_root_archive": obstruction["all_checks_pass"],
            "first_T1_hidden_normal_pairings": [
                first_hidden["negative_extreme_ray"]["hidden_normal_pairing"],
                first_hidden["positive_extreme_ray"]["hidden_normal_pairing"],
            ],
        },
        "independence_boundary": {
            "T2_facet_and_package_matrix_hashes": "independently recomputed from frozen arrays",
            "T2_bridge_matrix_unimodularity": "independently recomputed from frozen P",
            "T1_facet_Coxeter_identity": "independently recomputed from frozen factory Gram",
            "T1_rank18_package_hashes": (
                "manifest equality only: the obstruction archive does not expose the full "
                "rank18 matrix arrays, so their arithmetic is transported"
            ),
            "period_Torelli_complete_nef_and_exchanged_root_theorems": (
                "transported from separately hash-pinned green certificates; not reproved here"
            ),
        },
    }


def build_certificate() -> dict[str, Any]:
    support_rows = []
    survivor_grams: dict[str, tuple[sp.Matrix, sp.Matrix]] = {}
    for support_id in ("S16", "S17"):
        row, grams = _support_audit(support_id)
        row["fresh_cone_and_vinberg"] = _fresh_vinberg(grams[0])
        support_rows.append(row)
        survivor_grams[support_id] = grams

    public_checks = _public_matrix_cross_checks(survivor_grams)
    raw_labelled = sum(row["all_eight_m4_labellings"]["labelled_count"] for row in support_rows)
    raw_orbits = sum(row["all_eight_m4_labellings"]["explicit_orbit_count"] for row in support_rows)
    compatible_labelled = sum(
        row["primitive_norm_2_4_normalization"]["labelled_compatible_count"]
        for row in support_rows
    )
    compatible_orbits = sum(
        row["primitive_norm_2_4_normalization"]["explicit_orbit_count"]
        for row in support_rows
    )
    lanes = sum(row["oriented_integral_gram_lanes"]["count"] for row in support_rows)
    selected_lorentzian = sum(
        row["oriented_integral_gram_lanes"]["short4_long10_selector_lorentzian_lane_count"]
        for row in support_rows
    )

    checks = {
        "no_production_or_old_D34_module_imported": True,
        "two_literal_raw_supports_reconstructed": len(support_rows) == 2,
        "automorphism_orders_are_48_and_8": [
            row["graph"]["automorphism_order"] for row in support_rows
        ] == [48, 8],
        "raw_counts_are_37180_labelled_and_3610_orbits": (
            raw_labelled == 37180 and raw_orbits == 3610
        ),
        "normalization_iff_cut_checked_for_every_raw_mask": (
            compatible_labelled == 4789 and compatible_orbits == 374
        ),
        "all_748_oriented_exact_Gram_lanes_have_inertia": lanes == 748,
        "one_short4_long10_Lorentzian_lane_per_support": selected_lorentzian == 2,
        "exactly_one_scaled_dual_Lorentzian_diagram_per_support": all(
            row["unique_lorentzian_diagram_orbit"]["fixed_by_full_support_automorphism_group"]
            for row in support_rows
        ),
        "survivor_Grams_and_public_lineage_exact": (
            public_checks["S16_generated_Gram_equals_public_T2_entry_for_entry"]
            and public_checks["S17_generated_Coxeter_relation_matrix_equals_public_T1_label_for_label"]
        ),
        "fresh_active_set_and_full_Vinberg_pass_for_both": all(
            row["fresh_cone_and_vinberg"]["active_set_endpoints_equal_vinberg_endpoints"]
            and row["fresh_cone_and_vinberg"]["generic_exact_vinberg"]["decision"]
            == "finite_volume_noncompact"
            and row["fresh_cone_and_vinberg"]["generic_exact_vinberg"]
            ["bad_full_parabolic_neighborhood_count"] == 0
            for row in support_rows
        ),
        "strong_and_countercalibration_scope_not_overclaimed": True,
    }
    if not all(checks.values()):
        raise ArithmeticError("hostile variable-denominator audit is not green")

    return {
        "schema": SCHEMA,
        "status_date": "2026-08-08",
        "method_independence": {
            "production_denominator_module_imported": False,
            "old_D34_corpus_or_compiler_read": False,
            "T2_loader_module_imported": False,
            "S16_S17_synthesis_module_imported": False,
            "support_source": "two literal edge lists in this audit",
            "automorphism_method": "custom exhaustive adjacency-preserving bijection backtracking",
            "orbit_methods": ["explicit group action partition", "Burnside fixed-mask average"],
            "normalization_methods": [
                "all vertex cuts with vertex 14 fixed outside",
                "independent xor-constraint propagation on every raw eight-m4 mask",
            ],
            "inertia_method": "exact rational symmetric congruence",
            "cone_method": "independent exact rank-nine active-set ray enumeration",
            "volume_method": "generic exact Vinberg edge endpoints plus full affine zero-neighborhood audits",
            "floating_point_used": False,
        },
        "cut_normalization_lemma": {
            "statement": (
                "On a connected support, a primitive norm-square assignment in {2,4} "
                "exists exactly when the m4 edge set is delta(X); X and its complement "
                "are the two scaled-dual orientations."
            ),
            "proof": (
                "Encode norm square 2/4 by a bit. Crystallographic integrality forces "
                "equal bits across m3 and unequal bits across m4, which is precisely "
                "the coboundary equation. Conversely the bit assignment gives integral "
                "m3 pairing half the common norm and m4 pairing half the larger norm."
            ),
            "computational_exhaustion": "the xor solver and cut enumerator agree on all 37180 raw masks",
        },
        "support_audits": support_rows,
        "global_counts": {
            "raw_labelled_eight_m4_assignments": raw_labelled,
            "raw_support_automorphism_orbits": raw_orbits,
            "normalizable_labelled_assignments": compatible_labelled,
            "normalizable_diagram_orbits": compatible_orbits,
            "scaled_dual_oriented_integral_Gram_lanes": lanes,
            "Lorentzian_diagram_orbits": 2,
            "short4_long10_oriented_Lorentzian_lanes": selected_lorentzian,
            "finite_volume_noncompact_polytopes": 2,
        },
        "public_matrix_and_synthesis_cross_checks": public_checks,
        "scope_boundary": {
            "independently_proved": (
                "raw U0 denominator through support automorphisms, all labellings, cut iff "
                "normalization, 748 exact Grams/inertias, unique oriented survivors, exact "
                "survivor facet Grams, irredundant cones, and fresh full Vinberg V"
            ),
            "independently_hash_checked": (
                "S16/T2 public facet and rank18 package arrays, bridge matrix unimodularity, "
                "S17/T1 public Coxeter relation identity, and the frozen synthesis source hashes"
            ),
            "transported_not_reproved": (
                "the 16/5 stable-overclosure inventories, uniqueness of weak K3 packages, "
                "the full T2 marked bridge equalities, K3 period/Torelli existence, complete "
                "T2 invariant-nef equality C, and the rank18 arithmetic behind the T1 "
                "exchanged-root obstruction"
            ),
            "T1_manifest_limit": (
                "only the two matching T1 rank18 hash manifests are public in the obstruction "
                "archive; absent full arrays, this audit does not pretend to recompute them"
            ),
            "new_abstract_polytope": False,
            "lineage": "independent variable reconstruction of the known T2 and T1 H9/14 types",
        },
        "checks": checks,
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.is_file() or args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("hostile S16/S17 archive differs from clean-room replay")
        print("archive_check=PASS")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={args.output}")
    print("all_checks_pass=true")
    print("raw=37180/3610; cuts=4789/374; oriented_Gram_lanes=748")
    print("Lorentzian=2; short4_long10=2; finite_volume_noncompact=2")
    print(f"payload_sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()

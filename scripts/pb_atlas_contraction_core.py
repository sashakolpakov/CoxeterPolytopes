#!/usr/bin/env python3
"""Arrays-only exact contraction core for PB-Atlas numerical Stage A.

The input is a normalized physical Coble branch given only by an integral
Picard Gram matrix, primitive oriented wall vectors, marked boundary indices
and weights, and the canonical class.  No surface name, Coxeter type label,
stored endpoint Gram, contraction witness, or output surface occurs in the
recognizer.

Every current (-1)-curve test is the exact Schur complement

    c^2 - (c.W) Gram(W)^(-1) (W.c) = -1,

with the simultaneous adjunction test K_current.c = -1.  To keep large wall
systems finite in practice, reachable subsets are enumerated modulo the full
automorphism group of the marked weighted wall Gram.  Stabilizer orbits of
the next wall give an exact quotient traversal; expanding the terminal
orbits recovers every labelled frame.

At codimension rho-1 the orthogonal complement is checked to be P2.  At
codimension rho-2 its two-dimensional nef cone is reconstructed from the
input inequalities and checked against a Hirzebruch F_e lattice.  Each
terminal orbit carries a complete endpoint/blow-up coordinate program for a
representative, including a sequential contraction witness, total
exceptional transforms, every pushed wall class, and the marked boundary.

This module deliberately stops before incidence equations, nonemptiness,
bad-factor generation, strong-Coble Hilbert conditions, hidden-wall
exclusion, same-P return, and any L4 claim.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import networkx as nx
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form


INPUT_SCHEMA = "coxeter.pb_atlas.physical_coble_branch.v1"
OUTPUT_SCHEMA = "coxeter.pb_atlas.numerical_contraction_stage_a.v1"


def _json_hash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _integer_matrix(value: Any, name: str) -> sp.Matrix:
    matrix = sp.Matrix(value)
    if any(entry.q != 1 for entry in matrix):
        raise ValueError(f"{name} is not integral")
    return matrix.applyfunc(int)


def _rows(matrix: sp.Matrix) -> list[list[int]]:
    if any(entry.q != 1 for entry in matrix):
        raise ArithmeticError("an integral matrix became nonintegral")
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _vector(vector: sp.Matrix) -> list[int]:
    if vector.cols != 1 or any(entry.q != 1 for entry in vector):
        raise ArithmeticError("an integral column vector became nonintegral")
    return [int(entry) for entry in vector]


def _primitive_entries(values: Sequence[int], *, orient: bool = False) -> tuple[int, ...]:
    divisor = math.gcd(*(abs(int(value)) for value in values))
    if divisor == 0:
        raise ArithmeticError("zero vector has no primitive normalization")
    answer = tuple(int(value) // divisor for value in values)
    if orient and next(value for value in answer if value) < 0:
        answer = tuple(-value for value in answer)
    return answer


def _exact_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Exact rational congruence elimination (shared project algorithm)."""

    work = sp.Matrix(matrix)
    positive = 0
    negative = 0
    while work.rows:
        pivot = next((index for index in range(work.rows) if work[index, index]), None)
        if pivot is None:
            pair = next(
                (
                    (left, right)
                    for left in range(work.rows)
                    for right in range(left + 1, work.rows)
                    if work[left, right]
                ),
                None,
            )
            if pair is None:
                return positive, negative, work.rows
            left, right = pair
            order = [left, right] + [
                index for index in range(work.rows) if index not in pair
            ]
            work = work.extract(order, order)
            change = sp.eye(work.rows)
            change[1, 0] = 1
            work = change.T * work * change
            pivot = 0
        if pivot:
            order = [pivot] + [index for index in range(work.rows) if index != pivot]
            work = work.extract(order, order)
        diagonal = sp.Rational(work[0, 0])
        if diagonal > 0:
            positive += 1
        else:
            negative += 1
        if work.rows == 1:
            return positive, negative, 0
        column = work[1:, 0]
        work = work[1:, 1:] - column * column.T / diagonal
    return positive, negative, 0


@dataclass(frozen=True)
class PhysicalBranch:
    gram: sp.Matrix
    walls: tuple[sp.Matrix, ...]
    boundary_indices: tuple[int, ...]
    boundary_weights: tuple[int, ...]
    canonical: sp.Matrix
    wall_gram: sp.Matrix
    canonical_pairings: tuple[int, ...]
    wall_sublattice_index: int
    input_payload: dict[str, Any]

    @property
    def rank(self) -> int:
        return self.gram.rows

    @property
    def wall_count(self) -> int:
        return len(self.walls)

    @property
    def allowed_indices(self) -> tuple[int, ...]:
        boundary = set(self.boundary_indices)
        return tuple(index for index in range(self.wall_count) if index not in boundary)

    @property
    def boundary_weight_by_index(self) -> tuple[int, ...]:
        result = [0] * self.wall_count
        for index, weight in zip(self.boundary_indices, self.boundary_weights):
            result[index] = weight
        return tuple(result)


def parse_physical_branch(payload: Mapping[str, Any]) -> PhysicalBranch:
    """Validate and normalize one arrays-only physical branch."""

    allowed_keys = {
        "schema",
        "picard_gram",
        "wall_vectors",
        "boundary_indices_zero_based",
        "boundary_weights",
        "canonical_class",
    }
    unknown = set(payload) - allowed_keys
    if unknown:
        raise ValueError(
            "the numerical core accepts arrays only; unknown input fields: "
            + ", ".join(sorted(unknown))
        )
    if payload.get("schema") not in (None, INPUT_SCHEMA):
        raise ValueError("unsupported physical-branch schema")

    gram = _integer_matrix(payload["picard_gram"], "picard_gram")
    if gram.rows != gram.cols or gram != gram.T:
        raise ValueError("picard_gram must be square symmetric")
    rank = gram.rows
    if abs(int(gram.det())) != 1:
        raise ValueError("a normalized smooth rational Picard lattice must be unimodular")
    inertia = _exact_inertia(gram)
    if inertia != (1, rank - 1, 0):
        raise ValueError(f"Picard signature is {inertia}, not (1,{rank - 1},0)")

    raw_walls = payload["wall_vectors"]
    if not isinstance(raw_walls, list) or not raw_walls:
        raise ValueError("wall_vectors must be a nonempty array")
    walls = tuple(_integer_matrix(row, "wall_vector") for row in raw_walls)
    if any(wall.shape != (rank, 1) for wall in walls):
        raise ValueError("every wall vector must have Picard-rank length")
    if any(math.gcd(*(abs(int(entry)) for entry in wall)) != 1 for wall in walls):
        raise ValueError("every wall vector must be primitive")
    wall_matrix = sp.Matrix.hstack(*walls)
    if wall_matrix.rank() != rank:
        raise ValueError("the primitive walls do not span the Picard space")
    wall_hnf = hermite_normal_form(wall_matrix)
    if wall_hnf.shape != (rank, rank):
        raise ValueError("the wall-column HNF does not have full Picard rank")
    wall_sublattice_index = abs(int(wall_hnf.det()))
    wall_gram = wall_matrix.T * gram * wall_matrix
    if any(
        wall_gram[left, right] < 0
        for left in range(len(walls))
        for right in range(left)
    ):
        raise ValueError("wall orientation is not the effective acute orientation")

    boundary_indices = tuple(int(value) for value in payload["boundary_indices_zero_based"])
    if not boundary_indices or len(set(boundary_indices)) != len(boundary_indices):
        raise ValueError("boundary indices must be a nonempty set")
    if any(index < 0 or index >= len(walls) for index in boundary_indices):
        raise ValueError("a boundary index lies outside the wall array")
    boundary_weights = tuple(int(value) for value in payload["boundary_weights"])
    if len(boundary_weights) != len(boundary_indices) or any(
        weight <= 0 for weight in boundary_weights
    ):
        raise ValueError("boundary weights must be positive and parallel the boundary indices")

    canonical = _integer_matrix(payload["canonical_class"], "canonical_class")
    if canonical.shape != (rank, 1):
        raise ValueError("canonical_class must have Picard-rank length")
    boundary_sum = sum(
        (weight * walls[index] for index, weight in zip(boundary_indices, boundary_weights)),
        sp.zeros(rank, 1),
    )
    if boundary_sum != -2 * canonical:
        raise ValueError("weighted boundary sum is not -2K")
    if int((canonical.T * gram * canonical)[0]) != 10 - rank:
        raise ValueError("K^2 is incompatible with a smooth rational surface of this rank")
    standard = sp.eye(rank)
    if not all(
        (
            int((canonical.T * gram * standard[:, index])[0])
            - int(gram[index, index])
        )
        % 2
        == 0
        for index in range(rank)
    ):
        raise ValueError("canonical_class is not characteristic")
    canonical_pairings = tuple(
        int((canonical.T * gram * wall)[0]) for wall in walls
    )

    normalized_payload = {
        "schema": INPUT_SCHEMA,
        "picard_gram": _rows(gram),
        "wall_vectors": [_vector(wall) for wall in walls],
        "boundary_indices_zero_based": list(boundary_indices),
        "boundary_weights": list(boundary_weights),
        "canonical_class": _vector(canonical),
    }
    return PhysicalBranch(
        gram=gram,
        walls=walls,
        boundary_indices=boundary_indices,
        boundary_weights=boundary_weights,
        canonical=canonical,
        wall_gram=wall_gram,
        canonical_pairings=canonical_pairings,
        wall_sublattice_index=wall_sublattice_index,
        input_payload=normalized_payload,
    )


def _permutation_lift(
    branch: PhysicalBranch, permutation: Sequence[int]
) -> sp.Matrix | None:
    """Recover and audit the unique Q-linear lift of a wall permutation."""

    wall_matrix = sp.Matrix.hstack(*branch.walls)
    pivot_columns = tuple(int(index) for index in wall_matrix.rref()[1])
    if len(pivot_columns) != branch.rank:
        raise ArithmeticError("the wall span lost rank while constructing a lift")
    source_basis = wall_matrix[:, pivot_columns]
    target_basis = sp.Matrix.hstack(
        *(branch.walls[int(permutation[index])] for index in pivot_columns)
    )
    lift = target_basis * source_basis.inv()
    if any(entry.q != 1 for entry in lift):
        return None
    lift = lift.applyfunc(int)
    if abs(int(lift.det())) != 1:
        return None
    permuted_walls = sp.Matrix.hstack(
        *(branch.walls[int(permutation[index])] for index in range(branch.wall_count))
    )
    if (
        lift * wall_matrix != permuted_walls
        or lift.T * branch.gram * lift != branch.gram
        or lift * branch.canonical != branch.canonical
    ):
        return None
    return lift


def _marked_wall_gram_permutations(
    branch: PhysicalBranch,
) -> tuple[tuple[int, ...], ...]:
    graph = nx.Graph()
    boundary_weights = branch.boundary_weight_by_index
    for index in range(branch.wall_count):
        graph.add_node(
            index,
            norm=int(branch.wall_gram[index, index]),
            canonical_pairing=branch.canonical_pairings[index],
            boundary_weight=boundary_weights[index],
        )
    for left in range(branch.wall_count):
        for right in range(left):
            pairing = int(branch.wall_gram[left, right])
            if pairing:
                graph.add_edge(left, right, pairing=pairing)
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph,
        graph,
        node_match=lambda left, right: left == right,
        edge_match=lambda left, right: left == right,
    )
    return tuple(
        sorted(
            {
                tuple(int(mapping[index]) for index in range(branch.wall_count))
                for mapping in matcher.isomorphisms_iter()
            }
        )
    )


def wall_gram_automorphisms(branch: PhysicalBranch) -> tuple[tuple[int, ...], ...]:
    """Return exactly the marked Gram permutations that lift integrally."""

    gram_permutations = _marked_wall_gram_permutations(branch)
    # Saturation makes every generator permutation integral.  For a proper
    # full-rank wall sublattice we retain the branch but explicitly filter the
    # unsafe symmetries; the identity always survives, so enumeration remains
    # exact (only less compressed).
    permutations = (
        gram_permutations
        if branch.wall_sublattice_index == 1
        else tuple(
            permutation
            for permutation in gram_permutations
            if _permutation_lift(branch, permutation) is not None
        )
    )
    if not permutations:
        raise ArithmeticError("the identity wall permutation failed its lattice lift")
    validate_automorphisms(branch, permutations)
    return permutations


def validate_automorphisms(
    branch: PhysicalBranch, permutations: Sequence[Sequence[int]]
) -> None:
    expected = tuple(range(branch.wall_count))
    boundary_weights = branch.boundary_weight_by_index
    normalized_rows = tuple(tuple(map(int, permutation)) for permutation in permutations)
    normalized = set(normalized_rows)
    if len(normalized) != len(normalized_rows):
        raise ValueError("the supplied automorphism list contains duplicates")
    if expected not in normalized:
        raise ValueError("the supplied automorphism list omits the identity")
    if any(tuple(sorted(permutation)) != expected for permutation in normalized):
        raise ValueError("a wall automorphism is not a permutation")

    # Canonical masks and stabilizer orbits are valid only for a group action,
    # not for an arbitrary collection of individually valid permutations.
    # Verify subgroup closure efficiently by greedily adding generators and
    # rebuilding their finite permutation subgroup.  Every generated word
    # must stay in the supplied set; after all supplied rows have been
    # considered, equality with that generated subgroup proves closure and
    # inverses as well.  This costs O(|G| * number_of_generators), rather than
    # the quadratic multiplication table of the usually much larger group.
    generators: list[tuple[int, ...]] = []
    generated = {expected}
    for candidate in sorted(normalized):
        if candidate in generated:
            continue
        generators.append(candidate)
        generated = {expected}
        frontier = [expected]
        while frontier:
            current = frontier.pop()
            for generator in generators:
                product = tuple(generator[current[index]] for index in expected)
                if product not in normalized:
                    raise ValueError(
                        "the supplied automorphism list is not closed under composition"
                    )
                if product not in generated:
                    generated.add(product)
                    frontier.append(product)
    if generated != normalized:
        raise ValueError("the supplied automorphism list is not a permutation group")

    for permutation in normalized:
        for left in range(branch.wall_count):
            if (
                boundary_weights[left] != boundary_weights[permutation[left]]
                or branch.canonical_pairings[left]
                != branch.canonical_pairings[permutation[left]]
            ):
                raise ValueError("a supplied permutation changes marked wall data")
            for right in range(branch.wall_count):
                if (
                    branch.wall_gram[left, right]
                    != branch.wall_gram[permutation[left], permutation[right]]
                ):
                    raise ValueError("a supplied permutation changes the wall Gram")

        if branch.wall_sublattice_index != 1 and _permutation_lift(
            branch, permutation
        ) is None:
            raise ValueError(
                "a supplied nonsaturated wall-Gram symmetry has no integral Picard lift"
            )

    # When the wall columns generate Picard integrally, every marked Gram
    # permutation defines an integral map on those generators.  Gram equality
    # preserves all Q-relations (the walls span and the form is nondegenerate),
    # so the map is well-defined; its inverse permutation makes it unimodular.
    # Canonical pairings then force K to be fixed.  Thus saturation is the
    # precise automatic lift gate.  In the nonsaturated case the explicit
    # lift test above replaces this shortcut.


def transform_mask(mask: int, permutation: Sequence[int]) -> int:
    answer = 0
    remaining = mask
    while remaining:
        bit = remaining & -remaining
        index = bit.bit_length() - 1
        answer |= 1 << int(permutation[index])
        remaining -= bit
    return answer


def mask_indices(mask: int) -> tuple[int, ...]:
    answer: list[int] = []
    remaining = mask
    while remaining:
        bit = remaining & -remaining
        answer.append(bit.bit_length() - 1)
        remaining -= bit
    return tuple(answer)


def mask_hex(mask: int, wall_count: int) -> str:
    width = max(1, (wall_count + 3) // 4)
    return f"0x{mask:0{width}x}"


@dataclass(frozen=True)
class FrameState:
    indices: tuple[int, ...]
    inverse: tuple[tuple[int, ...], ...]


def _inverse_for_indices(
    branch: PhysicalBranch, indices: Sequence[int]
) -> tuple[tuple[int, ...], ...]:
    if not indices:
        return tuple()
    matrix = branch.wall_gram.extract(indices, indices)
    inverse = matrix.inv()
    if any(entry.q != 1 for entry in inverse):
        raise ArithmeticError("a reachable exceptional Gram has nonintegral inverse")
    return tuple(
        tuple(int(inverse[left, right]) for right in range(inverse.cols))
        for left in range(inverse.rows)
    )


def _eligible_extension(
    branch: PhysicalBranch,
    state: FrameState,
    candidate: int,
) -> bool:
    size = len(state.indices)
    vector = [int(branch.wall_gram[index, candidate]) for index in state.indices]
    product = [
        sum(state.inverse[left][right] * vector[right] for right in range(size))
        for left in range(size)
    ]
    square = int(branch.wall_gram[candidate, candidate]) - sum(
        vector[index] * product[index] for index in range(size)
    )
    canonical_pairing = branch.canonical_pairings[candidate] - sum(
        branch.canonical_pairings[state.indices[index]] * product[index]
        for index in range(size)
    )
    return square == -1 and canonical_pairing == -1


def _candidate_orbit_representatives(
    branch: PhysicalBranch,
    mask: int,
    stabilizer: Sequence[Sequence[int]],
) -> tuple[int, ...]:
    unassigned = {
        index for index in branch.allowed_indices if not (mask & (1 << index))
    }
    representatives: list[int] = []
    while unassigned:
        representative = min(unassigned)
        orbit = {int(permutation[representative]) for permutation in stabilizer}
        if not orbit <= unassigned:
            # Stabilizers preserve both the mask and the marked boundary set.
            raise ArithmeticError("a stabilizer candidate orbit escaped the allowed complement")
        unassigned -= orbit
        representatives.append(representative)
    return tuple(representatives)


def enumerate_frame_orbits(
    branch: PhysicalBranch,
    *,
    automorphisms: Sequence[Sequence[int]] | None = None,
    maximum_depth: int | None = None,
) -> dict[str, Any]:
    """Exhaust all reachable contraction frames modulo exact wall symmetry."""

    permutations = (
        wall_gram_automorphisms(branch)
        if automorphisms is None
        else tuple(tuple(map(int, permutation)) for permutation in automorphisms)
    )
    validate_automorphisms(branch, permutations)
    depth_limit = branch.rank - 1 if maximum_depth is None else int(maximum_depth)
    if depth_limit < 0 or depth_limit > branch.rank - 1:
        raise ValueError("maximum_depth lies outside the surface-MMP range")

    canonical_cache: dict[int, int] = {}
    orbit_size: dict[int, int] = {}
    orbit_members: dict[int, tuple[int, ...]] = {}
    stabilizer_cache: dict[int, tuple[tuple[int, ...], ...]] = {}

    def canonical(mask: int) -> int:
        if mask in canonical_cache:
            return canonical_cache[mask]
        members = tuple(sorted({transform_mask(mask, permutation) for permutation in permutations}))
        representative = members[0]
        for member in members:
            canonical_cache[member] = representative
        orbit_size[representative] = len(members)
        orbit_members[representative] = members
        return representative

    def stabilizer(mask: int) -> tuple[tuple[int, ...], ...]:
        if mask not in stabilizer_cache:
            stabilizer_cache[mask] = tuple(
                permutation
                for permutation in permutations
                if transform_mask(mask, permutation) == mask
            )
        return stabilizer_cache[mask]

    current: dict[int, FrameState] = {0: FrameState(tuple(), tuple())}
    level_census = [
        {
            "depth": 0,
            "orbit_count": 1,
            "labelled_frame_count": 1,
            "transition_orbit_count": 0,
        }
    ]
    retained: dict[int, dict[int, FrameState]] = {}
    for depth in range(1, depth_limit + 1):
        following: dict[int, FrameState] = {}
        transition_orbits = 0
        for mask, state in current.items():
            for candidate in _candidate_orbit_representatives(
                branch, mask, stabilizer(mask)
            ):
                if not _eligible_extension(branch, state, candidate):
                    continue
                transition_orbits += 1
                representative = canonical(mask | (1 << candidate))
                if representative not in following:
                    indices = mask_indices(representative)
                    following[representative] = FrameState(
                        indices,
                        _inverse_for_indices(branch, indices),
                    )
        current = following
        labelled_count = sum(orbit_size[mask] for mask in current)
        level_census.append(
            {
                "depth": depth,
                "orbit_count": len(current),
                "labelled_frame_count": labelled_count,
                "transition_orbit_count": transition_orbits,
            }
        )
        if depth in (branch.rank - 2, branch.rank - 1):
            retained[depth] = dict(current)

    return {
        "automorphisms": permutations,
        "group_order": len(permutations),
        "level_census": level_census,
        "retained_states": retained,
        "orbit_members": orbit_members,
        "orbit_sizes": orbit_size,
        "canonical_cache_size": len(canonical_cache),
    }


def _orthogonal_data(
    branch: PhysicalBranch, mask: int
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]:
    indices = mask_indices(mask)
    contracted = sp.Matrix.hstack(*(branch.walls[index] for index in indices))
    exceptional_gram = contracted.T * branch.gram * contracted
    projector = (
        sp.eye(branch.rank)
        - contracted * exceptional_gram.inv() * contracted.T * branch.gram
    )
    if any(entry.q != 1 for entry in projector):
        raise ArithmeticError("the contraction projector is not integral")
    projector = projector.applyfunc(int)
    complement = hermite_normal_form(projector)
    expected_rank = branch.rank - len(indices)
    if complement.shape != (branch.rank, expected_rank):
        raise ArithmeticError("HNF did not return the full orthogonal complement lattice")
    complement_gram = complement.T * branch.gram * complement
    canonical_image = projector * branch.canonical
    canonical_coordinates = (
        complement_gram.inv() * complement.T * branch.gram * canonical_image
    )
    if any(entry.q != 1 for entry in canonical_coordinates):
        raise ArithmeticError("the endpoint canonical class is nonintegral")
    return projector, complement, complement_gram, canonical_coordinates


def _two_dimensional_nef_rays(
    branch: PhysicalBranch, complement: sp.Matrix
) -> tuple[tuple[int, int], tuple[int, int]]:
    inequalities = sp.Matrix.hstack(*branch.walls).T * branch.gram * complement
    rays: set[tuple[int, int]] = set()
    for row in range(inequalities.rows):
        left = int(inequalities[row, 0])
        right = int(inequalities[row, 1])
        if left == right == 0:
            continue
        for candidate in ((right, -left), (-right, left)):
            values = inequalities * sp.Matrix(candidate)
            if all(value >= 0 for value in values):
                rays.add(_primitive_entries(candidate))
    if len(rays) != 2:
        raise ArithmeticError(
            f"the contracted rank-two nef cone has {len(rays)} rays instead of two"
        )
    return tuple(sorted(rays))  # type: ignore[return-value]


def classify_endpoint(branch: PhysicalBranch, mask: int) -> dict[str, Any]:
    projector, complement, gram, canonical = _orthogonal_data(branch, mask)
    endpoint_rank = complement.cols
    if endpoint_rank == 1:
        generator = complement[:, 0]
        intersections = [
            int((wall.T * branch.gram * generator)[0]) for wall in branch.walls
        ]
        if all(value <= 0 for value in intersections):
            generator = -generator
            intersections = [-value for value in intersections]
        if not all(value >= 0 for value in intersections):
            raise ArithmeticError("the rank-one complement has no nef orientation")
        if (
            int((generator.T * branch.gram * generator)[0]) != 1
            or int((branch.canonical.T * branch.gram * generator)[0]) != -3
        ):
            raise ArithmeticError("a codimension rho-1 frame is not a P2 endpoint")
        return {
            "kind": "P2",
            "model_parameter_e": None,
            "projector": projector,
            "model_basis": generator,
            "nef_rays_in_model_basis": [(1,)],
        }
    if endpoint_rank != 2:
        raise ValueError("endpoint classification is only defined in ranks one and two")

    rays = _two_dimensional_nef_rays(branch, complement)
    candidates: list[tuple[tuple[Any, ...], int, sp.Matrix, sp.Matrix]] = []
    wall_matrix = sp.Matrix.hstack(*branch.walls)
    for first, second in ((rays[0], rays[1]), (rays[1], rays[0])):
        fiber_coordinates = sp.Matrix(first)
        positive_coordinates = sp.Matrix(second)
        e = int((positive_coordinates.T * gram * positive_coordinates)[0])
        if not (
            (fiber_coordinates.T * gram * fiber_coordinates)[0] == 0
            and (fiber_coordinates.T * gram * positive_coordinates)[0] == 1
            and (canonical.T * gram * fiber_coordinates)[0] == -2
            and e >= 0
            and (canonical.T * gram * positive_coordinates)[0] == -(e + 2)
        ):
            continue
        section_coordinates = positive_coordinates - e * fiber_coordinates
        model_basis = sp.Matrix.hstack(
            complement * section_coordinates,
            complement * fiber_coordinates,
        )
        model_gram = model_basis.T * branch.gram * model_basis
        if model_gram != sp.Matrix([[-e, 1], [1, 0]]):
            raise ArithmeticError("the proposed Hirzebruch basis has the wrong Gram")
        endpoint_coordinates = (
            model_gram.inv() * model_basis.T * branch.gram * projector * wall_matrix
        ).T
        if any(entry.q != 1 for entry in endpoint_coordinates):
            raise ArithmeticError("a wall image is nonintegral on the ruled endpoint")
        # This key is unchanged by Picard basis mutation and wall permutation.
        profile = tuple(
            sorted(
                (
                    branch.boundary_weight_by_index[index],
                    int(endpoint_coordinates[index, 0]),
                    int(endpoint_coordinates[index, 1]),
                )
                for index in range(branch.wall_count)
            )
        )
        candidates.append((profile, e, model_basis, endpoint_coordinates))
    if not candidates:
        raise ArithmeticError("a codimension rho-2 frame is not a Hirzebruch endpoint")
    profile, e, model_basis, _ = min(candidates, key=lambda value: value[0])
    return {
        "kind": "F_e",
        "model_parameter_e": e,
        "minimal": e != 1,
        "projector": projector,
        "model_basis": model_basis,
        "nef_rays_in_model_basis": [(0, 1), (1, e)],
        "canonical_wall_image_profile": profile,
    }


def _lexicographic_sequential_order(branch: PhysicalBranch, target_mask: int) -> tuple[int, ...]:
    memo: dict[int, tuple[int, ...] | None] = {0: tuple()}

    def visit(mask: int) -> tuple[int, ...] | None:
        if mask in memo:
            return memo[mask]
        for last in mask_indices(mask):
            previous = mask & ~(1 << last)
            prefix = visit(previous)
            if prefix is None:
                continue
            state = FrameState(mask_indices(previous), _inverse_for_indices(branch, mask_indices(previous)))
            if _eligible_extension(branch, state, last):
                memo[mask] = (*prefix, last)
                return memo[mask]
        memo[mask] = None
        return None

    answer = visit(target_mask)
    if answer is None:
        raise ArithmeticError("a retained frame has no sequential (-1)-contraction order")
    return answer


def _total_exceptional_transforms(
    branch: PhysicalBranch, order: Sequence[int]
) -> tuple[sp.Matrix, ...]:
    strict: list[sp.Matrix] = []
    total: list[sp.Matrix] = []
    for wall_index in order:
        curve = branch.walls[wall_index]
        if strict:
            strict_matrix = sp.Matrix.hstack(*strict)
            strict_gram = strict_matrix.T * branch.gram * strict_matrix
            pairings = strict_matrix.T * branch.gram * curve
            exceptional = curve - strict_matrix * strict_gram.inv() * pairings
        else:
            exceptional = curve
        if any(entry.q != 1 for entry in exceptional):
            raise ArithmeticError("a total exceptional transform is nonintegral")
        exceptional = exceptional.applyfunc(int)
        if (
            int((exceptional.T * branch.gram * exceptional)[0]) != -1
            or int((branch.canonical.T * branch.gram * exceptional)[0]) != -1
            or any((old.T * branch.gram * exceptional)[0] for old in total)
        ):
            raise ArithmeticError("the sequential exceptional basis failed its exact audit")
        strict.append(curve)
        total.append(exceptional)
    return tuple(total)


def endpoint_program(branch: PhysicalBranch, mask: int) -> dict[str, Any]:
    """Emit the complete numerical packet needed by the later incidence stage."""

    endpoint = classify_endpoint(branch, mask)
    order = _lexicographic_sequential_order(branch, mask)
    exceptionals = _total_exceptional_transforms(branch, order)
    model_basis: sp.Matrix = endpoint["model_basis"]
    basis = sp.Matrix.hstack(model_basis, *exceptionals)
    if basis.shape != (branch.rank, branch.rank) or abs(int(basis.det())) != 1:
        raise ArithmeticError("the endpoint plus total exceptionals is not a unimodular basis")
    model_rank = model_basis.cols
    e = endpoint["model_parameter_e"]
    expected_model_gram = (
        sp.Matrix([[1]]) if endpoint["kind"] == "P2" else sp.Matrix([[-e, 1], [1, 0]])
    )
    expected_gram = sp.diag(expected_model_gram, *([-1] * len(exceptionals)))
    actual_gram = basis.T * branch.gram * basis
    if actual_gram != expected_gram:
        raise ArithmeticError("the endpoint/blow-up basis Gram is incorrect")
    coordinates = basis.inv() * sp.Matrix.hstack(*branch.walls)
    canonical_coordinates = basis.inv() * branch.canonical
    if any(entry.q != 1 for entry in coordinates) or any(
        entry.q != 1 for entry in canonical_coordinates
    ):
        raise ArithmeticError("the endpoint program acquired nonintegral coordinates")
    expected_canonical = (
        [-3]
        if endpoint["kind"] == "P2"
        else [-2, -(int(e) + 2)]
    ) + [1] * len(exceptionals)
    if _vector(canonical_coordinates) != expected_canonical:
        raise ArithmeticError("the endpoint canonical-coordinate formula changed")

    coordinate_rows = [
        [int(coordinates[column, index]) for column in range(branch.rank)]
        for index in range(branch.wall_count)
    ]
    endpoint_rows = [row[:model_rank] for row in coordinate_rows]
    selected = set(mask_indices(mask))
    profile_records = sorted(
        (
            branch.boundary_weight_by_index[index],
            index in selected,
            int(branch.wall_gram[index, index]),
            branch.canonical_pairings[index],
            tuple(endpoint_rows[index]),
            tuple(
                sorted(int(branch.wall_gram[index, chosen]) for chosen in selected)
            ),
        )
        for index in range(branch.wall_count)
    )
    trace = []
    for position, (wall_index, exceptional) in enumerate(zip(order, exceptionals)):
        trace.append(
            {
                "step": position + 1,
                "source_wall_index_zero_based": wall_index,
                "current_total_exceptional_in_input_basis": _vector(exceptional),
                "current_square": int((exceptional.T * branch.gram * exceptional)[0]),
                "current_canonical_pairing": int(
                    (branch.canonical.T * branch.gram * exceptional)[0]
                ),
            }
        )

    return {
        "endpoint_kind": endpoint["kind"],
        "hirzebruch_e": e,
        "hirzebruch_minimal": endpoint.get("minimal"),
        "contracted_wall_indices_zero_based": list(mask_indices(mask)),
        "sequential_contraction_order_zero_based": list(order),
        "strict_contracted_wall_gram": _rows(
            branch.wall_gram.extract(mask_indices(mask), mask_indices(mask))
        ),
        "model_basis_in_input_picard_coordinates": _rows(model_basis.T),
        "total_exceptional_basis_in_input_picard_coordinates": [
            _vector(exceptional) for exceptional in exceptionals
        ],
        "endpoint_blowup_basis_gram": _rows(actual_gram),
        "canonical_coordinates_in_endpoint_blowup_basis": expected_canonical,
        "nef_rays_in_model_basis": [
            list(ray) for ray in endpoint["nef_rays_in_model_basis"]
        ],
        "all_wall_coordinates_in_endpoint_blowup_basis": coordinate_rows,
        "all_wall_endpoint_image_coordinates": endpoint_rows,
        "boundary": [
            {
                "wall_index_zero_based": index,
                "anti_bicanonical_weight": weight,
                "endpoint_blowup_coordinates": coordinate_rows[index],
                "endpoint_image_coordinates": endpoint_rows[index],
            }
            for index, weight in zip(branch.boundary_indices, branch.boundary_weights)
        ],
        "contraction_trace": trace,
        "permutation_invariant_profile_sha256": _json_hash(profile_records),
        "stage_B_payload_complete": True,
    }


def _endpoint_orbits(
    branch: PhysicalBranch,
    states: Mapping[int, FrameState],
    orbit_members: Mapping[int, Sequence[int]],
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    records: list[dict[str, Any]] = []
    histogram: dict[str, int] = {}
    for representative in sorted(states):
        program = endpoint_program(branch, representative)
        members = tuple(orbit_members[representative])
        key = (
            "P2"
            if program["endpoint_kind"] == "P2"
            else f"F_{program['hirzebruch_e']}"
        )
        histogram[key] = histogram.get(key, 0) + len(members)
        records.append(
            {
                "representative_mask": mask_hex(representative, branch.wall_count),
                "orbit_size": len(members),
                "member_masks": [mask_hex(mask, branch.wall_count) for mask in members],
                "program": program,
            }
        )
    return records, dict(sorted(histogram.items()))


def certify_branch(
    payload: Mapping[str, Any],
    *,
    automorphisms: Sequence[Sequence[int]] | None = None,
    retain_automorphisms: bool = True,
) -> dict[str, Any]:
    branch = parse_physical_branch(payload)
    enumeration = enumerate_frame_orbits(branch, automorphisms=automorphisms)
    ruled_states = enumeration["retained_states"].get(branch.rank - 2, {})
    plane_states = enumeration["retained_states"].get(branch.rank - 1, {})
    ruled, ruled_histogram = _endpoint_orbits(
        branch, ruled_states, enumeration["orbit_members"]
    )
    plane, plane_histogram = _endpoint_orbits(
        branch, plane_states, enumeration["orbit_members"]
    )
    if set(plane_histogram) - {"P2"}:
        raise ArithmeticError("a codimension rho-1 orbit is not a plane endpoint")
    if any(not key.startswith("F_") for key in ruled_histogram):
        raise ArithmeticError("a codimension rho-2 orbit is not a ruled endpoint")

    permutations = enumeration["automorphisms"]
    result = {
        "schema": OUTPUT_SCHEMA,
        "scope": (
            "PB-Atlas numerical Stage A only: exact all-contraction frames and endpoint "
            "coordinate programs; no incidence nonemptiness, bad-factor generation, "
            "strong-Coble Hilbert gate, hidden-wall exclusion, same-P return, or L4"
        ),
        "input": {
            "rank": branch.rank,
            "wall_count": branch.wall_count,
            "boundary_count": len(branch.boundary_indices),
            "boundary_indices_zero_based": list(branch.boundary_indices),
            "boundary_weights": list(branch.boundary_weights),
            "array_payload_sha256": _json_hash(branch.input_payload),
            "only_arrays_entered_recognizer": True,
        },
        "validation": {
            "picard_determinant": int(branch.gram.det()),
            "picard_signature": list(_exact_inertia(branch.gram)),
            "wall_span_rank": sp.Matrix.hstack(*branch.walls).rank(),
            "wall_sublattice_index": branch.wall_sublattice_index,
            "canonical_square": int(
                (branch.canonical.T * branch.gram * branch.canonical)[0]
            ),
            "weighted_boundary_sum_is_minus_two_K": True,
            "all_walls_primitive_and_effectively_oriented": True,
        },
        "symmetry_quotient": {
            "full_marked_wall_gram_automorphism_group_order": enumeration["group_order"],
            "automorphism_permutations_zero_based": (
                [list(permutation) for permutation in permutations]
                if retain_automorphisms
                else None
            ),
            "canonical_mask_cache_size": enumeration["canonical_cache_size"],
            "exactness": (
                "only integrally liftable marked-Gram permutations are retained; at each orbit "
                "representative, candidates are partitioned by its full stabilizer, so every "
                "labelled transition is conjugate to exactly one tested orbit"
            ),
        },
        "contraction_level_census": enumeration["level_census"],
        "hirzebruch_endpoint_orbits": ruled,
        "hirzebruch_labelled_histogram": ruled_histogram,
        "P2_endpoint_orbits": plane,
        "P2_labelled_count": plane_histogram.get("P2", 0),
        "stage_B_contract": {
            "available": [
                "one valid sequential blowdown order for every terminal orbit",
                "total exceptional transforms and exact endpoint/blow-up basis",
                "all wall and marked-boundary coordinates in that basis",
                "P2 degree or F_e bidegree/section-fiber data",
            ],
            "still_missing": [
                "symbolic incidence and proximity equation generation",
                "collision/Jacobian/irreducibility/Hilbert bad factors",
                "nonempty regular-open proof",
                "hidden-wall exclusion and literal same-P return",
            ],
        },
        "all_checks_pass": True,
    }
    return result


def _load_cli_payload(path: Path, case_index: int | None) -> Mapping[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if case_index is None:
        return payload
    cases = payload.get("cases")
    if not isinstance(cases, list):
        raise ValueError("--case-index requires a fixture object with a cases array")
    return cases[case_index]["input"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--case-index", type=int)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--omit-automorphisms", action="store_true")
    arguments = parser.parse_args()
    payload = _load_cli_payload(arguments.input, arguments.case_index)
    certificate = certify_branch(
        payload, retain_automorphisms=not arguments.omit_automorphisms
    )
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if arguments.output is None:
        print(rendered, end="")
    else:
        arguments.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

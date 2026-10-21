"""Exact finite enumeration from rational Cartan cycle packets.

This module implements the effective middle of the bounded-frame proof in
``results/ag_field_funnel_bounded_frame_native_host_closure.md`` for the
ground field Q and small rank.  It is deliberately mathematical rather than
a report-text checker:

* all connected labelled edge-product assignments in a finite box are
  enumerated;
* the complete simple-cycle packet is built and tested for Q-integrality;
* one Cartan matrix is reconstructed, exactly, modulo diagonal root gauge;
* an integral gauge is found by prime-valuation difference constraints;
* Lorentzian signature and the finite-volume simplex links are checked by
  exact congruence elimination; and
* a repository chamber certificate can be emitted and independently replayed.

The current backend is Q-only.  General number fields require a Minkowski
box enumerator, ideal-valued gauges, and exact sign oracles at all real
embeddings; the data interfaces here are designed so those are separable
extensions rather than changes to the cycle reconstruction theorem.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from math import isqrt
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.exact.matrix import Inertia, inertia
from coxeter_polytopes.verify.chamber import (
    ChamberVerification,
    verify_chamber_geometry,
)


Matrix = tuple[tuple[Fraction, ...], ...]
Cycle = tuple[int, ...]
Edge = tuple[int, int]


class PacketError(ValueError):
    """A cycle packet is incomplete or violates an exact compatibility."""


def _fraction_string(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _edge(i: int, j: int) -> Edge:
    if i == j:
        raise ValueError("a graph edge needs two distinct vertices")
    return (i, j) if i < j else (j, i)


def canonical_cycle(vertices: Iterable[int]) -> Cycle:
    """Canonicalize a simple cycle modulo rotation and reversal.

    Length-two backtracks are represented by the sorted edge.  For length at
    least three this quotient is valid for symmetrizable Cartan matrices,
    where forward and reverse cycle products agree.
    """

    cycle = tuple(vertices)
    if len(cycle) < 2:
        raise ValueError("a cycle must have length at least two")
    if len(set(cycle)) != len(cycle):
        raise ValueError("cycle vertices must be distinct (the closing vertex is implicit)")
    if len(cycle) == 2:
        return _edge(*cycle)
    rotations = []
    for oriented in (cycle, tuple(reversed(cycle))):
        rotations.extend(
            oriented[offset:] + oriented[:offset] for offset in range(len(oriented))
        )
    return min(rotations)


def _connected(rank: int, edges: Iterable[Edge]) -> bool:
    adjacency = [set() for _ in range(rank)]
    for i, j in edges:
        adjacency[i].add(j)
        adjacency[j].add(i)
    reached = {0}
    frontier = [0]
    while frontier:
        i = frontier.pop()
        for j in adjacency[i] - reached:
            reached.add(j)
            frontier.append(j)
    return len(reached) == rank


def simple_cycles(rank: int, edges: Iterable[Edge]) -> tuple[Cycle, ...]:
    """List every unoriented simple graph cycle, including edge backtracks."""

    edge_set = frozenset(_edge(*edge) for edge in edges)
    cycles: set[Cycle] = set(edge_set)
    for length in range(3, rank + 1):
        for vertices in combinations(range(rank), length):
            first = min(vertices)
            remainder = tuple(vertex for vertex in vertices if vertex != first)
            for order in permutations(remainder):
                cycle = (first, *order)
                # Fix the remaining reversal ambiguity before the edge test.
                if cycle[1] > cycle[-1]:
                    continue
                if all(
                    _edge(cycle[position], cycle[(position + 1) % length])
                    in edge_set
                    for position in range(length)
                ):
                    cycles.add(canonical_cycle(cycle))
    return tuple(sorted(cycles, key=lambda cycle: (len(cycle), cycle)))


def rational_square_root(value: Fraction) -> Fraction | None:
    """Return the nonnegative rational square root, or ``None``."""

    value = Fraction(value)
    if value < 0:
        return None
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator * numerator != value.numerator:
        return None
    if denominator * denominator != value.denominator:
        return None
    return Fraction(numerator, denominator)


@dataclass(frozen=True, slots=True)
class CyclePacket:
    """A complete tuple of simple closed Cartan products over Q."""

    rank: int
    values: tuple[tuple[Cycle, Fraction], ...]

    def __post_init__(self) -> None:
        if self.rank < 2:
            raise PacketError("a connected packet must have rank at least two")
        normalized: list[tuple[Cycle, Fraction]] = []
        seen: set[Cycle] = set()
        for raw_cycle, raw_value in self.values:
            cycle = canonical_cycle(raw_cycle)
            if cycle != tuple(raw_cycle):
                raise PacketError(f"cycle {raw_cycle} is not canonically represented")
            if any(vertex < 0 or vertex >= self.rank for vertex in cycle):
                raise PacketError(f"cycle {cycle} has a vertex outside the packet rank")
            if cycle in seen:
                raise PacketError(f"cycle {cycle} occurs more than once")
            seen.add(cycle)
            normalized.append((cycle, Fraction(raw_value)))
        if tuple(normalized) != tuple(sorted(normalized, key=lambda item: (len(item[0]), item[0]))):
            raise PacketError("packet cycles must be sorted by length and lexicographically")
        object.__setattr__(self, "values", tuple(normalized))

    @classmethod
    def from_edge_products(
        cls,
        rank: int,
        edge_products: Mapping[Edge, Fraction | int],
        *,
        require_integral: bool = True,
    ) -> "CyclePacket":
        """Build the unique acute symmetrizable packet over Q.

        The length-r cycle has sign ``(-1)^r`` and square equal to the
        product of its length-two edge invariants.  Failure of that product
        to be a square in Q means that the minimal cycle field is not Q.
        """

        normalized: dict[Edge, Fraction] = {}
        for raw_edge, raw_value in edge_products.items():
            edge = _edge(*raw_edge)
            value = Fraction(raw_value)
            if value <= 0:
                raise PacketError("listed graph edges must have positive products")
            if edge in normalized and normalized[edge] != value:
                raise PacketError(f"edge {edge} has conflicting products")
            normalized[edge] = value
        if not _connected(rank, normalized):
            raise PacketError("the nonorthogonality graph is disconnected")

        values: list[tuple[Cycle, Fraction]] = []
        for cycle in simple_cycles(rank, normalized):
            if len(cycle) == 2:
                value = normalized[_edge(*cycle)]
            else:
                squared = Fraction(1)
                for position, i in enumerate(cycle):
                    squared *= normalized[_edge(i, cycle[(position + 1) % len(cycle)])]
                root = rational_square_root(squared)
                if root is None:
                    raise PacketError(
                        f"cycle {cycle} has nonsquare edge-product {squared}; its field is not Q"
                    )
                value = -root if len(cycle) % 2 else root
            if require_integral and value.denominator != 1:
                raise PacketError(f"cycle {cycle} has nonintegral value {value}")
            values.append((cycle, value))
        return cls(rank, tuple(values))

    @property
    def dictionary(self) -> dict[Cycle, Fraction]:
        return dict(self.values)

    @property
    def edge_products(self) -> dict[Edge, Fraction]:
        return {
            _edge(*cycle): value
            for cycle, value in self.values
            if len(cycle) == 2
        }

    @property
    def edges(self) -> tuple[Edge, ...]:
        return tuple(self.edge_products)

    def validate_complete(self) -> None:
        if not _connected(self.rank, self.edges):
            raise PacketError("the packet graph is disconnected")
        expected = simple_cycles(self.rank, self.edges)
        actual = tuple(cycle for cycle, _ in self.values)
        if actual != expected:
            missing = sorted(set(expected) - set(actual))
            extra = sorted(set(actual) - set(expected))
            raise PacketError(f"packet is not complete: missing={missing}, extra={extra}")
        for edge, value in self.edge_products.items():
            if value <= 0:
                raise PacketError(f"edge product {edge} must be positive")

    def relabel(self, old_to_new: Sequence[int]) -> "CyclePacket":
        if sorted(old_to_new) != list(range(self.rank)):
            raise ValueError("a relabelling must be a permutation of the vertices")
        relabelled = []
        for cycle, value in self.values:
            relabelled.append(
                (canonical_cycle(old_to_new[vertex] for vertex in cycle), value)
            )
        return CyclePacket(
            self.rank,
            tuple(sorted(relabelled, key=lambda item: (len(item[0]), item[0]))),
        )

    def key(self) -> tuple[tuple[Cycle, Fraction], ...]:
        return self.values

    def canonicalized(self) -> "CyclePacket":
        """Return the lexicographically least packet under facet permutation."""

        if self.rank > 8:
            raise ValueError("factorial packet canonicalization is restricted to rank at most 8")
        return min(
            (self.relabel(permutation) for permutation in permutations(range(self.rank))),
            key=CyclePacket.key,
        )

    def to_json(self) -> dict[str, object]:
        return {
            "rank": self.rank,
            "cycles": [
                {"vertices": list(cycle), "value": _fraction_string(value)}
                for cycle, value in self.values
            ],
        }


def _tree_data(rank: int, edges: Sequence[Edge]) -> tuple[list[int | None], set[Edge]]:
    adjacency = [set() for _ in range(rank)]
    for i, j in edges:
        adjacency[i].add(j)
        adjacency[j].add(i)
    parent: list[int | None] = [None] * rank
    parent[0] = 0
    queue = [0]
    for i in queue:
        for j in sorted(adjacency[i]):
            if parent[j] is None:
                parent[j] = i
                queue.append(j)
    if any(value is None for value in parent):
        raise PacketError("the packet graph is disconnected")
    tree = {_edge(i, value) for i, value in enumerate(parent) if i and value is not None}
    return parent, tree


def _tree_path(rank: int, tree: set[Edge], start: int, finish: int) -> tuple[int, ...]:
    adjacency = [set() for _ in range(rank)]
    for i, j in tree:
        adjacency[i].add(j)
        adjacency[j].add(i)
    predecessor: dict[int, int | None] = {start: None}
    queue = [start]
    for i in queue:
        if i == finish:
            break
        for j in sorted(adjacency[i]):
            if j not in predecessor:
                predecessor[j] = i
                queue.append(j)
    if finish not in predecessor:
        raise PacketError("the proposed spanning tree is disconnected")
    reverse = [finish]
    while reverse[-1] != start:
        previous = predecessor[reverse[-1]]
        assert previous is not None
        reverse.append(previous)
    return tuple(reversed(reverse))


def directed_cycle_product(cartan: Matrix, cycle: Sequence[int]) -> Fraction:
    value = Fraction(1)
    for position, i in enumerate(cycle):
        value *= cartan[i][cycle[(position + 1) % len(cycle)]]
    return value


def cartan_cycle_packet(cartan: Matrix) -> CyclePacket:
    rank = len(cartan)
    if any(len(row) != rank for row in cartan):
        raise ValueError("Cartan matrix must be square")
    edges = {
        (i, j)
        for i in range(rank)
        for j in range(i + 1, rank)
        if cartan[i][j] or cartan[j][i]
    }
    if any(bool(cartan[i][j]) != bool(cartan[j][i]) for i, j in edges):
        raise PacketError("Cartan zero pattern is not symmetric")
    values = tuple(
        (cycle, directed_cycle_product(cartan, cycle))
        for cycle in simple_cycles(rank, edges)
    )
    return CyclePacket(rank, values)


def reconstruct_cartan(packet: CyclePacket) -> Matrix:
    """Reconstruct the unique Cartan matrix class in a spanning-tree gauge."""

    packet.validate_complete()
    rank = packet.rank
    values = packet.dictionary
    edge_products = packet.edge_products
    parent, tree = _tree_data(rank, packet.edges)
    cartan = [[Fraction(0) for _ in range(rank)] for _ in range(rank)]
    for i in range(rank):
        cartan[i][i] = Fraction(2)

    # The gauge a_parent,child=-1 fixes all diagonal freedom except a common scalar.
    for child in range(1, rank):
        p = parent[child]
        assert p is not None
        c = edge_products[_edge(p, child)]
        cartan[p][child] = Fraction(-1)
        cartan[child][p] = -c

    # A fundamental cycle fixes each remaining directed entry.
    for u, v in packet.edges:
        if (u, v) in tree:
            continue
        path = _tree_path(rank, tree, v, u)
        cycle = (u, *path[:-1])
        tree_product = Fraction(1)
        for position in range(len(path) - 1):
            tree_product *= cartan[path[position]][path[position + 1]]
        cycle_value = values[canonical_cycle(cycle)]
        cartan[u][v] = cycle_value / tree_product
        cartan[v][u] = edge_products[(u, v)] / cartan[u][v]

    result = tuple(tuple(row) for row in cartan)
    if any(result[i][j] >= 0 for i, j in packet.edges):
        raise PacketError("the packet does not reconstruct an acute Cartan matrix")
    rebuilt = cartan_cycle_packet(result)
    if rebuilt != packet:
        raise PacketError("the complete cycle packet fails reconstruction compatibility")
    symmetrizer(result)  # Also checks forward/reverse compatibility on all cycles.
    return result


def diagonal_rescale(cartan: Matrix, scales: Sequence[Fraction | int]) -> Matrix:
    rank = len(cartan)
    values = tuple(Fraction(value) for value in scales)
    if len(values) != rank or any(not value for value in values):
        raise ValueError("diagonal gauge needs one nonzero scale per vertex")
    return tuple(
        tuple(values[j] * cartan[i][j] / values[i] for j in range(rank))
        for i in range(rank)
    )


def symmetrizer(cartan: Matrix) -> tuple[Fraction, ...]:
    """Return positive d with d_i a_ij=d_j a_ji, normalized by d_0=1."""

    rank = len(cartan)
    if any(len(row) != rank for row in cartan):
        raise PacketError("Cartan matrix is not square")
    if any(cartan[i][i] != 2 for i in range(rank)):
        raise PacketError("Cartan diagonal must be 2")
    adjacency = [set() for _ in range(rank)]
    for i in range(rank):
        for j in range(i + 1, rank):
            if bool(cartan[i][j]) != bool(cartan[j][i]):
                raise PacketError("Cartan zero pattern is not symmetric")
            if cartan[i][j]:
                if cartan[i][j] >= 0 or cartan[j][i] >= 0:
                    raise PacketError("off-diagonal Cartan entries must be negative")
                adjacency[i].add(j)
                adjacency[j].add(i)
    values: list[Fraction | None] = [None] * rank
    values[0] = Fraction(1)
    queue = [0]
    for i in queue:
        assert values[i] is not None
        for j in sorted(adjacency[i]):
            proposed = values[i] * cartan[i][j] / cartan[j][i]
            if proposed <= 0:
                raise PacketError("Cartan matrix has no positive symmetrizer")
            if values[j] is None:
                values[j] = proposed
                queue.append(j)
            elif values[j] != proposed:
                raise PacketError("forward and reverse cycle products disagree")
    if any(value is None for value in values):
        raise PacketError("Cartan graph is disconnected")
    return tuple(value for value in values if value is not None)


def cartan_gram(cartan: Matrix) -> Matrix:
    diagonal = symmetrizer(cartan)
    gram = tuple(
        tuple(diagonal[i] * cartan[i][j] for j in range(len(cartan)))
        for i in range(len(cartan))
    )
    if any(gram[i][j] != gram[j][i] for i in range(len(gram)) for j in range(i)):
        raise PacketError("internal error: the symmetrized Cartan matrix is not symmetric")
    return gram


def _prime_factors(value: int) -> set[int]:
    value = abs(value)
    factors: set[int] = set()
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            factors.add(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if value > 1:
        factors.add(value)
    return factors


def _valuation_integer(value: int, prime: int) -> int:
    value = abs(value)
    exponent = 0
    while value and value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def _valuation(value: Fraction, prime: int) -> int:
    return _valuation_integer(value.numerator, prime) - _valuation_integer(
        value.denominator, prime
    )


def integral_cartan_gauge(cartan: Matrix) -> tuple[Matrix, tuple[Fraction, ...]]:
    """Find a diagonal Q-gauge in which every Cartan entry is integral.

    For each prime p, integrality asks for potentials x_i=v_p(lambda_i)

        x_j - x_i >= -v_p(a_ij).

    The obstruction on a directed cycle is precisely nonintegrality of its
    closed Cartan product.  Bellman--Ford-style longest-path relaxation
    therefore gives a constructive version of the cyclic integrality lemma.
    """

    rank = len(cartan)
    packet = cartan_cycle_packet(cartan)
    if any(value.denominator != 1 for _, value in packet.values):
        raise PacketError("nonintegral closed cycle obstructs an integral Cartan gauge")
    directed_edges = [
        (i, j)
        for i in range(rank)
        for j in range(rank)
        if i != j and cartan[i][j]
    ]
    primes: set[int] = set()
    for i, j in directed_edges:
        value = cartan[i][j]
        primes |= _prime_factors(value.numerator)
        primes |= _prime_factors(value.denominator)

    exponents = [[0 for _ in primes] for _ in range(rank)]
    ordered_primes = sorted(primes)
    for prime_index, prime in enumerate(ordered_primes):
        potential = [0] * rank
        for pass_index in range(rank):
            changed = False
            for i, j in directed_edges:
                required = potential[i] - _valuation(cartan[i][j], prime)
                if potential[j] < required:
                    potential[j] = required
                    changed = True
            if not changed:
                break
            if pass_index == rank - 1:
                raise PacketError(
                    f"positive p-adic constraint cycle at p={prime}; packet is not integral"
                )
        shift = min(potential)
        for vertex in range(rank):
            exponents[vertex][prime_index] = potential[vertex] - shift

    scales = []
    for vertex in range(rank):
        scale = 1
        for prime_index, prime in enumerate(ordered_primes):
            scale *= prime ** exponents[vertex][prime_index]
        scales.append(Fraction(scale))
    gauged = diagonal_rescale(cartan, scales)
    if any(entry.denominator != 1 for row in gauged for entry in row):
        raise PacketError("internal error: valuation gauge did not integralize the Cartan matrix")
    return gauged, tuple(scales)


def _submatrix(matrix: Matrix, omitted: int) -> Matrix:
    keep = [index for index in range(len(matrix)) if index != omitted]
    return tuple(tuple(matrix[i][j] for j in keep) for i in keep)


@dataclass(frozen=True, slots=True)
class SimplexAnalysis:
    full_inertia: Inertia
    vertex_inertias: tuple[Inertia, ...]
    vertex_kinds: tuple[str, ...]
    lorentzian: bool
    finite_volume: bool
    compact: bool

    def to_json(self) -> dict[str, object]:
        def inertia_json(value: Inertia) -> list[int]:
            return [value.positive, value.negative, value.zero]

        return {
            "full_inertia": inertia_json(self.full_inertia),
            "vertex_inertias": [inertia_json(value) for value in self.vertex_inertias],
            "vertex_kinds": list(self.vertex_kinds),
            "lorentzian": self.lorentzian,
            "finite_volume": self.finite_volume,
            "compact": self.compact,
        }


def analyze_simplex_gram(gram: Matrix) -> SimplexAnalysis:
    rank = len(gram)
    full = inertia(gram)
    expected = Inertia(rank - 1, 1, 0)
    vertex_inertias = tuple(inertia(_submatrix(gram, omitted)) for omitted in range(rank))
    kinds = []
    for value in vertex_inertias:
        if value == Inertia(rank - 1, 0, 0):
            kinds.append("finite")
        elif value == Inertia(rank - 2, 0, 1):
            kinds.append("ideal")
        else:
            kinds.append("hyperideal")
    lorentzian = full == expected
    finite_volume = lorentzian and all(kind != "hyperideal" for kind in kinds)
    return SimplexAnalysis(
        full,
        vertex_inertias,
        tuple(kinds),
        lorentzian,
        finite_volume,
        finite_volume and all(kind == "finite" for kind in kinds),
    )


@dataclass(frozen=True, slots=True)
class CartanCandidate:
    packet: CyclePacket
    cartan: Matrix
    integral_cartan: Matrix
    integral_gauge: tuple[Fraction, ...]
    gram: Matrix
    analysis: SimplexAnalysis

    @property
    def candidate_id(self) -> str:
        encoded = json.dumps(self.packet.to_json(), sort_keys=True, separators=(",", ":"))
        return sha256(encoded.encode("utf-8")).hexdigest()[:16]

    def to_json(self, *, include_matrices: bool = True) -> dict[str, object]:
        result: dict[str, object] = {
            "candidate_id": self.candidate_id,
            "packet": self.packet.to_json(),
            "analysis": self.analysis.to_json(),
        }
        if include_matrices:
            matrix_json = lambda matrix: [
                [_fraction_string(value) for value in row] for row in matrix
            ]
            result.update(
                {
                    "cartan_tree_gauge": matrix_json(self.cartan),
                    "integral_gauge_scales": [
                        _fraction_string(value) for value in self.integral_gauge
                    ],
                    "integral_cartan": matrix_json(self.integral_cartan),
                    "symmetrized_gram": matrix_json(self.gram),
                }
            )
        return result


def reconstruct_candidate(packet: CyclePacket) -> CartanCandidate:
    cartan = reconstruct_cartan(packet)
    integral_cartan, gauge = integral_cartan_gauge(cartan)
    gram = cartan_gram(integral_cartan)
    return CartanCandidate(
        packet,
        cartan,
        integral_cartan,
        gauge,
        gram,
        analyze_simplex_gram(gram),
    )


RATIONAL_COXETER_LABELS: dict[Fraction, int | str] = {
    Fraction(0): 2,
    Fraction(1): 3,
    Fraction(2): 4,
    Fraction(3): 6,
    Fraction(4): "parallel",
}


def _solve(matrix: Matrix, right: Sequence[Fraction]) -> tuple[Fraction, ...]:
    size = len(matrix)
    if len(right) != size:
        raise ValueError("linear-system dimensions disagree")
    work = [list(matrix[i]) + [Fraction(right[i])] for i in range(size)]
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            raise PacketError("candidate Gram matrix is singular")
        work[column], work[pivot] = work[pivot], work[column]
        pivot_value = work[column][column]
        work[column] = [value / pivot_value for value in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            multiple = work[row][column]
            work[row] = [
                work[row][index] - multiple * work[column][index]
                for index in range(size + 1)
            ]
    return tuple(work[i][-1] for i in range(size))


def _golden_json(value: Fraction | int) -> dict[str, str]:
    return Golden(Fraction(value)).to_json()


def build_q_simplex_certificate(candidate: CartanCandidate) -> dict[str, object]:
    """Serialize a finite-volume rational simplex for independent replay."""

    if not candidate.analysis.finite_volume:
        raise PacketError("only finite-volume simplex candidates have chamber certificates")
    rank = candidate.packet.rank
    if rank < 3:
        raise PacketError("the chamber verifier starts in hyperbolic dimension two")
    edge_products = candidate.packet.edge_products
    relations: list[list[int | str]] = [[1 for _ in range(rank)] for _ in range(rank)]
    for i in range(rank):
        for j in range(i):
            value = edge_products.get((j, i), Fraction(0))
            if value not in RATIONAL_COXETER_LABELS:
                raise PacketError(
                    f"edge product {value} is not a recognized rational Coxeter label"
                )
            relations[i][j] = relations[j][i] = RATIONAL_COXETER_LABELS[value]

    identity_roots = [
        [_golden_json(int(i == j)) for j in range(rank)] for i in range(rank)
    ]
    gram = [[_golden_json(value) for value in row] for row in candidate.gram]
    interior = _solve(candidate.gram, [Fraction(-1)] * rank)
    vertices = []
    cusps = []
    for omitted, kind in enumerate(candidate.analysis.vertex_kinds):
        active = [index for index in range(rank) if index != omitted]
        if kind == "finite":
            vertices.append(active)
        elif kind == "ideal":
            cusps.append(active)

    return {
        "schema_version": "0.1.0",
        "certificate_kind": (
            "compact_reflective" if candidate.analysis.compact else "cofinite_noncompact"
        ),
        "field": "Q",
        "dimension": rank - 1,
        "lattice_gram": gram,
        "simple_roots": identity_roots,
        "wall_gram": gram,
        "coxeter_relations": relations,
        "completeness": {
            "method": "finite_volume_coxeter_simplex_principal_subgram_criterion",
            "evidence": {
                "candidate_id": candidate.candidate_id,
                "integral_closed_cycle_packet": True,
                "integral_cartan_gauge": True,
                "all_vertex_links_spherical_or_affine": True,
            },
        },
        "geometry": {
            "nonempty": True,
            "interior_point": [_golden_json(value) for value in interior],
            "vertices": vertices,
            "cusps": cusps,
        },
    }


def verify_q_simplex_candidate(candidate: CartanCandidate) -> ChamberVerification:
    return verify_chamber_geometry(build_q_simplex_certificate(candidate))


@dataclass(frozen=True, slots=True)
class Census:
    rank: int
    edge_products: tuple[Fraction, ...]
    assignment_count: int
    connected_assignment_count: int
    integral_q_packet_assignment_count: int
    unlabelled_packet_count: int
    lorentzian_packet_count: int
    finite_volume_candidates: tuple[CartanCandidate, ...]
    independently_verified_count: int

    @property
    def compact_count(self) -> int:
        return sum(candidate.analysis.compact for candidate in self.finite_volume_candidates)

    @property
    def noncompact_count(self) -> int:
        return len(self.finite_volume_candidates) - self.compact_count

    def to_json(self, *, include_matrices: bool = True) -> dict[str, object]:
        return {
            "schema": "rational_cartan_cycle_census/v1",
            "scope": {
                "field": "Q",
                "rank": self.rank,
                "dimension": self.rank - 1,
                "allowed_nonzero_edge_products": [
                    _fraction_string(value) for value in self.edge_products
                ],
                "connected_frames_only": True,
                "acute_symmetrizable_cycle_packets": True,
                "quotient": "facet permutations and diagonal root gauge",
                "chamber_type": "finite-volume Coxeter simplices",
            },
            "counts": {
                "labelled_pair_assignments": self.assignment_count,
                "labelled_connected_assignments": self.connected_assignment_count,
                "labelled_integral_Q_cycle_packets": self.integral_q_packet_assignment_count,
                "unlabelled_integral_Q_cycle_packets": self.unlabelled_packet_count,
                "unlabelled_lorentzian_packets": self.lorentzian_packet_count,
                "finite_volume_simplex_candidates": len(self.finite_volume_candidates),
                "compact": self.compact_count,
                "noncompact": self.noncompact_count,
                "independent_exact_geometry_replays": self.independently_verified_count,
            },
            "candidates": [
                candidate.to_json(include_matrices=include_matrices)
                for candidate in self.finite_volume_candidates
            ],
        }


def enumerate_q_simplex_candidates(
    rank: int,
    edge_products: Iterable[Fraction | int] = (1, 2, 3, 4),
    *,
    max_assignments: int = 2_000_000,
    independent_replay: bool = True,
) -> Census:
    """Exhaust a finite rational cycle box and return all simplex candidates.

    Zero (orthogonality, Coxeter label 2) is automatically included for each
    pair.  The supplied values are the allowed *nonzero* length-two Cartan
    products.  Results are quotiented exactly by all facet permutations.
    """

    if rank < 3:
        raise ValueError("simplex enumeration requires rank at least three")
    if rank > 8:
        raise ValueError("the exact factorial isomorphism quotient is limited to rank eight")
    allowed = tuple(sorted(set(Fraction(value) for value in edge_products)))
    if not allowed or any(value <= 0 for value in allowed):
        raise ValueError("edge products must be a nonempty set of positive rationals")
    pairs = tuple(combinations(range(rank), 2))
    assignment_count = (len(allowed) + 1) ** len(pairs)
    if assignment_count > max_assignments:
        raise ValueError(
            f"search has {assignment_count} labelled assignments, above cap {max_assignments}; "
            "partition the edge box or raise max_assignments explicitly"
        )

    connected_count = 0
    integral_packet_count = 0
    representatives: dict[tuple[tuple[Cycle, Fraction], ...], CyclePacket] = {}
    for values in product((Fraction(0), *allowed), repeat=len(pairs)):
        assignment = {
            pair: value for pair, value in zip(pairs, values, strict=True) if value
        }
        if not _connected(rank, assignment):
            continue
        connected_count += 1
        try:
            packet = CyclePacket.from_edge_products(rank, assignment)
        except PacketError:
            continue
        integral_packet_count += 1
        canonical = packet.canonicalized()
        representatives.setdefault(canonical.key(), canonical)

    candidates = [reconstruct_candidate(packet) for packet in representatives.values()]
    candidates.sort(key=lambda candidate: candidate.packet.key())
    lorentzian_count = sum(candidate.analysis.lorentzian for candidate in candidates)
    finite_volume = tuple(
        candidate for candidate in candidates if candidate.analysis.finite_volume
    )
    verified = 0
    if independent_replay:
        for candidate in finite_volume:
            result = verify_q_simplex_candidate(candidate)
            if result.dimension != rank - 1 or result.facet_count != rank:
                raise PacketError("independent chamber replay returned inconsistent dimensions")
            verified += 1

    return Census(
        rank,
        allowed,
        assignment_count,
        connected_count,
        integral_packet_count,
        len(representatives),
        lorentzian_count,
        finite_volume,
        verified,
    )

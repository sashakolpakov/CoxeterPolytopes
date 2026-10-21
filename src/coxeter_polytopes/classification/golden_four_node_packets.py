"""Exact four-node golden AGCat packets with at most one dotted edge.

This module moves the golden decision queue beyond isolated triangles.  It
enumerates, modulo the full four-vertex symmetry, every connected four-node
packet whose edges are all finite/orthogonal, and every connected packet
with exactly one dotted edge.  In the latter tranche the dotted value ranges
over the complete 6,869-element golden frame alphabet; the other five slots
range over all six finite/orthogonal Coxeter products.

Cycle squareclasses are the exact symmetrizability gate.  For a connected
edge-labelled graph, the product on every graph cycle must be a square in
``Q(sqrt(5))``.  Conversely those equations let one propagate a positive
symmetrizer along a spanning tree, and every chord equation then has a
golden square root.  Since the edge products are integral and the golden
integers are integrally closed, all resulting directed cycle products are
integral.

The normalized four-node Gram determinant is evaluated from pair products
and the compatible triangle/four-cycle roots.  At the nonselected real
place, positive definiteness is checked by exact Sylvester signs.  At the
selected place, the finite-only tranche retains positive-definite or
Lorentzian restrictions; the one-dotted tranche must be Lorentzian.  This
is a complete four-node *local* tranche, not a rank-23--30 graph census or a
high-dimensional nonexistence theorem.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
from itertools import permutations, product

from coxeter_polytopes.classification.golden_agcat_manifest import (
    golden_coxeter_frame_pair_alphabet,
    golden_finite_pair_values,
)
from coxeter_polytopes.classification.golden_triangle_packets import (
    GoldenElementSquareclass,
    _pair_mul,
    _pair_sign,
    _pair_square_root,
    _squareclass_xor,
    golden_element_squareclass,
)


Pair = tuple[int, int]

# Pair slots are kept in this order throughout the certificate.
EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
EDGE_INDEX = {edge: index for index, edge in enumerate(EDGES)}
TRIANGLES = (
    (0, 1, 3),  # 0-1-2
    (0, 2, 4),  # 0-1-3
    (1, 2, 5),  # 0-2-3
    (3, 4, 5),  # 1-2-3
)
FOUR_CYCLES = (
    (0, 3, 5, 2),  # 0-1-2-3-0
    (0, 4, 5, 1),  # 0-1-3-2-0
    (1, 3, 4, 2),  # 0-2-1-3-0
)
ALL_CYCLES = TRIANGLES + FOUR_CYCLES
OPPOSITE_EDGE_PAIRS = ((0, 5), (1, 4), (2, 3))


def _pair_add(*values: Pair) -> Pair:
    return sum(value[0] for value in values), sum(value[1] for value in values)


def _pair_scale(multiplier: int, value: Pair) -> Pair:
    return multiplier * value[0], multiplier * value[1]


def _all_nonzero(values: tuple[Pair, ...]) -> bool:
    return all(value != (0, 0) for value in values)


@lru_cache(maxsize=None)
def _cycle_root(sorted_edges: tuple[Pair, ...]) -> Pair:
    cycle_product = (1, 0)
    for edge in sorted_edges:
        cycle_product = _pair_mul(cycle_product, edge)
    return _pair_square_root(cycle_product)


def golden_triangle_determinant(edges: tuple[Pair, ...], slots: tuple[int, int, int]) -> Pair:
    """Return the normalized determinant on one compatible vertex triple."""

    triangle_edges = tuple(edges[slot] for slot in slots)
    determinant = _pair_add((8, 0), _pair_scale(-2, _pair_add(*triangle_edges)))
    if _all_nonzero(triangle_edges):
        # An acute odd cycle is the negative identity-positive square root.
        determinant = _pair_add(
            determinant,
            _pair_scale(-2, _cycle_root(tuple(sorted(triangle_edges)))),
        )
    return determinant


def golden_four_node_determinant(edges: tuple[Pair, ...]) -> Pair:
    """Return the exact normalized four-node Gram determinant.

    Expanding the determinant by permutation cycle type gives the identity,
    six transpositions, three disjoint-transposition pairs, four triangles,
    and three unoriented Hamilton cycles.  In an acute realization triangle
    roots enter negatively and four-cycle roots positively.
    """

    if len(edges) != 6:
        raise ValueError("a four-node complete pair-slot record has six entries")
    determinant = _pair_add((16, 0), _pair_scale(-4, _pair_add(*edges)))
    for first, second in OPPOSITE_EDGE_PAIRS:
        determinant = _pair_add(determinant, _pair_mul(edges[first], edges[second]))
    for cycle in TRIANGLES:
        cycle_edges = tuple(edges[slot] for slot in cycle)
        if _all_nonzero(cycle_edges):
            determinant = _pair_add(
                determinant,
                _pair_scale(-4, _cycle_root(tuple(sorted(cycle_edges)))),
            )
    for cycle in FOUR_CYCLES:
        cycle_edges = tuple(edges[slot] for slot in cycle)
        if _all_nonzero(cycle_edges):
            determinant = _pair_add(
                determinant,
                _pair_scale(-2, _cycle_root(tuple(sorted(cycle_edges)))),
            )
    return determinant


def _edge_slot(first: int, second: int) -> int:
    return EDGE_INDEX[tuple(sorted((first, second)))]


def _permute_pattern(pattern: tuple[int, ...], permutation: tuple[int, ...]) -> tuple[int, ...]:
    output = [0] * 6
    for slot, (first, second) in enumerate(EDGES):
        output[_edge_slot(permutation[first], permutation[second])] = pattern[slot]
    return tuple(output)


def _connected(pattern: tuple[int, ...]) -> bool:
    reached = {0}
    for _ in range(4):
        for slot, (first, second) in enumerate(EDGES):
            if pattern[slot] and (first in reached or second in reached):
                reached.update((first, second))
    return len(reached) == 4


@lru_cache(maxsize=1)
def finite_only_canonical_patterns() -> tuple[tuple[int, ...], ...]:
    """Connected finite/orthogonal patterns modulo all of ``S_4``."""

    vertex_permutations = tuple(permutations(range(4)))
    output = []
    for pattern in product(range(6), repeat=6):
        if not _connected(pattern):
            continue
        if pattern == min(_permute_pattern(pattern, p) for p in vertex_permutations):
            output.append(pattern)
    return tuple(output)


@lru_cache(maxsize=1)
def one_dotted_canonical_patterns() -> tuple[tuple[int, ...], ...]:
    """Finite complements modulo the stabilizer of the unique dotted edge.

    The marker ``6`` occupies slot ``01``.  Since the dotted edge is unique,
    every ``S_4`` orbit has a representative in this slice, and its residual
    symmetry is exactly the endpoint/complement group ``S_2 x S_2``.
    """

    stabilizer = (
        (0, 1, 2, 3),
        (1, 0, 2, 3),
        (0, 1, 3, 2),
        (1, 0, 3, 2),
    )
    output = []
    for complement in product(range(6), repeat=5):
        pattern = (6,) + complement
        if not _connected(pattern):
            continue
        if pattern == min(_permute_pattern(pattern, p) for p in stabilizer):
            output.append(pattern)
    return tuple(output)


def _xor_all(classes: tuple[GoldenElementSquareclass, ...]) -> GoldenElementSquareclass:
    output = classes[0]
    for squareclass in classes[1:]:
        output = _squareclass_xor(output, squareclass)
    return output


def _topology(pattern: tuple[int, ...]) -> tuple[str, int, tuple[int, ...], int]:
    degrees = [0, 0, 0, 0]
    edge_count = 0
    for slot, (first, second) in enumerate(EDGES):
        if pattern[slot]:
            edge_count += 1
            degrees[first] += 1
            degrees[second] += 1
    degree_sequence = tuple(sorted(degrees))
    names = {
        (3, (1, 1, 1, 3)): "star_K1_3",
        (3, (1, 1, 2, 2)): "path_P4",
        (4, (1, 2, 2, 3)): "triangle_with_pendant",
        (4, (2, 2, 2, 2)): "cycle_C4",
        (5, (2, 2, 3, 3)): "diamond_K4_minus_edge",
        (6, (3, 3, 3, 3)): "complete_K4",
    }
    return names[(edge_count, degree_sequence)], edge_count, degree_sequence, edge_count - 3


@dataclass(frozen=True, slots=True)
class GoldenFourNodeStratum:
    tranche: str
    topology: str
    nonzero_edge_count: int
    degree_sequence: tuple[int, ...]
    cycle_rank: int
    selected_leading_triple_sign: int
    selected_full_determinant_sign: int
    conjugate_leading_triple_sign: int
    conjugate_full_determinant_sign: int
    count: int


@dataclass(frozen=True, slots=True)
class GoldenFourNodePacketCensus:
    finite_value_count: int
    dotted_value_count: int
    finite_only_connected_orbit_count: int
    one_dotted_complement_orbit_count: int
    one_dotted_raw_orbit_count: int
    finite_only_cycle_compatible_count: int
    one_dotted_cycle_compatible_count: int
    finite_only_conjugate_positive_count: int
    one_dotted_conjugate_positive_count: int
    finite_only_retained_count: int
    finite_only_retained_selected_positive_count: int
    finite_only_retained_selected_lorentzian_count: int
    one_dotted_retained_count: int
    one_dotted_retained_selected_lorentzian_count: int
    one_dotted_selected_nonlorentzian_count: int
    finite_only_selected_leading_triple_degenerate_count: int
    finite_only_selected_full_degenerate_count: int
    one_dotted_selected_leading_triple_degenerate_count: int
    one_dotted_selected_full_degenerate_count: int
    finite_only_conjugate_leading_triple_degenerate_count: int
    finite_only_conjugate_full_degenerate_count: int
    one_dotted_conjugate_leading_triple_degenerate_count: int
    one_dotted_conjugate_full_degenerate_count: int
    strata: tuple[GoldenFourNodeStratum, ...]
    finite_only_compatible_stream_sha256: str
    finite_only_retained_stream_sha256: str
    one_dotted_compatible_stream_sha256: str
    one_dotted_retained_stream_sha256: str


@dataclass(frozen=True, slots=True)
class GoldenOneDottedPathRankPacket:
    native_rank: int
    hyperbolic_dimension: int
    admissible_dotted_value_count: int
    admissible_dotted_stream_sha256: str
    sample_dotted_product: Pair
    sample_full_determinant: Pair


def golden_path_leading_determinant(rank: int, dotted_product: Pair) -> Pair:
    """Determinant of the rank-``rank`` one-dotted label-3 path.

    The first edge has product ``c`` and every later edge has product one.
    The tridiagonal recurrence ``D_j=2D_{j-1}-D_{j-2}``, starting with
    ``D_1=2`` and ``D_2=4-c``, gives ``D_j=2j-(j-1)c``.
    """

    if rank < 1:
        raise ValueError("path rank must be positive")
    return 2 * rank - (rank - 1) * dotted_product[0], -(rank - 1) * dotted_product[1]


@lru_cache(maxsize=1)
def golden_one_dotted_path_rank_packets() -> tuple[GoldenOneDottedPathRankPacket, ...]:
    """Exhaust the one-dotted, otherwise label-3 path family in ranks 23--30.

    At the selected place every dotted ``c>4`` makes ``D_j<0`` for all
    ``j>=2``, hence the path has exactly one negative direction.  At the
    conjugate place Sylvester positivity through rank ``r`` is equivalent to
    ``c' < 2r/(r-1)`` because these thresholds decrease with ``j``.
    """

    dotted_values = tuple(
        value
        for value in golden_coxeter_frame_pair_alphabet()
        if (value - 4).sign(0) > 0
    )
    sample = (2, 2)  # 2+2*phi; its conjugate is 4-2*phi in (0,1).
    output = []
    for rank in range(23, 31):
        admissible = tuple(
            value
            for value in dotted_values
            if ((rank - 1) * value - 2 * rank).sign(1) < 0
        )
        digest = sha256()
        for value in admissible:
            digest.update(f"{int(value.a)},{int(value.b)}\n".encode())
        if sample not in tuple((int(value.a), int(value.b)) for value in admissible):
            raise AssertionError("the uniform golden path witness left the admissible set")
        output.append(
            GoldenOneDottedPathRankPacket(
                native_rank=rank,
                hyperbolic_dimension=rank - 1,
                admissible_dotted_value_count=len(admissible),
                admissible_dotted_stream_sha256=digest.hexdigest(),
                sample_dotted_product=sample,
                sample_full_determinant=golden_path_leading_determinant(rank, sample),
            )
        )
    return tuple(output)


def _record(edges: tuple[Pair, ...], leading: Pair, determinant: Pair) -> bytes:
    return (
        ";".join(f"{value[0]},{value[1]}" for value in (*edges, leading, determinant))
        + "\n"
    ).encode()


@lru_cache(maxsize=1)
def golden_four_node_packet_census() -> GoldenFourNodePacketCensus:
    """Enumerate the exact connected four-node at-most-one-dotted tranche."""

    finite_values = tuple(row.product for row in golden_finite_pair_values())
    finite_pairs = tuple((int(value.a), int(value.b)) for value in finite_values)
    finite_classes = tuple(
        None if not value else golden_element_squareclass(value)
        for value in finite_values
    )
    trivial_class = golden_element_squareclass(finite_values[1] * finite_values[1])

    dotted_values = tuple(
        value
        for value in golden_coxeter_frame_pair_alphabet()
        if (value - 4).sign(0) > 0
    )
    dotted_by_class: dict[GoldenElementSquareclass, list[Pair]] = defaultdict(list)
    for value in dotted_values:
        dotted_by_class[golden_element_squareclass(value)].append(
            (int(value.a), int(value.b))
        )
    for values in dotted_by_class.values():
        values.sort()
    all_dotted_pairs = tuple(
        sorted((int(value.a), int(value.b)) for value in dotted_values)
    )

    strata: Counter[tuple[object, ...]] = Counter()
    finite_compatible_digest = sha256()
    finite_retained_digest = sha256()
    dotted_compatible_digest = sha256()
    dotted_retained_digest = sha256()
    counts: Counter[str] = Counter()

    # Every finite-only cycle must have trivial element squareclass.
    for pattern in finite_only_canonical_patterns():
        compatible = True
        for cycle in ALL_CYCLES:
            if all(pattern[slot] for slot in cycle):
                cycle_class = _xor_all(
                    tuple(finite_classes[pattern[slot]] for slot in cycle)  # type: ignore[misc]
                )
                if cycle_class != trivial_class:
                    compatible = False
                    break
        if not compatible:
            continue

        edges = tuple(finite_pairs[index] for index in pattern)
        leading = golden_triangle_determinant(edges, TRIANGLES[0])
        determinant = golden_four_node_determinant(edges)
        signs = (
            _pair_sign(leading, 0),
            _pair_sign(determinant, 0),
            _pair_sign(leading, 1),
            _pair_sign(determinant, 1),
        )
        selected_leading, selected_full, conjugate_leading, conjugate_full = signs
        conjugate_positive = conjugate_leading > 0 and conjugate_full > 0
        selected_positive = selected_leading > 0 and selected_full > 0
        # All selected two-node minors are positive in this tranche.  Hence a
        # negative full determinant means exactly one negative eigenvalue;
        # positive leading 3x3 and 4x4 determinants give positive definiteness.
        selected_lorentzian = selected_full < 0
        retained = conjugate_positive and (selected_positive or selected_lorentzian)

        counts["finite_compatible"] += 1
        counts["finite_conjugate_positive"] += conjugate_positive
        counts["finite_retained"] += retained
        counts["finite_retained_positive"] += retained and selected_positive
        counts["finite_retained_lorentzian"] += retained and selected_lorentzian
        counts["finite_selected_leading_degenerate"] += selected_leading == 0
        counts["finite_selected_full_degenerate"] += selected_full == 0
        counts["finite_conjugate_leading_degenerate"] += conjugate_leading == 0
        counts["finite_conjugate_full_degenerate"] += conjugate_full == 0

        topology, edge_count, degrees, cycle_rank = _topology(pattern)
        strata[("finite_only", topology, edge_count, degrees, cycle_rank, *signs)] += 1
        encoded = _record(edges, leading, determinant)
        finite_compatible_digest.update(encoded)
        if retained:
            finite_retained_digest.update(encoded)

    # With one unique dotted edge fixed at 01, cycles containing 01 all
    # prescribe the same required dotted squareclass.  Cycles avoiding it
    # must already have trivial finite squareclass.  If 01 is a bridge there
    # is no requirement, so all 6,869 dotted values occur.
    for pattern in one_dotted_canonical_patterns():
        required_classes: list[GoldenElementSquareclass] = []
        compatible = True
        for cycle in ALL_CYCLES:
            if not all(pattern[slot] for slot in cycle):
                continue
            finite_cycle_classes = tuple(
                finite_classes[pattern[slot]]
                for slot in cycle
                if slot != 0
            )
            cycle_class = _xor_all(finite_cycle_classes)  # type: ignore[arg-type]
            if 0 in cycle:
                required_classes.append(cycle_class)
            elif cycle_class != trivial_class:
                compatible = False
                break
        if required_classes and any(
            squareclass != required_classes[0]
            for squareclass in required_classes[1:]
        ):
            compatible = False
        if not compatible:
            dotted_choices: tuple[Pair, ...] = ()
        elif required_classes:
            dotted_choices = tuple(dotted_by_class.get(required_classes[0], ()))
        else:
            dotted_choices = all_dotted_pairs

        topology, edge_count, degrees, cycle_rank = _topology(pattern)
        for dotted in dotted_choices:
            edges = (dotted,) + tuple(finite_pairs[index] for index in pattern[1:])
            leading = golden_triangle_determinant(edges, TRIANGLES[0])
            determinant = golden_four_node_determinant(edges)
            signs = (
                _pair_sign(leading, 0),
                _pair_sign(determinant, 0),
                _pair_sign(leading, 1),
                _pair_sign(determinant, 1),
            )
            selected_leading, selected_full, conjugate_leading, conjugate_full = signs
            conjugate_positive = conjugate_leading > 0 and conjugate_full > 0
            selected_lorentzian = selected_leading < 0 and selected_full < 0
            retained = conjugate_positive and selected_lorentzian

            counts["dotted_compatible"] += 1
            counts["dotted_conjugate_positive"] += conjugate_positive
            counts["dotted_retained"] += retained
            counts["dotted_retained_lorentzian"] += retained and selected_lorentzian
            counts["dotted_selected_nonlorentzian"] += not selected_lorentzian
            counts["dotted_selected_leading_degenerate"] += selected_leading == 0
            counts["dotted_selected_full_degenerate"] += selected_full == 0
            counts["dotted_conjugate_leading_degenerate"] += conjugate_leading == 0
            counts["dotted_conjugate_full_degenerate"] += conjugate_full == 0

            strata[("exactly_one_dotted", topology, edge_count, degrees, cycle_rank, *signs)] += 1
            encoded = _record(edges, leading, determinant)
            dotted_compatible_digest.update(encoded)
            if retained:
                dotted_retained_digest.update(encoded)

    strata_rows = tuple(
        GoldenFourNodeStratum(*key, count)
        for key, count in sorted(strata.items())
    )
    return GoldenFourNodePacketCensus(
        finite_value_count=len(finite_values),
        dotted_value_count=len(dotted_values),
        finite_only_connected_orbit_count=len(finite_only_canonical_patterns()),
        one_dotted_complement_orbit_count=len(one_dotted_canonical_patterns()),
        one_dotted_raw_orbit_count=len(one_dotted_canonical_patterns()) * len(dotted_values),
        finite_only_cycle_compatible_count=counts["finite_compatible"],
        one_dotted_cycle_compatible_count=counts["dotted_compatible"],
        finite_only_conjugate_positive_count=counts["finite_conjugate_positive"],
        one_dotted_conjugate_positive_count=counts["dotted_conjugate_positive"],
        finite_only_retained_count=counts["finite_retained"],
        finite_only_retained_selected_positive_count=counts["finite_retained_positive"],
        finite_only_retained_selected_lorentzian_count=counts["finite_retained_lorentzian"],
        one_dotted_retained_count=counts["dotted_retained"],
        one_dotted_retained_selected_lorentzian_count=counts["dotted_retained_lorentzian"],
        one_dotted_selected_nonlorentzian_count=counts["dotted_selected_nonlorentzian"],
        finite_only_selected_leading_triple_degenerate_count=counts[
            "finite_selected_leading_degenerate"
        ],
        finite_only_selected_full_degenerate_count=counts[
            "finite_selected_full_degenerate"
        ],
        one_dotted_selected_leading_triple_degenerate_count=counts[
            "dotted_selected_leading_degenerate"
        ],
        one_dotted_selected_full_degenerate_count=counts[
            "dotted_selected_full_degenerate"
        ],
        finite_only_conjugate_leading_triple_degenerate_count=counts[
            "finite_conjugate_leading_degenerate"
        ],
        finite_only_conjugate_full_degenerate_count=counts[
            "finite_conjugate_full_degenerate"
        ],
        one_dotted_conjugate_leading_triple_degenerate_count=counts[
            "dotted_conjugate_leading_degenerate"
        ],
        one_dotted_conjugate_full_degenerate_count=counts[
            "dotted_conjugate_full_degenerate"
        ],
        strata=strata_rows,
        finite_only_compatible_stream_sha256=finite_compatible_digest.hexdigest(),
        finite_only_retained_stream_sha256=finite_retained_digest.hexdigest(),
        one_dotted_compatible_stream_sha256=dotted_compatible_digest.hexdigest(),
        one_dotted_retained_stream_sha256=dotted_retained_digest.hexdigest(),
    )


__all__ = [name for name in globals() if not name.startswith("_")]

"""Exact abstract residual census for graph-shadow record 12.

The graph shadow is ``K_{2,2} + 10I``.  This module does not assert convex
or Coxeter realizability.  It reconstructs the selected facet complements
from the two balanced branch matrices and performs exact incidence audits.

The corrected exhaustive search never applies the invalid one-branch purity
condition.  It uses only the following safe one-sided consequence: a minimal
transversal of one branch that meets a branch special is automatically a
minimal transversal of the coupled clutter, so its order must be four or
five.  Isolate-only branch transversals are retained until the full coupled
blocker is formed.  The resulting 371 objects are abstract incidence and
mod-two homology-sphere candidates; PL, convex-polytopal, and Coxeter
realizability are separate questions.
"""

from __future__ import annotations

from collections import Counter, deque
from hashlib import sha256
from itertools import combinations
import json
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


TARGET_RECORD_INDEX = 12
TARGET_GRAPH6 = "M???????????o?o??"
TARGET_VERTEX_COUNT = 14
TARGET_DIMENSION = 9
FULL_MASK = (1 << TARGET_VERTEX_COUNT) - 1

ACTIVE_A = (0, 1)
ACTIVE_B = (12, 13)
ISOLATES = tuple(range(2, 12))
LOCAL_ISOLATES = tuple(range(10))
TARGET_EDGES = ((0, 12), (1, 12), (0, 13), (1, 13))
TARGET_EDGE_MASKS = tuple(sum(1 << vertex for vertex in edge) for edge in TARGET_EDGES)
BALANCED_CUT_MASK = sum(1 << vertex for vertex in range(5, 10))
MATRIX_ALL = (1 << 25) - 1

PAIR_INDEX = {
    pair: index for index, pair in enumerate(combinations(LOCAL_ISOLATES, 2))
}

EXPECTED_X_GRAPH_COUNT = 12_172
EXPECTED_NORMALIZED_CUT_STATE_COUNT = 56_310
EXPECTED_LOCAL_PAIR_FEASIBLE_STATE_COUNT = 1_391_382
EXPECTED_SAFE_SPECIAL_BRANCH_STATE_COUNT = 6_184
EXPECTED_COLOURED_BRANCH_ORBIT_COUNT = 1_308
EXPECTED_PROFILE_FACET_TYPE_PAIR_COUNT = 22_858
EXPECTED_RELATIVE_EMBEDDING_COUNT = 2_725_719
EXPECTED_FULL_BLOCKER_PASS_COUNT = 42_000
EXPECTED_SURVIVING_ORDERED_TYPE_COUNT = 100
EXPECTED_BALANCED_BRANCH_TYPE_COUNT = 10
EXPECTED_FULL_ORBIT_COUNT = 371
EXPECTED_SCHEDULE_ORBIT_COUNTS = {0: 48, 1: 95, 2: 144, 3: 63, 4: 21}
EXPECTED_UNFILTERED_LOCAL_STREAM_SHA256 = (
    "ed12edfbff40c2fab730ed3f1740bcf6e5272554f85a807468c7cb8e515961ca"
)
EXPECTED_SAFE_SPECIAL_STREAM_SHA256 = (
    "c91021fe69bc047c8f2eea5f43bf708d9025492c8e3cac3552205b857731e4ee"
)
EXPECTED_COLOURED_BRANCH_STREAM_SHA256 = (
    "640217770888c33fe60fc2d1ffd6dbb05c499c9307c40d4ec9642a482fca1905"
)
EXPECTED_COUPLED_PAIR_STREAM_SHA256 = (
    "a189faa3d1faeb1f4db21630d2fbb816f8d3d4f1421d04eb251112d37945f336"
)
EXPECTED_FULL_ORBIT_STREAM_SHA256 = (
    "4d978594b6be052328b4056d6b05c6249e80e8b156a09d2eba8057f3b7daa23c"
)


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def vertices(value: int, order: int = TARGET_VERTEX_COUNT) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(order) if value & (1 << vertex))


def adjacency_from_graph6(word: str) -> tuple[int, ...]:
    order, edges = decode_graph6(word)
    if order != 10:
        raise ValueError("record-12 residual X graph must have order ten")
    result = [0] * order
    for first, second in edges:
        result[first] |= 1 << second
        result[second] |= 1 << first
    return tuple(result)


def pair_mask_from_adjacency(adjacency: Sequence[int]) -> int:
    if len(adjacency) != 10:
        raise ValueError("record-12 branch adjacency must have order ten")
    return sum(
        1 << PAIR_INDEX[first, second]
        for first, second in combinations(LOCAL_ISOLATES, 2)
        if adjacency[first] & (1 << second)
    )


def matrix_to_pair_mask(matrix: int) -> int:
    if matrix < 0 or matrix & ~MATRIX_ALL:
        raise ValueError("record-12 branch matrix is not 5 by 5")
    return sum(
        1 << PAIR_INDEX[first, second + 5]
        for first in range(5)
        for second in range(5)
        if matrix & (1 << (5 * first + second))
    )


def pair_mask_to_matrix(value: int) -> int:
    result = 0
    for first, second in combinations(LOCAL_ISOLATES, 2):
        if not value & (1 << PAIR_INDEX[first, second]):
            continue
        if not (first < 5 <= second):
            raise ValueError("record-12 balanced graph contains a within-side edge")
        result |= 1 << (5 * first + second - 5)
    return result


def _pair_edge(value: int, first: int, second: int) -> bool:
    if first > second:
        first, second = second, first
    return bool(value & (1 << PAIR_INDEX[first, second]))


def _cut_edge(cut_mask: int, first: int, second: int) -> bool:
    return bool(cut_mask & (1 << first)) != bool(cut_mask & (1 << second))


def switched_triangle_free(x_mask: int, cut_mask: int) -> bool:
    return all(
        not all(
            _pair_edge(x_mask, first, second)
            != _cut_edge(cut_mask, first, second)
            for first, second in combinations(triple, 2)
        )
        for triple in combinations(LOCAL_ISOLATES, 3)
    )


def branch_triples(x_mask: int, cut_mask: int, y_mask: int = 0) -> tuple[int, ...]:
    """Return the branch triple clutter on ten isolates plus two specials."""

    result: set[int] = set()
    for isolate in LOCAL_ISOLATES:
        if y_mask & (1 << isolate):
            result.add(mask((10, 11, isolate)))
    for first, second in combinations(LOCAL_ISOLATES, 2):
        x = _pair_edge(x_mask, first, second)
        delta = _cut_edge(cut_mask, first, second)
        if x:
            result.add(mask((11, first, second)))
        if x != delta:
            result.add(mask((10, first, second)))
    for triple in combinations(LOCAL_ISOLATES, 3):
        if sum(_pair_edge(x_mask, *pair) for pair in combinations(triple, 2)) % 2:
            result.add(mask(triple))
    return tuple(sorted(result))


def minimal_transversals(edges: Iterable[int]) -> tuple[int, ...]:
    """Dualize a small finite clutter without a SAT/SMT backend."""

    blocker = [0]
    for edge in sorted(set(edges)):
        candidates = []
        for row in blocker:
            if row & edge:
                candidates.append(row)
                continue
            bits = edge
            while bits:
                bit = bits & -bits
                candidates.append(row | bit)
                bits ^= bit
        candidates = sorted(set(candidates), key=lambda row: (row.bit_count(), row))
        blocker = []
        for row in candidates:
            if not any(prior & row == prior for prior in blocker):
                blocker.append(row)
    return tuple(blocker)


def add_y_to_branch_blocker(
    base_blocker: Sequence[int], y_mask: int
) -> tuple[int, ...]:
    """Dualize the added triples ``{s,t,i}`` from the y support exactly."""

    specials = mask((10, 11))
    candidate_set = set()
    for row in base_blocker:
        if row & specials:
            candidate_set.add(row)
        else:
            candidate_set.update((row | y_mask, row | (1 << 10), row | (1 << 11)))
    candidates = sorted(candidate_set, key=lambda row: (row.bit_count(), row))
    result = []
    for row in candidates:
        if not any(prior & row == prior for prior in result):
            result.append(row)
    return tuple(result)


def allowed_z_values(
    x: int, delta_p: int, delta_q: int, y_first: int, y_second: int
) -> tuple[int, ...]:
    return tuple(
        z
        for z in (0, 1)
        if y_first
        + y_second
        + x
        + (x ^ delta_p)
        + z
        + (z ^ delta_q)
        in (0, 2)
    )


def branch_has_local_pair_extension(x_mask: int, p_mask: int, y_mask: int) -> bool:
    q_mask = p_mask ^ y_mask
    return all(
        allowed_z_values(
            int(_pair_edge(x_mask, first, second)),
            int(_cut_edge(p_mask, first, second)),
            int(_cut_edge(q_mask, first, second)),
            int(bool(y_mask & (1 << first))),
            int(bool(y_mask & (1 << second))),
        )
        for first, second in combinations(LOCAL_ISOLATES, 2)
    )


def selected_complements_from_pair_masks(
    x_mask: int,
    z_mask: int,
    p_mask: int = BALANCED_CUT_MASK,
    q_mask: int = BALANCED_CUT_MASK,
    y_mask: int = 0,
) -> tuple[int, ...]:
    p_triples = branch_triples(x_mask, p_mask, y_mask)
    q_triples = branch_triples(z_mask, q_mask, y_mask)
    result: set[int] = set()
    for triple in p_triples:
        global_triple = 0
        for local in LOCAL_ISOLATES:
            if triple & (1 << local):
                global_triple |= 1 << (local + 2)
        if triple & (1 << 10):
            global_triple |= 1 << 12
        if triple & (1 << 11):
            global_triple |= 1 << 13
        result.add(mask(ACTIVE_A) | global_triple)
    for triple in q_triples:
        global_triple = 0
        for local in LOCAL_ISOLATES:
            if triple & (1 << local):
                global_triple |= 1 << (local + 2)
        if triple & (1 << 10):
            global_triple |= 1 << 0
        if triple & (1 << 11):
            global_triple |= 1 << 1
        result.add(mask(ACTIVE_B) | global_triple)
    return tuple(sorted(result))


def selected_complements_from_matrices(
    x_matrix: int, z_matrix: int
) -> tuple[int, ...]:
    return selected_complements_from_pair_masks(
        matrix_to_pair_mask(x_matrix), matrix_to_pair_mask(z_matrix)
    )


def face_table_from_complements(complements: Sequence[int]) -> tuple[bool, ...]:
    table = [False] * (1 << TARGET_VERTEX_COUNT)
    for complement in complements:
        facet = FULL_MASK ^ complement
        subset = facet
        while True:
            table[subset] = True
            if subset == 0:
                break
            subset = (subset - 1) & facet
    return tuple(table)


def minimal_nonfaces_from_table(
    face_table: Sequence[bool],
) -> dict[int, tuple[int, ...]]:
    result: dict[int, list[int]] = {}
    for value in range(1, 1 << TARGET_VERTEX_COUNT):
        if face_table[value]:
            continue
        bits = value
        minimal = True
        while bits:
            bit = bits & -bits
            if not face_table[value ^ bit]:
                minimal = False
                break
            bits ^= bit
        if minimal:
            result.setdefault(value.bit_count(), []).append(value)
    return {size: tuple(rows) for size, rows in sorted(result.items())}


def _ridge_histogram(complements: Sequence[int]) -> Counter[int]:
    groups: Counter[int] = Counter()
    for complement in complements:
        outside = FULL_MASK ^ complement
        while outside:
            bit = outside & -outside
            groups[complement | bit] += 1
            outside ^= bit
    return Counter(groups.values())


def _dual_connected(complements: Sequence[int]) -> bool:
    values = set(complements)
    if not values:
        return False
    root = min(values)
    seen = {root}
    queue = deque([root])
    while queue:
        value = queue.popleft()
        outside = FULL_MASK ^ value
        while outside:
            added = outside & -outside
            inside = value
            while inside:
                removed = inside & -inside
                neighbor = (value | added) ^ removed
                if neighbor in values and neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
                inside ^= removed
            outside ^= added
    return len(seen) == len(values)


def incidence_audit(x_matrix: int, z_matrix: int) -> dict[str, object]:
    complements = selected_complements_from_matrices(x_matrix, z_matrix)
    table = face_table_from_complements(complements)
    nonfaces = minimal_nonfaces_from_table(table)
    histogram = {size: len(rows) for size, rows in nonfaces.items()}
    a = histogram.get(4, 0)
    face_counts = Counter(
        value.bit_count() for value, is_face in enumerate(table) if is_face
    )
    vector = tuple(face_counts[size] for size in range(1, TARGET_DIMENSION + 1))
    ridge = _ridge_histogram(complements)
    verified = (
        set(nonfaces) <= {2, 4, 5}
        and set(nonfaces.get(2, ())) == set(TARGET_EDGE_MASKS)
        and a in EXPECTED_SCHEDULE_ORBIT_COUNTS
        and len(complements) == 130 - 2 * a
        and ridge == Counter({2: vector[-2]})
        and _dual_connected(complements)
    )
    return {
        "missing_four_count": a,
        "facet_count": len(complements),
        "face_vector": vector,
        "minimal_nonface_histogram": histogram,
        "ridge_degree_histogram": dict(ridge),
        "dual_connected": _dual_connected(complements),
        "selected_complement_sha256": sha256(
            b"".join(value.to_bytes(2, "little") for value in complements)
        ).hexdigest(),
        "verified": verified,
    }


def canonical_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def verify_archive(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record12_exact_371_orbit_whole_shadow_residual":
        raise ValueError("wrong record-12 residual archive kind")
    if (
        document.get("record_index") != TARGET_RECORD_INDEX
        or document.get("graph6") != TARGET_GRAPH6
        or document.get("theorem_status")
        != "certified_exact_371_orbit_abstract_residual"
    ):
        raise ValueError("record-12 residual target changed")
    census = document.get("census")
    if not isinstance(census, Mapping):
        raise ValueError("record-12 residual census is absent")
    expected_counts = {
        "x_graph_count": EXPECTED_X_GRAPH_COUNT,
        "normalized_cut_states": EXPECTED_NORMALIZED_CUT_STATE_COUNT,
        "locally_pair_feasible_states": EXPECTED_LOCAL_PAIR_FEASIBLE_STATE_COUNT,
        "safe_special_branch_states": EXPECTED_SAFE_SPECIAL_BRANCH_STATE_COUNT,
        "coloured_branch_orbits": EXPECTED_COLOURED_BRANCH_ORBIT_COUNT,
        "profile_and_facet_compatible_ordered_type_pairs": (
            EXPECTED_PROFILE_FACET_TYPE_PAIR_COUNT
        ),
        "relative_embedded_graphs": EXPECTED_RELATIVE_EMBEDDING_COUNT,
        "full_blocker_passing_relative_embeddings": EXPECTED_FULL_BLOCKER_PASS_COUNT,
        "surviving_ordered_type_pairs": EXPECTED_SURVIVING_ORDERED_TYPE_COUNT,
        "surviving_balanced_branch_types": EXPECTED_BALANCED_BRANCH_TYPE_COUNT,
        "full_graph_shadow_orbits": EXPECTED_FULL_ORBIT_COUNT,
    }
    if any(census.get(key) != value for key, value in expected_counts.items()):
        raise ValueError("record-12 residual stage count changed")
    streams = document.get("generation_streams")
    if not isinstance(streams, Mapping):
        raise ValueError("record-12 generation streams are absent")
    expected_streams = {
        "unfiltered_local_ridge_states": (
            EXPECTED_LOCAL_PAIR_FEASIBLE_STATE_COUNT,
            EXPECTED_UNFILTERED_LOCAL_STREAM_SHA256,
        ),
        "safe_special_branch_states": (
            EXPECTED_SAFE_SPECIAL_BRANCH_STATE_COUNT,
            EXPECTED_SAFE_SPECIAL_STREAM_SHA256,
        ),
        "coloured_branch_orbits": (
            EXPECTED_COLOURED_BRANCH_ORBIT_COUNT,
            EXPECTED_COLOURED_BRANCH_STREAM_SHA256,
        ),
        "coupled_full_blocker_passes": (
            EXPECTED_FULL_BLOCKER_PASS_COUNT,
            EXPECTED_COUPLED_PAIR_STREAM_SHA256,
        ),
        "full_automorphism_orbits": (
            EXPECTED_FULL_ORBIT_COUNT,
            EXPECTED_FULL_ORBIT_STREAM_SHA256,
        ),
    }
    for name, (line_count, raw_digest) in expected_streams.items():
        row = streams.get(name)
        if (
            not isinstance(row, Mapping)
            or row.get("line_count") != line_count
            or row.get("raw_stream_sha256") != raw_digest
        ):
            raise ValueError(f"record-12 {name} stream changed")
    branches = document.get("coloured_branch_orbits")
    orbits = document.get("residual_orbits")
    if not isinstance(branches, list) or len(branches) != EXPECTED_COLOURED_BRANCH_ORBIT_COUNT:
        raise ValueError("record-12 coloured branch list changed")
    if not isinstance(orbits, list) or len(orbits) != EXPECTED_FULL_ORBIT_COUNT:
        raise ValueError("record-12 residual orbit list changed")
    raw_multiplicity_sum = 0
    for index, branch in enumerate(branches):
        if not isinstance(branch, Mapping) or branch.get("index") != index:
            raise ValueError("malformed record-12 coloured branch orbit")
        word = str(branch["x_graph6"])
        p_mask = int(branch["p_mask"])
        y_mask = int(branch["y_mask"])
        x_mask = pair_mask_from_adjacency(adjacency_from_graph6(word))
        if (
            not switched_triangle_free(x_mask, 0)
            or not switched_triangle_free(x_mask, p_mask)
            or not branch_has_local_pair_extension(x_mask, p_mask, y_mask)
        ):
            raise ValueError("archived record-12 branch violates a ridge equation")
        blocker = minimal_transversals(branch_triples(x_mask, p_mask, y_mask))
        specials = mask((10, 11))
        if any(
            row & specials and row.bit_count() not in (4, 5)
            for row in blocker
        ):
            raise ValueError("archived record-12 branch violates the safe filter")
        raw_multiplicity_sum += int(branch["raw_representation_count"])
    if raw_multiplicity_sum != EXPECTED_SAFE_SPECIAL_BRANCH_STATE_COUNT:
        raise ValueError("record-12 coloured branch multiplicities changed")

    observed_schedule: Counter[int] = Counter()
    representatives = set()
    orbit_multiplicity_sum = 0
    for orbit in orbits:
        if not isinstance(orbit, Mapping):
            raise ValueError("malformed record-12 residual orbit")
        x_matrix = int(orbit["x_matrix"])
        z_matrix = int(orbit["z_matrix"])
        key = (x_matrix, z_matrix)
        if key in representatives:
            raise ValueError("duplicate record-12 residual representative")
        representatives.add(key)
        audit = incidence_audit(x_matrix, z_matrix)
        if not audit["verified"]:
            raise AssertionError("record-12 residual incidence audit failed")
        if (
            orbit.get("missing_four_count") != audit["missing_four_count"]
            or orbit.get("facet_count") != audit["facet_count"]
            or orbit.get("selected_complement_sha256")
            != audit["selected_complement_sha256"]
        ):
            raise ValueError("record-12 residual orbit payload changed")
        observed_schedule[int(audit["missing_four_count"])] += 1
        orbit_multiplicity_sum += int(orbit["relative_embedding_multiplicity"])
    if orbit_multiplicity_sum != EXPECTED_FULL_BLOCKER_PASS_COUNT:
        raise ValueError("record-12 quotient multiplicities changed")
    if dict(sorted(observed_schedule.items())) != EXPECTED_SCHEDULE_ORBIT_COUNTS:
        raise ValueError("record-12 residual schedule histogram changed")
    if document.get("schedule_orbit_counts") != {
        str(key): value for key, value in EXPECTED_SCHEDULE_ORBIT_COUNTS.items()
    }:
        raise ValueError("record-12 archived schedule histogram changed")
    digest = canonical_digest(document)
    if document.get("certificate_sha256") != digest:
        raise ValueError("record-12 residual archive digest changed")
    return {
        "record_index": TARGET_RECORD_INDEX,
        "graph6": TARGET_GRAPH6,
        "coloured_branch_orbits": len(branches),
        "residual_orbits": len(orbits),
        "schedule_orbit_counts": dict(sorted(observed_schedule.items())),
        "relative_embedding_multiplicity_sum": orbit_multiplicity_sum,
        "certificate_sha256": digest,
    }

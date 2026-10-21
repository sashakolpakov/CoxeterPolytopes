"""Solver-independent exclusion certificate for H9/14 record 54.

The proof uses only finite integer enumeration.  A vertex-link ``g4`` double
count removes the schedules with two or three missing four-faces.  The
remaining ``a=0`` and ``a=1`` branches are decided by blocker duality.

For ``a=0``, local ``g4`` bounds reduce four missing five-faces to 235
possible multisets of active parts.  The eight isolated graph vertices are
then encoded by their four-bit row-membership columns.  Enumerating column
*multiplicities* rather than labelled vertices gives a complete 1,422-state
calculation.  For ``a=1``, the 718 possible missing four-faces form six exact
graph-automorphism types, of which only two survive local ``g4``.  A direct
three-face blocker enumeration finishes those two types.

No SAT/SMT backend is imported or trusted here.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, combinations_with_replacement, permutations
from math import comb, factorial
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record92 import (
    graph_vertex_cover_masks,
    independent_masks,
    minimal_graph_vertex_covers,
    summed_vertex_link_g4,
)
from coxeter_polytopes.combinatorics.missing_edges import (
    TARGET_VERTEX_COUNT,
    decode_graph6,
    macaulay_raise,
)


TARGET_RECORD_INDEX = 54
TARGET_GRAPH6 = "M?????????O?W?W??"
TARGET_MISSING_EDGES = (
    (0, 11),
    (1, 12),
    (2, 12),
    (1, 13),
    (2, 13),
)
ACTIVE_VERTICES = (0, 1, 2, 11, 12, 13)
ISOLATED_VERTICES = tuple(range(3, 11))
MINIMAL_GRAPH_COVERS = (
    (0, 1, 2),
    (1, 2, 11),
    (0, 12, 13),
    (11, 12, 13),
)

EXPECTED_A0_HISTOGRAM = {
    60: 25_200,
    62: 26_880,
    64: 75_600,
    65: 13_440,
    66: 174_720,
    67: 16_800,
    68: 414_120,
    69: 107_520,
    70: 449_400,
    72: 287_420,
    74: 248_640,
    76: 270_340,
    82: 3_360,
    84: 59_360,
    86: 25_760,
    88: 20_160,
    90: 4_200,
    92: 1_120,
    95: 224,
    98: 4_900,
}
EXPECTED_A0_CANONICAL_ORBIT_HISTOGRAM = {
    60: 2,
    62: 2,
    64: 8,
    65: 1,
    66: 6,
    67: 2,
    68: 29,
    69: 2,
    70: 21,
    72: 21,
    74: 9,
    76: 27,
    82: 1,
    84: 12,
    86: 3,
    88: 2,
    90: 6,
    92: 1,
    95: 1,
    98: 3,
}
EXPECTED_A1_AGGREGATE_HISTOGRAM = {
    60: 3_360,
    62: 23_520,
    64: 20_160,
    65: 26_880,
    66: 107_520,
    67: 33_600,
}


def _mask(vertices: Iterable[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def _vertices(mask: int) -> tuple[int, ...]:
    return tuple(
        vertex for vertex in range(TARGET_VERTEX_COUNT) if mask & (1 << vertex)
    )


def _edge_masks(graph6: str) -> tuple[int, ...]:
    vertex_count, edges = decode_graph6(graph6)
    if vertex_count != TARGET_VERTEX_COUNT:
        raise ValueError("record-54 certificate requires fourteen vertices")
    return tuple(_mask(edge) for edge in edges)


def _minimal_transversals(
    hyperedges: Sequence[int], maximum_size: int
) -> tuple[int, ...]:
    values: list[int] = []
    for size in range(1, maximum_size + 1):
        for vertices in combinations(range(TARGET_VERTEX_COUNT), size):
            candidate = _mask(vertices)
            if not all(candidate & edge for edge in hyperedges):
                continue
            if all(
                not all((candidate ^ (1 << vertex)) & edge for edge in hyperedges)
                for vertex in vertices
            ):
                values.append(candidate)
    return tuple(values)


def _base_link_g3_g4(graph6: str) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Return link ``g3,g4`` when every independent four/five-set is a face."""

    independent = {size: independent_masks(graph6, size) for size in range(2, 6)}
    g3_values: list[int] = []
    g4_values: list[int] = []
    for vertex in range(TARGET_VERTEX_COUNT):
        counts = {
            size: sum(bool(mask & (1 << vertex)) for mask in independent[size])
            for size in range(2, 6)
        }
        f0, f1, f2, f3 = (counts[size] for size in range(2, 6))
        g3_values.append(f2 - 7 * f1 + 28 * f0 - 84)
        g4_values.append(f3 - 6 * f2 + 21 * f1 - 56 * f0 + 126)
    return tuple(g3_values), tuple(g4_values)


def _active_automorphisms() -> tuple[dict[int, int], ...]:
    """Enumerate the full automorphism group on the six active vertices."""

    edge_set = {frozenset(edge) for edge in TARGET_MISSING_EDGES}
    values: list[dict[int, int]] = []
    for image in permutations(ACTIVE_VERTICES):
        mapping = dict(zip(ACTIVE_VERTICES, image))
        mapped = {
            frozenset((mapping[first], mapping[second]))
            for first, second in TARGET_MISSING_EDGES
        }
        if mapped == edge_set:
            values.append(mapping)
    return tuple(values)


def _active_orbit_key(mask: int, automorphisms: Sequence[Mapping[int, int]]):
    active = tuple(vertex for vertex in ACTIVE_VERTICES if mask & (1 << vertex))
    canonical = min(
        tuple(sorted(automorphism[vertex] for vertex in active))
        for automorphism in automorphisms
    )
    isolate_count = sum(bool(mask & (1 << vertex)) for vertex in ISOLATED_VERTICES)
    return canonical, isolate_count


@dataclass(frozen=True)
class MissingFourOrbit:
    active_signature: tuple[int, ...]
    isolate_count: int
    multiplicity: int
    representative_mask: int
    extension_count: int
    forced_missing_five_count: int
    survives_strict_local_g4: bool

    def to_record(self) -> dict[str, object]:
        return {
            "active_signature": list(self.active_signature),
            "isolate_count": self.isolate_count,
            "multiplicity": self.multiplicity,
            "representative_mask": self.representative_mask,
            "representative_vertices": list(_vertices(self.representative_mask)),
            "extension_count": self.extension_count,
            "forced_missing_five_count": self.forced_missing_five_count,
            "survives_strict_local_g4": self.survives_strict_local_g4,
        }


def missing_four_orbits(
    graph6: str,
) -> tuple[MissingFourOrbit, ...]:
    """Classify all independent four-sets under ``Aut(H) x S8``."""

    automorphisms = _active_automorphisms()
    if len(automorphisms) != 16:
        raise AssertionError("record-54 active graph should have 16 automorphisms")
    independent_fives = set(independent_masks(graph6, 5))
    base_g3, base_g4 = _base_link_g3_g4(graph6)
    grouped: dict[tuple[tuple[int, ...], int], list[int]] = defaultdict(list)
    for face in independent_masks(graph6, 4):
        grouped[_active_orbit_key(face, automorphisms)].append(face)

    result: list[MissingFourOrbit] = []
    for (active_signature, isolate_count), faces in sorted(grouped.items()):
        extension_counts = {
            sum(
                (face | (1 << vertex)) in independent_fives
                for vertex in range(TARGET_VERTEX_COUNT)
                if not face & (1 << vertex)
            )
            for face in faces
        }
        if len(extension_counts) != 1:
            raise AssertionError("extension count is not constant on an orbit")
        extension_count = extension_counts.pop()
        representative = min(faces)
        extension_vertices = {
            vertex
            for vertex in range(TARGET_VERTEX_COUNT)
            if not representative & (1 << vertex)
            and (representative | (1 << vertex)) in independent_fives
        }
        no_a5_g4 = tuple(
            base_g4[vertex]
            + (6 - extension_count if representative & (1 << vertex) else 0)
            - (1 if vertex in extension_vertices else 0)
            for vertex in range(TARGET_VERTEX_COUNT)
        )
        # Missing five-faces only decrease link g4, so a nonpositive entry is
        # an unconditional contradiction.  The two c=7 types have admissible
        # Macaulay intervals, checked again in the a=1 enumeration.
        survives = all(value >= 1 for value in no_a5_g4)
        result.append(
            MissingFourOrbit(
                active_signature=active_signature,
                isolate_count=isolate_count,
                multiplicity=len(faces),
                representative_mask=representative,
                extension_count=extension_count,
                forced_missing_five_count=10 - extension_count,
                survives_strict_local_g4=survives,
            )
        )

    if sum(orbit.multiplicity for orbit in result) != 718:
        raise AssertionError("the six missing-four orbits do not cover I4(H)")
    # Silence a false-positive unused warning while documenting that the base
    # g3 values are part of the independently reconstructed data.
    if base_g3 != (2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2):
        raise AssertionError("record-54 base vertex-link g3 values changed")
    return tuple(result)


def _a5_incidence_intervals(
    graph6: str, missing_four: int
) -> tuple[tuple[int, int], ...]:
    independent_fives = set(independent_masks(graph6, 5))
    extension_vertices = {
        vertex
        for vertex in range(TARGET_VERTEX_COUNT)
        if not missing_four & (1 << vertex)
        and (missing_four | (1 << vertex)) in independent_fives
    }
    extension_count = len(extension_vertices)
    base_g3, base_g4 = _base_link_g3_g4(graph6)
    intervals: list[tuple[int, int]] = []
    for vertex in range(TARGET_VERTEX_COUNT):
        in_four = bool(missing_four & (1 << vertex))
        g3 = base_g3[vertex] - int(in_four)
        maximum_g4 = macaulay_raise(g3, 3)
        before_a5 = (
            base_g4[vertex]
            + (6 - extension_count if in_four else 0)
            - int(vertex in extension_vertices)
        )
        intervals.append((max(0, before_a5 - maximum_g4), before_a5 - 1))
    return tuple(intervals)


def _small_transversal_from_base_cover(
    missing_fives: Sequence[int], base_minimal_covers: Sequence[int]
) -> bool:
    """Whether a base cover extends to a transversal of size at most four."""

    for cover in base_minimal_covers:
        budget = 4 - cover.bit_count()
        missed = tuple(face for face in missing_fives if not face & cover)
        if not missed:
            return True
        if budget == 0:
            continue
        if budget != 1:
            raise AssertionError("record-54 base covers should have size 3 or 4")
        common = missed[0]
        for face in missed[1:]:
            common &= face
        if common:
            return True
    return False


@dataclass(frozen=True)
class A1TypeAudit:
    active_signature: tuple[int, ...]
    missing_four_orbit_size: int
    representative_mask: int
    candidate_missing_five_count: int
    local_incidence_feasible_triples: int
    no_small_transversal_triples: int
    active_free_no_small_transversal_triples: int
    five_transversal_histogram: tuple[tuple[int, int], ...]

    def to_record(self) -> dict[str, object]:
        return {
            "active_signature": list(self.active_signature),
            "missing_four_orbit_size": self.missing_four_orbit_size,
            "representative_mask": self.representative_mask,
            "representative_vertices": list(_vertices(self.representative_mask)),
            "candidate_missing_five_count": self.candidate_missing_five_count,
            "local_incidence_feasible_triples": self.local_incidence_feasible_triples,
            "no_small_transversal_triples": self.no_small_transversal_triples,
            "active_free_no_small_transversal_triples": (
                self.active_free_no_small_transversal_triples
            ),
            "five_transversal_histogram": {
                str(count): multiplicity
                for count, multiplicity in self.five_transversal_histogram
            },
        }


def audit_a1(
    graph6: str, orbits: Sequence[MissingFourOrbit]
) -> tuple[tuple[A1TypeAudit, ...], tuple[tuple[int, int], ...]]:
    """Exhaust the two locally feasible ``a=1`` orbit representatives."""

    survivors = tuple(orbit for orbit in orbits if orbit.survives_strict_local_g4)
    if len(survivors) != 2 or {orbit.extension_count for orbit in survivors} != {7}:
        raise AssertionError("record-54 should have exactly two c=7 A4 types")

    independent_fives = independent_masks(graph6, 5)
    edge_masks = _edge_masks(graph6)
    prepared = []
    for orbit in survivors:
        missing_four = orbit.representative_mask
        intervals = _a5_incidence_intervals(graph6, missing_four)
        if any(lower != 0 or upper < 0 for lower, upper in intervals):
            raise AssertionError("surviving A4 type has an impossible A5 interval")
        candidates = tuple(
            face
            for face in independent_fives
            if face & missing_four != missing_four
            and all(
                not face & (1 << vertex) or intervals[vertex][1] >= 1
                for vertex in range(TARGET_VERTEX_COUNT)
            )
        )
        if len(candidates) != 462:
            raise AssertionError("a=1 A5 candidate ground should have C(11,5) sets")
        base_covers = _minimal_transversals(edge_masks + (missing_four,), 4)
        if Counter(cover.bit_count() for cover in base_covers) != {3: 3, 4: 3}:
            raise AssertionError("unexpected base blocker for an a=1 type")
        cover_fives = tuple(
            transversal
            for transversal in combinations(range(TARGET_VERTEX_COUNT), 5)
            for transversal in (_mask(transversal),)
            if all(transversal & edge for edge in edge_masks + (missing_four,))
        )
        if len(cover_fives) != 177:
            raise AssertionError("unexpected a=1 base five-cover count")
        prepared.append((orbit, intervals, candidates, base_covers, cover_fives))

    if prepared[0][2] != prepared[1][2]:
        raise AssertionError("the two a=1 types should have the same A5 ground")
    candidates = prepared[0][2]
    cap_one_mask = _mask((11, 12, 13))
    isolate_mask = _mask(ISOLATED_VERTICES)
    local_count = 0
    no_small_counts = [0, 0]
    active_free_counts = [0, 0]
    histograms = [Counter(), Counter()]
    for first_index in range(len(candidates) - 2):
        first = candidates[first_index]
        for second_index in range(first_index + 1, len(candidates) - 1):
            second = candidates[second_index]
            if first & second & cap_one_mask:
                continue
            first_second = first & second
            for third_index in range(second_index + 1, len(candidates)):
                third = candidates[third_index]
                if (first | second) & third & cap_one_mask:
                    continue
                if first_second & third & isolate_mask:
                    continue
                local_count += 1
                missing_fives = (first, second, third)
                for type_index, (_, _, _, base_covers, cover_fives) in enumerate(
                    prepared
                ):
                    if _small_transversal_from_base_cover(
                        missing_fives, base_covers
                    ):
                        continue
                    no_small_counts[type_index] += 1
                    if not any(face & cap_one_mask for face in missing_fives):
                        active_free_counts[type_index] += 1
                    facet_count = sum(
                        all(transversal & face for face in missing_fives)
                        for transversal in cover_fives
                    )
                    histograms[type_index][facet_count] += 1

    audits = tuple(
        A1TypeAudit(
            active_signature=orbit.active_signature,
            missing_four_orbit_size=orbit.multiplicity,
            representative_mask=orbit.representative_mask,
            candidate_missing_five_count=len(candidates),
            local_incidence_feasible_triples=local_count,
            no_small_transversal_triples=no_small_counts[index],
            active_free_no_small_transversal_triples=active_free_counts[index],
            five_transversal_histogram=tuple(sorted(histograms[index].items())),
        )
        for index, (orbit, _, _, _, _) in enumerate(prepared)
    )
    aggregate: Counter[int] = Counter()
    for audit in audits:
        for facet_count, multiplicity in audit.five_transversal_histogram:
            aggregate[facet_count] += audit.missing_four_orbit_size * multiplicity
    return audits, tuple(sorted(aggregate.items()))


def _independent_active_masks() -> tuple[int, ...]:
    edge_set = {frozenset(edge) for edge in TARGET_MISSING_EDGES}
    values: list[int] = []
    for size in range(4):
        for vertices in combinations(ACTIVE_VERTICES, size):
            if all(
                frozenset(edge) not in edge_set for edge in combinations(vertices, 2)
            ):
                values.append(_mask(vertices))
    return tuple(values)


def _isolate_column_distributions(
    allowed_patterns: tuple[int, ...], row_sums: tuple[int, int, int, int]
) -> tuple[tuple[int, ...], ...]:
    """All multiplicity vectors for eight labelled-isolate column patterns."""

    values: list[tuple[int, ...]] = []
    counts = [0] * len(allowed_patterns)

    def recurse(position: int, columns_left: int, sums_left: tuple[int, ...]):
        if position == len(allowed_patterns):
            if columns_left == 0 and not any(sums_left):
                values.append(tuple(counts))
            return
        pattern = allowed_patterns[position]
        maximum = columns_left
        for row in range(4):
            if pattern & (1 << row):
                maximum = min(maximum, sums_left[row])
        for multiplicity in range(maximum + 1):
            next_columns = columns_left - multiplicity
            next_sums = tuple(
                sums_left[row]
                - multiplicity * int(bool(pattern & (1 << row)))
                for row in range(4)
            )
            # Every remaining labelled isolate contributes at most one to a
            # row.  This is a necessary, safe pruning inequality.
            if all(0 <= value <= next_columns for value in next_sums):
                counts[position] = multiplicity
                recurse(position + 1, next_columns, next_sums)
        counts[position] = 0

    recurse(0, len(ISOLATED_VERTICES), row_sums)
    return tuple(values)


@dataclass(frozen=True)
class A0Audit:
    active_mask_count: int
    active_part_multiset_count: int
    column_distribution_signature_count: int
    compressed_state_count: int
    canonical_orbit_count: int
    canonical_orbit_histogram: tuple[tuple[int, int], ...]
    no_small_transversal_family_count: int
    five_transversal_histogram: tuple[tuple[int, int], ...]

    def to_record(self) -> dict[str, object]:
        return {
            "active_mask_count": self.active_mask_count,
            "active_part_multiset_count": self.active_part_multiset_count,
            "column_distribution_signature_count": (
                self.column_distribution_signature_count
            ),
            "compressed_state_count": self.compressed_state_count,
            "canonical_orbit_count": self.canonical_orbit_count,
            "canonical_orbit_histogram": {
                str(count): multiplicity
                for count, multiplicity in self.canonical_orbit_histogram
            },
            "no_small_transversal_family_count": (
                self.no_small_transversal_family_count
            ),
            "five_transversal_histogram": {
                str(count): multiplicity
                for count, multiplicity in self.five_transversal_histogram
            },
        }


def audit_a0(graph6: str, minimal_covers: Sequence[int]) -> A0Audit:
    """Complete compressed blocker enumeration for four missing five-faces."""

    active_masks = _independent_active_masks()
    if len(active_masks) != 21:
        raise AssertionError("record-54 active graph should have 21 independent sets")

    def pattern(active_mask: int) -> int:
        return sum(
            bool(active_mask & cover) << index
            for index, cover in enumerate(minimal_covers)
        )

    active_configs: list[tuple[int, int, int, int]] = []
    for config in combinations_with_replacement(active_masks, 4):
        # Active vertices have base g4=2, so strict positivity permits each
        # active vertex in at most one missing five-face.
        if any(config[first] & config[second] for first in range(4) for second in range(first)):
            continue
        patterns = tuple(pattern(value) for value in config)
        # A minimal graph cover must miss at least two rows; otherwise it is
        # itself, or extends by one point to, a transversal of size <=4.
        if any(
            sum(bool(value & (1 << bit)) for value in patterns) > 2
            for bit in range(len(minimal_covers))
        ):
            continue
        active_configs.append(config)
    if len(active_configs) != 235:
        raise AssertionError("record-54 should have 235 active-part multisets")

    graph_cover_fives = graph_vertex_cover_masks(graph6, 5)
    if len(graph_cover_fives) != 198:
        raise AssertionError("record-54 graph should have 198 five-covers")
    distribution_cache: dict[
        tuple[tuple[int, ...], tuple[int, int, int, int]],
        tuple[tuple[int, ...], ...],
    ] = {}
    histogram: Counter[int] = Counter()
    compressed_states = 0
    canonical_states: dict[
        tuple[tuple[int, ...], tuple[int, ...]], int
    ] = {}
    active_automorphisms = _active_automorphisms()
    row_permutations = tuple(permutations(range(4)))
    active_transforms = tuple(
        {
            active_mask: _mask(
                automorphism[vertex]
                for vertex in ACTIVE_VERTICES
                if active_mask & (1 << vertex)
            )
            for active_mask in active_masks
        }
        for automorphism in active_automorphisms
    )
    column_transforms = {
        row_permutation: tuple(
            sum(
                bool(column & (1 << row_permutation[row])) << row
                for row in range(4)
            )
            for column in range(1 << 4)
        )
        for row_permutation in row_permutations
    }

    def canonical_state(
        config: tuple[int, int, int, int], columns: tuple[int, ...]
    ) -> tuple[tuple[int, ...], tuple[int, ...]]:
        best: tuple[tuple[int, ...], tuple[int, ...]] | None = None
        for transform in active_transforms:
            transformed_active = tuple(
                transform[active_mask] for active_mask in config
            )
            for row_permutation in row_permutations:
                permuted_active = tuple(
                    transformed_active[row_permutation[row]] for row in range(4)
                )
                column_transform = column_transforms[row_permutation]
                permuted_columns = tuple(
                    sorted(column_transform[column] for column in columns)
                )
                candidate = (permuted_active, permuted_columns)
                if best is None or candidate < best:
                    best = candidate
        if best is None:
            raise AssertionError("empty canonical-state orbit")
        return best

    for config in active_configs:
        patterns = tuple(pattern(value) for value in config)
        missed_rows = tuple(
            sum(
                1 << row
                for row, value in enumerate(patterns)
                if not value & (1 << cover_index)
            )
            for cover_index in range(len(minimal_covers))
        )
        # One isolate has a four-bit column recording the missing five-faces
        # containing it.  Strict local g4 forbids 1111.  The blocker condition
        # forbids every column containing all rows missed by one graph cover.
        allowed = tuple(
            column
            for column in range(1 << 4)
            if column != 0b1111
            and all(column & rows != rows for rows in missed_rows)
        )
        row_sums = tuple(5 - value.bit_count() for value in config)
        signature = (allowed, row_sums)
        distributions = distribution_cache.setdefault(
            signature, _isolate_column_distributions(allowed, row_sums)
        )
        ordered_histogram: Counter[int] = Counter()
        for counts in distributions:
            # Equal active parts can occur only for the empty active mask.
            # Their isolate rows must differ, or two missing faces coincide.
            if any(
                config[first] == config[second]
                and not any(
                    multiplicity
                    and bool(column & (1 << first))
                    != bool(column & (1 << second))
                    for column, multiplicity in zip(allowed, counts)
                )
                for first in range(4)
                for second in range(first)
            ):
                continue
            columns = tuple(
                column
                for column, multiplicity in zip(allowed, counts)
                for _ in range(multiplicity)
            )
            if len(columns) != len(ISOLATED_VERTICES):
                raise AssertionError("isolate column count changed")
            missing_fives = list(config)
            for isolate_index, column in enumerate(columns):
                isolate_bit = 1 << ISOLATED_VERTICES[isolate_index]
                for row in range(4):
                    if column & (1 << row):
                        missing_fives[row] |= isolate_bit
            if any(face.bit_count() != 5 for face in missing_fives):
                raise AssertionError("column state does not produce five-faces")
            facet_count = sum(
                all(transversal & face for face in missing_fives)
                for transversal in graph_cover_fives
            )
            canonical_key = canonical_state(config, columns)
            previous_facet_count = canonical_states.setdefault(
                canonical_key, facet_count
            )
            if previous_facet_count != facet_count:
                raise AssertionError(
                    "facet count is not invariant on a canonical orbit"
                )
            labelled_isolate_multiplicity = factorial(8)
            for multiplicity in counts:
                labelled_isolate_multiplicity //= factorial(multiplicity)
            ordered_histogram[facet_count] += labelled_isolate_multiplicity
            compressed_states += 1

        # Rows with equal active part were temporarily labelled.  Since the
        # resulting five-faces are distinct, their permutation action is free.
        row_symmetry = 1
        for multiplicity in Counter(config).values():
            row_symmetry *= factorial(multiplicity)
        for facet_count, multiplicity in ordered_histogram.items():
            if multiplicity % row_symmetry:
                raise AssertionError("row-symmetry quotient is nonintegral")
            histogram[facet_count] += multiplicity // row_symmetry

    return A0Audit(
        active_mask_count=len(active_masks),
        active_part_multiset_count=len(active_configs),
        column_distribution_signature_count=len(distribution_cache),
        compressed_state_count=compressed_states,
        canonical_orbit_count=len(canonical_states),
        canonical_orbit_histogram=tuple(
            sorted(Counter(canonical_states.values()).items())
        ),
        no_small_transversal_family_count=sum(histogram.values()),
        five_transversal_histogram=tuple(sorted(histogram.items())),
    )


@dataclass(frozen=True)
class Record54Audit:
    graph6: str
    missing_edges: tuple[tuple[int, int], ...]
    minimal_graph_covers: tuple[int, ...]
    schedule_link_g4_sums: tuple[tuple[int, int], ...]
    excluded_high_a_branches: tuple[tuple[int, int, int], ...]
    missing_four_orbits: tuple[MissingFourOrbit, ...]
    a0: A0Audit
    a1_types: tuple[A1TypeAudit, ...]
    a1_aggregate_histogram: tuple[tuple[int, int], ...]

    @property
    def verified(self) -> bool:
        return (
            self.missing_edges == TARGET_MISSING_EDGES
            and self.minimal_graph_covers
            == tuple(_mask(cover) for cover in MINIMAL_GRAPH_COVERS)
            and self.schedule_link_g4_sums
            == ((0, 24), (1, 18), (2, 12), (3, 6))
            and all(total < minimum for _, total, minimum in self.excluded_high_a_branches)
            and len(self.missing_four_orbits) == 6
            and [orbit.multiplicity for orbit in self.missing_four_orbits]
            == [70, 112, 224, 32, 224, 56]
            and [orbit.extension_count for orbit in self.missing_four_orbits]
            == [10, 9, 7, 7, 8, 8]
            and [orbit.survives_strict_local_g4 for orbit in self.missing_four_orbits]
            == [False, False, True, True, False, False]
            and self.a0.active_mask_count == 21
            and self.a0.active_part_multiset_count == 235
            and self.a0.column_distribution_signature_count == 35
            and self.a0.compressed_state_count == 1_422
            and self.a0.canonical_orbit_count == 159
            and dict(self.a0.canonical_orbit_histogram)
            == EXPECTED_A0_CANONICAL_ORBIT_HISTOGRAM
            and sum(dict(self.a0.canonical_orbit_histogram).values()) == 159
            and self.a0.no_small_transversal_family_count == 2_229_164
            and dict(self.a0.five_transversal_histogram) == EXPECTED_A0_HISTOGRAM
            and sum(dict(self.a0.five_transversal_histogram).values())
            == 2_229_164
            and 96 not in dict(self.a0.five_transversal_histogram)
            and len(self.a1_types) == 2
            and [audit.missing_four_orbit_size for audit in self.a1_types] == [224, 32]
            and all(audit.candidate_missing_five_count == 462 for audit in self.a1_types)
            and all(audit.local_incidence_feasible_triples == 492_030 for audit in self.a1_types)
            and all(audit.no_small_transversal_triples == 840 for audit in self.a1_types)
            and all(
                audit.active_free_no_small_transversal_triples == 840
                for audit in self.a1_types
            )
            and [dict(audit.five_transversal_histogram) for audit in self.a1_types]
            == [{64: 90, 65: 120, 66: 480, 67: 150}, {60: 105, 62: 735}]
            and dict(self.a1_aggregate_histogram) == EXPECTED_A1_AGGREGATE_HISTOGRAM
            and sum(dict(self.a1_aggregate_histogram).values()) == 256 * 840
            and 94 not in dict(self.a1_aggregate_histogram)
        )

    def to_record(self) -> dict[str, object]:
        return {
            "graph6": self.graph6,
            "missing_edges": [list(edge) for edge in self.missing_edges],
            "minimal_graph_cover_masks": list(self.minimal_graph_covers),
            "minimal_graph_cover_vertices": [
                list(_vertices(mask)) for mask in self.minimal_graph_covers
            ],
            "schedule_link_g4_sums": [
                {"missing_four_count": a, "summed_vertex_link_g4": total}
                for a, total in self.schedule_link_g4_sums
            ],
            "strict_vertex_link_g4_required_sum": 14,
            "excluded_high_a_branches": [
                {
                    "missing_four_count": a,
                    "summed_vertex_link_g4": total,
                    "required_minimum": minimum,
                }
                for a, total, minimum in self.excluded_high_a_branches
            ],
            "a0": {
                "forced_missing_five_count": 4,
                "required_facet_count": 96,
                **self.a0.to_record(),
                "required_count_absent": 96
                not in dict(self.a0.five_transversal_histogram),
            },
            "a1": {
                "independent_missing_four_count": 718,
                "orbit_count": len(self.missing_four_orbits),
                "orbits": [orbit.to_record() for orbit in self.missing_four_orbits],
                "locally_feasible_orbit_count": len(self.a1_types),
                "forced_missing_five_count_on_surviving_types": 3,
                "required_facet_count": 94,
                "types": [audit.to_record() for audit in self.a1_types],
                "aggregate_five_transversal_histogram": {
                    str(count): multiplicity
                    for count, multiplicity in self.a1_aggregate_histogram
                },
                "required_count_absent": 94
                not in dict(self.a1_aggregate_histogram),
            },
            "blocker_conclusion": (
                "If a wrong-sized minimal transversal exists, purity fails. "
                "Otherwise every five-transversal is minimal, and its count "
                "must equal the scheduled facet count.  The exhaustive a=0 "
                "and a=1 histograms omit 96 and 94 respectively."
            ),
            "verified": self.verified,
        }


def audit_record54(
    graph6: str, schedules: Sequence[Mapping[str, object]]
) -> Record54Audit:
    """Recompute the complete solver-independent record-54 exclusion."""

    if graph6 != TARGET_GRAPH6:
        raise ValueError("certificate is scoped to H9/14 record 54")
    _, decoded_edges = decode_graph6(graph6)
    if tuple(decoded_edges) != TARGET_MISSING_EDGES:
        raise ValueError("record-54 missing-edge graph changed")
    parsed_schedules: dict[int, tuple[int, ...]] = {}
    for schedule in schedules:
        missing_four_count = schedule.get("missing_four_count")
        f_vector = schedule.get("f_vector")
        if not isinstance(missing_four_count, int) or not isinstance(
            f_vector, (list, tuple)
        ):
            raise ValueError("malformed record-54 Dehn--Sommerville schedule")
        parsed_schedules[missing_four_count] = tuple(int(value) for value in f_vector)
    if tuple(sorted(parsed_schedules)) != (0, 1, 2, 3):
        raise ValueError("record-54 schedules must be exactly a=0,1,2,3")
    link_sums = tuple(
        (a, summed_vertex_link_g4(parsed_schedules[a]))
        for a in sorted(parsed_schedules)
    )
    high_a = tuple((a, total, 14) for a, total in link_sums if a >= 2)

    minimal_covers = minimal_graph_vertex_covers(graph6)
    expected_covers = tuple(_mask(cover) for cover in MINIMAL_GRAPH_COVERS)
    if minimal_covers != expected_covers:
        raise AssertionError("record-54 minimal graph covers changed")
    orbits = missing_four_orbits(graph6)
    a0 = audit_a0(graph6, minimal_covers)
    a1_types, a1_aggregate = audit_a1(graph6, orbits)
    audit = Record54Audit(
        graph6=graph6,
        missing_edges=tuple(decoded_edges),
        minimal_graph_covers=minimal_covers,
        schedule_link_g4_sums=link_sums,
        excluded_high_a_branches=high_a,
        missing_four_orbits=orbits,
        a0=a0,
        a1_types=a1_types,
        a1_aggregate_histogram=a1_aggregate,
    )
    if not audit.verified:
        raise AssertionError("record-54 arithmetic certificate did not verify")
    return audit

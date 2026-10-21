#!/usr/bin/env python3
"""Enumerate record-30 a=0 target clutters with an unordered Boolean model.

The search variables are the 1,372 possible independent A5 minimal nonfaces
and the 326 five-covers of the fixed 3K2.  Facet variables are linked exactly
to the blocker of the chosen eight A5s.  Dynamic lex-leader constraints use a
full-action lex canonicalizer: noncanonical models teach the solver a valid
global symmetry breaker, while canonical models are emitted once and blocked.

Z3 is only the search backend.  A frozen theorem certificate must replay the
emitted canonical model set and the final exhaustion independently.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path
import subprocess
import tempfile
from time import monotonic

import z3


N = 14
TARGET_ROWS = 8
TARGET_FACETS = 132
ACTIVE = frozenset((0, 1, 2, 11, 12, 13))
EDGES = ((0, 11), (1, 12), (2, 13))
EDGE_MASKS = tuple((1 << left) | (1 << right) for left, right in EDGES)
ROOT = Path(__file__).resolve().parents[1]
CANONICALIZER_SOURCE = ROOT / "scripts/h9_14_record30_a0_lex_canonicalizer.cpp"
ORBIT_FORMS_SOURCE = ROOT / "scripts/h9_14_record30_a0_orbit_forms.cpp"


def mask(vertices) -> int:
    return sum(1 << vertex for vertex in vertices)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(N) if value & (1 << vertex))


def independent(value: int) -> bool:
    return all(value & edge != edge for edge in EDGE_MASKS)


def hits_edges(value: int) -> bool:
    return all(value & edge for edge in EDGE_MASKS)


def local_cap(value: int) -> int:
    size = value.bit_count()
    active = sum(value >> vertex & 1 for vertex in ACTIVE)
    if size == 1:
        return {0: 7, 1: 3}[active]
    if size == 2:
        return {0: 7, 1: 3, 2: 1}[active]
    return {0: 6, 1: 3, 2: 1, 3: 0}[active]


def transform(value: int, permutation: tuple[int, ...]) -> int:
    return sum(1 << permutation[vertex] for vertex in vertices(value))


def base_generators() -> tuple[tuple[int, ...], ...]:
    result = []

    def moved(*cycles: tuple[int, ...]) -> tuple[int, ...]:
        permutation = list(range(N))
        for cycle in cycles:
            for source, target in zip(cycle, cycle[1:] + cycle[:1], strict=True):
                permutation[source] = target
        return tuple(permutation)

    # Wreath-product generators for Aut(3K2).
    result.append(moved((0, 1), (11, 12)))
    result.append(moved((0, 1, 2), (11, 12, 13)))
    result.append(moved((0, 11)))
    # Adjacent transpositions generate S8 on the isolates.
    result.extend(moved((left, left + 1)) for left in range(3, 10))
    return tuple(result)


class Model:
    def __init__(self, timeout_ms: int, initial_generators: bool) -> None:
        self.rows = tuple(sorted(
            value for value in range(1 << N)
            if value.bit_count() == 5 and independent(value)
        ))
        self.facets = tuple(sorted(
            value for value in range(1 << N)
            if value.bit_count() == 5 and hits_edges(value)
        ))
        if (len(self.rows), len(self.facets)) != (1372, 326):
            raise AssertionError("record-30 universe count changed")
        self.row_index = {value: index for index, value in enumerate(self.rows)}
        self.selected = tuple(z3.Bool(f"a5_{index}") for index in range(len(self.rows)))
        self.facet = tuple(z3.Bool(f"facet_{index}") for index in range(len(self.facets)))
        self.solver = z3.Solver()
        self.solver.set(timeout=timeout_ms)
        self.solver.add(z3.PbEq([(item, 1) for item in self.selected], TARGET_ROWS))
        self.solver.add(z3.PbEq([(item, 1) for item in self.facet], TARGET_FACETS))

        for size in (1, 2, 3):
            for chosen in combinations(range(N), size):
                value = mask(chosen)
                if not independent(value):
                    continue
                self.solver.add(z3.PbLe([
                    (self.selected[index], 1)
                    for index, row in enumerate(self.rows) if row & value == value
                ], local_cap(value)))

        for size in (3, 4):
            for chosen in combinations(range(N), size):
                value = mask(chosen)
                if hits_edges(value):
                    self.solver.add(z3.Or(*(
                        self.selected[index]
                        for index, row in enumerate(self.rows) if not row & value
                    )))

        for facet_index, complement in enumerate(self.facets):
            disjoint_rows = tuple(
                self.selected[index]
                for index, row in enumerate(self.rows) if not row & complement
            )
            for row in disjoint_rows:
                self.solver.add(z3.Or(z3.Not(self.facet[facet_index]), z3.Not(row)))
            self.solver.add(z3.Or(self.facet[facet_index], *disjoint_rows))

        facet_index = {value: index for index, value in enumerate(self.facets)}
        for chosen in combinations(range(N), 6):
            six = mask(chosen)
            members = tuple(
                facet_index[six ^ (1 << vertex)]
                for vertex in chosen if six ^ (1 << vertex) in facet_index
            )
            if members:
                degree = z3.Sum(*(z3.If(self.facet[index], 1, 0) for index in members))
                self.solver.add(z3.Or(degree == 0, degree == 2))

        self.lex_actions: set[tuple[int, ...]] = set()
        self.lex_cuts: set[tuple[tuple[int, ...], int]] = set()
        self.last_lex_cut: tuple[tuple[int, ...], int] | None = None
        self.lex_count = 0
        if initial_generators:
            for permutation in base_generators():
                self.add_lex_leader(permutation)

    def enable_double_lex(self) -> None:
        """Select the lossless row-sorted/isolate-column-sorted orbit slice."""
        positions = tuple(tuple(
            z3.Bool(f"position_{rank}_row_{index}")
            for index in range(len(self.rows))
        ) for rank in range(TARGET_ROWS))
        ranks = []
        for rank, choices in enumerate(positions):
            self.solver.add(z3.PbEq([(choice, 1) for choice in choices], 1))
            ranks.append(z3.Sum(*(
                z3.If(choice, index, 0) for index, choice in enumerate(choices)
            )))
            for index, choice in enumerate(choices):
                self.solver.add(z3.Or(z3.Not(choice), self.selected[index]))
        for left, right in zip(ranks, ranks[1:]):
            self.solver.add(left < right)
        columns = []
        for vertex in range(3, 11):
            bits = []
            for rank, choices in enumerate(positions):
                bit = z3.Bool(f"position_{rank}_vertex_{vertex}")
                self.solver.add(bit == z3.Or(*(
                    choices[index] for index, row in enumerate(self.rows)
                    if row & (1 << vertex)
                )))
                bits.append(bit)
            columns.append(z3.Sum(*(
                z3.If(bit, 1 << rank, 0) for rank, bit in enumerate(bits)
            )))
        for left, right in zip(columns, columns[1:]):
            self.solver.add(left <= right)

    def add_lex_leader(self, permutation: tuple[int, ...]) -> bool:
        """Require characteristic(selected) >= its image (tuple-min leader)."""
        if permutation in self.lex_actions:
            return False
        if sorted(permutation) != list(range(N)):
            raise AssertionError("invalid vertex permutation")
        image: list[z3.BoolRef | None] = [None] * len(self.rows)
        for index, row in enumerate(self.rows):
            image[self.row_index[transform(row, permutation)]] = self.selected[index]
        if any(item is None for item in image):
            raise AssertionError("row action is not a permutation")
        tag = self.lex_count
        self.lex_count += 1
        prefix: z3.BoolRef = z3.Bool(f"lex_{tag}_prefix_0")
        self.solver.add(prefix)
        for index, (left, right) in enumerate(zip(self.selected, image, strict=True)):
            assert right is not None
            # With 0 < 1, reject the first unequal coordinate 0,1.
            self.solver.add(z3.Or(z3.Not(prefix), left, z3.Not(right)))
            if index + 1 < len(self.rows):
                following = z3.Bool(f"lex_{tag}_prefix_{index + 1}")
                self.solver.add(following == z3.And(prefix, left == right))
                prefix = following
        self.lex_actions.add(permutation)
        return True

    def add_violated_lex_cut(
        self,
        permutation: tuple[int, ...],
        family: tuple[int, ...],
        canonical: tuple[int, ...],
    ) -> bool:
        """Add the one violated prefix clause from characteristic >= image."""
        image: list[z3.BoolRef | None] = [None] * len(self.rows)
        image_index: list[int | None] = [None] * len(self.rows)
        for index, row in enumerate(self.rows):
            target = self.row_index[transform(row, permutation)]
            image[target] = self.selected[index]
            image_index[target] = index
        if any(item is None for item in image) or any(item is None for item in image_index):
            raise AssertionError("row action is not a permutation")
        family_set = set(family)
        canonical_set = set(canonical)
        first = next(
            index for index, row in enumerate(self.rows)
            if (row in family_set) != (row in canonical_set)
        )
        if self.rows[first] in family_set or self.rows[first] not in canonical_set:
            raise AssertionError("canonical characteristic vector is not larger")
        key = (permutation, first)
        if key in self.lex_cuts:
            return False
        prefix_differences = []
        for index in range(first):
            assert image[index] is not None
            prefix_differences.append(self.selected[index] != image[index])
        assert image[first] is not None
        self.solver.add(z3.Or(
            *prefix_differences, self.selected[first], z3.Not(image[first])
        ))
        self.lex_cuts.add(key)
        self.last_lex_cut = key
        return True

    def family(self, model: z3.ModelRef) -> tuple[int, ...]:
        return tuple(
            self.rows[index] for index, item in enumerate(self.selected)
            if z3.is_true(model.eval(item, model_completion=True))
        )

    def block_canonical(self, family: tuple[int, ...]) -> None:
        if len(family) != TARGET_ROWS:
            raise AssertionError("terminal family size changed")
        self.solver.add(z3.Or(*(z3.Not(self.selected[self.row_index[row]]) for row in family)))

    def block_forms(self, forms: tuple[tuple[int, ...], ...]) -> None:
        for family in forms:
            self.block_canonical(family)


def canonicalize(executable: Path, family: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    document = json.loads(subprocess.run(
        [str(executable), *map(str, family)], cwd=ROOT, check=True,
        text=True, capture_output=True,
    ).stdout)
    canonical = tuple(map(int, document["canonical"]))
    permutation = tuple(map(int, document["permutation"]))
    if tuple(sorted(transform(row, permutation) for row in family)) != canonical:
        raise AssertionError("canonical action replay failed")
    return canonical, permutation


def characteristic_greater(left: tuple[int, ...], right: tuple[int, ...]) -> bool:
    """For equal-size sets, tuple-min masks mean characteristic-vector max."""
    left_set = set(left)
    right_set = set(right)
    first = min(left_set ^ right_set)
    return first in left_set


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-orbits", type=int, default=10000)
    parser.add_argument("--timeout-ms", type=int, default=300_000)
    parser.add_argument("--dynamic-lex", action="store_true")
    args = parser.parse_args()
    started = monotonic()
    search = Model(args.timeout_ms, False)
    if not args.dynamic_lex:
        search.enable_double_lex()
    canonical_models: list[tuple[int, ...]] = []
    noncanonical_models = 0
    status = "max_orbits"
    blocked_forms = 0
    with tempfile.TemporaryDirectory(prefix="record30_a0_unordered_") as temporary:
        executable = Path(temporary) / (
            "lex_canonicalizer" if args.dynamic_lex else "orbit_forms"
        )
        subprocess.run([
            "clang++", "-std=c++20", "-O3", "-Wall", "-Wextra", "-Werror",
            str(CANONICALIZER_SOURCE if args.dynamic_lex else ORBIT_FORMS_SOURCE),
            "-o", str(executable),
        ], cwd=ROOT, check=True)
        while len(canonical_models) < args.max_orbits:
            result = search.solver.check()
            if result == z3.unsat:
                status = "exhausted"
                break
            if result != z3.sat:
                status = str(result)
                break
            family = search.family(search.solver.model())
            if not args.dynamic_lex:
                orbit = json.loads(subprocess.run(
                    [str(executable), *map(str, family)], cwd=ROOT, check=True,
                    text=True, capture_output=True,
                ).stdout)
                forms = tuple(tuple(map(int, form)) for form in orbit["forms"])
                if not forms or family not in forms:
                    raise AssertionError("double-lex model absent from its orbit slice")
                representative = min(forms)
                if representative in canonical_models:
                    raise AssertionError("a fully blocked double-lex orbit repeated")
                canonical_models.append(representative)
                search.block_forms(forms)
                blocked_forms += len(forms)
                print(json.dumps({
                    "event": "canonical_orbit",
                    "orbit": len(canonical_models),
                    "representative": representative,
                    "double_lex_form_count": len(forms),
                    "blocked_form_count": blocked_forms,
                    "elapsed_seconds": round(monotonic() - started, 3),
                }), flush=True)
                continue
            canonical, permutation = canonicalize(executable, family)
            if family != canonical:
                if not characteristic_greater(canonical, family):
                    raise AssertionError("canonical tuple is not the characteristic-vector leader")
                if not search.add_violated_lex_cut(
                    permutation, family, canonical
                ):
                    raise AssertionError("a learned lex cut failed to exclude its model")
                noncanonical_models += 1
                print(json.dumps({
                    "event": "learned_symmetry",
                    "noncanonical_models": noncanonical_models,
                    "lex_cut_count": len(search.lex_cuts),
                    "family": family,
                    "canonical": canonical,
                    "cut_first_index": search.last_lex_cut[1] if search.last_lex_cut else None,
                    "cut_permutation": permutation,
                    "elapsed_seconds": round(monotonic() - started, 3),
                }), flush=True)
                continue
            if family in canonical_models:
                raise AssertionError("canonical model repeated")
            canonical_models.append(family)
            search.block_canonical(family)
            print(json.dumps({
                "event": "canonical_orbit",
                "orbit": len(canonical_models),
                "family": family,
                "noncanonical_models": noncanonical_models,
                "lex_action_count": len(search.lex_actions),
                "elapsed_seconds": round(monotonic() - started, 3),
            }), flush=True)
    print(json.dumps({
        "kind": "exploratory_h9_14_record30_a0_unordered_orbit_stream",
        "status": status,
        "canonical_orbit_count": len(canonical_models),
        "canonical_models": canonical_models,
        "noncanonical_models_rejected": noncanonical_models,
        "lex_action_count": len(search.lex_actions),
        "lex_cut_count": len(search.lex_cuts),
        "initial_generator_count": 0,
        "double_lex_quotient": not args.dynamic_lex,
        "blocked_double_lex_form_count": blocked_forms,
        "elapsed_seconds": round(monotonic() - started, 6),
        "solver_independent_certificate": False,
        "scope": "Compact a=0 incidence targets only; no finite-volume claim.",
    }), flush=True)


if __name__ == "__main__":
    main()

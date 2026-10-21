#!/usr/bin/env python3
"""Build a proof-producing SAT replay of the Record-30 ``a=0`` census.

This is deliberately dual to the production connected-ridge search.  Its
primary variables choose eight independent five-element minimal nonfaces.
The 326 possible five-element facet complements are then defined as their
literal blocker.  The formula imposes the frozen local caps, destroys every
three/four-element transversal, requires exactly 132 blocker elements, and
requires every nonempty ridge group to have degree two.

For a finite completeness check we fix the type-(3 active + 2 isolate)
facet-complement seed ``{0,1,2,3,4}``.  Every 132-facet target has such a
complement and the full group Aut(3K2) x S8 is transitive on them.  Every
labelled image of each supplied target orbit containing that seed is blocked
literally.  UNSAT therefore says that the supplied orbits exhaust a *larger*
necessary-incidence class than the connected production census: connectivity
and the remaining face-vector entries are not needed by this formula.

The emitted DIMACS file is intended for a proof-producing solver such as
CaDiCaL and an independent DRAT/LRAT checker.  Merely obtaining ``unsat`` from
an SMT/SAT process is not treated as a frozen certificate.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path
import sys
from collections import Counter
from typing import Sequence


N = 14
ROW_COUNT = 8
FACET_COUNT = 132
FULL = (1 << N) - 1
ACTIVE = (0, 1, 2, 11, 12, 13)
ISOLATES = tuple(range(3, 11))
EDGES = ((0, 11), (1, 12), (2, 13))
EDGE_MASKS = tuple((1 << left) | (1 << right) for left, right in EDGES)
SEED = sum(1 << vertex for vertex in (0, 1, 2, 3, 4))
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CENSUS = (
    ROOT / "data/enumerations/h9_14_record30_a0_connected_ridge_census.json"
)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(N) if value >> vertex & 1)


def independent(value: int) -> bool:
    return all(value & edge != edge for edge in EDGE_MASKS)


def hits_edges(value: int) -> bool:
    return all(value & edge for edge in EDGE_MASKS)


def local_cap(value: int) -> int:
    size = value.bit_count()
    active = sum(value >> vertex & 1 for vertex in ACTIVE)
    if size == 1:
        return (7, 3)[active]
    if size == 2:
        return (7, 3, 1)[active]
    if size == 3:
        return (6, 3, 1, 0)[active]
    raise ValueError("local caps are frozen only in sizes one through three")


def canonical_sha(value: object) -> str:
    return sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode()
    ).hexdigest()


class CNF:
    """Small auditable CNF builder with exact unary threshold counters."""

    def __init__(self) -> None:
        self.variable_count = 0
        self.clauses: list[tuple[int, ...]] = []

    def variable(self) -> int:
        self.variable_count += 1
        return self.variable_count

    def clause(self, *literals: int) -> None:
        if not literals:
            self.clauses.append(())
            return
        if any(literal == 0 for literal in literals):
            raise ValueError("zero is not a CNF literal")
        values = set(literals)
        if any(-literal in values for literal in values):
            return
        self.clauses.append(tuple(sorted(values, key=lambda item: (abs(item), item))))

    def _equiv_or_and(
        self,
        output: int,
        left: int | bool,
        lower: int | bool,
        item: int,
    ) -> None:
        """Encode ``output <-> left or (lower and item)`` exactly."""
        if left is True:
            self.clause(output)
            return
        if left is False:
            if lower is False:
                self.clause(-output)
                return
            if lower is True:
                self.clause(-output, item)
                self.clause(output, -item)
                return
            self.clause(-output, lower)
            self.clause(-output, item)
            self.clause(output, -lower, -item)
            return

        self.clause(-left, output)
        if lower is True:
            self.clause(-item, output)
        elif lower is not False:
            self.clause(-lower, -item, output)
        if lower is not True:
            if lower is False:
                self.clause(-output, left)
            else:
                self.clause(-output, left, lower)
        self.clause(-output, left, item)

    def thresholds(
        self, variables: Sequence[int], limit: int
    ) -> tuple[int | bool, ...]:
        """Return exact predicates ``sum(variables) >= j`` for ``j<=limit``."""
        if limit < 0:
            raise ValueError("negative threshold limit")
        previous: list[int | bool] = [True] + [False] * limit
        for position, item in enumerate(variables, start=1):
            current: list[int | bool] = [True] + [False] * limit
            for threshold in range(1, min(position, limit) + 1):
                output = self.variable()
                self._equiv_or_and(
                    output,
                    previous[threshold],
                    previous[threshold - 1],
                    item,
                )
                current[threshold] = output
            previous = current
        return tuple(previous)

    def at_most(self, variables: Sequence[int], bound: int) -> None:
        if bound < 0:
            self.clause()
            return
        if len(variables) <= bound:
            return
        threshold = self.thresholds(variables, bound + 1)[bound + 1]
        if threshold is True:
            self.clause()
        elif threshold is not False:
            self.clause(-threshold)

    def exactly(self, variables: Sequence[int], target: int) -> None:
        if target < 0 or target > len(variables):
            self.clause()
            return
        thresholds = self.thresholds(variables, target + 1)
        lower = thresholds[target]
        upper = thresholds[target + 1]
        if lower is False:
            self.clause()
        elif lower is not True:
            self.clause(lower)
        if upper is True:
            self.clause()
        elif upper is not False:
            self.clause(-upper)

    def lex_greater_equal(
        self, left: Sequence[int], right: Sequence[int]
    ) -> None:
        """Encode ``left >=lex right`` for Boolean vectors with ``False<True``."""
        if len(left) != len(right):
            raise ValueError("lex vectors have different lengths")
        unequal_positions = [
            (a, b) for a, b in zip(left, right, strict=True) if a != b
        ]
        prefix: int | bool = True
        for position, (a, b) in enumerate(unequal_positions):
            # A first difference False,True is forbidden.
            if prefix is True:
                self.clause(a, -b)
            else:
                self.clause(-prefix, a, -b)
            if position + 1 == len(unequal_positions):
                break
            following = self.variable()
            # following <-> prefix and (a == b).
            if prefix is not True:
                self.clause(-following, prefix)
            self.clause(-following, -a, b)
            self.clause(-following, a, -b)
            if prefix is True:
                self.clause(-a, -b, following)
                self.clause(a, b, following)
            else:
                self.clause(-prefix, -a, -b, following)
                self.clause(-prefix, a, b, following)
            prefix = following

    def write(self, path: Path) -> str:
        digest = sha256()
        with path.open("w", encoding="ascii") as stream:
            header = f"p cnf {self.variable_count} {len(self.clauses)}\n"
            stream.write(header)
            digest.update(header.encode("ascii"))
            for clause in self.clauses:
                line = " ".join(map(str, clause)) + " 0\n"
                stream.write(line)
                digest.update(line.encode("ascii"))
        return digest.hexdigest()


def normalize_targets(document: object) -> tuple[tuple[int, ...], ...]:
    if not isinstance(document, dict):
        raise ValueError("census JSON must be an object")
    raw_targets = document.get("target_orbit_representatives")
    if not isinstance(raw_targets, list):
        raise ValueError("census has no target representative list")
    declared = document.get("target_orbit_count")
    if declared != len(raw_targets):
        raise ValueError("declared target orbit count does not match list")
    targets = []
    for target_index, raw_family in enumerate(raw_targets):
        if not isinstance(raw_family, list) or len(raw_family) != ROW_COUNT:
            raise ValueError(f"target {target_index} does not have eight A5s")
        family = []
        for raw in raw_family:
            if isinstance(raw, int):
                value = raw
            elif isinstance(raw, list) and all(isinstance(v, int) for v in raw):
                if len(set(raw)) != len(raw):
                    raise ValueError(f"target {target_index} repeats a vertex")
                value = sum(1 << vertex for vertex in raw)
            else:
                raise ValueError(f"target {target_index} has a malformed A5")
            if value.bit_count() != 5 or not independent(value):
                raise ValueError(f"target {target_index} has an invalid A5")
            family.append(value)
        normalized = tuple(sorted(family))
        if len(set(normalized)) != ROW_COUNT:
            raise ValueError(f"target {target_index} repeats an A5")
        targets.append(normalized)
    if len(set(targets)) != len(targets):
        raise ValueError("target list repeats a family literally")
    return tuple(targets)


def active_actions() -> tuple[tuple[int, ...], ...]:
    actions = []
    for pair_permutation in permutations(range(3)):
        for flips in product((0, 1), repeat=3):
            action = list(range(N))
            for source, (left, right) in enumerate(EDGES):
                moved_left, moved_right = EDGES[pair_permutation[source]]
                if flips[source]:
                    moved_left, moved_right = moved_right, moved_left
                action[left] = moved_left
                action[right] = moved_right
            actions.append(tuple(action))
    if len(set(actions)) != 48:
        raise AssertionError("Aut(3K2) active action count changed")
    return tuple(actions)


def seed_stabilizer_generators() -> tuple[tuple[int, ...], ...]:
    """Generators of the order-8640 stabilizer of the canonical seed."""
    result = []

    def moved(*cycles: Sequence[int]) -> tuple[int, ...]:
        action = list(range(N))
        for cycle in cycles:
            for source, target in zip(
                cycle, cycle[1:] + cycle[:1], strict=True
            ):
                action[source] = target
        if transform(SEED, action) != SEED:
            raise AssertionError("claimed seed-stabilizer generator moves seed")
        return tuple(action)

    result.append(moved((0, 1), (11, 12)))
    result.append(moved((0, 1, 2), (11, 12, 13)))
    result.append(moved((3, 4)))
    result.extend(moved((left, left + 1)) for left in range(5, 10))
    if len(set(result)) != 8:
        raise AssertionError("seed-stabilizer generator count changed")
    return tuple(result)


def transform(value: int, action: Sequence[int]) -> int:
    return sum(1 << action[vertex] for vertex in vertices(value))


def transform_family(
    family: Sequence[int], action: Sequence[int]
) -> tuple[int, ...]:
    return tuple(sorted(transform(row, action) for row in family))


def characteristic_lex_greater_equal(
    left: Sequence[int], right: Sequence[int]
) -> bool:
    """Compare equal-cardinality set characteristic vectors with False<True."""
    left_set = set(left)
    right_set = set(right)
    if len(left_set) != len(left) or len(right_set) != len(right):
        raise ValueError("characteristic family repeats an element")
    difference = left_set ^ right_set
    return not difference or min(difference) in left_set


def seed_stabilizer_lex_admissible(family: Sequence[int]) -> bool:
    return all(
        characteristic_lex_greater_equal(
            family, transform_family(family, generator)
        )
        for generator in seed_stabilizer_generators()
    )


def seed_forms(
    family: Sequence[int], actions: Sequence[Sequence[int]]
) -> set[tuple[int, ...]]:
    """Enumerate all full-group images whose blocker contains ``SEED``."""
    active_versions = tuple(
        tuple(transform(row & sum(1 << v for v in ACTIVE), action) for row in family)
        for action in actions
    )
    isolate_patterns = tuple(
        sum(((row >> vertex) & 1) << index for index, vertex in enumerate(ISOLATES))
        for row in family
    )
    result: set[tuple[int, ...]] = set()
    for isolate_permutation in permutations(ISOLATES):
        isolate_parts = tuple(
            sum(
                ((pattern >> index) & 1) << isolate_permutation[index]
                for index in range(8)
            )
            for pattern in isolate_patterns
        )
        for active_parts in active_versions:
            moved = tuple(
                sorted(
                    active | isolate
                    for active, isolate in zip(
                        active_parts, isolate_parts, strict=True
                    )
                )
            )
            if all(SEED & row for row in moved):
                result.add(moved)
    return result


def semantic_target_audit(family: Sequence[int]) -> dict[str, object]:
    """Check directly that a supplied orbit representative satisfies the model."""
    normalized = tuple(sorted(family))
    if len(normalized) != ROW_COUNT or len(set(normalized)) != ROW_COUNT:
        raise ValueError("a semantic target must contain eight distinct A5s")
    if any(value.bit_count() != 5 or not independent(value) for value in normalized):
        raise ValueError("a semantic target contains an invalid A5")

    local_slacks = []
    for size in (1, 2, 3):
        for chosen in combinations(range(N), size):
            subset = sum(1 << vertex for vertex in chosen)
            if not independent(subset):
                continue
            incidence = sum(row & subset == subset for row in normalized)
            local_slacks.append(local_cap(subset) - incidence)
    if min(local_slacks) < 0:
        raise ValueError("a semantic target violates a frozen local cap")

    surviving_small = []
    for size in (3, 4):
        for chosen in combinations(range(N), size):
            transversal = sum(1 << vertex for vertex in chosen)
            if hits_edges(transversal) and all(
                row & transversal for row in normalized
            ):
                surviving_small.append(transversal)
    if surviving_small:
        raise ValueError("a semantic target has a size-three/four transversal")

    blocker = tuple(
        complement
        for complement in range(1 << N)
        if complement.bit_count() == 5
        and hits_edges(complement)
        and all(complement & row for row in normalized)
    )
    if len(blocker) != FACET_COUNT:
        raise ValueError("a semantic target does not have a 132-element blocker")
    blocker_set = set(blocker)
    ridge_degrees = []
    for six_vertices in combinations(range(N), 6):
        six = sum(1 << vertex for vertex in six_vertices)
        degree = sum(
            six ^ (1 << vertex) in blocker_set for vertex in six_vertices
        )
        if degree:
            ridge_degrees.append(degree)
    histogram = Counter(ridge_degrees)
    if histogram != {2: 594}:
        raise ValueError("a semantic target violates the ridge-degree gate")
    return {
        "A5_masks": list(normalized),
        "A5_family_sha256": canonical_sha(normalized),
        "minimum_local_cap_slack": min(local_slacks),
        "surviving_small_transversal_count": 0,
        "blocker_count": len(blocker),
        "blocker_sha256": sha256(
            b"".join(value.to_bytes(2, "little") for value in blocker)
        ).hexdigest(),
        "ridge_degree_histogram": {"2": 594},
    }


def build(
    targets: Sequence[Sequence[int]], *, seed_stabilizer_lex: bool = False,
    exact_face_vector: bool = False,
) -> tuple[CNF, dict[str, object]]:
    rows = tuple(
        value
        for value in range(1 << N)
        if value.bit_count() == 5 and independent(value)
    )
    facets = tuple(
        value
        for value in range(1 << N)
        if value.bit_count() == 5 and hits_edges(value)
    )
    if (len(rows), len(facets)) != (1372, 326):
        raise AssertionError("Record-30 universes changed")
    row_index = {value: index for index, value in enumerate(rows)}
    facet_index = {value: index for index, value in enumerate(facets)}

    target_audits = [semantic_target_audit(target) for target in targets]

    formula = CNF()
    row_variables = tuple(formula.variable() for _ in rows)
    facet_variables = tuple(formula.variable() for _ in facets)
    formula.exactly(row_variables, ROW_COUNT)

    local_constraint_count = 0
    for size in (1, 2, 3):
        for chosen in combinations(range(N), size):
            subset = sum(1 << vertex for vertex in chosen)
            if not independent(subset):
                continue
            containing = tuple(
                row_variables[index]
                for index, row in enumerate(rows)
                if row & subset == subset
            )
            formula.at_most(containing, local_cap(subset))
            local_constraint_count += 1

    small_transversal_count = 0
    for size in (3, 4):
        for chosen in combinations(range(N), size):
            transversal = sum(1 << vertex for vertex in chosen)
            if not hits_edges(transversal):
                continue
            disjoint = tuple(
                row_variables[index]
                for index, row in enumerate(rows)
                if not row & transversal
            )
            formula.clause(*disjoint)
            small_transversal_count += 1

    blocker_binary_clause_count = 0
    for index, complement in enumerate(facets):
        facet = facet_variables[index]
        disjoint = tuple(
            row_variables[row]
            for row, support in enumerate(rows)
            if not support & complement
        )
        for selected in disjoint:
            formula.clause(-facet, -selected)
            blocker_binary_clause_count += 1
        formula.clause(facet, *disjoint)
    formula.exactly(facet_variables, FACET_COUNT)

    face_vector_indicator_counts: dict[str, int] = {}
    if exact_face_vector:
        # An independent six/seven-set is absent precisely when it contains
        # at least one selected A5.  The target f-vector therefore requires
        # respectively 56 and 160 absent independent sets.
        for size, missing_target in ((6, 56), (7, 160)):
            missing_variables = []
            for chosen in combinations(range(N), size):
                support = sum(1 << vertex for vertex in chosen)
                if not independent(support):
                    continue
                contained = tuple(
                    row_variables[index]
                    for index, row in enumerate(rows)
                    if row & support == row
                )
                missing = formula.variable()
                for selected in contained:
                    formula.clause(-selected, missing)
                formula.clause(-missing, *contained)
                missing_variables.append(missing)
            formula.exactly(missing_variables, missing_target)
            face_vector_indicator_counts[str(size)] = len(missing_variables)

    ridge_group_count = 0
    for six_vertices in combinations(range(N), 6):
        six = sum(1 << vertex for vertex in six_vertices)
        members = tuple(
            facet_variables[facet_index[six ^ (1 << vertex)]]
            for vertex in six_vertices
            if six ^ (1 << vertex) in facet_index
        )
        if not members:
            continue
        ridge_group_count += 1
        for triple in combinations(members, 3):
            formula.clause(*(-item for item in triple))
        for item in members:
            formula.clause(-item, *(other for other in members if other != item))

    formula.clause(facet_variables[facet_index[SEED]])

    stabilizer_generators = seed_stabilizer_generators()
    if seed_stabilizer_lex:
        for action in stabilizer_generators:
            image: list[int | None] = [None] * len(rows)
            for source, row in enumerate(rows):
                target = row_index[transform(row, action)]
                image[target] = row_variables[source]
            if any(item is None for item in image):
                raise AssertionError("seed stabilizer does not permute A5 rows")
            formula.lex_greater_equal(
                row_variables, tuple(int(item) for item in image)
            )

    actions = active_actions()
    blocked: set[tuple[int, ...]] = set()
    per_target = []
    for target_index, target in enumerate(targets):
        forms = seed_forms(target, actions)
        overlap = blocked & forms
        if overlap:
            raise ValueError(
                f"target {target_index} overlaps an earlier full-group orbit"
            )
        blocked.update(forms)
        per_target.append(len(forms))
    for family in sorted(blocked):
        formula.clause(*(-row_variables[row_index[row]] for row in family))

    metadata = {
        "kind": "h9_14_record30_a0_dual_blocked_cnf",
        "schema_version": "1.0.0",
        "row_universe_count": len(rows),
        "facet_complement_universe_count": len(facets),
        "selected_A5_count": ROW_COUNT,
        "selected_facet_complement_count": FACET_COUNT,
        "local_constraint_count": local_constraint_count,
        "small_transversal_constraint_count": small_transversal_count,
        "blocker_binary_clause_count": blocker_binary_clause_count,
        "ridge_group_count": ridge_group_count,
        "canonical_seed": list(vertices(SEED)),
        "active_group_order": len(actions),
        "isolate_group_order": 40320,
        "full_group_order": len(actions) * 40320,
        "seed_stabilizer_order": 8640,
        "seed_stabilizer_generator_count": len(stabilizer_generators),
        "seed_stabilizer_lex_encoded": seed_stabilizer_lex,
        "input_target_orbit_count": len(targets),
        "input_target_semantic_audits": target_audits,
        "seed_form_counts": per_target,
        "blocked_seed_form_count": len(blocked),
        "blocked_seed_form_families_sha256": canonical_sha(sorted(blocked)),
        "variable_count": formula.variable_count,
        "clause_count": len(formula.clauses),
        "connectivity_encoded": False,
        "face_vector_encoded": exact_face_vector,
        "face_vector_missing_indicator_counts": face_vector_indicator_counts,
        "stronger_than_needed_if_unsat": (
            "The formula omits connectivity; UNSAT excludes every target in "
            "this larger necessary class."
            if exact_face_vector
            else "The formula omits connectivity and the unforced face-vector "
            "entries; UNSAT excludes every target in this larger necessary class."
        ),
        "finite_volume_noncompact_statement": False,
        "scope": (
            "Dual compact necessary-incidence exhaustion only. An UNSAT solver "
            "status is theorem-grade only together with a checked proof object."
        ),
    }
    return formula, metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--census", type=Path, default=DEFAULT_CENSUS)
    parser.add_argument("--cnf", type=Path, required=True)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--seed-stabilizer-lex", action="store_true")
    parser.add_argument("--exact-face-vector", action="store_true")
    args = parser.parse_args()
    targets = normalize_targets(json.loads(args.census.read_text(encoding="utf-8")))
    formula, metadata = build(
        targets,
        seed_stabilizer_lex=args.seed_stabilizer_lex,
        exact_face_vector=args.exact_face_vector,
    )
    metadata["cnf_sha256"] = formula.write(args.cnf)
    rendered = json.dumps(metadata, indent=2, sort_keys=True) + "\n"
    if args.metadata is not None:
        args.metadata.write_text(rendered, encoding="utf-8")
    sys.stdout.write(rendered)


if __name__ == "__main__":
    main()

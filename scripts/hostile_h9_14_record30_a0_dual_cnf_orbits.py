#!/usr/bin/env python3
"""No-project-import orbit audit for the Record-30 ``a=0`` dual CNF.

This checker intentionally does not import the CNF generator.  It constructs
``Aut(3K2) x S8`` from elementary permutation generators, closes the full
orbit of every supplied eight-row family, and filters the orbit by literal
incidence with the canonical type-three seed.  It then checks that the tail
of a supplied DIMACS file consists *exactly* of the corresponding negative
eight-row blocking clauses.

The independent orbit-stabilizer calculation is also recorded.  It provides
a second route to every seed-form count: if an orbit contains ``o`` families
and a representative has ``t`` type-three blocker elements, then the number
of orbit members containing the fixed seed is ``o*t/224``.
"""

from __future__ import annotations

import argparse
from collections import deque
from hashlib import sha256
import json
from pathlib import Path
from typing import Iterable, Sequence


N = 14
FULL_GROUP_ORDER = 1_935_360
TYPE_THREE_ORBIT_SIZE = 224
SEED = sum(1 << vertex for vertex in (0, 1, 2, 3, 4))
EDGES = ((0, 11), (1, 12), (2, 13))
ACTIVE = frozenset(vertex for edge in EDGES for vertex in edge)
ISOLATES = tuple(range(3, 11))
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CENSUS = (
    ROOT / "data/enumerations/h9_14_record30_a0_connected_ridge_census.json"
)


Family = tuple[int, ...]
Action = tuple[int, ...]


def canonical_sha(value: object) -> str:
    rendered = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    )
    return sha256(rendered.encode("ascii")).hexdigest()


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def cycle_action(*cycles: Sequence[int]) -> Action:
    action = list(range(N))
    for cycle in cycles:
        if len(set(cycle)) != len(cycle):
            raise ValueError("a permutation cycle repeats a vertex")
        for source, target in zip(cycle, cycle[1:] + cycle[:1], strict=True):
            action[source] = target
    if len(set(action)) != N:
        raise ValueError("cycles do not define a permutation")
    return tuple(action)


def compose(left: Action, right: Action) -> Action:
    """Return ``left after right``."""
    return tuple(left[right[vertex]] for vertex in range(N))


def action_closure(generators: Sequence[Action]) -> set[Action]:
    identity = tuple(range(N))
    reached = {identity}
    queue = [identity]
    for current in queue:
        for generator in generators:
            moved = compose(generator, current)
            if moved not in reached:
                reached.add(moved)
                queue.append(moved)
    return reached


def full_group_generators() -> tuple[Action, ...]:
    # One matching-edge flip and the two standard generators of S3 on the
    # matching edges generate C2 wr S3.  Adjacent isolate swaps generate S8.
    active = (
        cycle_action((0, 11)),
        cycle_action((0, 1), (11, 12)),
        cycle_action((0, 1, 2), (11, 12, 13)),
    )
    isolates = tuple(cycle_action((left, left + 1)) for left in range(3, 10))
    if len(action_closure(active)) != 48:
        raise AssertionError("active generators do not generate Aut(3K2)")
    if len(action_closure(isolates)) != 40_320:
        raise AssertionError("isolate generators do not generate S8")
    return active + isolates


def seed_stabilizer_generators() -> tuple[Action, ...]:
    # The seed chooses one endpoint from each matching edge and two isolates.
    # Its stabilizer is S3 x S2 x S6, of order 8640.
    return (
        cycle_action((0, 1), (11, 12)),
        cycle_action((0, 1, 2), (11, 12, 13)),
        cycle_action((3, 4)),
        *(cycle_action((left, left + 1)) for left in range(5, 10)),
    )


def transform_mask(mask: int, action: Action) -> int:
    return sum(
        1 << action[vertex] for vertex in range(N) if mask >> vertex & 1
    )


def transform_family(family: Sequence[int], action: Action) -> Family:
    return tuple(sorted(transform_mask(row, action) for row in family))


def family_orbit(family: Family, generators: Sequence[Action]) -> set[Family]:
    reached = {family}
    queue = [family]
    for current in queue:
        for generator in generators:
            moved = transform_family(current, generator)
            if moved not in reached:
                reached.add(moved)
                queue.append(moved)
    return reached


def hits_matching(mask: int) -> bool:
    return all(mask >> left & 1 or mask >> right & 1 for left, right in EDGES)


def independent(mask: int) -> bool:
    return all(not (mask >> left & 1 and mask >> right & 1) for left, right in EDGES)


def normalize_targets(document: object) -> tuple[Family, ...]:
    if not isinstance(document, dict):
        raise ValueError("census must be a JSON object")
    raw_targets = document.get("target_orbit_representatives")
    if not isinstance(raw_targets, list) or document.get("target_orbit_count") != len(raw_targets):
        raise ValueError("census target list/count is malformed")
    targets = []
    for target_index, raw_family in enumerate(raw_targets):
        if not isinstance(raw_family, list) or len(raw_family) != 8:
            raise ValueError(f"target {target_index} does not have eight rows")
        rows = []
        for raw_row in raw_family:
            if not isinstance(raw_row, list) or len(raw_row) != 5:
                raise ValueError(f"target {target_index} has a malformed row")
            if len(set(raw_row)) != 5 or any(
                not isinstance(vertex, int) or not 0 <= vertex < N
                for vertex in raw_row
            ):
                raise ValueError(f"target {target_index} has an invalid row")
            mask = sum(1 << vertex for vertex in raw_row)
            if not independent(mask):
                raise ValueError(f"target {target_index} row contains a missing edge")
            rows.append(mask)
        family = tuple(sorted(rows))
        if len(set(family)) != 8:
            raise ValueError(f"target {target_index} repeats a row")
        targets.append(family)
    if len(set(targets)) != len(targets):
        raise ValueError("census repeats a target literally")
    return tuple(targets)


def type_three_blocker_count(family: Family) -> int:
    count = 0
    for left in EDGES[0]:
        for middle in EDGES[1]:
            for right in EDGES[2]:
                for first_position in range(len(ISOLATES)):
                    for second_position in range(first_position + 1, len(ISOLATES)):
                        complement = sum(
                            1 << vertex
                            for vertex in (
                                left,
                                middle,
                                right,
                                ISOLATES[first_position],
                                ISOLATES[second_position],
                            )
                        )
                        if all(complement & row for row in family):
                            count += 1
    return count


def characteristic_lex_ge(left: Family, right: Family) -> bool:
    difference = set(left) ^ set(right)
    return not difference or min(difference) in set(left)


def lex_maximum(families: Iterable[Family]) -> Family:
    # For False < True, the first characteristic-vector difference favours
    # the set containing the smallest differing row index.
    iterator = iter(families)
    best = next(iterator)
    for family in iterator:
        if characteristic_lex_ge(family, best):
            best = family
    return best


def parse_clause(line: str) -> tuple[int, ...]:
    raw = tuple(map(int, line.split()))
    if not raw or raw[-1] != 0 or any(literal == 0 for literal in raw[:-1]):
        raise ValueError("malformed DIMACS clause")
    return raw[:-1]


def verify_cnf_tail(
    path: Path, expected: set[tuple[int, ...]]
) -> tuple[int, int, str]:
    with path.open("r", encoding="ascii") as stream:
        header = stream.readline().split()
        if len(header) != 4 or header[:2] != ["p", "cnf"]:
            raise ValueError("malformed DIMACS header")
        variable_count, declared_clause_count = map(int, header[2:])
        tail: deque[tuple[int, ...]] = deque(maxlen=len(expected))
        actual_clause_count = 0
        for line in stream:
            if not line.strip() or line.startswith("c"):
                continue
            tail.append(parse_clause(line))
            actual_clause_count += 1
    if actual_clause_count != declared_clause_count:
        raise ValueError("DIMACS clause count does not match its header")
    if len(tail) != len(expected) or len(set(tail)) != len(tail):
        raise ValueError("DIMACS blocking tail is truncated or duplicated")
    if set(tail) != expected:
        raise ValueError(
            "DIMACS blocking tail differs from the independently closed orbits"
        )
    return variable_count, declared_clause_count, file_sha256(path)


def audit(census: Path, cnf: Path, metadata: Path | None) -> dict[str, object]:
    targets = normalize_targets(json.loads(census.read_text(encoding="utf-8")))
    generators = full_group_generators()
    stabilizer_generators = seed_stabilizer_generators()
    if len(action_closure(stabilizer_generators)) != 8_640:
        raise AssertionError("seed-stabilizer generators have the wrong order")
    if any(transform_mask(SEED, action) != SEED for action in stabilizer_generators):
        raise AssertionError("a seed-stabilizer generator moves the seed")

    all_orbits: set[Family] = set()
    all_seed_forms: set[Family] = set()
    target_audits = []
    for target_index, target in enumerate(targets):
        orbit = family_orbit(target, generators)
        overlap = all_orbits & orbit
        if overlap:
            raise ValueError(f"target {target_index} overlaps an earlier orbit")
        all_orbits.update(orbit)
        seed_forms = {
            family for family in orbit if all(SEED & row for row in family)
        }
        if all_seed_forms & seed_forms:
            raise ValueError(f"target {target_index} repeats an earlier seed form")
        all_seed_forms.update(seed_forms)

        type_three = type_three_blocker_count(target)
        if FULL_GROUP_ORDER % len(orbit):
            raise ArithmeticError("family orbit size does not divide group order")
        stabilizer_order = FULL_GROUP_ORDER // len(orbit)
        numerator = len(orbit) * type_three
        if numerator % TYPE_THREE_ORBIT_SIZE:
            raise ArithmeticError("orbit-incidence double count is not integral")
        double_count = numerator // TYPE_THREE_ORBIT_SIZE
        if double_count != len(seed_forms):
            raise ArithmeticError("direct and double-counted seed forms disagree")

        maximum = lex_maximum(seed_forms)
        if not all(
            characteristic_lex_ge(
                maximum, transform_family(maximum, generator)
            )
            for generator in stabilizer_generators
        ):
            raise AssertionError("lex maximum fails a seed-stabilizer generator")
        target_audits.append(
            {
                "target_index": target_index,
                "full_family_orbit_size": len(orbit),
                "family_stabilizer_order": stabilizer_order,
                "type_three_blocker_count": type_three,
                "direct_seed_form_count": len(seed_forms),
                "orbit_incidence_seed_form_count": double_count,
                "lex_admissible_witness_A5_masks": list(maximum),
            }
        )

    rows = tuple(
        mask
        for mask in range(1 << N)
        if mask.bit_count() == 5 and independent(mask)
    )
    if len(rows) != 1_372:
        raise AssertionError("independent five-set universe changed")
    row_index = {row: index + 1 for index, row in enumerate(rows)}
    expected_clauses = {
        tuple(-row_index[row] for row in family) for family in all_seed_forms
    }
    if len(expected_clauses) != len(all_seed_forms):
        raise AssertionError("two seed forms map to the same DIMACS clause")
    variable_count, clause_count, cnf_sha = verify_cnf_tail(cnf, expected_clauses)

    family_hash = canonical_sha(sorted(all_seed_forms))
    clause_hash = canonical_sha(sorted(expected_clauses))
    if metadata is not None:
        declared = json.loads(metadata.read_text(encoding="utf-8"))
        expected_counts = [audit["direct_seed_form_count"] for audit in target_audits]
        if declared.get("seed_form_counts") != expected_counts:
            raise ValueError("metadata seed-form counts disagree")
        if declared.get("blocked_seed_form_count") != len(all_seed_forms):
            raise ValueError("metadata blocked-family count disagrees")
        if declared.get("blocked_seed_forms_sha256") != family_hash:
            raise ValueError("metadata blocked-family digest disagrees")
        if declared.get("cnf_sha256") != cnf_sha:
            raise ValueError("metadata CNF digest disagrees")
        if declared.get("variable_count") != variable_count:
            raise ValueError("metadata variable count disagrees")
        if declared.get("clause_count") != clause_count:
            raise ValueError("metadata clause count disagrees")

    result = {
        "kind": "h9_14_record30_a0_dual_cnf_no_import_orbit_audit",
        "schema_version": "1.0.0",
        "census_sha256": file_sha256(census),
        "cnf_sha256": cnf_sha,
        "cnf_variable_count": variable_count,
        "cnf_clause_count": clause_count,
        "full_group_order": FULL_GROUP_ORDER,
        "type_three_complement_orbit_size": TYPE_THREE_ORBIT_SIZE,
        "input_target_orbit_count": len(targets),
        "target_orbit_audits": target_audits,
        "blocked_seed_form_count": len(all_seed_forms),
        "blocked_seed_families_sha256": family_hash,
        "blocking_clauses_sha256": clause_hash,
        "all_target_orbits_pairwise_disjoint": True,
        "every_seed_count_matches_orbit_incidence_double_count": True,
        "every_orbit_has_a_seed_stabilizer_lex_witness": True,
        "cnf_tail_is_exactly_the_independent_blocking_clause_set": True,
        "scope": (
            "Structurally independent orbit/block-clause replay only; this "
            "script does not validate the SAT solver proof."
        ),
    }
    result["summary_sha256"] = canonical_sha(result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--census", type=Path, default=DEFAULT_CENSUS)
    parser.add_argument("--cnf", type=Path, required=True)
    parser.add_argument("--metadata", type=Path)
    args = parser.parse_args()
    print(
        json.dumps(
            audit(args.census, args.cnf, args.metadata),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()

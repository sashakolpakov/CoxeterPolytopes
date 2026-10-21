#!/usr/bin/env python3
"""Hostile replay for the Record-30 ``a=0`` dual CNF and proof.

The verifier independently reconstructs all seed-containing images of the
input target orbits, checks that the final DIMACS clauses are *exactly* their
eight-negative-literal blockers, removes precisely those clauses, and asks an
external SAT solver for a positive witness from every known orbit.  Each
witness is decoded back to its eight A5s and 132-element blocker and every CNF
clause is evaluated literally.  Optionally it also invokes an external DRAT
checker on the original blocked CNF/proof pair.
"""

from __future__ import annotations

import argparse
from collections import deque
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import Sequence

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.independent_h9_14_record30_a0_dual_cnf import (
    EDGE_MASKS,
    DEFAULT_CENSUS,
    ISOLATES,
    N,
    ROW_COUNT,
    canonical_sha,
    hits_edges,
    independent,
    normalize_targets,
    semantic_target_audit,
)

ACTIVE_LABELS = (0, 1, 2, 11, 12, 13)
SEED = sum(1 << vertex for vertex in (0, 1, 2, 3, 4))


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def parse_clause(line: str) -> tuple[int, ...]:
    values = tuple(map(int, line.split()))
    if not values or values[-1] != 0 or any(value == 0 for value in values[:-1]):
        raise ValueError("malformed DIMACS clause")
    return values[:-1]


def universes() -> tuple[tuple[int, ...], tuple[int, ...]]:
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
    return rows, facets


def _transform(value: int, action: Sequence[int]) -> int:
    return sum(
        1 << action[vertex] for vertex in range(N) if value >> vertex & 1
    )


def _brute_active_actions() -> tuple[tuple[int, ...], ...]:
    """Independently scan all 6! active permutations preserving the matching."""
    expected_edges = set(EDGE_MASKS)
    result = []
    for images in permutations(ACTIVE_LABELS):
        action = list(range(N))
        for source, target in zip(ACTIVE_LABELS, images, strict=True):
            action[source] = target
        moved_edges = {_transform(edge, action) for edge in EDGE_MASKS}
        if moved_edges == expected_edges:
            result.append(tuple(action))
    if len(result) != 48 or len(set(result)) != 48:
        raise AssertionError("independent active-group scan did not find 48 actions")
    return tuple(result)


def _stabilizer_generators() -> tuple[tuple[int, ...], ...]:
    def moved(*cycles: tuple[int, ...]) -> tuple[int, ...]:
        action = list(range(N))
        for cycle in cycles:
            for source, target in zip(
                cycle, cycle[1:] + cycle[:1], strict=True
            ):
                action[source] = target
        if _transform(SEED, action) != SEED:
            raise AssertionError("independent stabilizer generator moves seed")
        return tuple(action)

    return (
        moved((0, 1), (11, 12)),
        moved((0, 1, 2), (11, 12, 13)),
        moved((3, 4)),
        *(moved((left, left + 1)) for left in range(5, 10)),
    )


def _lex_admissible(family: Sequence[int]) -> bool:
    selected = set(family)
    for action in _stabilizer_generators():
        moved = {_transform(row, action) for row in family}
        difference = selected ^ moved
        if difference and min(difference) not in selected:
            return False
    return True


def _independent_orbit_forms(
    family: Sequence[int], active_actions: Sequence[Sequence[int]]
) -> tuple[set[tuple[int, ...]], set[tuple[int, ...]]]:
    """Literal full-group traversal, separate from the CNF generator helpers."""
    active_mask = sum(1 << vertex for vertex in ACTIVE_LABELS)
    active_versions = tuple(
        tuple(_transform(row & active_mask, action) for row in family)
        for action in active_actions
    )
    isolate_patterns = tuple(
        sum(
            ((row >> vertex) & 1) << position
            for position, vertex in enumerate(ISOLATES)
        )
        for row in family
    )
    full_orbit: set[tuple[int, ...]] = set()
    with_seed: set[tuple[int, ...]] = set()
    for isolate_images in permutations(ISOLATES):
        isolate_parts = tuple(
            sum(
                ((pattern >> position) & 1) << isolate_images[position]
                for position in range(8)
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
            full_orbit.add(moved)
            if all(SEED & row for row in moved):
                with_seed.add(moved)
    return full_orbit, with_seed


def expected_blocks(
    targets: Sequence[Sequence[int]], rows: Sequence[int]
) -> tuple[
    set[tuple[int, ...]],
    set[tuple[int, ...]],
    tuple[int, ...],
    tuple[tuple[int, ...], ...],
    tuple[dict[str, int], ...],
]:
    row_index = {value: index for index, value in enumerate(rows)}
    actions = _brute_active_actions()
    all_forms: set[tuple[int, ...]] = set()
    counts = []
    positive_forms = []
    orbit_audits = []
    for target_index, target in enumerate(targets):
        full_orbit, forms = _independent_orbit_forms(target, actions)
        if all_forms & forms:
            raise ValueError(
                f"target {target_index} overlaps a preceding full-group orbit"
            )
        all_forms.update(forms)
        counts.append(len(forms))
        admissible = sorted(
            family for family in forms if _lex_admissible(family)
        )
        if not admissible:
            raise ValueError(
                f"target {target_index} has no seed-stabilizer lex image"
            )
        positive_forms.append(admissible[0])
        blocker = tuple(
            complement
            for complement in range(1 << N)
            if complement.bit_count() == 5
            and hits_edges(complement)
            and all(complement & row for row in target)
        )
        type_three = sum(
            (complement & sum(1 << v for v in ACTIVE_LABELS)).bit_count() == 3
            for complement in blocker
        )
        orbit_size = len(full_orbit)
        if 1935360 % orbit_size:
            raise AssertionError("target orbit size does not divide full group order")
        if len(forms) * 224 != orbit_size * type_three:
            raise AssertionError("seed-form double-count identity failed")
        orbit_audits.append(
            {
                "full_orbit_size": orbit_size,
                "automorphism_order": 1935360 // orbit_size,
                "type_three_blocker_count": type_three,
                "seed_form_count": len(forms),
                "seed_double_count": len(forms) * 224,
            }
        )
    clauses = {
        tuple(-(row_index[row] + 1) for row in family) for family in all_forms
    }
    if len(clauses) != len(all_forms):
        raise AssertionError("two orbit images produced the same blocking clause")
    return (
        clauses,
        all_forms,
        tuple(counts),
        tuple(positive_forms),
        tuple(orbit_audits),
    )


def split_and_verify_blocks(
    cnf: Path,
    expected: set[tuple[int, ...]],
    base_path: Path,
) -> tuple[int, int]:
    with cnf.open("r", encoding="ascii") as source:
        header = source.readline().split()
        if len(header) != 4 or header[:2] != ["p", "cnf"]:
            raise ValueError("malformed DIMACS header")
        variable_count, clause_count = map(int, header[2:])
        base_count = clause_count - len(expected)
        if base_count <= 0:
            raise ValueError("blocking clauses exceed the formula")
        tail: deque[tuple[int, ...]] = deque(maxlen=len(expected))
        parsed = 0
        with base_path.open("w", encoding="ascii") as base:
            base.write(f"p cnf {variable_count} {base_count}\n")
            for line in source:
                if not line.strip() or line.startswith("c"):
                    continue
                clause = parse_clause(line)
                if parsed < base_count:
                    base.write(" ".join(map(str, clause)) + " 0\n")
                tail.append(clause)
                parsed += 1
    if parsed != clause_count:
        raise ValueError("DIMACS header clause count is false")
    if len(tail) != len(expected) or set(tail) != expected:
        missing = expected - set(tail)
        extra = set(tail) - expected
        raise ValueError(
            f"blocking-clause mismatch: missing={len(missing)} extra={len(extra)}"
        )
    if len(set(tail)) != len(tail):
        raise ValueError("DIMACS repeats a target blocking clause")
    if any(len(clause) != ROW_COUNT or any(literal >= 0 for literal in clause) for clause in tail):
        raise ValueError("a target blocker is not eight negative literals")
    return variable_count, base_count


def write_forced_formula(
    base: Path,
    output: Path,
    variable_count: int,
    base_count: int,
    unit_variables: Sequence[int],
) -> None:
    with base.open("r", encoding="ascii") as source, output.open(
        "w", encoding="ascii"
    ) as target:
        header = source.readline().split()
        if header != ["p", "cnf", str(variable_count), str(base_count)]:
            raise ValueError("temporary base header changed")
        target.write(
            f"p cnf {variable_count} {base_count + len(unit_variables)}\n"
        )
        for line in source:
            target.write(line)
        for variable in unit_variables:
            target.write(f"{variable} 0\n")


def parse_witness(output: str, variable_count: int) -> tuple[bool, ...]:
    literals = []
    for line in output.splitlines():
        if line.startswith("v "):
            literals.extend(map(int, line.split()[1:]))
    literals = [literal for literal in literals if literal]
    assignment: list[bool | None] = [None] * (variable_count + 1)
    for literal in literals:
        variable = abs(literal)
        if variable > variable_count:
            raise ValueError("SAT witness assigns an out-of-range variable")
        value = literal > 0
        if assignment[variable] is not None and assignment[variable] != value:
            raise ValueError("SAT witness contradicts itself")
        assignment[variable] = value
    if any(value is None for value in assignment[1:]):
        missing = sum(value is None for value in assignment[1:])
        raise ValueError(f"SAT witness omits {missing} variables")
    return tuple(bool(value) for value in assignment)


def evaluate_dimacs(path: Path, assignment: Sequence[bool]) -> int:
    with path.open("r", encoding="ascii") as stream:
        header = stream.readline().split()
        clause_count = int(header[3])
        checked = 0
        for line in stream:
            if not line.strip() or line.startswith("c"):
                continue
            clause = parse_clause(line)
            if not any(
                assignment[abs(literal)] == (literal > 0)
                for literal in clause
            ):
                raise ValueError(f"decoded SAT witness falsifies clause {checked}")
            checked += 1
    if checked != clause_count:
        raise ValueError("literal witness replay saw the wrong clause count")
    return checked


def positive_witness_audits(
    solver: Path,
    base: Path,
    variable_count: int,
    base_count: int,
    rows: Sequence[int],
    facets: Sequence[int],
    families: Sequence[Sequence[int]],
    temporary: Path,
) -> list[dict[str, object]]:
    row_index = {value: index for index, value in enumerate(rows)}
    audits = []
    for target_index, family in enumerate(families):
        forced = temporary / f"target_{target_index}.cnf"
        unit_variables = tuple(row_index[row] + 1 for row in family)
        write_forced_formula(
            base, forced, variable_count, base_count, unit_variables
        )
        process = subprocess.run(
            [str(solver), "--no-colors", str(forced)],
            text=True,
            capture_output=True,
            check=False,
        )
        if process.returncode != 10 or "s SATISFIABLE" not in process.stdout:
            raise ValueError(
                f"positive target {target_index} did not solve SAT: "
                f"returncode={process.returncode}"
            )
        assignment = parse_witness(process.stdout, variable_count)
        checked = evaluate_dimacs(forced, assignment)
        selected_rows = tuple(
            row for index, row in enumerate(rows) if assignment[index + 1]
        )
        if selected_rows != tuple(family):
            raise ValueError(f"positive target {target_index} decoded incorrectly")
        facet_offset = len(rows) + 1
        selected_facets = tuple(
            complement
            for index, complement in enumerate(facets)
            if assignment[facet_offset + index]
        )
        expected_facets = tuple(
            complement
            for complement in facets
            if all(complement & row for row in family)
        )
        if selected_facets != expected_facets or len(selected_facets) != 132:
            raise ValueError(
                f"positive target {target_index} blocker decoded incorrectly"
            )
        audits.append(
            {
                "target_index": target_index,
                "forced_A5_masks": list(family),
                "forced_A5_family_sha256": canonical_sha(family),
                "decoded_selected_A5_count": len(selected_rows),
                "decoded_selected_facet_count": len(selected_facets),
                "decoded_blocker_sha256": sha256(
                    b"".join(
                        value.to_bytes(2, "little") for value in selected_facets
                    )
                ).hexdigest(),
                "literal_clause_count_checked": checked,
                "solver_returncode": process.returncode,
            }
        )
    return audits


def run_drat_checker(checker: Path, cnf: Path, proof: Path) -> dict[str, object]:
    process = subprocess.run(
        [str(checker), str(cnf), str(proof), "-i"],
        text=True,
        capture_output=True,
        check=False,
    )
    transcript = process.stdout + process.stderr
    if process.returncode != 0 or "VERIFIED" not in transcript:
        raise ValueError(
            f"DRAT verification failed with return code {process.returncode}"
        )
    return {
        "checker_returncode": process.returncode,
        "checker_transcript_sha256": sha256(transcript.encode()).hexdigest(),
        "checker_reported_verified": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--census", type=Path, default=DEFAULT_CENSUS)
    parser.add_argument("--cnf", type=Path, required=True)
    parser.add_argument("--proof", type=Path)
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--drat-trim", type=Path)
    args = parser.parse_args()

    document = json.loads(args.census.read_text(encoding="utf-8"))
    targets = normalize_targets(document)
    semantic_audits = [semantic_target_audit(target) for target in targets]
    rows, facets = universes()
    expected, form_families, counts, positive_forms, orbit_audits = (
        expected_blocks(targets, rows)
    )
    with tempfile.TemporaryDirectory(prefix="record30_a0_dual_hostile_") as raw:
        temporary = Path(raw)
        base = temporary / "unblocked.cnf"
        variable_count, base_count = split_and_verify_blocks(
            args.cnf, expected, base
        )
        positives = positive_witness_audits(
            args.solver,
            base,
            variable_count,
            base_count,
            rows,
            facets,
            positive_forms,
            temporary,
        )

    proof_audit = None
    if args.proof is not None:
        if args.drat_trim is None:
            raise ValueError("--proof requires --drat-trim")
        proof_audit = run_drat_checker(args.drat_trim, args.cnf, args.proof)
        proof_audit["proof_sha256"] = file_sha256(args.proof)
        proof_audit["proof_size_bytes"] = args.proof.stat().st_size

    result = {
        "kind": "h9_14_record30_a0_dual_cnf_hostile_replay",
        "schema_version": "1.0.0",
        "cnf_sha256": file_sha256(args.cnf),
        "cnf_size_bytes": args.cnf.stat().st_size,
        "variable_count": variable_count,
        "clause_count": base_count + len(expected),
        "base_unblocked_clause_count": base_count,
        "target_orbit_count": len(targets),
        "target_semantic_audits": semantic_audits,
        "seed_form_counts": list(counts),
        "blocked_seed_form_count": len(expected),
        "blocked_seed_form_families_sha256": canonical_sha(
            sorted(form_families)
        ),
        "blocking_clauses_sha256": canonical_sha(sorted(expected)),
        "all_blocking_clauses_exactly_reconstructed": True,
        "independent_orbit_stabilizer_audits": list(orbit_audits),
        "imports_generator_group_or_seed_form_helpers": False,
        "positive_target_witness_count": len(positives),
        "positive_target_witness_audits": positives,
        "all_positive_witnesses_literally_satisfy_unblocked_cnf": True,
        "drat_proof_audit": proof_audit,
        "finite_volume_noncompact_statement": False,
        "scope": (
            "Hostile proof replay for compact necessary-incidence exhaustion; "
            "no ideal-face or finite-volume noncompact inference."
        ),
    }
    result["summary_sha256"] = canonical_sha(result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

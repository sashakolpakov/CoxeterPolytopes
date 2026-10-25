#!/usr/bin/env python3
"""Prepared structural regressions for the supplemental finite57 tail guard.

SOURCE ONLY until coordinating and different-author review, then admitted
execution. This driver is not executed by being prepared. It imports the
new key guard in ONE process, never a legacy mathematical worker. It uses
one owned temporary fixture directory, restores each changed input and
removes the directory on success or failure. No enumeration, certificate
arithmetic, KKT check, optimizer, Lean, TeX or child process is invoked.

A positive key-only result does not validate numeric field associations.
Permuting whole valid packing rows is allowed; unexpected or repeated
ordered budget KEYS must fail. Source clearance binds this driver and the
guard source; this driver additionally enforces the seven fixed input pins.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import tempfile
import types

GUARD_NAME = "verify_finite57_tail_key_coverage_supplement_20261024.py"
TWO = "results/fresh_nikulin_two_budget_frontier_20261015.json"
HEIGHT = "results/fresh_nikulin_height_charge_20261015.json"
UNIT = "results/fresh_nikulin_unit_charge_20261015.json"
MAX_FILE_BYTES = 4 * 1024 * 1024
MAX_FIXTURE_BYTES = 8 * 1024 * 1024
FIXED_INPUTS = {
    TWO: "1a32a2aa37b44b73ab03a322c136e846cea5c37028718e5241cfdad779edfdc9",
    HEIGHT: "6427dadcaa1577d0416b31f48da7b67425eb647d4661ab464c2e7f55f7098b6b",
    "results/fresh_nikulin_height_charge_20261015_review.json": "bebad068dac2e18fc70830adbb3fc9613487414ad2f1b3e0061dd399e9bfd449",
    "results/fresh_nikulin_defect_frontier_20261015.json": "3c8897e52d781b97c9b5de222e276234c84b8a95770f7039d2d2b99880614d04",
    "results/fresh_nikulin_defect_frontier_20261015_review.json": "a45b3764e1ba9fdc71e3c28a03d36368e1b4df0ded8f26f4dd585a26e9c63da2",
    UNIT: "c6b85ef3c5abbcd53c9f90347f7e129ecd52285fac4672cea0b67563a67bd888",
    "results/fresh_nikulin_unit_charge_20261015_review.json": "7b1d281d52b2d2e171b2ffe1b4155ca55708e333889c4b2525af579de5f88c14",
}
DEFECT = "results/fresh_nikulin_defect_frontier_20261015.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def bounded_read(path: Path) -> bytes:
    with path.open("rb") as stream:
        raw = stream.read(MAX_FILE_BYTES + 1)
    require(len(raw) <= MAX_FILE_BYTES, "regression input/source exceeds four MiB")
    return raw


def duplicate_e7(data: dict) -> None:
    rows = data["same_row_E7_exclusions"]
    rows[:] = [deepcopy(rows[0]) for _ in range(4)]


def missing_e7(data: dict) -> None:
    data["same_row_E7_exclusions"].pop()


def alter_e7_side(data: dict) -> None:
    data["same_row_E7_exclusions"][0]["side"] = 1


def alter_e7_flag(data: dict) -> None:
    data["same_row_E7_exclusions"][0]["own_omitted"] = True


def duplicate_cycle_wrapper(data: dict) -> None:
    rows = data["classical_A"]["small_cycles"]
    rows[1] = deepcopy(rows[0])


def replace_widened_key(data: dict) -> None:
    # Keep 580 rows and a valid cycle-2 node range, but replace the required
    # first key (2,(0,),(1,)) by the unlisted key (2,(1,),(0,)).
    row = data["widened_165_checks"]["cycle_cases"][0]
    require(row["cycle_size"] == 2, "expected the frozen first cycle-2 row")
    row["left"], row["right"] = [1], [0]


def duplicate_exceptional_descriptor(data: dict) -> None:
    rows = data["exceptional_cases"]
    rows[1] = deepcopy(rows[0])


def wrong_unit_e7_nodes(data: dict) -> None:
    rows = [row for row in data["exceptional_cases"]
            if row["type"] == "E7" and row["forbidden"] is None]
    require(len(rows) == 1, "expected one ordinary E7 descriptor")
    rows[0]["two_contact_case"]["nodes"] = [[0, 4], [1, 6]]


def unexpected_reversed_budget(data: dict) -> None:
    rows = data["ordered_budget_packing"]
    require(rows[0]["endpoint_budgets"] == ["10/9", "4/3"], "expected frozen first budget")
    rows[0]["endpoint_budgets"] = ["4/3", "10/9"]


def repeated_budget(data: dict) -> None:
    rows = data["ordered_budget_packing"]
    rows[1]["endpoint_budgets"] = deepcopy(rows[0]["endpoint_budgets"])


def missing_tail(data: dict) -> None:
    data["all_rank_tails"]["BCD"].pop()


def permute_valid_packing_rows(data: dict) -> None:
    data["ordered_budget_packing"].reverse()


def run(root: Path, fixture_parent: Path) -> dict:
    # The only imported repository module is the new structural guard.
    # No mathematical project script or subprocess is loaded.
    guard_path = Path(__file__).resolve().with_name(GUARD_NAME)
    guard_bytes = bounded_read(guard_path)
    guard_digest = hashlib.sha256(guard_bytes).hexdigest()
    guard = types.ModuleType("finite57_tail_key_guard_supplement")
    guard.__file__ = str(guard_path)
    # Execute exactly the bounded bytes just hashed. A normal source loader
    # would reread the path outside this cap and could create a local cache.
    exec(compile(guard_bytes, str(guard_path), "exec"), guard.__dict__)

    originals = {}
    for relative, expected in FIXED_INPUTS.items():
        raw = bounded_read(root / relative)
        require(hashlib.sha256(raw).hexdigest() == expected, "fixed input pin mismatch: " + relative)
        originals[relative] = raw
    original_total = sum(len(raw) for raw in originals.values())
    require(original_total <= MAX_FIXTURE_BYTES, "fixed fixture exceeds eight MiB")
    records = []
    fixture_parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="finite57-tail-keys-", dir=fixture_parent) as temporary:
        fixture_root = Path(temporary)

        def write_fixture(relative: str, raw: bytes) -> None:
            require(len(raw) <= MAX_FILE_BYTES, "mutated fixture file exceeds four MiB")
            require(original_total - len(originals[relative]) + len(raw) <= MAX_FIXTURE_BYTES,
                    "mutated fixture exceeds eight MiB")
            path = fixture_root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("wb") as stream:
                stream.write(raw)

        for relative, raw in originals.items():
            write_fixture(relative, raw)

        positive = guard.check(fixture_root)
        require(positive["status"] == "PASS_FINITE57_TAIL_KEY_COVERAGE_ONLY", "unexpected guard status")
        require(all(positive["sources"][path] == digest for path, digest in FIXED_INPUTS.items()),
                "positive fixture receipt does not bind all seven fixed inputs")
        require(positive["sources"]["scripts/" + GUARD_NAME] == guard_digest,
                "guard source changed after its executed bounded bytes were read")
        records.append({"case": "fixed_inputs", "result": "accepted_key_coverage_only"})

        cases = [
            ("four_duplicate_e7_rows", TWO, duplicate_e7, "four same-row E7 supports: duplicate key"),
            ("missing_e7_row", TWO, missing_e7, "four same-row E7 supports: key mismatch"),
            ("altered_e7_side", TWO, alter_e7_side, "four same-row E7 supports: key mismatch"),
            ("altered_e7_own_flag", TWO, alter_e7_flag, "four same-row E7 supports: key mismatch"),
            ("duplicate_cycle_wrapper", HEIGHT, duplicate_cycle_wrapper, "height producer cycles: duplicate key"),
            ("replaced_widened580_key", DEFECT, replace_widened_key, "all widened cycle scalar receipts: key mismatch"),
            ("duplicate_exceptional_descriptor", HEIGHT, duplicate_exceptional_descriptor, "height producer exceptional: duplicate key"),
            ("wrong_unit_e7_nodes", UNIT, wrong_unit_e7_nodes, "ordinary E7 two-contact support must be"),
            ("unexpected_reversed_ordered_budget", TWO, unexpected_reversed_budget, "ordered packing receipts: key mismatch"),
            ("repeated_ordered_budget", TWO, repeated_budget, "ordered packing receipts: duplicate key"),
            ("missing_BCD_tail_key", TWO, missing_tail, "five BCD analytic-tail receipts: key mismatch"),
        ]
        for name, relative, mutate, expected_message in cases:
            data = json.loads(originals[relative], object_pairs_hook=guard.reject_duplicate_object_keys)
            mutate(data)
            changed = json.dumps(data, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
            write_fixture(relative, changed)
            try:
                try:
                    guard.check(fixture_root)
                except guard.CoverageError as error:
                    message = str(error)
                    require(expected_message in message, name + ": rejection came from the wrong guard")
                    records.append({"case": name, "result": "required_rejection", "message": message})
                else:
                    raise ValueError(name + ": invalid key fixture was accepted")
            finally:
                # Restore the only changed file before the next case, including
                # when the expected rejection was absent or of the wrong kind.
                write_fixture(relative, originals[relative])

        # Whole-row permutations preserve all required keys and their values.
        # Reassigning numeric associations while preserving the key set is also
        # outside a key-only guard; that remains a legacy arithmetic premise.
        data = json.loads(originals[TWO], object_pairs_hook=guard.reject_duplicate_object_keys)
        permute_valid_packing_rows(data)
        changed = json.dumps(data, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
        write_fixture(TWO, changed)
        try:
            result = guard.check(fixture_root)
            require(result["status"] == "PASS_FINITE57_TAIL_KEY_COVERAGE_ONLY", "valid row permutation rejected")
            records.append({"case": "valid_whole_packing_row_permutation", "result": "accepted_key_coverage_only"})
        finally:
            write_fixture(TWO, originals[TWO])
    # TemporaryDirectory removes only the directory owned by this invocation.
    # A return is reached only after its cleanup context completed.
    source_path = Path(__file__).resolve()
    return {
        "status": "PASS_FINITE57_TAIL_KEY_GUARD_STRUCTURAL_REGRESSIONS_ONLY",
        "scope": "One process; two positive key-only cases and eleven required structural rejections. No certificate arithmetic, geometric enumeration, packing optimization, Lean, TeX or full-manuscript verification.",
        "fixed_input_pins": FIXED_INPUTS,
        "guard_sha256": guard_digest,
        "driver_sha256": hashlib.sha256(bounded_read(source_path)).hexdigest(),
        "fixture_limits_bytes": {"per_file": MAX_FILE_BYTES, "total": MAX_FIXTURE_BYTES},
        "owned_fixture_cleanup": "TemporaryDirectory context exited before receipt construction",
        "cases": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--fixture-parent", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    fixture_parent = args.fixture_parent or root / "formal" / ".proof-batches"
    print(json.dumps(run(root, fixture_parent.resolve()), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()


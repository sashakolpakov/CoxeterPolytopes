#!/usr/bin/env python3
"""Supplemental key coverage for the frozen finite57 unbounded-tail chain.

SOURCE PREPARED 2026-10-24; not executed or accepted by this file's creation.
Run only through the coordinating agent's admitted exclusive bounded runner.

This guard performs JSON structure/key checks only. It does not import any
project worker, reconstruct Gram matrices, enumerate geometric supports,
optimize, verify KKT arithmetic, or prove the finite57 theorem. Expected
small-A and inner exceptional keys are cross-bound to outputs whose complete
mathematical domain is established by the separately reviewed legacy loops.
Literal descriptor, same-row E7, packing, tail and contact-class keys are
constructed here independently of their supplied receipt rows.

Legacy producers/reviewers are intentionally unchanged. A successful run
means KEY COVERAGE ONLY and must be combined with their numerical checks.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Callable, Hashable

STEM = "fresh_nikulin_"
SUFFIX = "_20261015"
MAX_INPUT_BYTES = 32 * 1024 * 1024


class CoverageError(ValueError):
    """A required certificate key domain was not matched exactly."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CoverageError(message)


def integer(value: object) -> int:
    require(type(value) is int, "expected an integer, excluding bool")
    return value


def flag(value: object) -> bool:
    require(type(value) is bool, "expected a boolean flag")
    return value


def text_value(value: object) -> str:
    require(type(value) is str, "expected a string")
    return value


def sequence(value: object) -> list:
    require(type(value) is list, "expected a JSON array")
    return value


def nodes(value: object) -> tuple[int, ...]:
    out = tuple(integer(z) for z in sequence(value))
    require(len(set(out)) == len(out), "duplicate node inside a support")
    return out


def contact_key(row: dict) -> tuple[int, int]:
    side, node = integer(row["side"]), integer(row["node"])
    require(side in (0, 1), "contact side must be 0 or 1")
    return side, node


def support_key(row: dict) -> tuple[tuple[int, int], ...]:
    out = []
    for raw in sequence(row["nodes"]):
        pair = sequence(raw)
        require(len(pair) == 2, "a contact key must have two entries")
        side, node = (integer(z) for z in pair)
        require(side in (0, 1), "contact side must be 0 or 1")
        out.append((side, node))
    require(len(set(out)) == len(out), "duplicate contact in a support")
    return tuple(out)


def descriptor(row: dict) -> tuple[str, int | None]:
    name, forbidden = text_value(row["type"]), row["forbidden"]
    if forbidden is not None:
        forbidden = integer(forbidden)
    return name, forbidden


def index_rows(
    value: object,
    key: Callable[[dict], Hashable],
    label: str,
    expected: set[Hashable] | None = None,
) -> dict:
    result = {}
    for row in sequence(value):
        require(type(row) is dict, label + ": expected object rows")
        item = key(row)
        require(item not in result, label + ": duplicate key " + repr(item))
        result[item] = row
    if expected is not None:
        actual = set(result)
        require(
            actual == expected,
            label + ": key mismatch; missing=" + repr(sorted(expected - actual, key=repr))
            + "; extra=" + repr(sorted(actual - expected, key=repr)),
        )
    return result


def reject_duplicate_object_keys(pairs: list) -> dict:
    result = {}
    for key, value in pairs:
        require(key not in result, "duplicate JSON object key " + repr(key))
        result[key] = value
    return result


def check(root: Path) -> dict:
    sources = {}
    checks = []

    def load(name: str, review: bool = False) -> dict:
        relative = "results/" + STEM + name + SUFFIX + ("_review" if review else "") + ".json"
        path = root / relative
        require(path.stat().st_size <= MAX_INPUT_BYTES, relative + ": oversized input")
        raw = path.read_bytes()
        sources[relative] = hashlib.sha256(raw).hexdigest()
        data = json.loads(raw, object_pairs_hook=reject_duplicate_object_keys)
        require(type(data) is dict, relative + ": expected top-level object")
        return data

    two = load("two_budget_frontier")
    old = load("height_charge")
    old_review = load("height_charge", review=True)
    defect = load("defect_frontier")
    defect_review = load("defect_frontier", review=True)
    unit = load("unit_charge")
    unit_review = load("unit_charge", review=True)

    # Independently prescribed four full supports, including both node states.
    same_expected = {
        (0, False, (4, 7), (6,)),
        (0, True, (4, 7, 3), (6,)),
        (1, False, (3,), (4, 7)),
        (1, True, (3,), (4, 7, 6)),
    }
    index_rows(
        two["same_row_E7_exclusions"],
        lambda r: (integer(r["side"]), flag(r["own_omitted"]), nodes(r["left"]), nodes(r["right"])),
        "four same-row E7 supports", same_expected,
    )
    checks.append({"domain": "independent four same-row E7 full supports", "keys": 4})

    # The inner legacy numerical verifiers reconstruct these support domains.
    # Here we strictly bind all outer wrappers and widened scalar receipts to
    # those already-domain-verified witness keys; this is not a new enumeration.
    cycle_sizes = set(range(2, 18))
    cycle_tables = [
        ("height producer cycles", old["classical_A"]["small_cycles"]),
        ("height review cycles", old_review["complete_cycle_support_cases"]),
        ("defect review cycles", defect_review["widened_cycle_KKT"]["cycles"]),
    ]
    cycle_indexes = [index_rows(rows, lambda r: integer(r["cycle_size"]), label, cycle_sizes)
                     for label, rows in cycle_tables]

    def full_support(row: dict, left: str, right: str, m: int) -> tuple:
        U, V = nodes(row[left]), nodes(row[right])
        require(bool(U) and bool(V), "full supports must both be nonempty")
        require(all(0 <= k < m for k in U + V), "cycle support node outside its cycle")
        return U, V

    flat_expected = set()
    for m in sorted(cycle_sizes):
        reference = index_rows(
            cycle_indexes[0][m]["cases"],
            lambda r, m=m: full_support(r, "left_support", "right_support", m),
            "height producer inner cycle " + str(m),
        )
        expected = set(reference)
        for (label, _), outer in zip(cycle_tables[1:], cycle_indexes[1:]):
            index_rows(outer[m]["cases"], lambda r, m=m: full_support(r, "left", "right", m),
                       label + " inner " + str(m), expected)
        flat_expected.update((m, U, V) for U, V in expected)
    require(len(flat_expected) == 580, "legacy cycle reference must contain 580 distinct keys")

    def widened_cycle_key(row: dict) -> tuple:
        m = integer(row["cycle_size"])
        U, V = full_support(row, "left", "right", m)
        return m, U, V

    index_rows(defect["widened_165_checks"]["cycle_cases"], widened_cycle_key,
               "all widened cycle scalar receipts", flat_expected)
    checks.append({"domain": "strict 16 cycle wrappers and cross-bound inner/widened supports", "keys": 580})

    # Literal exceptional descriptor coverage, independent of receipt contents.
    ranks = {"A1": 1, "E6": 6, "E7": 7, "E8": 8, "F4": 4, "G2": 2}
    exceptional_expected = {(name, f) for name, r in ranks.items() for f in [None] + list(range(r + 1))}
    exceptional_tables = [
        ("height producer exceptional", old["exceptional_cases"]),
        ("height review exceptional", old_review["exceptional_cases"]),
        ("defect widened exceptional", defect["widened_165_checks"]["exceptional_records"]),
        ("defect residual exceptional", defect["residual_exceptional_cases"]),
        ("defect physical exceptional", defect_review["physical_exceptional"]["records"]),
        ("unit producer exceptional", unit["exceptional_cases"]),
        ("unit review exceptional", unit_review["exceptional_cases"]),
    ]
    exceptional_indexes = [index_rows(rows, descriptor, label, exceptional_expected)
                           for label, rows in exceptional_tables]
    height, href, widened, residual, physical, uprod, ureview = exceptional_indexes
    single_total = multiple_total = unit_total = 0
    for desc in sorted(exceptional_expected, key=repr):
        # Expected INNER support keys come from independently generated review
        # outputs, whose geometric completeness is an explicit inherited input.
        single_keys = set(index_rows(href[desc]["single_cases"], contact_key, "height reference singletons"))
        multiple_keys = set(index_rows(href[desc]["multiple_cases"], support_key, "height reference multiples"))
        for label, rows in (("height singleton receipts", height[desc]["single_cases"]),
                            ("widened singleton receipts", widened[desc]["singletons"]),
                            ("physical singleton output", physical[desc]["singletons"])):
            index_rows(rows, contact_key, label + repr(desc), single_keys)
        for label, rows in (("height multiple receipts", height[desc]["multiple_cases"]),
                            ("physical multiple output", physical[desc]["multiple"])):
            index_rows(rows, support_key, label + repr(desc), multiple_keys)
        subunit_keys = set(index_rows(physical[desc]["subunit_coordinates"], contact_key, "physical subunit reference"))
        index_rows(residual[desc]["all_diagonals_below_one"], contact_key,
                   "residual subunit receipts" + repr(desc), subunit_keys)
        unit_keys = set(index_rows(ureview[desc]["single_contacts"], contact_key, "unit singleton reference"))
        index_rows(uprod[desc]["single_contact_cases"], contact_key,
                   "unit singleton receipts" + repr(desc), unit_keys)
        single_total += len(single_keys)
        multiple_total += len(multiple_keys)
        unit_total += len(unit_keys)
    require((single_total, multiple_total, unit_total) == (71, 16, 46), "exceptional reference key counts changed")
    require({desc for desc in exceptional_expected if uprod[desc]["two_contact_case"] is not None}
            == {("E7", None)}, "unit multiple-support descriptor domain changed")
    checks.append({"domain": "independent 40 exceptional descriptors; cross-bound inner support receipts",
                   "original_singletons": 71, "original_multiples": 16, "unit_singletons": 46})

    def failure_key(row: dict) -> tuple:
        return descriptor(row) + (support_key(row),)

    failure_reference = index_rows(old_review["possible_defects"], failure_key, "height 15 failure reference")
    require(len(failure_reference) == 15, "failure reference must have 15 distinct occurrences")
    index_rows(defect_review["physical_exceptional"]["defect_occurrences"], failure_key,
               "physical 15 failure output", set(failure_reference))
    removed = {("E7", None, ((0, 4), (0, 7))), ("E7", None, ((1, 4), (1, 7)))}
    require(removed <= set(failure_reference), "both prescribed same-row failure occurrences are required")
    index_rows(two["remaining_failure_floors"], failure_key, "remaining 13 failure floors",
               set(failure_reference) - removed)
    checks.append({"domain": "cross-bound 15-to-13 exceptional occurrence keys", "remaining": 13})

    def budget_key(row: dict) -> tuple:
        return tuple(text_value(z) for z in sequence(row["endpoint_budgets"]))

    index_rows(two["ordered_budget_packing"], budget_key, "ordered packing receipts",
               {("10/9", "4/3"), ("1", "3/2"), ("3/2", "1")})
    index_rows(two["all_rank_tails"]["BCD"],
               lambda r: (text_value(r["type"]), flag(r["modified"]), integer(r["rank_start"])),
               "five BCD analytic-tail receipts",
               {("B", False, 8), ("B", True, 8), ("C", False, 8), ("C", True, 8), ("D", False, 18)})
    index_rows(two["all_rank_tails"]["A"],
               lambda r: (integer(r["parity"]), tuple(text_value(z) for z in sequence(r["cycle_size"]))),
               "two A analytic-tail receipts", {(0, ("18", "2")), (1, ("19", "2"))})
    checks.append({"domain": "independent ordered packing and analytic-tail receipt keys", "packing": 3, "BCD": 5, "A": 2})

    # The 24 prescribed coordinate classes, with integer sides distinguished
    # from the symbolic sides 'r' and 'r-1'. No polynomial is evaluated here.
    coordinate_expected = set()
    for typ, mod, side, classes in (
        ("B", False, 0, ("1", "2<=k<=r-1", "r")),
        ("B", False, "r", ("0", "1", "2<=k<=r-1")),
        ("B", True, "r-1", ("0", "1", "2<=k<=r-2", "r")),
        ("B", True, 0, ("1", "2<=k<=r-1", "r")),
        ("C", False, 0, ("1<=k<=r-1", "r")),
        ("C", True, 1, ("0", "2<=k<=r-1", "r")),
        ("C", True, "r", ("0", "1<=k<=r-1")),
        ("D", False, 0, ("1", "2<=k<=r-2", "r-1", "r")),
    ):
        coordinate_expected.update((typ, mod, side, cls) for cls in classes)

    def coordinate_key(row: dict) -> tuple:
        side = row["side"]
        require(type(side) is int or type(side) is str, "invalid formal side type")
        return text_value(row["type"]), flag(row["modified"]), side, text_value(row["node_class"])

    index_rows(unit["classical_all_rank_checks"]["prefix_coordinate_cases"], coordinate_key,
               "unit producer 24 coordinate-class receipts", coordinate_expected)
    index_rows(defect_review["all_rank_BCD"]["classes"], coordinate_key,
               "defect review 24 coordinate-class receipts", coordinate_expected)
    cycle_charge_expected = {(str(m), (1,)) for m in (7, 8, 9)} | {
        ("2*t + 10", (1,)), ("2*t + 10", (2,)),
        ("2*t + 11", (1,)), ("2*t + 11", (2,)),
    }
    index_rows(unit["classical_all_rank_checks"]["cycle_cases"],
               lambda r: (text_value(r["m"]), nodes(r["allowed_distance"])),
               "unit producer cycle-class receipts", cycle_charge_expected)
    checks.append({"domain": "independent classical coordinate/cycle class receipt keys", "coordinate_classes": 24, "cycle_classes": 7})

    source_path = Path(__file__).resolve()
    sources["scripts/" + source_path.name] = hashlib.sha256(source_path.read_bytes()).hexdigest()
    return {
        "status": "PASS_FINITE57_TAIL_KEY_COVERAGE_ONLY",
        "scope": "Structure/key validation and explicit receipt cross-binding only; no geometric support enumeration, KKT/certificate arithmetic, packing optimization, Lean or full-manuscript verification.",
        "independent_literal_domains": ["same-row E7 full supports", "16 small-cycle wrappers", "40 exceptional descriptors", "three ordered packing keys", "five BCD and two A tail keys", "24 coordinate and seven cycle class keys"],
        "inherited_inner_domain_premises": ["complete original-budget 580 support reconstruction in the legacy height/defect reviewers", "complete literal exceptional subset reconstruction and resulting review-output keys", "legacy residual finite census and all numerical inequalities remain separately required"],
        "checks": checks,
        "sources": sources,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    result = check(args.root.resolve())
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()


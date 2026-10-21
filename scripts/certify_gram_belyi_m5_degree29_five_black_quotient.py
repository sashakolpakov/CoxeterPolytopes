#!/usr/bin/env python3
"""Exact marked-passport/master-system quotient for five black vertices.

This is a combinatorial certificate only.  It enumerates the degree-29
black and white partitions, distinguishes the marked lambda part and ordered
source endpoints, and then forgets the exact white ramification index while
retaining the binary simple/critical endpoint status used by the four-black
coefficient systems.  It launches no coefficient resultant and proves no
five-black no-go.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any, Iterator


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/gram_belyi_m5_degree29_five_black_quotient.json"
SCHEMA = "coxeter.gram_belyi.m5_degree29_five_black_quotient.v1"
DEGREE = 29
BLACK_VERTICES = 5
WHITE_VERTICES = DEGREE + 1 - BLACK_VERTICES


def fixed_length_partitions(
    total: int, length: int, maximum: int | None = None
) -> Iterator[tuple[int, ...]]:
    """Yield decreasing positive partitions of ``total`` with fixed length."""

    if length == 0:
        if total == 0:
            yield ()
        return
    if total < length:
        return
    largest = min(total - length + 1, maximum or total)
    for first in range(largest, 0, -1):
        for tail in fixed_length_partitions(total - first, length - 1, first):
            yield (first,) + tail


def ordered_choice_counts(multiplicities: Counter[int]) -> tuple[int, int]:
    """Choices of one part and two ordered parts, quotienting equal parts."""

    distinct = sum(value > 0 for value in multiplicities.values())
    repeatable = sum(value >= 2 for value in multiplicities.values())
    return distinct, distinct * (distinct - 1) + repeatable


def black_counts_by_partition_marking() -> tuple[int, int, int]:
    """Mark lambda, then zero, one, or two ordered endpoint parts."""

    marked = marked_one = marked_two = 0
    for partition in fixed_length_partitions(DEGREE, BLACK_VERTICES):
        multiplicities = Counter(partition)
        for lambda_degree in multiplicities:
            available = multiplicities.copy()
            available[lambda_degree] -= 1
            one, two = ordered_choice_counts(available)
            marked += 1
            marked_one += one
            marked_two += two
    return marked, marked_one, marked_two


def canonical_black_records(endpoint_count: int) -> list[tuple[int, ...]]:
    """Direct coefficient records with 0, 1, or 2 black endpoints fixed.

    The first entry is the marked lambda degree.  The following
    ``endpoint_count`` entries retain endpoint order.  The remaining black
    degrees are sorted increasingly because their roots are unlabelled.
    """

    if endpoint_count not in (0, 1, 2):
        raise ValueError(endpoint_count)

    records: set[tuple[int, ...]] = set()

    def extend(prefix: tuple[int, ...], slots: int, remaining: int) -> None:
        if slots == 0:
            tail_length = BLACK_VERTICES - len(prefix)
            for tail in fixed_length_partitions(remaining, tail_length):
                records.add(prefix + tuple(reversed(tail)))
            return
        for value in range(1, remaining):
            extend(prefix + (value,), slots - 1, remaining - value)

    for lambda_degree in range(1, DEGREE):
        extend((lambda_degree,), endpoint_count, DEGREE - lambda_degree)
    return sorted(records)


def record_digest(records: list[tuple[int, ...]]) -> str:
    encoded = json.dumps(records, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def partition_label(partition: tuple[int, ...]) -> str:
    multiplicities = Counter(partition)
    terms = []
    for degree in sorted(multiplicities, reverse=True):
        count = multiplicities[degree]
        terms.append(str(degree) if count == 1 else f"{degree}^{count}")
    return "(" + ",".join(terms) + ")"


def white_partition_rows() -> list[dict[str, Any]]:
    rows = []
    for partition in fixed_length_partitions(DEGREE, WHITE_VERTICES):
        multiplicities = Counter(partition)
        critical = {
            degree: count for degree, count in multiplicities.items() if degree > 1
        }
        critical_types = len(critical)
        simple_count = multiplicities[1]
        single = {
            "S": int(simple_count >= 1),
            "C": critical_types,
        }
        ordered = {
            "SS": int(simple_count >= 2),
            "SC": int(simple_count >= 1) * critical_types,
            "CS": int(simple_count >= 1) * critical_types,
            "CC": critical_types * (critical_types - 1)
            + sum(count >= 2 for count in critical.values()),
        }
        one, two = ordered_choice_counts(multiplicities)
        if sum(single.values()) != one or sum(ordered.values()) != two:
            raise ArithmeticError((partition, single, ordered, one, two))
        rows.append(
            {
                "partition": partition_label(partition),
                "part_multiplicities": {
                    str(degree): multiplicities[degree]
                    for degree in sorted(multiplicities, reverse=True)
                },
                "one_white_endpoint_choices": single,
                "two_ordered_white_endpoint_choices": ordered,
            }
        )
    return rows


def build_archive() -> dict[str, Any]:
    partition_counts = black_counts_by_partition_marking()
    records = [canonical_black_records(endpoint_count) for endpoint_count in range(3)]
    direct_counts = tuple(len(value) for value in records)
    if partition_counts != direct_counts:
        raise ArithmeticError((partition_counts, direct_counts))
    marked, marked_one, marked_two = direct_counts

    white_rows = white_partition_rows()
    white_partition_count = len(white_rows)
    single = {
        status: sum(row["one_white_endpoint_choices"][status] for row in white_rows)
        for status in ("S", "C")
    }
    ordered = {
        status: sum(
            row["two_ordered_white_endpoint_choices"][status] for row in white_rows
        )
        for status in ("SS", "SC", "CS", "CC")
    }

    status_rows = [
        ("BB", marked_two, white_partition_count),
        ("BW-S", marked_one, single["S"]),
        ("BW-C", marked_one, single["C"]),
        ("WB-S", marked_one, single["S"]),
        ("WB-C", marked_one, single["C"]),
        ("WW-SS", marked, ordered["SS"]),
        ("WW-SC", marked, ordered["SC"]),
        ("WW-CS", marked, ordered["CS"]),
        ("WW-CC", marked, ordered["CC"]),
    ]
    refined = [
        {
            "endpoint_block": block,
            "binary_status_master_systems": master_count,
            "marked_passport_multiplier": multiplier,
            "raw_marked_passports": master_count * multiplier,
        }
        for block, master_count, multiplier in status_rows
    ]

    raw_by_color = {
        "BB": marked_two * white_partition_count,
        "BW": marked_one * sum(single.values()),
        "WB": marked_one * sum(single.values()),
        "WW": marked * sum(ordered.values()),
    }
    masters_by_color = {
        "BB": marked_two,
        "BW": marked_one * len(single),
        "WB": marked_one * len(single),
        "WW": marked * len(ordered),
    }

    return {
        "schema": SCHEMA,
        "status": "EXACT_FIVE_BLACK_PASSPORT_MASTER_QUOTIENT",
        "degree": DEGREE,
        "black_vertex_count": BLACK_VERTICES,
        "white_vertex_count": WHITE_VERTICES,
        "white_excess_ramification": DEGREE - WHITE_VERTICES,
        "marking_convention": (
            "B(lambda)=0; B/W retain the ordered colours of x=0,x=1; "
            "equal-valency cycle identities are quotiented"
        ),
        "black_degree_record_counts": {
            "M0_marked_lambda_remaining_four_unordered": marked,
            "M1_plus_one_black_endpoint_remaining_three_unordered": marked_one,
            "M2_plus_two_ordered_black_endpoints_remaining_two_unordered": marked_two,
        },
        "black_record_list_sha256": {
            "M0": record_digest(records[0]),
            "M1": record_digest(records[1]),
            "M2": record_digest(records[2]),
        },
        "black_count_crosscheck": {
            "partition_marking": list(partition_counts),
            "direct_canonical_records": list(direct_counts),
            "agree": True,
        },
        "white_partitions": white_rows,
        "white_choice_totals": {
            "partition_count": white_partition_count,
            "one_white_endpoint": single,
            "two_ordered_white_endpoints": ordered,
            "one_white_endpoint_exact_degree_choices": sum(single.values()),
            "two_ordered_white_endpoint_exact_degree_choices": sum(ordered.values()),
        },
        "refined_endpoint_blocks": refined,
        "raw_marked_passports_by_endpoint_colors": {
            **raw_by_color,
            "total": sum(raw_by_color.values()),
        },
        "binary_status_master_systems_by_endpoint_colors": {
            **masters_by_color,
            "total": sum(masters_by_color.values()),
        },
        "raw_marked_passport_count": sum(raw_by_color.values()),
        "binary_status_master_system_count": sum(masters_by_color.values()),
        "quotient_scope": (
            "necessary coefficient systems for no-go proofs only; exact white "
            "ramification indices impose additional multiplicity and inequation data"
        ),
        "coefficient_resultants_run": False,
        "five_black_no_go_claimed": False,
        "global_degree_29_status": "OPEN",
    }


def verify_archive(archive: dict[str, Any]) -> None:
    if archive["black_degree_record_counts"] != {
        "M0_marked_lambda_remaining_four_unordered": 1292,
        "M1_plus_one_black_endpoint_remaining_three_unordered": 4225,
        "M2_plus_two_ordered_black_endpoints_remaining_two_unordered": 11011,
    }:
        raise ArithmeticError("black record quotient changed")
    if archive["white_choice_totals"]["one_white_endpoint"] != {"S": 5, "C": 7}:
        raise ArithmeticError("one-white status count changed")
    if archive["white_choice_totals"]["two_ordered_white_endpoints"] != {
        "SS": 5,
        "SC": 7,
        "CS": 7,
        "CC": 7,
    }:
        raise ArithmeticError("two-white status count changed")
    if archive["raw_marked_passports_by_endpoint_colors"] != {
        "BB": 55055,
        "BW": 50700,
        "WB": 50700,
        "WW": 33592,
        "total": 190047,
    }:
        raise ArithmeticError("raw five-black census changed")
    if archive["binary_status_master_systems_by_endpoint_colors"] != {
        "BB": 11011,
        "BW": 8450,
        "WB": 8450,
        "WW": 5168,
        "total": 33079,
    }:
        raise ArithmeticError("binary-status master quotient changed")
    if archive["coefficient_resultants_run"] or archive["five_black_no_go_claimed"]:
        raise ArithmeticError("quotient overclaimed")


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    archive = build_archive()
    verify_archive(archive)
    encoded = canonical_bytes(archive)
    digest = hashlib.sha256(encoded).hexdigest()
    if arguments.write:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_bytes(encoded)
    if arguments.check:
        if not arguments.output.exists() or arguments.output.read_bytes() != encoded:
            raise SystemExit("five-black quotient archive is stale")
    print(
        json.dumps(
            {
                "status": archive["status"],
                "raw_marked_passports": archive["raw_marked_passport_count"],
                "binary_status_master_systems": archive[
                    "binary_status_master_system_count"
                ],
                "sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()

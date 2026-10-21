#!/usr/bin/env python3
"""Clean-room audit of the golden exactly-two-dotted four-node tree census.

This file deliberately imports no ``coxeter_polytopes`` module.  It rebuilds
the relevant part of ``Z[phi]`` from integer pairs, enumerates the dotted
alphabet, checks the three tree orbits under ``S_4``, and applies the two real
Gram gates with exact integer arithmetic.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations, permutations, product
from math import isqrt
from pathlib import Path


FRAME_BOUND = 62
ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data/enumerations/"
    "golden_agcat_two_dotted_four_node_trees_independent_audit.json"
)
OWNER_OUTPUT = (
    ROOT / "data/enumerations/golden_agcat_two_dotted_four_node_trees.json"
)
REPORT = (
    ROOT
    / "results/golden_agcat_two_dotted_four_node_trees_independent_hostile_audit.md"
)
EXPECTED_DOTTED_ORDER_SHA256 = (
    "96a3040813293c87465a98eeb6bda0cff6ab0c9c5d6a973fe9020b8e0139ec16"
)


@dataclass(frozen=True, slots=True)
class ZPhi:
    """The golden integer ``a+b*phi``, where ``phi**2=phi+1``."""

    a: int
    b: int = 0

    def __add__(self, other: object) -> ZPhi:
        other = _coerce(other)
        return ZPhi(self.a + other.a, self.b + other.b)

    def __radd__(self, other: object) -> ZPhi:
        return self + other

    def __sub__(self, other: object) -> ZPhi:
        other = _coerce(other)
        return ZPhi(self.a - other.a, self.b - other.b)

    def __rsub__(self, other: object) -> ZPhi:
        return _coerce(other) - self

    def __mul__(self, other: object) -> ZPhi:
        other = _coerce(other)
        return ZPhi(
            self.a * other.a + self.b * other.b,
            self.a * other.b + self.b * other.a + self.b * other.b,
        )

    def __rmul__(self, other: object) -> ZPhi:
        return self * other

    def conjugate(self) -> ZPhi:
        return ZPhi(self.a + self.b, -self.b)

    def sign(self, place: int = 0) -> int:
        value = self if place == 0 else self.conjugate()
        # 2*(a+b*phi)=(2*a+b)+b*sqrt(5).  The comparison below is
        # exact; equality in the irrational branch cannot occur for b!=0.
        rational = 2 * value.a + value.b
        radical = value.b
        if radical == 0:
            return (rational > 0) - (rational < 0)
        if radical > 0:
            if rational >= 0:
                return 1
            return (5 * radical * radical > rational * rational) * 2 - 1
        if rational <= 0:
            return -1
        return (rational * rational > 5 * radical * radical) * 2 - 1


def _coerce(value: object) -> ZPhi:
    if isinstance(value, ZPhi):
        return value
    if isinstance(value, int):
        return ZPhi(value)
    return NotImplemented


def _floor(value: ZPhi, place: int = 0) -> int:
    if place == 1:
        value = value.conjugate()
    if value.b == 0:
        return value.a
    magnitude = abs(value.b)
    floor_m_phi = (magnitude + isqrt(5 * magnitude * magnitude)) // 2
    if value.b > 0:
        return value.a + floor_m_phi
    return value.a - floor_m_phi - 1


def _ceil(value: ZPhi, place: int = 0) -> int:
    floor = _floor(value, place)
    coefficient = value.b if place == 0 else -value.b
    return floor if coefficient == 0 else floor + 1


def dotted_values() -> tuple[ZPhi, ...]:
    """Rebuild the complete ``4<x<62^2, 0<x'<4`` alphabet."""

    selected_upper = FRAME_BOUND * FRAME_BOUND
    # From x-x'=b*sqrt(5), sqrt(5)>2, 0<=x<3844 and 0<=x'<4.
    b_min = -1
    b_max = (selected_upper - 1) // 2
    output: list[ZPhi] = []
    for b in range(b_min, b_max + 1):
        minus_b_phi = ZPhi(0, -b)
        lower = max(_ceil(minus_b_phi, 0), _ceil(minus_b_phi, 1))
        upper = min(
            _ceil(ZPhi(selected_upper, -b), 0) - 1,
            _ceil(ZPhi(4, -b), 1) - 1,
        )
        for a in range(lower, upper + 1):
            value = ZPhi(a, b)
            if (
                value.sign(0) >= 0
                and (value - selected_upper).sign(0) < 0
                and value.sign(1) >= 0
                and (value - 4).sign(1) < 0
                and (value - 4).sign(0) > 0
            ):
                output.append(value)
    if len(output) != len(set(output)):
        raise AssertionError("clean-room dotted enumeration has duplicates")
    return tuple(output)


def _insertion_sort(values: tuple[ZPhi, ...], place: int) -> tuple[ZPhi, ...]:
    """Exact sort without floats or an imported golden comparison helper."""

    output: list[ZPhi] = []
    for value in values:
        low = 0
        high = len(output)
        while low < high:
            middle = (low + high) // 2
            if (value - output[middle]).sign(place) < 0:
                high = middle
            else:
                low = middle + 1
        output.insert(low, value)
    return tuple(output)


def _canonical_graph(
    edges: tuple[tuple[int, int, str], ...], permutation: tuple[int, ...]
) -> tuple[tuple[int, int, str], ...]:
    transformed = []
    for first, second, label in edges:
        first, second = sorted((permutation[first], permutation[second]))
        transformed.append((first, second, label))
    return tuple(sorted(transformed))


def _tree_orbit_counts(dotted_label_count: int) -> dict[str, int]:
    """Brute-force labelled trees and quotient by all 24 vertex permutations."""

    vertices = range(4)
    possible_edges = tuple(combinations(vertices, 2))
    vertex_permutations = tuple(permutations(vertices))
    representatives = set()
    for edge_set in combinations(possible_edges, 3):
        adjacency = {vertex: set() for vertex in vertices}
        for first, second in edge_set:
            adjacency[first].add(second)
            adjacency[second].add(first)
        seen = set()
        frontier = [0]
        while frontier:
            vertex = frontier.pop()
            if vertex in seen:
                continue
            seen.add(vertex)
            frontier.extend(adjacency[vertex] - seen)
        if len(seen) != 4:
            continue
        for finite_edge_index in range(3):
            dotted_edges = [index for index in range(3) if index != finite_edge_index]
            for dotted_labels in product(range(dotted_label_count), repeat=2):
                labels = ["F"] * 3
                for index, dotted_label in zip(dotted_edges, dotted_labels):
                    labels[index] = f"D{dotted_label}"
                decorated = tuple(
                    (first, second, label)
                    for (first, second), label in zip(edge_set, labels)
                )
                representatives.add(
                    min(
                        _canonical_graph(decorated, permutation)
                        for permutation in vertex_permutations
                    )
                )
    counts = {"star": 0, "disjoint_path": 0, "adjacent_path": 0}
    for representative in representatives:
        degrees = [0, 0, 0, 0]
        dotted_edges = []
        for first, second, label in representative:
            degrees[first] += 1
            degrees[second] += 1
            if label.startswith("D"):
                dotted_edges.append({first, second})
        if sorted(degrees) == [1, 1, 1, 3]:
            counts["star"] += 1
        elif dotted_edges[0].isdisjoint(dotted_edges[1]):
            counts["disjoint_path"] += 1
        else:
            counts["adjacent_path"] += 1
    return counts


def count_star(values: tuple[ZPhi, ...], finite: ZPhi) -> int:
    """Unordered dotted pair; exact conjugate Schur-complement gate."""

    right = len(values) - 1
    total = 0
    for left, first in enumerate(values):
        while right >= left and (first + values[right] + finite - 4).sign(1) >= 0:
            right -= 1
        if right < left:
            break
        total += right - left + 1
    return total


def count_adjacent_path(values: tuple[ZPhi, ...], finite: ZPhi) -> int:
    """Ordered remote/adjacent dotted pair; exact conjugate determinant gate."""

    right = len(values) - 1
    total = 0
    for remote in values:
        determinant = (4 - remote) * (4 - finite) - 4 * values[right]
        while right >= 0 and determinant.sign(1) <= 0:
            right -= 1
            if right >= 0:
                determinant = (4 - remote) * (4 - finite) - 4 * values[right]
        total += right + 1
    return total


def count_disjoint_path(
    selected_sorted: tuple[ZPhi, ...], finite: ZPhi
) -> tuple[int, ZPhi | None, ZPhi | None]:
    """Unordered pairs passing selected index-one and conjugate-PD gates."""

    total = 0
    witness: tuple[ZPhi | None, ZPhi | None] = (None, None)
    four_finite = 4 * finite
    for left, first in enumerate(selected_sorted):
        first_gap = first - 4
        for second in selected_sorted[left:]:
            selected_determinant = first_gap * (second - 4) - four_finite
            if selected_determinant.sign(0) >= 0:
                break
            conjugate_determinant = (4 - first) * (4 - second) - four_finite
            if conjugate_determinant.sign(1) > 0:
                total += 1
                if witness == (None, None):
                    witness = first, second
    return total, *witness


def _pair_json(value: ZPhi | None) -> dict[str, int] | None:
    if value is None:
        return None
    return {"a": value.a, "b": value.b}


def build_payload() -> dict[str, object]:
    values = dotted_values()
    conjugate_sorted = _insertion_sort(values, 1)
    selected_sorted = _insertion_sort(values, 0)
    digest = sha256(
        "".join(f"{value.a},{value.b}\n" for value in conjugate_sorted).encode()
    ).hexdigest()
    finite_rows = (
        ("3", ZPhi(1)),
        ("4", ZPhi(2)),
        ("5", ZPhi(1, 1)),
        ("6", ZPhi(3)),
        ("10", ZPhi(2, 1)),
    )

    if len(values) != 6_869:
        raise AssertionError(f"expected 6869 dotted values, got {len(values)}")
    if digest != EXPECTED_DOTTED_ORDER_SHA256:
        raise AssertionError(f"unexpected clean-room dotted order digest: {digest}")
    orbit_checks = {}
    for label_count in (1, 2, 3):
        unordered_small = label_count * (label_count + 1) // 2
        expected = {
            "star": unordered_small,
            "disjoint_path": unordered_small,
            "adjacent_path": label_count * label_count,
        }
        observed = _tree_orbit_counts(label_count)
        if observed != expected:
            raise AssertionError(
                f"S4 orbit count at N={label_count}: {observed} != {expected}"
            )
        orbit_checks[str(label_count)] = observed

    unordered = len(values) * (len(values) + 1) // 2
    ordered = len(values) * len(values)
    raw = len(finite_rows) * (2 * unordered + ordered)
    rows: list[dict[str, object]] = []
    for label, finite in finite_rows:
        disjoint, first, second = count_disjoint_path(selected_sorted, finite)
        rows.append(
            {
                "finite_label": label,
                "finite_product": _pair_json(finite),
                "star_retained": count_star(conjugate_sorted, finite),
                "disjoint_path_retained": disjoint,
                "adjacent_path_retained": count_adjacent_path(
                    conjugate_sorted, finite
                ),
                "first_disjoint_witness": {
                    "left": _pair_json(first),
                    "right": _pair_json(second),
                },
            }
        )
    observed_rows = tuple(
        (
            row["finite_label"],
            row["star_retained"],
            row["disjoint_path_retained"],
            row["adjacent_path_retained"],
        )
        for row in rows
    )
    expected_rows = (
        ("3", 6_634_798, 7, 17_689_960),
        ("4", 2_948_447, 5, 11_792_207),
        ("5", 9_650_285, 44, 21_334_861),
        ("6", 737_201, 2, 5_895_415),
        ("10", 5_052_550, 23, 15_436_911),
    )
    if observed_rows != expected_rows:
        raise AssertionError(f"clean-room rows changed: {observed_rows!r}")
    star = sum(int(row["star_retained"]) for row in rows)
    disjoint = sum(int(row["disjoint_path_retained"]) for row in rows)
    adjacent = sum(int(row["adjacent_path_retained"]) for row in rows)
    retained = star + disjoint + adjacent
    if (raw, star, disjoint, adjacent, retained) != (
        471_865_955,
        25_023_281,
        81,
        72_149_354,
        97_172_716,
    ):
        raise AssertionError("clean-room aggregate counts changed")

    witness_x = ZPhi(3, 1)
    witness_z = ZPhi(2, 2)
    witness_determinant = (witness_x - 4) * (witness_z - 4) - 4
    if witness_determinant != ZPhi(0, -2):
        raise AssertionError("the explicit disjoint-path determinant changed")
    if witness_determinant.sign(0) >= 0 or witness_determinant.sign(1) <= 0:
        raise AssertionError("the explicit disjoint-path witness lost its signs")

    return {
        "schema": (
            "golden-agcat-two-dotted-four-node-trees-independent-hostile-"
            "audit-0.1.0"
        ),
        "independence": {
            "imports_project_modules": False,
            "arithmetic": (
                "integer pairs (a,b) for a+b*phi; signs compare integer "
                "squares against 5b^2"
            ),
            "counting": (
                "independent exact two-pointer joins plus a direct selected-"
                "ordered scan for the thin disjoint locus"
            ),
        },
        "alphabet": {
            "dotted_value_count": len(values),
            "conjugate_order_sha256": digest,
        },
        "s4_orbit_checks_for_small_alphabets": orbit_checks,
        "canonical_count_formula": "5*(2*C(N+1,2)+N^2), N=6869",
        "canonical_before_gram": raw,
        "retained_by_topology": {
            "star": star,
            "disjoint_path": disjoint,
            "adjacent_path": adjacent,
        },
        "retained_total": retained,
        "rows": rows,
        "explicit_counterexample_to_old_inertia_claim": {
            "finite_label": "3",
            "left_dotted": _pair_json(witness_x),
            "right_dotted": _pair_json(witness_z),
            "determinant": _pair_json(witness_determinant),
            "selected_sign": witness_determinant.sign(0),
            "conjugate_sign": witness_determinant.sign(1),
            "selected_inertia": [3, 1, 0],
            "conjugate_inertia": [4, 0, 0],
        },
        "integral_cartan_audit": (
            "PASS: on a tree orient each product c as (-1,-c); starting "
            "with the product of all c_e and dividing along each oriented "
            "root path gives a totally positive integral diagonal "
            "symmetrizer"
        ),
        "scope": (
            "exactly two dotted edges on connected four-node trees only; "
            "no cyclic packets are enumerated"
        ),
        "verdict": "PASS after correcting the disjoint-path locus",
    }


def report() -> str:
    return """# Independent hostile audit: golden two-dotted four-node trees

This is a clean-room audit: its certifier imports no project module and
rebuilds `Z[phi]`, the `4 < x < 62^2`, `0 < x' < 4` dotted alphabet, exact
real-place signs, and all counts from integer pairs `(a,b)`.

The original blanket rejection of the disjoint-dotted path was false.  Its
two indefinite edge blocks are coupled by the finite middle edge, so their
negative indices do not add.  For example, with finite label `3`, endpoint
products

```text
x = 3 + phi,   z = 2 + 2 phi
```

give full determinant `-2 phi`: it is negative at the selected embedding
and positive at the conjugate embedding.  Exact leading-principal-minor
signs give selected inertia `(3,1,0)` and conjugate inertia `(4,0,0)`.

After that correction, the audit passes.  It independently obtains 6,869
dotted values and the same conjugate-order SHA-256, and brute `S_4` orbit
enumeration on small alphabets proves the roles used by the closed count:
stars and disjoint paths use unordered dotted pairs, while adjacent paths
use ordered remote/adjacent pairs.  Thus the raw count is

```text
5 * (2*C(6869+1,2) + 6869^2) = 471,865,955.
```

The exact retained rows `(star, disjoint, adjacent)` for finite labels
`3,4,5,6,10` are respectively

```text
(6,634,798,  7, 17,689,960)
(2,948,447,  5, 11,792,207)
(9,650,285, 44, 21,334,861)
(  737,201,  2,  5,895,415)
(5,052,550, 23, 15,436,911).
```

Therefore exactly **97,172,716** packets survive: 25,023,281 stars, 81
disjoint-dotted paths, and 72,149,354 adjacent-dotted paths.  All comparisons
are strict and exact.  The conjugate-PD formulas follow from Schur complements
and Sylvester's criterion; the selected disjoint path has index one exactly
when its full determinant is negative.

The integral Cartan/symmetrizer claim also passes.  Orient each tree edge
product `c` as `(-1,-c)`.  With the product of all edge products at a root,
divide recursively by `c` along oriented edges.  This gives an integral,
totally positive diagonal symmetrizer, with no cycle-consistency condition.

Strict scope: this audit covers connected four-node trees with exactly two
dotted edges.  It does not enumerate cyclic two-dotted packets.
"""


def verify_owner_artifact(audit: dict[str, object]) -> None:
    """Compare the frozen owner JSON to the independently rebuilt census."""

    owner = json.loads(OWNER_OUTPUT.read_text())
    if owner["alphabet"]["dotted_value_count"] != audit["alphabet"][
        "dotted_value_count"
    ]:
        raise SystemExit("owner dotted alphabet count differs from clean-room audit")
    if owner["alphabet"]["dotted_values_sorted_by_conjugate_sha256"] != audit[
        "alphabet"
    ]["conjugate_order_sha256"]:
        raise SystemExit("owner dotted alphabet order differs from clean-room audit")
    if owner["canonical_tree_packets_before_gram_gate"] != audit[
        "canonical_before_gram"
    ]:
        raise SystemExit("owner raw tree count differs from clean-room audit")
    if owner["two_place_retained_count"] != audit["retained_total"]:
        raise SystemExit("owner retained total differs from clean-room audit")
    expected_topologies = {
        "star": audit["retained_by_topology"]["star"],
        "disjoint_dotted_path": audit["retained_by_topology"]["disjoint_path"],
        "adjacent_dotted_path": audit["retained_by_topology"]["adjacent_path"],
    }
    if owner["retained_by_topology"] != expected_topologies:
        raise SystemExit("owner topology totals differ from clean-room audit")
    owner_rows = tuple(
        (
            row["finite_label"],
            row["star_retained"],
            row["disjoint_path_retained"],
            row["adjacent_path_retained"],
        )
        for row in owner["rows"]
    )
    audit_rows = tuple(
        (
            row["finite_label"],
            row["star_retained"],
            row["disjoint_path_retained"],
            row["adjacent_path_retained"],
        )
        for row in audit["rows"]
    )
    if owner_rows != audit_rows:
        raise SystemExit("owner per-label rows differ from clean-room audit")
    scope = owner["scope"]
    if (
        scope["connected_graph"] != "tree"
        or scope["dotted_edge_count"] != 2
        or "cycles" not in scope["not_covered"]
    ):
        raise SystemExit("owner artifact no longer states the audited tree scope")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    payload_text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    report_text = report()
    if args.verify:
        if OUTPUT.read_text() != payload_text:
            raise SystemExit("frozen independent audit JSON differs")
        if REPORT.read_text() != report_text:
            raise SystemExit("frozen independent audit report differs")
        verify_owner_artifact(payload)
        print("PASS")
        return
    OUTPUT.write_text(payload_text)
    REPORT.write_text(report_text)
    print(OUTPUT)
    print(REPORT)


if __name__ == "__main__":
    main()

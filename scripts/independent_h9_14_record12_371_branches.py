#!/usr/bin/env python3
"""No-import replay of the corrected record-12 branch quotient.

The script regenerates the 1,391,382 locally feasible ``(X,p,y)`` rows from
the frozen order-ten triangle-free graph stream, applies only the proved
special-meeting blocker filter, and independently quotients the 6,184 rows by
coloured graph isomorphism.  It imports no record-12 census/residual module
and does not invoke nauty for the quotient replay.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
X_ARCHIVE = ROOT / "data/enumerations/h9_14_record12_x_graphs.json"
LOCAL_ALL = Path("/tmp/record12_local_all.txt")
SAFE_RAW = Path("/tmp/record12_branch_safe_special_correct.txt")
COLOURED = Path("/tmp/record12_branch_safe_special_correct_orbits.txt")

EXPECTED_X_CERTIFICATE = (
    "84e80de3be27e5bec86791e4dc5721093cf6160692a327820040625fd0cc9056"
)
EXPECTED_LOCAL_ALL_SHA256 = (
    "ed12edfbff40c2fab730ed3f1740bcf6e5272554f85a807468c7cb8e515961ca"
)
EXPECTED_SAFE_RAW_SHA256 = (
    "c91021fe69bc047c8f2eea5f43bf708d9025492c8e3cac3552205b857731e4ee"
)
EXPECTED_COLOURED_SHA256 = (
    "640217770888c33fe60fc2d1ffd6dbb05c499c9307c40d4ec9642a482fca1905"
)

N = 10
FULL10 = (1 << N) - 1
SPECIALS = (1 << 10) | (1 << 11)
PAIR_INDEX = {pair: index for index, pair in enumerate(combinations(range(N), 2))}


def canonical_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def bits(value: int):
    while value:
        bit = value & -value
        yield bit
        value ^= bit


def decode(word: str) -> tuple[int, ...]:
    codes = [ord(character) - 63 for character in word]
    if not codes or codes[0] != N:
        raise ValueError("expected an order-ten short graph6 word")
    payload = []
    for code in codes[1:]:
        payload.extend((code >> shift) & 1 for shift in range(5, -1, -1))
    result = [0] * N
    at = 0
    for second in range(1, N):
        for first in range(second):
            if payload[at]:
                result[first] |= 1 << second
                result[second] |= 1 << first
            at += 1
    if any(payload[at:]):
        raise ValueError("nonzero graph6 padding")
    return tuple(result)


def edge(rows: Sequence[int], first: int, second: int) -> bool:
    return bool(rows[first] & (1 << second))


def cut_edge(cut: int, first: int, second: int) -> bool:
    return bool(cut & (1 << first)) != bool(cut & (1 << second))


def branch_triples(rows: Sequence[int], cut: int, y: int) -> tuple[int, ...]:
    result = set()
    for vertex in range(N):
        if y & (1 << vertex):
            result.add(SPECIALS | (1 << vertex))
    for first, second in combinations(range(N), 2):
        x = edge(rows, first, second)
        delta = cut_edge(cut, first, second)
        if x:
            result.add((1 << 11) | (1 << first) | (1 << second))
        if x != delta:
            result.add((1 << 10) | (1 << first) | (1 << second))
    for triple in combinations(range(N), 3):
        if sum(edge(rows, *pair) for pair in combinations(triple, 2)) % 2:
            result.add(sum(1 << vertex for vertex in triple))
    return tuple(sorted(result))


def minimalize(rows: Iterable[int]) -> tuple[int, ...]:
    ordered = sorted(set(rows), key=lambda row: (row.bit_count(), row))
    result = []
    for row in ordered:
        if not any(prior & row == prior for prior in result):
            result.append(row)
    return tuple(result)


def blocker(clutter: Iterable[int]) -> tuple[int, ...]:
    """Incremental blocker, including ``b(empty)={empty}``."""

    rows = (0,)
    for hyperedge in sorted(set(clutter)):
        candidates = []
        for row in rows:
            if row & hyperedge:
                candidates.append(row)
            else:
                candidates.extend(row | bit for bit in bits(hyperedge))
        rows = minimalize(candidates)
    return rows


def add_y_rows(base_blocker: Sequence[int], y: int) -> tuple[int, ...]:
    """Blocker of ``T union {{s,t,i}:i in y}`` from ``b(T)``.

    A base blocker row meeting ``s`` or ``t`` already hits every added row.
    A row avoiding both must be enlarged minimally by ``s``, by ``t``, or by
    all of ``y``.  Global minimalization gives the displayed formula.  For an
    empty base clutter ``base_blocker=(0,)``; when ``y=0`` this correctly
    returns ``(0,)``, and otherwise it returns the minimal rows ``s,t,y``.
    """

    candidates = []
    for row in base_blocker:
        if row & SPECIALS:
            candidates.append(row)
        else:
            candidates.extend((row | y, row | (1 << 10), row | (1 << 11)))
    return minimalize(candidates)


def safe_special_filter(blocker_rows: Sequence[int]) -> bool:
    return all(
        not row & SPECIALS or row.bit_count() in (4, 5)
        for row in blocker_rows
    )


def forbidden_cut_tables() -> dict[tuple[int, int, int], int]:
    result = {}
    for first, second, isolate in (
        (first, second, isolate)
        for first, second in combinations(range(N), 2)
        for isolate in range(N)
        if isolate not in (first, second)
    ):
        invalid = 0
        for compressed in range(1 << (N - 1)):
            cut = compressed << 1  # normalized by 0 outside the cut.
            endpoints_equal = bool(cut & (1 << first)) == bool(
                cut & (1 << second)
            )
            isolate_differs = bool(cut & (1 << isolate)) != bool(
                cut & (1 << first)
            )
            if endpoints_equal and isolate_differs:
                invalid |= 1 << compressed
        result[first, second, isolate] = invalid
    return result


FORBIDDEN_CUTS = forbidden_cut_tables()


def _bit_indices(value: int):
    while value:
        bit = value & -value
        yield bit.bit_length() - 1
        value ^= bit


def valid_switched_cuts(rows: Sequence[int]) -> tuple[int, ...]:
    """All normalized cuts for which switching preserves triangle-freeness.

    Since ``X`` is triangle-free, switching creates a triangle precisely on
    an induced ``K2+I`` whose isolated vertex is separated from the edge by
    the cut.  A 512-bit mask enforces all such constraints at once.
    """

    invalid = 0
    for triple in combinations(range(N), 3):
        present = [pair for pair in combinations(triple, 2) if edge(rows, *pair)]
        if len(present) == 1:
            first, second = present[0]
            isolate = next(vertex for vertex in triple if vertex not in present[0])
            invalid |= FORBIDDEN_CUTS[first, second, isolate]
    all_compressed = (1 << (1 << (N - 1))) - 1
    return tuple(2 * index for index in _bit_indices(all_compressed ^ invalid))


def submasks(value: int) -> tuple[int, ...]:
    result = []
    subset = value
    while True:
        result.append(subset)
        if subset == 0:
            break
        subset = (subset - 1) & value
    return tuple(result)


def locally_feasible_y(rows: Sequence[int], cut: int) -> tuple[int, ...]:
    """Closed form of the six-term zero-or-two local ridge equation.

    There is never a forced-1 value of the opposite edge.  The support of
    ``y`` lies in one cut side, and a chosen vertex may have no ``X``-neighbor
    in that same side.
    """

    first_side = FULL10 ^ cut
    second_side = cut
    safe_first = sum(
        1 << vertex
        for vertex in range(N)
        if first_side & (1 << vertex) and not rows[vertex] & first_side
    )
    safe_second = sum(
        1 << vertex
        for vertex in range(N)
        if second_side & (1 << vertex) and not rows[vertex] & second_side
    )
    return tuple(sorted(set(submasks(safe_first)) | set(submasks(safe_second))))


def colours(cut: int, y: int) -> tuple[int, ...]:
    return tuple(
        2 * int(bool(cut & (1 << vertex)))
        + int(bool(y & (1 << vertex)))
        for vertex in range(N)
    )


def refined_labels(rows: Sequence[int], initial: Sequence[int]) -> tuple[int, ...]:
    labels = tuple(initial)
    while True:
        signatures = []
        for vertex in range(N):
            counts = Counter(
                labels[other]
                for other in range(N)
                if rows[vertex] & (1 << other)
            )
            signatures.append((labels[vertex], tuple(sorted(counts.items()))))
        alphabet = {signature: index for index, signature in enumerate(sorted(set(signatures)))}
        updated = tuple(alphabet[signature] for signature in signatures)
        if updated == labels:
            return labels
        labels = updated


def isomorphism_key(rows: Sequence[int], colour_values: Sequence[int]) -> tuple[object, ...]:
    labels = refined_labels(rows, colour_values)
    cell_count = max(labels, default=-1) + 1
    cell_sizes = tuple(labels.count(cell) for cell in range(cell_count))
    edge_counts = Counter()
    for first, second in combinations(range(N), 2):
        if edge(rows, first, second):
            edge_counts[tuple(sorted((labels[first], labels[second])))] += 1
    vertex_data = tuple(
        sorted(
            (
                colour_values[vertex],
                labels[vertex],
                rows[vertex].bit_count(),
                tuple(
                    sum(
                        bool(rows[vertex] & (1 << other))
                        and labels[other] == cell
                        for other in range(N)
                    )
                    for cell in range(cell_count)
                ),
            )
            for vertex in range(N)
        )
    )
    return (
        tuple(Counter(colour_values)[value] for value in range(4)),
        sum(row.bit_count() for row in rows) // 2,
        cell_sizes,
        tuple(sorted(edge_counts.items())),
        vertex_data,
    )


def coloured_isomorphic(
    first: Sequence[int], first_colours: Sequence[int],
    second: Sequence[int], second_colours: Sequence[int],
) -> bool:
    first_labels = refined_labels(first, first_colours)
    second_labels = refined_labels(second, second_colours)
    if Counter(first_labels) != Counter(second_labels):
        return False
    mapping: dict[int, int] = {}
    used = set()

    def search() -> bool:
        if len(mapping) == N:
            return True
        choices = []
        for vertex in range(N):
            if vertex in mapping:
                continue
            images = [
                image
                for image in range(N)
                if image not in used
                and first_colours[vertex] == second_colours[image]
                and first_labels[vertex] == second_labels[image]
                and all(
                    bool(first[vertex] & (1 << prior))
                    == bool(second[image] & (1 << mapped))
                    for prior, mapped in mapping.items()
                )
            ]
            if not images:
                return False
            choices.append((len(images), vertex, images))
        _, vertex, images = min(choices)
        for image in images:
            mapping[vertex] = image
            used.add(image)
            if search():
                return True
            used.remove(image)
            del mapping[vertex]
        return False

    return search()


def replay() -> dict[str, object]:
    x_document = json.loads(X_ARCHIVE.read_text(encoding="utf-8"))
    if (
        x_document.get("certificate_sha256") != EXPECTED_X_CERTIFICATE
        or canonical_digest(x_document) != EXPECTED_X_CERTIFICATE
        or x_document.get("graph_count") != 12_172
    ):
        raise AssertionError("hostile record-12 X archive changed")
    for path, expected in (
        (LOCAL_ALL, EXPECTED_LOCAL_ALL_SHA256),
        (SAFE_RAW, EXPECTED_SAFE_RAW_SHA256),
        (COLOURED, EXPECTED_COLOURED_SHA256),
    ):
        if file_sha256(path) != expected:
            raise AssertionError(f"hostile intermediate hash changed: {path}")

    actual_all = LOCAL_ALL.open("r", encoding="ascii")
    actual_safe = SAFE_RAW.open("r", encoding="ascii")
    safe_rows: list[tuple[str, int, int, tuple[int, ...]]] = []
    all_count = safe_count = empty_base_count = 0
    try:
        for word in x_document["graph6"]:
            rows = decode(word)
            for cut in valid_switched_cuts(rows):
                base = blocker(branch_triples(rows, cut, 0))
                empty_base_count += int(base == (0,))
                for y in locally_feasible_y(rows, cut):
                    expected_line = f"{word} {cut} {y}\n"
                    if actual_all.readline() != expected_line:
                        raise AssertionError("hostile local-feasible stream mismatch")
                    all_count += 1
                    augmented = add_y_rows(base, y)
                    if not safe_special_filter(augmented):
                        continue
                    if actual_safe.readline() != expected_line:
                        raise AssertionError("hostile safe-special stream mismatch")
                    # Direct replay on every retained row separately checks the
                    # blocker lemma implementation, including empty clutters.
                    if blocker(branch_triples(rows, cut, y)) != augmented:
                        raise AssertionError("hostile y-row blocker lemma replay failed")
                    safe_rows.append((word, cut, y, rows))
                    safe_count += 1
        if actual_all.readline() or actual_safe.readline():
            raise AssertionError("hostile branch streams have trailing rows")
    finally:
        actual_all.close()
        actual_safe.close()
    if (all_count, safe_count) != (1_391_382, 6_184):
        raise AssertionError("hostile branch counts changed")

    representatives = []
    buckets: dict[tuple[object, ...], list[int]] = defaultdict(list)
    expected_multiplicities = []
    for index, line in enumerate(COLOURED.read_text(encoding="ascii").splitlines()):
        word, raw_cut, raw_y, raw_multiplicity = line.split()
        cut, y = int(raw_cut), int(raw_y)
        rows = decode(word)
        representative = (word, cut, y, rows, colours(cut, y))
        representatives.append(representative)
        expected_multiplicities.append(int(raw_multiplicity))
        buckets[isomorphism_key(rows, representative[4])].append(index)
    if len(representatives) != 1_308:
        raise AssertionError("hostile coloured representative count changed")

    observed = Counter()
    for word, cut, y, rows in safe_rows:
        colour_values = colours(cut, y)
        candidates = buckets.get(isomorphism_key(rows, colour_values), [])
        matches = [
            index
            for index in candidates
            if coloured_isomorphic(
                rows, colour_values, representatives[index][3], representatives[index][4]
            )
        ]
        if len(matches) != 1:
            raise AssertionError(
                f"hostile coloured quotient has {len(matches)} matches for {(word, cut, y)}"
            )
        observed[matches[0]] += 1
    if [observed[index] for index in range(len(representatives))] != expected_multiplicities:
        raise AssertionError("hostile coloured quotient multiplicities changed")

    summary: dict[str, object] = {
        "x_graphs": len(x_document["graph6"]),
        "locally_feasible_rows": all_count,
        "safe_special_rows": safe_count,
        "coloured_branch_orbits": len(representatives),
        "coloured_multiplicity_sum": sum(expected_multiplicities),
        "base_empty_clutter_count": empty_base_count,
        "empty_clutter_blocker": list(blocker(())),
        "all_retained_y_blockers_match_direct_recomputation": True,
        "quotient_replayed_without_nauty": True,
        "input_hashes": {
            "local_all": EXPECTED_LOCAL_ALL_SHA256,
            "safe_raw": EXPECTED_SAFE_RAW_SHA256,
            "coloured": EXPECTED_COLOURED_SHA256,
        },
    }
    payload = json.dumps(summary, sort_keys=True, separators=(",", ":")).encode()
    summary["summary_sha256"] = sha256(payload).hexdigest()
    return summary


def main() -> None:
    print(json.dumps(replay(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

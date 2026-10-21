#!/usr/bin/env python3
"""Quotient the safe, unfiltered record-12 branch stream by coloured S_10.

Input rows have the form ``graph6 p_mask y_mask``.  The four vertex colours
are ``(p_i,y_i)``.  We first relabel every graph so equal colours are
contiguous, then invoke nauty ``labelg`` once per colour multiplicity.  The
output is an exact list of coloured branch orbits; no blocker-purity
condition is used here.

This is an exploratory enumerator.  Its output is not a theorem certificate.
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
import shutil
import subprocess
from tempfile import TemporaryDirectory

from coxeter_polytopes.combinatorics.missing_edges import (
    decode_graph6,
    encode_graph6,
)


N = 10


def standardized_graph(word: str, p_mask: int, y_mask: int) -> tuple[str, str, int, int]:
    order, edges = decode_graph6(word)
    if order != N:
        raise ValueError("expected an order-ten graph6 word")
    colours = tuple(
        2 * int(bool(p_mask & (1 << vertex)))
        + int(bool(y_mask & (1 << vertex)))
        for vertex in range(N)
    )
    old_vertices = tuple(sorted(range(N), key=lambda vertex: (colours[vertex], vertex)))
    old_to_new = {old: new for new, old in enumerate(old_vertices)}
    relabelled_edges = (
        tuple(sorted((old_to_new[first], old_to_new[second])))
        for first, second in edges
    )
    sorted_colours = tuple(colours[old] for old in old_vertices)
    colour_string = "".join(chr(ord("a") + colour) for colour in sorted_colours)
    canonical_p = sum((colour // 2) << vertex for vertex, colour in enumerate(sorted_colours))
    canonical_y = sum((colour % 2) << vertex for vertex, colour in enumerate(sorted_colours))
    return encode_graph6(N, relabelled_edges), colour_string, canonical_p, canonical_y


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=Path("/tmp/record12_local_all.txt"),
    )
    parser.add_argument(
        "output",
        nargs="?",
        type=Path,
        default=Path("/tmp/record12_local_orbits.txt"),
    )
    args = parser.parse_args()

    labelg = shutil.which("labelg")
    if labelg is None:
        raise RuntimeError("nauty labelg is required")

    with TemporaryDirectory(prefix="record12-colour-") as raw_directory:
        directory = Path(raw_directory)
        handles: dict[str, object] = {}
        row_count = 0
        try:
            with args.input.open("r", encoding="ascii") as source:
                for line in source:
                    word, raw_p, raw_y = line.split()
                    standardized, colour_string, _, _ = standardized_graph(
                        word, int(raw_p), int(raw_y)
                    )
                    handle = handles.get(colour_string)
                    if handle is None:
                        handle = (directory / f"group-{colour_string}.g6").open(
                            "w", encoding="ascii"
                        )
                        handles[colour_string] = handle
                    handle.write(standardized + "\n")
                    row_count += 1
        finally:
            for handle in handles.values():
                handle.close()

        orbit_rows: list[tuple[str, int, int, int]] = []
        total_multiplicity = 0
        for colour_string in sorted(handles):
            input_path = directory / f"group-{colour_string}.g6"
            output_path = directory / f"canonical-{colour_string}.g6"
            subprocess.run(
                [labelg, "-q", "-f" + colour_string, str(input_path), str(output_path)],
                check=True,
            )
            multiplicities = Counter(output_path.read_text(encoding="ascii").splitlines())
            sorted_colours = tuple(ord(character) - ord("a") for character in colour_string)
            p_mask = sum((colour // 2) << vertex for vertex, colour in enumerate(sorted_colours))
            y_mask = sum((colour % 2) << vertex for vertex, colour in enumerate(sorted_colours))
            for canonical_word, multiplicity in sorted(multiplicities.items()):
                orbit_rows.append((canonical_word, p_mask, y_mask, multiplicity))
                total_multiplicity += multiplicity

        if total_multiplicity != row_count:
            raise AssertionError("coloured quotient lost branch rows")
        orbit_rows.sort(key=lambda row: (row[1].bit_count(), row[2].bit_count(), row[1], row[2], row[0]))
        with args.output.open("w", encoding="ascii") as output:
            for word, p_mask, y_mask, multiplicity in orbit_rows:
                output.write(f"{word} {p_mask} {y_mask} {multiplicity}\n")
        print(
            f"rows={row_count} colour_partitions={len(handles)} "
            f"coloured_orbits={len(orbit_rows)} multiplicity={total_multiplicity}"
        )


if __name__ == "__main__":
    main()

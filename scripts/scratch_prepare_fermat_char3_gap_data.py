"""Prepare GAP literals for the archived characteristic-3 Fermat line graph.

This is scratch support for ``scratch_fermat_char3_full_stabilizer.g``.  It
extracts only the exact integral Neron--Severi Gram matrix and the 112 line
classes from Kondo--Shimada's official data, computes the line-incidence
graph, and writes GAP-readable literals to a caller-selected scratch path.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import sympy as sp

from scripts.search_ss_fermat_quartic_char3_folds import _parse_array
DEFAULT_SOURCE = (
    ROOT
    / "data"
    / "external"
    / "kondo_shimada_fermat_quartic_char3"
    / "compdataAutFQChar3.txt"
)


def gap_literal(value: object) -> str:
    if isinstance(value, list):
        return "[" + ",".join(gap_literal(item) for item in value) + "]"
    return str(int(value))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    args = parser.parse_args()

    text = args.source.read_text(encoding="utf-8")
    gram = sp.Matrix(_parse_array(text, "FQGramS"))
    lines = sp.Matrix(_parse_array(text, "FQlineclasses"))
    intersections = lines * gram * lines.T
    if intersections.shape != (112, 112):
        raise ArithmeticError(f"unexpected line Gram shape {intersections.shape}")
    if any(intersections[i, i] != -2 for i in range(112)):
        raise ArithmeticError("a published line does not have square -2")
    adjacency = [
        [j + 1 for j in range(112) if i != j and intersections[i, j] == 1]
        for i in range(112)
    ]
    if any(len(row) != 30 for row in adjacency):
        raise ArithmeticError("the Hermitian line graph is not 30-regular")

    rows = [[int(value) for value in lines.row(i)] for i in range(112)]
    payload = (
        "archivedAdjacency := "
        + gap_literal(adjacency)
        + ";\narchivedLineClasses := "
        + gap_literal(rows)
        + ";\narchivedNSBasis := "
        + gap_literal(_parse_array(text, "FQNSbasis"))
        + ";\narchivedGramS := "
        + gap_literal(_parse_array(text, "FQGramS"))
        + ";\narchivedFrobenius := "
        + gap_literal(_parse_array(text, "FQAF"))
        + ";\n"
    )
    args.output.write_text(payload, encoding="utf-8")


if __name__ == "__main__":
    main()

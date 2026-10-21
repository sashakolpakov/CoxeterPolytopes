"""Export the frozen Jacobian--Kummer H12 facet Gram in CoxIter syntax.

This exporter deliberately uses only the exact pair relations.  A CoxIter
edge labelled ``1`` records an unspecified ultraparallel pair; its precise
distance is irrelevant to the elliptic/parabolic subdiagram test and f-vector.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "jacobian_kummer_h12_interim_facet_gram.json"
)


def _fraction(value: str) -> Fraction:
    return Fraction(value)


def _coxiter_label(left: Fraction, cross: Fraction, right: Fraction) -> int | None:
    if cross == 0:
        return None
    ratio = cross * cross / (left * right)
    if ratio == Fraction(1, 4):
        return 3
    if ratio == Fraction(1, 2):
        return 4
    if ratio == Fraction(3, 4):
        return 6
    if ratio == 1:
        return 0
    if ratio > 1:
        return 1
    raise ArithmeticError(f"unsupported Coxeter ratio {ratio}")


def render(path: Path) -> str:
    document = json.loads(path.read_text(encoding="utf-8"))
    status = document["status"]
    exact_pair_gate = status.get(
        "exact_coxeter_pair_gate",
        status.get("exact_full_coxeter_pair_gate", False),
    )
    if not (status["exact_irredundant_facets"] and exact_pair_gate):
        raise ArithmeticError("the interim artifact has not passed the exact pair gates")
    gram = [[_fraction(value) for value in row] for row in document["positive_wall_gram"]]
    facet_count = document["facet_count"]
    dimension = document["group_class"]["hyperbolic_dimension"]
    if len(gram) != facet_count or any(len(row) != facet_count for row in gram):
        raise ArithmeticError("the positive wall Gram has the wrong size")

    lines = [f"{facet_count} {dimension}"]
    for left in range(facet_count):
        if gram[left][left] <= 0:
            raise ArithmeticError("a positive wall Gram diagonal is not positive")
        for right in range(left + 1, facet_count):
            if gram[left][right] != gram[right][left] or gram[left][right] > 0:
                raise ArithmeticError("the positive wall Gram is not acute symmetric")
            label = _coxiter_label(
                gram[left][left], gram[left][right], gram[right][right]
            )
            if label is not None:
                lines.append(f"{left + 1} {right + 1} {label}")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    arguments = parser.parse_args()
    print(render(arguments.input), end="")


if __name__ == "__main__":
    main()

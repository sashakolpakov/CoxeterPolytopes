"""Export the H10 V4 Fermat cast into the generic exact Vinberg schema."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from scripts.search_ss_fermat_quartic_char3_folds import parse_source


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/enumerations/ss_fermat_quartic_char3_v4_fold_screen.json"
OUTPUT = ROOT / "data/enumerations/_scratch_fermat_char3_v4_h10_volume_input.json"


def text(value: sp.Rational) -> str:
    value = sp.Rational(value)
    return str(value.p) if value.q == 1 else f"{value.p}/{value.q}"


def main() -> None:
    artifact = json.loads(INPUT.read_text(encoding="utf-8"))
    selected = [
        cast
        for cast in artifact["casts"]
        if cast["subgroup"]["hyperbolic_dimension"] == 10
    ]
    if len(selected) != 1:
        raise ArithmeticError("expected one H10 V4 cast")
    cast = selected[0]
    if not (
        cast["strict_facet_gate"]["bad_pair_count"] == 0
        and cast["global_facet_reduction"]["exact_complete"] is True
    ):
        raise ArithmeticError("H10 V4 cast has not passed the complete facet gate")
    source = parse_source()
    inverse = source["gram"].inv()
    normals = [
        sp.Matrix(1, 22, row)
        for row in cast["walls"]["primitive_restricted_normal_covectors"]
    ]
    indices = [
        certificate["restricted_ray_index"] - 1
        for certificate in cast["strict_facet_gate"]["certificates"]
    ]
    facets = [normals[index] for index in indices]
    raw = sp.Matrix(
        [[(left * inverse * right.T)[0] for right in facets] for left in facets]
    )
    positive = -raw
    if not (
        len(facets) == 123
        and raw.rank() == 11
        and all(positive[i, i] > 0 for i in range(len(facets)))
        and all(
            positive[i, j] <= 0
            for i in range(len(facets))
            for j in range(i)
        )
    ):
        raise ArithmeticError("unexpected H10 V4 facet Gram")
    payload = {
        "status": {"all_facet_pairs_are_coxeter": True},
        "group_class": {"hyperbolic_dimension": 10},
        "facet_count": len(facets),
        "facet_restricted_ray_indices": [index + 1 for index in indices],
        "positive_wall_gram": [
            [text(positive[row, column]) for column in range(positive.cols)]
            for row in range(positive.rows)
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(OUTPUT)


if __name__ == "__main__":
    main()

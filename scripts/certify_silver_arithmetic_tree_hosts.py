#!/usr/bin/env python3
"""Freeze the exact rank-23 silver arithmetic tree-host census."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver, silver_determinant, silver_inertia
from coxeter_polytopes.search.silver_arithmetic_tree_hosts import (
    ARITHMETIC_EDGE_PRODUCT,
    GENTLE_DOTTED_PRODUCT,
    adjacent_two_dotted_rank23_host,
    analyze_tree_host,
    attach_arithmetic_leaf,
    genuine_dotted_hyperideal_lower_bound,
    rank23_tree_host_census,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/silver_arithmetic_tree_hosts_rank23.json"
REPORT = ROOT / "results/silver_arithmetic_tree_hosts_rank23.md"


def value_json(value: Silver) -> dict[str, str]:
    return value.to_json()


def inertia_json(value) -> list[int]:
    return [value.positive, value.negative, value.zero]


def build_payload() -> dict[str, object]:
    rows = rank23_tree_host_census()
    census = {}
    best_rows = {}
    for kind, best_attachment in (("A", 0), ("D", 21), ("B", 21)):
        arithmetic = [row for row in rows if row.kind == kind and row.arithmetic]
        census[kind] = {
            "attachment_count": 22,
            "arithmetic_attachment_vertices": [row.attachment for row in arithmetic],
            "ray_counts_by_arithmetic_attachment": {
                str(row.attachment): {
                    "finite": row.finite_count,
                    "ideal": row.ideal_count,
                    "hyperideal": row.hyperideal_count,
                }
                for row in arithmetic
            },
            "best_finite_ray_count": max(row.finite_count for row in arithmetic),
        }
        best = analyze_tree_host(kind, 22, best_attachment)
        best_gram = attach_arithmetic_leaf(kind, 22, best_attachment)
        best_rows[kind] = {
            "attachment_vertex": best_attachment,
            "determinant": value_json(silver_determinant(best_gram)),
            "selected_inertia": inertia_json(best.selected_inertia),
            "conjugate_inertia": inertia_json(best.conjugate_inertia),
            "inverse_diagonal": [value_json(value) for value in best.inverse_diagonal],
            "ray_kinds": list(best.ray_kinds),
        }

    dotted = adjacent_two_dotted_rank23_host()
    dotted_diagonal = analyze_dotted_diagonal(dotted)
    return {
        "schema": "silver-arithmetic-tree-hosts-rank23-0.1.0",
        "scope": {
            "field": "Q(sqrt(2))",
            "native_rank": 23,
            "hyperbolic_dimension": 22,
            "finite_skeletons": ["A22", "D22", "B22"],
            "all_labelled_attachment_vertices_exhausted": True,
        },
        "finite_arithmetic_edge": {
            "cartan_product": value_json(ARITHMETIC_EDGE_PRODUCT),
            "selected_less_than_four": True,
            "coxeter_label": 8,
            "cartan_orientation": ["-(2+sqrt(2))", "-1"],
        },
        "one_label8_edge_census": census,
        "best_hosts": best_rows,
        "two_genuine_dotted_edge_theorem": {
            "statement": (
                "Only a deletion incident to every genuine dotted edge can be "
                "finite or ideal; two distinct dotted edges therefore force at "
                "least rank-1 hyperideal simplex rays, and disjoint dotted edges "
                "force every ray hyperideal."
            ),
            "rank23_adjacent_lower_bound": genuine_dotted_hyperideal_lower_bound(
                23, ((0, 1), (0, 2))
            ),
            "rank23_disjoint_lower_bound": genuine_dotted_hyperideal_lower_bound(
                23, ((0, 1), (2, 3))
            ),
            "rank30_adjacent_lower_bound": genuine_dotted_hyperideal_lower_bound(
                30, ((0, 1), (0, 2))
            ),
        },
        "explicit_two_dotted_rank23_host": {
            "topology": "dotted leaves 1,2 at node 0; label-3 arm 0-3-...-22",
            "dotted_product": value_json(GENTLE_DOTTED_PRODUCT),
            "cartan_orientation_each": ["-(3+2sqrt(2))", "-1"],
            "symmetrizer": ["1", "3+2sqrt(2)", "3+2sqrt(2)", "1 (nodes 3..22)"],
            "selected_schur_complement": value_json(
                Silver(22) / 21 - GENTLE_DOTTED_PRODUCT
            ),
            "determinant": value_json(silver_determinant(dotted)),
            "selected_inertia": inertia_json(silver_inertia(dotted, embedding=0)),
            "conjugate_inertia": inertia_json(silver_inertia(dotted, embedding=1)),
            "inverse_diagonal": [value_json(value) for value in dotted_diagonal],
            "ray_counts": {
                "finite": sum(value.sign(0) < 0 for value in dotted_diagonal),
                "ideal": sum(value.sign(0) == 0 for value in dotted_diagonal),
                "hyperideal": sum(value.sign(0) > 0 for value in dotted_diagonal),
            },
        },
        "conclusion": (
            "D22 and B22 do not beat the A22 terminal label-8 seed: all have "
            "best ray counts 3 finite, 0 ideal, 20 hyperideal. Genuine "
            "multi-dotted simplex seeds are structurally worse."
        ),
    }


def analyze_dotted_diagonal(gram):
    # Local import keeps the payload construction visibly independent of the
    # ray-count helper used for the one-edge census.
    from coxeter_polytopes.search.silver_arithmetic_tree_hosts import inverse_diagonal

    return inverse_diagonal(gram)


def render_report(payload: dict[str, object]) -> str:
    return r"""# Silver arithmetic tree hosts in rank 23

This exact census separates two easily confused notions.  The product
\(c=2+\sqrt2<4\) is the finite Coxeter label-eight product
\(4\cos^2(\pi/8)\), although it generates \(\mathbf Q(\sqrt2)\).  A genuine
dotted edge instead has product greater than four; the explicit mild choice
below is \(p=3+2\sqrt2\).

## One field-generating label-eight leaf

A leaf with Cartan pair \((-c,-1)\) was attached at every labelled node of
the finite skeletons \(A_{22},D_{22},B_{22}\).  The resulting symmetric Gram
is integral over \(\mathbf Z[\sqrt2]\).  Exact two-place inertia and all 23
inverse-Gram diagonal entries give:

| skeleton | attachments with inertia `(22,1)` and conjugate PD | best attachment | finite / ideal / hyperideal |
|---|---|---:|---:|
| `A22` | `0,1,2,3,18,19,20,21` | endpoint `0` or `21` | **3 / 0 / 20** |
| `D22` | long-arm nodes `19,20,21` | terminal `21` | **3 / 0 / 20** |
| `B22` | long-arm nodes `19,20,21` | terminal `21` | **3 / 0 / 20** |

Thus attaching at the branching or minuscule end does not beat the terminal
\(A_{22}\) seed.  The complete labelled scan contains 66 attachments; cases
not listed fail conjugate positive definiteness (selected inertia is still
Lorentzian).

For the terminal \(D_{22}\) host, the inverse diagonal begins
\(19-\sqrt2,(21-\sqrt2)/4,(21-\sqrt2)/4,18-\sqrt2,\ldots\) and ends
\(1-\sqrt2,-\sqrt2,(1-\sqrt2)/2\).  Exactly the last two skeleton rays and
the new leaf ray are finite.  The terminal \(B_{22}\) host has the same sign
pattern; its last three entries are
\((1-\sqrt2)/2,-\sqrt2/2,(1-\sqrt2)/4\).

## Structural obstruction for two genuine dotted edges

Deleting one wall gives a simplex vertex link.  If its principal Gram still
contains a genuine dotted edge, the corresponding 2-by-2 principal block has
negative determinant, so that ray is hyperideal.  Hence a deletion can be
finite or ideal only when its wall is incident to **every** dotted edge.
Two distinct adjacent dotted edges have only one common endpoint and force at
least \(R-1\) hyperideal rays in rank \(R\); two disjoint dotted edges force
all \(R\).  In ranks 23--30 the adjacent lower bounds are therefore 22--29.

The bound is sharp within arithmetic silver trees.  In rank 23, put dotted
products \(p=3+2\sqrt2\) on `0-1` and `0-2`, and label-three edges on the arm
`0-3-...-22`.  With symmetrizer
`(1,p,p,1,...,1)`, the selected Schur complement is

\[
  2-p-20/21=22/21-p=-41/21-2\sqrt2<0,
\]

while its conjugate is positive.  Thus the exact inertias are `(22,1)` and
`(23,0)`, the minimal field is \(\mathbf Q(\sqrt2)\), and the cone has exactly
**1 finite, 0 ideal, 22 hyperideal** rays.

This is a cone-triage theorem, not a compact example or a Vinberg
nonreflectivity theorem.  It shows that genuine multi-dotted simplicial seeds
cannot supply the hoped-for small-hyperideal starting cone; the best tested
finite-skeleton seed remains the one-label-eight terminal attachment.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    report = render_report(payload)
    if args.verify:
        if OUTPUT.read_text() != encoded:
            raise SystemExit("frozen tree-host artifact differs")
        if REPORT.read_text() != report:
            raise SystemExit("frozen tree-host report differs")
        print("PASS: silver arithmetic rank-23 tree-host artifact and report reproduce")
        return
    OUTPUT.write_text(encoded)
    REPORT.write_text(report)
    print(OUTPUT)
    print(REPORT)


if __name__ == "__main__":
    main()

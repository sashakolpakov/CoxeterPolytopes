#!/usr/bin/env python3
"""Freeze the exact golden two-dotted four-node tree census."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from coxeter_polytopes.classification.golden_two_dotted_tree_packets import (
    golden_two_dotted_tree_census,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/golden_agcat_two_dotted_four_node_trees.json"
REPORT = ROOT / "results/golden_agcat_two_dotted_four_node_trees.md"


def _golden_json(value) -> dict[str, str]:
    return value.to_json()


def build_payload() -> dict[str, object]:
    census = golden_two_dotted_tree_census()
    return {
        "schema": "golden-agcat-two-dotted-four-node-trees-0.1.0",
        "scope": {
            "field": "Q(sqrt(5))",
            "vertices": 4,
            "connected_graph": "tree",
            "dotted_edge_count": 2,
            "finite_nonorthogonal_edge_count": 1,
            "not_covered": (
                "two-dotted four-node packets with cycles; packets with at least "
                "three dotted edges; full rank-23--30 compatibility; orders, "
                "conductors, stable lattices, and Vinberg hosts"
            ),
        },
        "alphabet": {
            "dotted_value_count": census.dotted_value_count,
            "finite_nonzero_value_count": census.finite_nonzero_value_count,
            "dotted_values_sorted_by_conjugate_sha256": (
                census.dotted_conjugate_order_sha256
            ),
        },
        "proof": {
            "cycle_gate": "vacuous because every packet is a tree",
            "integral_cartan_gauge": (
                "orient every tree edge as (-1,-c); propagation gives an integral "
                "Cartan matrix and a positive symmetrizer"
            ),
            "star": (
                "selected inertia is Lorentzian automatically; conjugate PD iff "
                "x'+y'+z'<4; x,y are unordered"
            ),
            "disjoint_dotted_path": (
                "for endpoint dotted products x,z and middle finite y, selected "
                "Lorentzian inertia is equivalent to (4-x)(4-z)-4y<0; "
                "conjugate PD is equivalent to the conjugate determinant >0"
            ),
            "adjacent_dotted_path": (
                "with remote dotted x, adjacent dotted y, finite z, selected "
                "inertia is Lorentzian and conjugate PD iff "
                "(4-x')(4-z')-4y'>0; x,y have ordered structural roles"
            ),
            "counting": (
                "exact conjugate ordering, monotone two-pointer inequalities, and "
                "an exact two-dimensional Fenwick dominance join; "
                "no floating-point comparisons"
            ),
        },
        "canonical_tree_packets_before_gram_gate": census.canonical_before_gram,
        "two_place_retained_count": census.retained_count,
        "retained_by_topology": {
            "star": census.star_retained_count,
            "disjoint_dotted_path": census.disjoint_path_retained_count,
            "adjacent_dotted_path": census.adjacent_path_retained_count,
        },
        "rows": [
            {
                "finite_label": row.finite_label,
                "finite_product": _golden_json(row.finite_product),
                "star_canonical_before_gram": row.star_canonical_before_gram,
                "star_retained": row.star_retained,
                "disjoint_path_canonical_before_gram": (
                    row.disjoint_path_canonical_before_gram
                ),
                "disjoint_path_retained": row.disjoint_path_retained,
                "adjacent_path_canonical_before_gram": (
                    row.adjacent_path_canonical_before_gram
                ),
                "adjacent_path_retained": row.adjacent_path_retained,
            }
            for row in census.rows
        ],
        "conclusion": (
            "Exactly-two-dotted tree packets are abundant, not empty.  This is a "
            "complete local four-node tranche and neither a compact example nor a "
            "high-dimensional nonexistence theorem."
        ),
        "discovery": None,
    }


def report() -> str:
    return """# Golden AGCat: exactly-two-dotted four-node trees

This certificate closes the complete connected four-node **tree** tranche
with exactly two dotted edges.  It uses the full 6,869-value dotted golden
alphabet and all five nonorthogonal finite Coxeter products.

Modulo graph automorphisms there are 471,865,955 canonical assignments before
the two-place Gram gate.  Exact conjugate-order counting retains 97,172,716:

- 25,023,281 stars;
- 81 paths with disjoint dotted edges (a thin determinant-sign locus);
- 72,149,354 paths with adjacent dotted edges.

The cycle-squareclass condition is vacuous on a tree, and an integral Cartan
gauge is obtained by orienting every edge product `c` as `(-1,-c)`.  At the
selected embedding the retained star and adjacent-path cases are automatically
Lorentzian.  At the conjugate embedding their positive-definiteness conditions
reduce exactly to

```text
star:          x' + y' + z' < 4,
adjacent path: (4-x')(4-z') - 4y' > 0.
```

The counts use exact arithmetic, monotone two-pointer joins, and an exact
two-dimensional dominance join, with no floating point.  Their abundance
shows that local AGCat compatibility does
not itself prove golden nonexistence.  Cyclic two-dotted packets, packets with
more dotted edges, global ranks 23--30, orders/conductors/stable lattices, and
complete Vinberg hosts remain open.  This certificate is neither a compact
example nor a global nonexistence theorem.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    payload_text = json.dumps(build_payload(), indent=2, sort_keys=True) + "\n"
    if args.verify:
        if OUTPUT.read_text() != payload_text:
            raise SystemExit("frozen two-dotted tree artifact differs")
        if REPORT.read_text() != report():
            raise SystemExit("frozen two-dotted tree report differs")
        print("PASS")
        return
    OUTPUT.write_text(payload_text)
    REPORT.write_text(report())
    print(OUTPUT)
    print(REPORT)


if __name__ == "__main__":
    main()

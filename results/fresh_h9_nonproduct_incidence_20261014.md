# Compact-H9 incidence search beyond the product proposals

No compact H9 example was found. The sixteen earlier 216-vertex proposals are now excluded for **every** Coxeter labeling and dotted distance by their combinatorial type. A new constructive attempt used globally polytopal, nonproduct incidences with sixteen facets and overlapping missing faces of orders four and five. All 51 relevant multiwedges over a nine-gon fail an exact necessary Coxeter-label test. These are restrictions on the specified combinatorial families, not a general obstruction in dimension nine.

The earlier [incidence-repair checkpoint](fresh_h9_incidence_repair_20261014.md) correctly treated its sixteen models as combinatorial proposals and its four metric trials as a bounded scout. The new type identification makes further metric work on those same models unnecessary.

| Family | Exact scope | Outcome |
|---|---|---|
| Sixteen previous models | Every saved maximal face, not just the 216 vertex count | Each is `Delta2 × Delta2 × Delta2 × I × I × I`; the whole type is impossible |
| Overlapping Lannér-four pilot | Two 16-facet label scaffolds, with a prescribed finite vertex and free right/dotted incidences | Solver-reported UNSAT; exploratory only |
| Sixteen-facet multiwedges over a nine-gon | All 6,435 positive weight compositions; 51 dihedral types remain after the missing-face order bound | Every remaining type has zero compatible mandatory Lannér/vertex states |

## Complete identification of the old proposals

The [stdlib checker](../scripts/fresh_h9_nonproduct_incidence_20261014_product_verify.py) finds three missing pairs and three missing triples in each model. These six sets are disjoint and partition the fifteen facets. It then checks that the saved maximal faces are **exactly** all sets formed by omitting one element of each of those six sets. Equality of the maximal-face lists proves equality of the entire downward-closed complexes and excludes hidden larger missing faces. The [certificate](fresh_h9_nonproduct_incidence_20261014_product_verify.json) records the literal factors of all sixteen models.

Alexandrov's Theorem B states that a compact hyperbolic Coxeter product of simplices must have at most three factors or be a cube. Our type has six factors, three of them triangles, so it is excluded. Equivalently, its disjoint minimal nonfaces give the forbidden product of Lannér diagrams in Theorem A. This conclusion is unrestricted in angle labels and dotted distances. [Alexandrov, *Lannér diagrams and combinatorial properties of compact hyperbolic Coxeter polytopes*, Theorems A–B](https://arxiv.org/html/2203.07248v2).

## A changed, globally polytopal starting family

The positive construction attempt uses the rational convex nine-gon with vertices `(i,i²)`, `0 <= i <= 8`. Its inequalities are

```text
(2i+1)x - y <= i(i+1),  i=0,...,7;
y - 8x <= 0.
```

A wedge over `a·x <= b` replaces that inequality by `a·x + z <= b` and `a·x - z <= b`, leaving the other inequalities independent of `z`. A bounded simple polytope remains bounded and simple: above a point of the old polytope the new coordinate satisfies `|z| <= b-a·x`. Repeating this elementary operation seven times gives an actual simple Euclidean nine-polytope with sixteen facets. The construction supplies global incidence closure before any Coxeter metric is sought.

Let `w_i >= 1` be the number of new facets over the i-th original side. Then `sum(w_i)=16`. In the dual nerve, replace original side i by a block `B_i` of size `w_i`. Its minimal nonfaces are precisely `B_i ∪ B_j` for nonadjacent sides of the nine-gon. A maximal face contains both complete blocks above an adjacent pair of sides and all but one element of every other block. Every maximal face has nine elements.

This is a different family from the prior fourteen-facet stellar/multiwedge schedules: both the base polygon and the facet budget have changed. It is also not a product of simplices: the minimal nonfaces overlap rather than partitioning the facet set.

In a compact Coxeter polytope every minimal nonface carries a Lannér diagram, whose order is at most five. Hence every nonadjacent weight sum must satisfy `w_i+w_j <= 5`. The checker enumerates all `binomial(15,8)=6435` positive compositions, imposes this necessary bound, and reduces by the dihedral symmetries of the nine-gon. Exactly 51 types remain; all have weights at most three. Their distribution is 19 types with two size-three blocks, 28 with one, and four with none.

For every type, the [exact script](../scripts/fresh_h9_nonproduct_incidence_20261014_multiwedge.py) constructs all maximal faces, cross-checks them against all nine-element sets avoiding the minimal nonfaces, and verifies that each eight-element ridge lies in exactly two maximal faces. There are between 272 and 448 actual Euclidean vertices per type. The [complete archive](fresh_h9_nonproduct_incidence_20261014_multiwedge.json) contains the weights, blocks, minimal nonfaces, and vertex incidences.

## Exact Coxeter-label attempt

Each missing four- or five-set must receive one of the complete classical Lannér diagrams of that order. The script regenerates all 108 labeled order-four placements and all 420 labeled order-five placements using the previously audited [classical-table implementation](../scripts/independent_h9_14_record30_a1_labels.py), whose hash is bound in the new archive.

The join enforces agreement on shared wall pairs. At every stage it tests the restrictions of every actual nine-wall vertex to the currently labeled walls. Unassigned meeting pairs are set to label two, giving the weakest possible nonnegative adjacency matrix. For normalized Gram `G=I-A`, with `A` nonnegative, positive definiteness is equivalent to `rho(A)<1`. Increasing an unassigned Coxeter label can only increase `A` entrywise and cannot repair a failed test. This justifies rejecting a partial assignment without enumerating the remaining labels.

All finite tests are exact finite-diagram classifications. Missing triples, arbitrary labels on pairs outside the selected Lannér diagrams, and dotted distances are relaxed, not silently fixed to a bounded alphabet. Consequently an empty join excludes every compact Coxeter realization of that combinatorial type.

All 51 joins become empty. Forty-eight types fail after three mandatory missing-face placements and three after four. The largest intermediate state set has only 420 elements, so the 10,000-state safety limit is never reached. No numerical predicate, SAT solver, Gram-rank fit, or arithmetic-lattice assumption enters this family exclusion. No metric candidate survives to a signature or compactness check.

The separate [two-recipe overlapping-cycle pilot](fresh_h9_nonproduct_incidence_20261014.json) is a weaker exploratory result. It prescribed two overlapping `(3,3,3,4)` Lannér cycles, including a new sixteenth facet, and a literal nine-node finite vertex, then allowed additional right pairs to become dotted. Both bounded solver calls returned UNSAT. No exported solver proof or broader exclusion is asserted for that pilot.

## Reproduction and limits

```sh
python scripts/fresh_h9_nonproduct_incidence_20261014_product_verify.py --verify
python scripts/fresh_h9_nonproduct_incidence_20261014_multiwedge.py --limit 51 --max-states 10000 --verify
```

Both saved replays passed. A separate [source and proof review](fresh_h9_nonproduct_incidence_20261014_multiwedge_review.md) checked the multiwedge formula, enumeration scope, source manifests, and all-label relaxation; it did not repeat the label joins. The product exclusion depends on the cited published theorem; the multiwedge label exclusion depends on the classical complete Lannér tables and the finite-diagram classification implemented in the bound helper. Neither is a proof that all nonproduct sixteen-facet nerves are impossible. Nonpolygon starting polytopes, other facet budgets, and different constructions remain outside the tested scope. The result supplies a concrete reason to change the starting incidences again, not a claim of a new polytope or a new universal dimension bound.

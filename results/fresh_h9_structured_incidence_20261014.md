# Structured product sections for compact H9

No compact H9 example was found. This round replaces random Gale proposals by explicit rational hyperplane sections of ten-dimensional products of simplices. It tests **9,846 parameter recipes**. All **7,784 recipes having at least fourteen actual facets** are excluded by exact necessary Coxeter conditions. The remaining 2,062 recipes have fewer facets and are explicitly outside this search's scope, rather than being excluded here.

The [producer](../scripts/fresh_h9_structured_incidence_20261014.py), [complete ledger](fresh_h9_structured_incidence_20261014.json), and [separate geometry audit](fresh_h9_structured_incidence_20261014_independent.json) preserve the literal recipes and all counterexamples. No earlier archive was changed. These bounded families do not exhaust product sections or compact H9 incidences.

## Correction to the proposed nonuniform-Gale direction

Nonuniform Gale data do not categorically provide new simplicial polytope types. A simplicial polytope admits a sufficiently small generic perturbation of its vertices preserving its face lattice. Every old facet has an affinely independent vertex set and strictly positive slack at every other vertex; these properties persist. Conversely, if arbitrarily small perturbations created a new facet, a convergent subsequence of its normalized supporting hyperplanes would support an old face containing its fixed set of nine vertex labels. Simpliciality forces that old face to be the old facet with exactly those labels. Thus no new facet can appear sufficiently close. The finitely many nonzero determinant conditions for generic position can then be met inside this open neighborhood.

For a nine-polytope with fifteen vertices, this gives a uniform rank-five Gale realization of the same type. A nonuniform parametrization might still be a useful search device, but it is not a separate incidence class omitted by uniform proposals. The prior bounded Gale searches remain nonexhaustive for the much simpler reason that only finitely many configurations were sampled.

## The explicit families

Use barycentric coordinates in each simplex factor. A factor with coefficient list `c_f` has nonnegative coordinates `lambda_fj` summing to one. Intersect the product with

`sum c_fj*lambda_fj = h/2`,

where `h` is positive and odd, and `h <= sum max(c_f)`. All coefficients are integers, so the hyperplane avoids every source vertex. Every section is a bounded, full-dimensional nine-polytope. Its vertices are precisely the intersections with crossing source edges. The source product is simple; each crossing edge has exactly nine incident source facets, making the section simple too.

| Family | Coefficient recipes | Parameter recipes | Large missing face | Small-face recipes / types |
|---|---|---:|---:|---:|
| Ten-cube | Sorted ten weights in `{1,2,3,4}`, lists `(0,w)` | 3,645 | 3,025 | 22 / 6 |
| Five triangles | Five lists `(0,a,b)`, `1<=a<b<=4`, modulo factor permutations | 2,158 | 1,664 | 0 / 0 |
| Mixed simplex products | The four dimension partitions below, arithmetic-progression lists with weights in `{1,2,3}` | 4,043 | 3,052 | 21 / 10 |

The four mixed dimension partitions are `(1,1,1,3,4)`, `(1,1,2,2,4)`, `(1,1,2,3,3)`, and `(1,2,2,2,3)`. Equal-dimensional factors are permuted safely; the weight assigned to a dimension-d factor gives coefficients `(0,w,...,dw)`. Reversing the entire functional identifies the opposite half of the level range. For triangles it also reverses each coefficient list to `(0,b-a,b)`, which stays in the same finite family.

The large-face column counts only recipes with at least fourteen actual facets. The respective numbers outside that facet scope are 598, 494, and 970.

A set of source facets means setting the corresponding barycentric coordinates to zero. If any factor loses every coordinate, the intersection is empty. Otherwise its functional range is exactly the sum of the factorwise minimum and maximum remaining coefficients. This gives an exact face-feasibility oracle. Each archived large missing face is empty, while deleting any one of its constraints makes it feasible. Its order exceeds five, contradicting the necessary Lannér bound for a compact Coxeter polytope.

## Label gates and the structural diagnosis

The 43 recipes passing the missing-face bound reduce to sixteen nerve types, with every duplicate accompanied by an explicit vertex-label permutation. Six are excluded by a literal induced obstruction from the previously audited dictionary. The other ten undergo complete joins of mandatory Lannér-four/five diagrams with all actual nine-facet vertex restrictions. Every join reaches zero states without truncation; the largest stored intermediate layer has 420 states. No necessary full-label survivor reaches a metric fit.

The sixteen-type count uses the existing colored-clutter canonicalizer. The separate audit checks every displayed duplicate isomorphism, but does not implement a second graph-isomorphism canonicalizer. No claim that these sixteen types were absent from every earlier catalogue is needed or made.

The joins leave every missing pair's dotted distance unspecified. They impose no finite label on a disjoint facet pair. Unassigned meeting pairs are relaxed to right angles only for a necessary positive-definiteness test; Perron–Frobenius monotonicity makes this safe for exclusion. The imported complete Lannér tables, finite-diagram classifier, and prior local-obstruction proofs remain explicit dependencies. No fifteen-facet-specific label restriction is applied to the twenty-facet cube sections.

All sixteen small-missing-face nerves have a nontrivial product decomposition: their minimal-nonface clutters split into disjoint connected components. This is an incidence statement, not an assertion that any hypothetical Coxeter metric makes those factors orthogonal. The actual exclusions use the recorded diagram gates. It explains the bounded search's failure: the simple coefficient patterns only retained small missing faces when they left a combinatorial product factor.

## Controlled moves through large intermediate missing faces

The [two-flip diagnostic](fresh_h9_structured_incidence_20261014_twoflip.json) tries to escape that product structure. For every saved type it exhausts first bistellar moves of sizes `2/8`, `3/7`, or `4/6` on the ten-facet carriers of source-product vertices. Exact star equalities certify 1,804 such moves. Their temporary minimal nonface has order eight, seven, or six respectively.

To return to the compact missing-face bound in one more flip, that temporary minimal nonface must cease to be missing. The second flip must therefore use it as its old missing set. The diagnostic checks every distinct source-vertex carrier containing it: **23,410 possibilities**. None has the required star. Every rejection includes a literal missing or unexpected facet witnessing the failed equality. The [separate stdlib replay](fresh_h9_structured_incidence_20261014_twoflip_independent.json) verifies all first moves, temporary minimal nonfaces, and second-carrier witnesses.

This closes only the displayed two-step repair mechanism. It says nothing about three or more flips, arbitrary ten-vertex carriers, or other realizations. No new abstract sphere from this diagnostic is counted as a realized polytope; there is no successful second move to realize.

## Verification and remaining scope

The [independent checker](../scripts/fresh_h9_structured_incidence_20261014_independent.py) imports neither the producer nor the label engine. It reconstructs the recipe list, checks 64,366 large-face feasibility assertions, enumerates 32,196 rational edge-intersection vertices across all 43 small-face recipes, recovers every minimal nonface from the complete downward face complex, and verifies all duplicate permutations and inherited-gate applications. It does not independently rederive the inherited diagram arithmetic. The primary saved replay repeats the ten complete label joins. Generation and saved-result replays pass for both the primary and the independent certificates, and for both two-flip checkers.

```sh
PYTHONPATH=src:. python scripts/fresh_h9_structured_incidence_20261014.py --verify
python scripts/fresh_h9_structured_incidence_20261014_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_structured_incidence_20261014_twoflip.py --verify
python scripts/fresh_h9_structured_incidence_20261014_twoflip_independent.py --verify
```

There is no unresolved saved label state, state-limit rejection, or numerical metric failure in this package. A next proposal should explicitly break the product-clutter decomposition and survive the order-five missing-face bound before metric work. More generic Gale samples, nonuniformity by itself, or another extension of this coefficient grid would not address the diagnosed obstacle. A realization-guided sequence with at least three nontrivial incidence moves is an untested possibility, not a claimed construction or completeness result.

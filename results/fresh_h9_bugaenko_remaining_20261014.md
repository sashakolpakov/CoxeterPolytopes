# Complete closure of the old surgery family and a changed operation order

No new compact H9 example was obtained. The **entire previous 14-facet surgery family is now excluded**: wedge any facet of the specified 11-facet Bugaenko compact H8 incidence, then truncate any two proper faces. The conclusion allows all finite Coxeter labels and arbitrary dotted distances, with no restriction to arithmetic metrics. A changed construction that wedges a newly created truncation facet gives a separate, still largely open family.

## Complete old-family theorem

The [frozen incidence census](fresh_h9_bugaenko_surgery_20261014_incidence.json) contains 73,513 safe ordered recipe representatives for two cuts of codimensions two through five. Of these, 67,527 have a minimal nonface of order greater than five and cannot be compact Coxeter polytopes. The remaining 5,986 recipes represent 2,953 distinct incidence types.

Every one of those 2,953 types now has an exact exclusion:

| Proof source | Distinct types |
|---|---:|
| New two-Lannér-block inertia theorem | 2,365 |
| Previously checked complete searches, outside that theorem | 112 |
| Newly replayed complete searches, outside that theorem | 476 |
| Total | **2,953** |

The [full-family ledger](fresh_h9_bugaenko_remaining_20261014_closure.json), [checker](../scripts/fresh_h9_bugaenko_remaining_20261014_closure.py), assigns each type a completed proof and checks the source hashes. There are **no unresolved types in this old catalogue**. This supersedes the partial counts in the [previous checkpoint](fresh_h9_bugaenko_closure_20261014.md), whose files remain unchanged.

## The structural improvement

Let A and B be disjoint minimal nonfaces of one of the order pairs

```text
(4,3), (5,3), (4,4), (5,4), (5,5).
```

Suppose their union contains no further minimal nonface. Then every union of proper subsets of A and B is an actual face and must have a spherical Coxeter diagram. The two full blocks themselves are Lannér diagrams. These conditions alone force their combined Gram matrix to have at least two negative eigenvalues, so it cannot occur in a hyperbolic Coxeter diagram.

The [exact table](fresh_h9_bugaenko_remaining_20261014_block_gate.json), [source](../scripts/fresh_h9_bugaenko_remaining_20261014_block_gate.py), proves this directly, without invoking a classification of whole products of simplices. Up to internal relabeling there are 44 triangle categories, nine order-four Lannér diagrams and five order-five diagrams. Their relevant pairs give 767 internal type pairs.

For each pair, assigning the least nonright cross label, 3, determines which cross pairs could possibly be nonright. Failure forces a cross pair to be right by monotonicity of the nonnegative cosine matrix. Every remaining cross-label category is then enumerated and tested against every maximal proper block union. Apart from the all-right choice, only **42 completions** survive: respectively 10, 27, 0, 2 and 3 for the five order pairs above. Every surviving cross label is 3.

The all-right case has two orthogonal negative blocks. For each other completion, use the first Lannér block as an invertible block A of negative index one. The exact Schur complement has strictly negative determinant, hence at least one negative eigenvalue of its own. Inertia additivity therefore gives at least two negative directions in the whole matrix.

All integer labels at least seven are covered, rather than cut off. A finite spherical diagram can contain such a label only in an isolated rank-two component, and the triangle Lannér test is constant throughout this category. Each of the 42 surviving completions contains at most one such internal label. Set `x=cos(pi/m)`. The elementary cubic for `cos(pi/7)` gives `9/10 < x < 1`; every relevant Schur determinant has strictly negative exact Bernstein coefficients on the larger closed interval `[9/10,1]`. This proves negativity for every `m>=7`.

The catalogue application checks the actual vertex list for every proper block union. It applies to 2,365 types, including all eight previously limited types. The [separate Fraction audit](fresh_h9_bugaenko_remaining_20261014_block_independent.json), [source](../scripts/fresh_h9_bugaenko_remaining_20261014_block_independent.py), rebuilds all cross-category tests with a separate finite-tree classifier, and all 42 Schur determinants with squarefree-radical arithmetic, Galois-conjugate inversion and permutation determinants. It imports neither the producer nor SymPy arithmetic. All signs and all 2,365 applications pass. A further [source/proof review and witness check](fresh_h9_bugaenko_remaining_20261014_block_gate_review.md) also passes. The classical complete Lannér tables remain an explicit shared input.

## Searches needed outside the structural theorem

The newly untouched ranges 650–1149 and 1150–1949 were tested at 500 evaluated states per type. They gave 499 and 799 exact exclusions, with cases 995 and 1506 respectively limited. The structural theorem closes 995 and skips most of the remaining untouched catalogue. Testing the 86 untouched types outside that theorem gives 85 exact exclusions and one limited type, 2421.

The last two types, 1506 and 2421, close at **218 and 421 states**, still using a 500-state limit. This uses a stronger exact relation step: discard Lannér table rows already incompatible with a spherical vertex; enforce consistency of overlapping table projections; and force any label common to every remaining table row. These are necessary consequences of the full join, not choices of an unassigned label. The [propagated search source](../scripts/fresh_h9_bugaenko_remaining_20261014_propagated.py) and the [1506](fresh_h9_bugaenko_remaining_20261014_propagated_1506_1507.json) and [2421](fresh_h9_bugaenko_remaining_20261014_propagated_2421_2422.json) certificates are retained.

The [new saved-search audit](fresh_h9_bugaenko_remaining_20261014_search_audit.json), [source](../scripts/fresh_h9_bugaenko_remaining_20261014_search_audit.py), reproduces all 476 newly necessary nonstructural records with separate finite-diagram and modular-determinant arithmetic. It reuses the DFS, canonicalizer and relation-propagation control flow; it is not an independently written exhaustive search. The other 112 nonstructural records retain their earlier frozen separate replay. No incomplete search is used as an exclusion, and no numerical metric failure is needed for the theorem.

## Why a larger first truncation adds no new case

A first cut at a face S creates the minimal nonface S. If its order exceeds five, a second cut must be at a proper subface T of S to remove it; otherwise S survives. Writing a and b for the new facets of the first and second cuts, there is a literal incidence identity:

```text
cut S, introducing a; then cut T, introducing b
=
cut T, introducing b; then cut {b} union (S minus T), introducing a.
```

At an old vertex v containing S, both sides give `v−t+a` for `t in T`, and `v−s−t+a+b` for `s in S minus T, t in T`. Vertices outside S give the same result directly. The new second face is indeed a face after the T-cut because S was an old face.

Both T and `{b} union (S minus T)` are minimal nonfaces in the final nerve. No old missing face can shorten either one, since all their old entries lie in the original face S; a missing face involving a cannot be contained in either support. Thus final admissibility forces both orders to be at most five. A second cut of order greater than five is itself a surviving forbidden missing face. Consequently every admissible two-cut construction, including an initial cut of codimension six through nine, is already represented in the closed small-cut family.

The [normalization certificate](fresh_h9_bugaenko_remaining_20261014_cut_normalization.json), [source](../scripts/fresh_h9_bugaenko_remaining_20261014_cut_normalization.py), verifies all 6,160 literal vertex/cardinality cases for the nested identity. The general minimal-nonface argument above is part of its mathematical scope. A separate [source/proof review and fresh replay](fresh_h9_bugaenko_remaining_20261014_cut_normalization_review.md) also passes.

## A genuinely changed constructive family

The next operation order is: **truncate an original H8 face, wedge the newly created truncation facet, then truncate again**. The wedge carrier is now one of the new facets, so this is not merely a larger parameter in the excluded wedge-first construction. Final nerves are explicitly checked against the entire old catalogue.

The [new incidence census](fresh_h9_bugaenko_remaining_20261014_newfacet.json), [source](../scripts/fresh_h9_bugaenko_remaining_20261014_newfacet.py), treats all 39 original-face orbits of codimensions two through five and every second cut of those codimensions. Safe reductions use the original block stabilizers split by first-face membership, together with the swap of the two new wedge facets. It produces:

| New-order census | Count |
|---|---:|
| Ordered recipe representatives | 11,729 |
| Recipes passing the missing-order test | 6,226 |
| Distinct passing 14-facet nerves | 5,536 |
| Types already in the excluded old catalogue | 1,040 |
| Additional types excluded by the block theorem | 1,220 |
| Genuinely new types outside both exclusions | 3,276 |

Each nerve is realized by ordinary sufficiently small face truncations and a wedge of a simple convex polytope. This certifies its incidence, not a Coxeter metric.

A [32-type pilot](fresh_h9_bugaenko_remaining_20261014_newfacet_search.json), [source](../scripts/fresh_h9_bugaenko_remaining_20261014_newfacet_search.py), prioritizes more dotted parameters and fewer large Lannér constraints. All 32 finish with exact label exclusions before reaching 500 states. No full necessary-label category survives, so no numerical fitting is run. The [separate pilot audit](fresh_h9_bugaenko_remaining_20261014_newfacet_audit.json), [source](../scripts/fresh_h9_bugaenko_remaining_20261014_newfacet_audit.py), rebuilds those 32 recipes directly from vertex sets, recovers their complete downward face complexes without the producer's clutter formula, and repeats their label exclusions with separate arithmetic. It does not independently enumerate the entire new family.

The changed family therefore still has **3,244 unsearched and unexcluded incidence types**, listed explicitly in the [frozen frontier checkpoint](fresh_h9_bugaenko_remaining_20261014_checkpoint.json). No conclusion about those types, about other surgery orders, or about compact H9 polytopes in general follows from the completed old-family theorem.

## Replay

```sh
python scripts/fresh_h9_bugaenko_remaining_20261014_block_gate.py --verify
python scripts/fresh_h9_bugaenko_remaining_20261014_block_independent.py --verify
python scripts/fresh_h9_bugaenko_remaining_20261014_cut_normalization.py --verify
python scripts/fresh_h9_bugaenko_remaining_20261014_closure.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h9_bugaenko_remaining_20261014_search_audit.py --verify
python scripts/fresh_h9_bugaenko_remaining_20261014_newfacet.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h9_bugaenko_remaining_20261014_newfacet_audit.py --verify
```

The complete old-family proof, the separate structural arithmetic, the 476 newly needed saved-search records, and the nested-cut identity pass. The new-order incidence census and 32-type pilot remain separately scoped; they supply the next concrete construction frontier, not a new example.

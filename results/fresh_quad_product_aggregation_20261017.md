# Canonical quadrilateral aggregation by the full wall star

Quadrilaterals admit a canonical aggregate even when their maximal product faces overlap: group them by their **entire set of n+2 incident walls**, and count every quadrilateral in that group once. A fixed star's quadrilaterals are exactly the K2,2 subgraphs of its positive-definite two-deletion graph. Dotted stars have the simpler structure of a unique maximal compact prism. All-finite stars need not have a unique maximal product under the same local hypotheses.

The [checker](../scripts/fresh_quad_product_aggregation_20261017.py) and [saved certificate](fresh_quad_product_aggregation_20261017.json) contain a complete deletion census for the 1,399 frozen dotted records and an exact overlapping-prism example from core 4202. This is an aggregation lemma and a local geometric diagnostic, not a feasible weighting rule, a complete arithmetic face catalogue, or a new dimension bound.

## The deletion graph

Let P be a full-dimensional acute-angled hyperbolic n-polytope with a strict timelike interior point H, and let Q be an actual compact quadrilateral two-face. Its wall star R consists of its n−2 containing walls and its four side walls. These n+2 roots span the ambient Lorentz space, so their Gram matrix has inertia `(n,1,1)`. To see the spanning statement, eliminate the positive containing-wall block: the four side normals of a nondegenerate compact hyperbolic polygon span its Lorentzian three-space. Failure to span would leave a common projective direction for all four side lines and could not bound a compact polygon.

Write the unique relation as

`sum(k_i n_i)=0`,

and partition its nonzero coefficients into positive indices P+ and negative indices P−; zero indices are retained separately. Define a graph E whose edges are the pairs `{i,j}` for which the Gram on `R \ {i,j}` is positive definite.

Each such deletion is an actual finite vertex, strict against every omitted wall and every other wall of the ambient polytope. Indeed, a positive definite Gram with nonpositive off-diagonal entries has nonnegative inverse. Projecting H off its normal span therefore adds a nonnegative combination of those normals. All other wall products remain strictly negative; the norm becomes more negative. This proves the finite vertex and its future orientation without assuming that the entire chamber is compact.

Dot the root relation with this vertex x. Only the two omitted terms remain:

`k_i (x·n_i) + k_j (x·n_j) = 0`.

Both displayed wall products are strictly negative. The two coefficients must consequently have opposite signs. Neither can vanish: if one vanished, both would, leaving a nonzero relation on the positive definite retained roots. Thus E is a bipartite graph on P+ and P−, and no zero-kernel root can be omitted at a finite vertex.

## Complete bipartite subgraphs give actual compact products

For nonempty subsets `I ⊂ P+`, `J ⊂ P−`, suppose every pair in `I × J` belongs to E. Put `W=R \ (I∪J)`. This is spherical because it is a principal subset of any of the corresponding finite-vertex Grams.

The full-rank root map

`x ↦ t,  t_i=−x·n_i`

identifies the ambient vector space with the hyperplane `sum(k_i t_i)=0`. On the common-W subspace, the cone defined by the R inequalities is therefore exactly

```
t_i >= 0 for i in I∪J,
sum(i in I) k_i t_i = sum(j in J) (-k_j) t_j.
```

Its extreme rays have precisely one positive coordinate from I and one from J. They are the finite vertex rays already obtained by deleting those two roots. Every such ray is future timelike. A positive section of their cone is consequently compact and **combinatorially** `Delta_(|I|−1) × Delta_(|J|−1)`. This is not a claim of a metric product decomposition.

All walls outside R are strict at every one of these vertices. They cannot cut their convex hull. Hence this is the entire actual face of P on W, with exactly `|I||J|` vertices. Inclusion-maximal complete bipartite subgraphs describe the maximal such product faces within the star.

## The canonical grouping does not require a unique maximal product

Taking `|I|=|J|=2` proves that every K2,2 in E gives an actual compact quadrilateral with full wall star exactly R. Conversely, every quadrilateral with wall star R has four finite corners; their two-root deletions give exactly one K2,2. No two different choices give the same face: their containing-wall sets W differ, while a face's complete containing-wall set is intrinsic.

Thus the actual compact quadrilaterals of P are partitioned by their literal full wall sets R. For each R, let `q(R)` be the number of K2,2 subgraphs in E. The valid aggregate requirement is

`sum of payments at all corners of those q(R) distinct quadrilaterals >= q(R)`.

Each quadrilateral occurs once globally in this sum. At a finite vertex indexed by an omitted pair `(i,j)`, its corner in a K2,2 is marked by the other positive/negative pair `(a,b)` in that K2,2. The checker records this distinction explicitly: omitted roots identify the vertex, while the other roots identify the unordered corner pair. There are exactly `4 q(R)` corner occurrences.

This permits compensation between forced quadrilaterals without selecting a maximal product. Overlapping maximal products must not have their demands simply added, because a shared quadrilateral would then be charged twice. Isomorphic stars from different catalogue records can be deduplicated separately; the geometric partition is by actual ambient facet sets, not by a chosen diagram label or record index.

## Dotted stars have a unique maximal prism

In a compact quadrilateral star, a dotted pair can only be one of the two opposite black pairs: all whites and every adjacent black pair occur together in a finite vertex. Call the dotted pair `{a,b}`. Its two roots have the same nonzero kernel sign, since the four finite corner deletions form the original K2,2.

Every positive definite two-root deletion must remove at least one of a,b. Since its omitted roots have opposite kernel signs, one sign side of E has no usable nodes besides a,b. Let J be their common-neighbor set on the other side. It contains at least two indices. Every K2,2 in E then uses `{a,b}` and two elements of J, and all of them lie in the single maximal compact prism

`Delta_1 × Delta_(|J|−1)`.

This prism has `2|J|` vertices and `binomial(|J|,2)` quadrilaterals. Its quadrilateral faces reuse the same full n+2-wall star and the same dotted base pair: in the simplex factor, every facet either contains an edge or meets it at an endpoint. In the segment factor, the two bases meet both endpoint fibers. Therefore every wall of the prism's star already occurs as a containing or side wall of each such quadrilateral. Applying the same deletion construction to any one of them recovers exactly this prism.

Consequently dotted quadrilaterals can also be partitioned globally by their unique maximal prism. Any purported larger product containing one of these quadrilaterals has all its defining and boundary walls in that same star, for the same edge-in-a-simplex incidence reason, and would give a larger biclique; none exists. If both opposite pairs are dotted, the prism is just the quadrilateral itself.

The saved census checks all two-root deletions of every one of the [1,399 exact dotted local records](fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017.md), recognizing positive definiteness through the complete finite Coxeter diagram catalogue for their allowed labels. Exact kernel, field, and arithmetic-validity checks are inherited from that frozen input.

| Prism dimension | Input records |
|---:|---:|
| 2 | 962 |
| 3 | 356 |
| 4 | 77 |
| 5 | 4 |

There are 8,985 positive definite two-root deletions. The per-record totals are 2,532 distinct K2,2 quadrilaterals and 10,128 corner occurrences. These totals are summed over the input records; no quotient by possible full-star isomorphisms across different records is claimed. The observed maximum prism dimension five is not a general dimension bound. Adding orthogonal white spectators changes neither this deletion graph nor the prism.

## Exact all-finite overlap in core 4202

The checker independently rebuilds the literal eight-root Gram of the [existing local arithmetic witness](fresh_sqrt5_adjacent_discount_witness_20261016.md). Put `phi=(1+sqrt(5))/2`. All root squares are two except root 6, whose square is four. Nonzero off-diagonal entries are

```
G05=G13=−phi,
G23=G25=G37=G45=−1,
G67=−2.
```

The physical inertia is `(6,1,1)` and the conjugate inertia is `(7,0,1)`. All 64 Cartan entries are integral. A kernel vector is proportional to

`[−(1+phi), phi, 1−phi, 2, −phi, −2phi, 1, 2]`.

Thus `P+={1,3,6,7}` and `P−={0,2,4,5}`. The exact checker evaluates all 28 two-root deletions, constructs every finite deletion vertex from a strict center, and verifies its norm, future orientation, and strict omitted products. There are precisely two maximal bicliques with at least two nodes on each side:

```
{1,3}   × {0,2,5},
{1,3,7} × {2,5}.
```

They give actual compact triangular-prism faces on white supports `{4,6,7}` and `{0,4,6}`. Their common quadrilateral has positive black roots `{1,3}`, negative black roots `{2,5}`, and white support `{0,4,6,7}`. All six vertices of each prism are future timelike, and their pairwise products are negative. The union cannot enlarge to `Delta_2 × Delta_2`: the missing deletion `{0,7}` has inertia `(5,1,0)` rather than being positive definite.

Each prism contains three quadrilaterals, but they share one. The canonical full-star aggregate therefore has **five quadrilaterals and twenty corner occurrences**, not six demands obtained by summing the two prism counts.

This is an exact counterexample to unique maximal-product grouping from the local acute, corank-one, and finite-corner hypotheses. The two product faces are genuinely compact, but the entire ambient H6 chamber is not asserted compact, and no compact completion containing it is supplied. Accordingly it does not by itself refute a stronger uniqueness claim with additional global compact-realizability hypotheses. The full-star grouping theorem above needs no such unproved hypothesis.

Validation: generation and fresh saved replay both passed. Syntax, direct source manifests, and report links passed. Source, saved data, and this report are final and frozen. No weighting coefficients or dimension claim are produced by this package.

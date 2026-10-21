# Dotted prisms force simplex bases: an exact incidence diagnostic

The saved combined dual violates **434 necessary prism/base incidence
inequalities**, even after taking conservative minima across every merged
quadrilateral provenance. This rejects that dual distribution as a possible
global incidence distribution. It is not an optimized bound for a successor
model and does not establish a dimension exclusion.

The [checker](../scripts/fresh_sqrt5_prism_simplex_base_incidence_20261017.py)
and [certificate](fresh_sqrt5_prism_simplex_base_incidence_20261017.json)
use the frozen [combined exact dual](fresh_sqrt5_rank13_combined_disjoint_ports_20261017.md),
[full-star quadrilateral aggregation](fresh_quad_product_aggregation_20261017.md),
and [disjoint simplex-port theorem](fresh_simplex_disjoint_ports_20261017.md).

## Geometric inequality

Let P be a simple compact polytope of dimension n. For k at least 3, a
k-dimensional prism Delta1 times Delta(k−1) has two simplex bases of dimension
k−1. Each base is contained in a unique maximal simplex face S, by the
disjoint-port theorem. This includes the edge case: every edge lies in a
unique maximal simplex, which may be the edge itself. For a square (k=2),
we count **all four edge bases**, avoiding an arbitrary choice of opposite
pair when both pairs of walls are dotted.

Suppose S has dimension d. Each of its (k−1)-simplex subfaces has exactly
n−k+1 incident k-faces in P. Precisely d−k+1 of these lie inside S. A prism
cannot be a face of a simplex, so at most n−d of the k-cofaces can be the
prisms being counted. There are binomial(d+1,k) possible base subfaces of S.
Thus, for each simplex type and each k,

`sum_Q m(Q,S,k) f_Q <= binomial(d+1,k) (n−d) f_S`.

Here m counts base incidences with multiplicity. Types can identify distinct
actual simplex faces; the inequality is summed over every face of that type.
The argument does not require arithmeticity. Construction_audit separately
checked the coface count, maximal-simplex uniqueness, square convention and
orientation of the deleted normal.

## Exact full-star map

For each certified dotted prism, delete one of its two base-facet normals.
The remaining n+1 normals form the full star of the **opposite** base.
Every positive-definite one-root deletion in this smaller star is a vertex
of its containing maximal simplex. The checker recounts all these deletions,
compares them with the complete two-root deletion graph of the original
star, and records the resulting simplex dimension. It never assumes a
dimension cap. The observed dimensions are 1 through 4.

The nonsingular Lorentz component of a base star is identified by exact
weighted-graph isomorphism, including the full squared product of any
remaining dotted bond. Its finite components and external spectators give
the remaining component multiset. This identifies the full normalized star.
The comparison uses the frozen 310-core / 12,420-star catalogue maps.

All 1,335 canonical prism cores and all 81,183 flexible full-star contexts
are processed. They yield 133,346 base-star types. Only 2,981 already occur
in the saved simplex catalogue and model. The missing types comprise 99,681
maximal-edge stars and 30,684 simplex stars of dimensions 2 through 4
(23,808 triangles, 4,208 tetrahedra and 2,668 four-dimensional simplexes).
Their masses are **unknown**, not zero.

## Conservative treatment of merged rows

Of the 1,582 positive quad rows in the exact dual, 243 have different base
coefficients among their merged provenances. For each base type separately,
the checker uses the minimum coefficient across every provenance of that
row. Therefore any reported excess persists under every nonnegative split
of the row mass among those provenances. All 81,183 maps are saved explicitly.
Passing these weakened inequalities would not prove that one simultaneous
split satisfies the actual stronger incidences.

There are 881 represented-base inequalities with positive conservative
load; 434 are violated. Another 2,676 positive requirements lead to missing
base variables and are recorded without comparison to zero.

For a compact five-dimensional prism, the two bases are four-simplexes and
the capacity is 13−4=9. Two exact violations are

- quad row 26302 versus simplex row 3406: excess `38911/625000000`;
- quad row 26319 versus simplex row 3572: excess `2030743/10000000000`.

These examples have nonzero simplex masses, so the violations are not
caused merely by a zero in the saved support. The highest-mass demand-10
quad row 26448 instead has bases at simplex rows 599 and 311; both satisfy
this coarse capacity test.

The saved minimal one-row witness is quad row 77718, source context 117644,
canonical prism core 691. It forces the triangle star 3768 / simplex row
2995, whose Lorentz core has bonds 3 and 10 and whose spectators are
`4 A1 + D7`. Its positive quad mass is `15699/1000000000`; the saved simplex
mass is zero. Capacity 33 cannot accommodate that positive load.

## Validation and scope

Generation and a fresh saved replay passed. No optimizer was called. This
is exact reconstruction and separation of one saved dual, with inherited
local arithmetic certificates and a separately reviewed combinatorial
argument. It is not an independent arithmetic recheck of the source cores,
a feasible signed-transfer rule, a complete face catalogue, or a new bound
on compact Coxeter dimensions.

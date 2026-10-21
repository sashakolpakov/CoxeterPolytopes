# A four-port facet coupling simplex blocks and selected dotted prisms

Let \(U\) be four edge ports at a vertex of a compact simple Coxeter polytope. Let \(C_S\) be the port block of each maximal simplex face through that vertex; include isolated edges as singleton blocks. Select any collection of the canonical **dotted** full-star prism aggregates. Define

\[
r(U)=\sum_S (|C_S\cap U|-1)_+,
\qquad
q_{\ge2}(U)=\#\{Q:Q\text{ occupies at least two square pairs inside }U\}.
\]

Then the following joint inequality holds:

\[
\boxed{r(U)+q_{\ge2}(U)\le3.}
\]

It gives a nonnegative charge with vertex cost three. It is a facet of the six-dimensional **necessary combinatorial domain** enumerated below. An exact fractional witness satisfies the complete simplex-only packing cone, all pair capacities, and all strengthened three-port inequalities \(2s+t+q\le3\), but violates this four-port inequality. No saved global dual has been scanned here, and no optimizer or dimension claim is involved.

## Why the selected quad footprints are stars

The [full-star quadrilateral theorem](fresh_quad_product_aggregation_20261017.md) identifies every selected dotted aggregate with a unique maximal compact prism

\[
\Delta_1\times\Delta_{k-1}.
\]

At one of its vertices, there is one segment direction \(h\) and \(k-1\) simplex-base directions \(C\). The quadrilateral two-faces at that vertex occupy exactly the pairs \(\{h,c\}\), \(c\in C\). Thus its footprint is a star. Its leaves lie in the port block of the unique maximal simplex containing that base. Restricting to \(U\) preserves both properties whenever the footprint is nonempty.

The [shared-base proof](fresh_prism_simplex_shared_capacity_20261017.md) proves the relevant three-dimensional completion: a triangle and two squares incident along its two edges at one vertex complete to a triangular prism. Therefore, if \(a,b\) lie in one simplex block and \(\{h,a\},\{h,b\}\) are selected square pairs, those two squares belong to the **same canonical aggregate**. Distinct aggregates cannot share an occupied pair, because a port pair determines one actual two-face and its full wall star is unique.

The checker recounts these star footprints on all 1,399 frozen dotted templates: 6,640 prism-vertex occurrences and 10,128 square-pair occurrences. At an omitted pair consisting of one dotted-base root and one opposite-factor root, the hub is the other dotted root and the leaves are the remaining opposite-factor roots. The calculation uses precisely the vertex union of all K2,2 quadrilaterals; it does **not** substitute every positive-definite two-deletion for a prism vertex. This is a literal template check, not a new Gram or arithmetic replay.

This restriction matters. A general non-dotted product aggregate, such as \(\Delta_2\times\Delta_2\), can have a K2,2 square footprint between two triangular blocks. Such aggregates are outside the selected \(Q\) family. Nothing here excludes them.

## Four-port proof

By the [disjoint maximal-simplex port theorem](fresh_simplex_disjoint_ports_20261017.md), the simplex blocks partition the four ports. Write only their nonempty intersections with \(U\).

| Partition sizes | \(r(U)\) | Maximum \(q_{\ge2}(U)\) |
|---|---:|---:|
| \(1+1+1+1\) | 0 | 0 |
| \(2+1+1\) | 1 | 2 |
| \(2+2\) | 2 | 1 |
| \(3+1\) | 2 | 1 |
| \(4\) | 3 | 0 |

Every footprint counted by \(q_{\ge2}\) needs a hub and at least two leaves in a common simplex block. Hence the first and last rows are immediate. For \(2+1+1\), only the two outside ports can serve as hubs, and the unique-square-pair property permits at most one counted aggregate at each hub. For \(3+1\), there is only one possible hub. Any two counted footprints would share a square pair because two subsets of at least two leaves in a three-element set intersect; they would be the same aggregate.

For \(2+2\), a counted footprint uses one entire two-port block as leaves and one port of the other block as hub. Consider two proposed counted aggregates. If their hubs lie on opposite sides, their two stars share a square pair, so they are the same aggregate; their combined three-edge path is not a star, a contradiction. If their hubs lie on the same side, the hubs form a triangle pair. Taking either common leaf, triangle-plus-two-squares completion forces the two aggregates to be the same. Their combined footprint is K2,2, again not a star. Therefore at most one exists. Every row proves the boxed inequality.

## Complete finite necessary-domain enumeration

The [standard-library exact checker](../scripts/fresh_sqrt5_joint_four_port_domain_20261017.py) enumerates all 15 set partitions of four labelled ports. Within each partition it enumerates every subset of nontriangular pairs as selected square pairs. This gives 337 patterns.

For each pattern, union-find imposes every forced equality of aggregates coming from a triangle pair and two incident squares. Each resulting group must be a star whose leaves lie in one simplex block. These conditions retain 322 patterns and reject 15. The rejected cases are exactly 12 three-edge paths and three four-cycles between two triangular pairs.

For an actual selected aggregate, this procedure recovers its entire restricted footprint: two or more of its edges have a common hub and leaves in one simplex block, so the closure joins them; a one-edge footprint is a singleton group. Conversely, every closure equality follows from the geometric completion theorem. The finite enumeration therefore contains every actual selected configuration. Retained configurations need not be geometrically realizable; enlarging the necessary domain is safe for the inequality.

Use the invariant coordinates

\[
(s_2,s_3,s_4,q_1,q_2,q_3),
\]

where \(s_i\) counts simplex blocks meeting \(U\) in exactly \(i\) ports, and \(q_j\) counts selected aggregates with exactly \(j\) square pairs inside \(U\). The enumeration has 25 distinct invariant vectors, of exact affine rank six. The inequality is

\[
s_2+2s_3+3s_4+q_2+q_3\le3.
\]

Its six tight vectors are

```
(0,0,1,0,0,0)  (0,1,0,0,0,1)  (0,1,0,0,1,0)
(1,0,0,0,2,0)  (1,0,0,1,2,0)  (2,0,0,0,1,0).
```

Their exact affine rank is five. Thus this is a facet of the convex hull of the retained **invariant** vectors, not a claim that the actual geometric configuration polytope has been completely characterized.

## Fractional separation beyond all three-port cuts

Give the two disjoint simplex blocks \(\{0,1\}\) and \(\{2,3\}\) mass one each. Their complete simplex-only packing cover is the single packing consisting of these two blocks, with mass one.

For each of the four ports, give mass \(1/3\) to the two-edge selected-prism star centered there, with both leaves in the opposite simplex block. Each triangular pair has total pair load one. Each cross-block square pair has load \(2/3\). Thus all pair capacities hold.

Every three-port set contains one of the two simplex pairs. It meets exactly three of the four prism stars. Consequently the strengthened joint inequality has

\[
2s+t+q=2\cdot1+0+3\cdot\tfrac13=3
\]

on every triple. Each older mixed-wedge load is only \(1/3\). Nevertheless the four-port load is

\[
r+q_{\ge2}=2+4\cdot\tfrac13=\tfrac{10}{3}>3.
\]

This is an exact fractional relaxation witness. It is not asserted to satisfy global coface flags, common arithmetic-form conditions, or actual geometric compatibility. Its purpose is to prove that the four-port inequality is not implied by the listed local simplex-only, pair, and three-port constraints.

## Use and validation

For a fixed literal four-port set, pay a simplex corner \((|C_S\cap U|-1)_+\) and a selected dotted-aggregate corner the indicator of at least two occupied pairs in \(U\). Average the **whole inequality** over the finite forest's automorphism orbit of \(U\). The vertex cost remains three. Keep one consistent mass for each actual quad provenance; do not mix independent representatives of merged rows.

Generation and a separate fresh saved replay both passed `PASS_EXACT_SELECTED_DOTTED_PRISM_FOUR_PORT_FACET`. All domain records, forced equalities, accepted and rejected patterns, exact ranks, template footprints and fractional inequalities are saved in the [certificate](fresh_sqrt5_joint_four_port_domain_20261017.json). There are no numerical or global optimizer calls. Applying the new charge to the saved dual is a separate future diagnostic.

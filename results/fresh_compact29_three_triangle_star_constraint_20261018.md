# A universal triangle constraint at an E7³+E8 compact vertex

Let P be a compact hyperbolic Coxeter polytope in H29 and let v be a
finite vertex whose wall diagram is E7⊔E7⊔E7⊔E8. In each E7 number the
simple edges as

\[
0-1-2-3-4-5,\qquad 2-6.
\]

Consider the three incident two-faces specified by the marked pairs
(2,4), (2,5), and (3,5) in each E7. At most one of the three E7 factors
can supply a triangle among these designated faces. Consequently at most
three of the nine designated two-faces are triangles, and at least six
have at least four sides.

This is a universal local incidence statement, without an arithmeticity
assumption. It does not prove a new face-charge saving, classify the
remaining polygonal faces, or improve the compact dimension bound.
The [producer](../scripts/fresh_compact29_three_triangle_star_constraint_20261018.py),
[exact certificate](fresh_compact29_three_triangle_star_constraint_20261018.json),
and [independent saved checker](../scripts/check_compact29_triangle_star_saved_20261018.py)
exhaust the necessary triangle contact domains. This strengthens the
[earlier fixed-first-root calculation](fresh_compact29_coupled_first_roots_20261018.md)
for the three stated triangle angles: it does not assume that an actual
polytope uses a previously selected local witness.

## Actual finite corners and untouched components

Normalize every wall root to square2. The 29 roots at v are independent
and span a positive-definite space T of type E7³+E8. Compact Coxeter
polytopes are simple, so the two-face corresponding to a pair a,b has
the other 27 vertex walls as its fixed supporting walls. If it is a
triangle, let c be its third side. The two other vertices have wall
sets obtained from the old 29 roots by replacing a, respectively b,
by c. Both sets must be positive definite. The root c is not a facet
through v.

We first exclude all contacts of c with the three untouched components.
A whole E8 component occurs in both new corners. It has no connected
extension in the finite Coxeter classification, so c is orthogonal to it.

Suppose instead that c contacts an untouched E7 component F. The finite
diagram F+c must be E8: it is a connected finite diagram of rank8
containing the induced E7. Because E8 has no further connected finite
extension, c must be orthogonal to every other old root appearing in
each of the two new corners. The union of those corners contains every
old root. Thus c is orthogonal to all old roots outside F. But then all
29 old roots together with c have positive-definite Gram matrix

\[
(F+c)\ \sqcup\ (T-F),
\]

of rank30. This is impossible in signature (29,1), whose positive index
is29. Therefore c is orthogonal to every untouched E7 component too.
All possible contacts are confined to the seven roots of the marked E7.
This reduction applies to arbitrary actual third walls, including ones
not present in the old finite witness catalogue.

## Why the finite label domain is complete

Each contact of c with an old E7 root occurs in at least one of the two
new finite corners. It is therefore right or a finite Coxeter bond;
parallel and dotted products cannot occur here.

In the finite Coxeter classification a bond of label at least6 is an
isolated I2(m) component. Suppose such a contact joins c to a white root
j, different from a,b. The root j occurs in both new corners. Its old
neighbors would have to lie entirely in the omitted singleton {a} and
entirely in {b}, which is impossible because E7 is connected. If the
high bond joins c to a marked root, that root would have to be a leaf
whose only old neighbor is the other mark, omitted at that corner.
The three designated pairs are nonadjacent, so this is impossible too.

Hence the complete contact label vector belongs to {2,3,4,5}⁷. This
finite domain follows from the actual duties; it is not a sampled label
cutoff. The producer tests all 4⁷ vectors for each marked pair. Its
finite corner classifier recognizes the exact finite Coxeter diagrams
combinatorially. The independent checker uses a separate rooted-tree
canonical encoding of the finite ADE, B, F, H and dihedral catalogue.

## Height and compactness exclusions

Write J for the inverse E7 Gram matrix, and let b be the seven products
of c with its roots. Since all other old components are orthogonal to c,
its projection p into T has square bJbᵀ. If e is the future unit vector
representing v, then

\[
c=p+h e,\qquad e^2=-1,\qquad
h=-c\cdot e>0,\qquad h^2=bJb^T-2.
\]

Strict positivity follows because v is strictly inside the halfspace of
the new third wall. A negative computed h² gives no such root; h²=0
would put the new wall through v. Both cases are excluded.

A compact acute Coxeter chamber has no positive-semidefinite singular
principal wall diagram. To see this directly, choose a connected
singular component and its positive kernel coefficients. Their linear
combination u of actual roots satisfies u²=0. For a strict interior
timelike point x, u·x<0, so u is nonzero and on the future sheet. Kernel
orthogonality holds on this component, and acuteness gives u·r≤0 at
every other wall. Thus the corresponding ideal point belongs to the
closed chamber, contradicting compactness. This uses actual wall
products; it is not an assumption that every principal diagram is a
face diagram.

For every finite-eligible vector with positive h², the producer tests
all principal subsets until it finds such a singular positive-semidefinite
subdiagram, or exhausts the domain. Exact arithmetic lies in
Q(√2,√5), as forced by the already proved finite label domain. This
choice of arithmetic field imposes no restriction on the polytope's
field of definition.

|Marked pair|All vectors|Finite-eligible|Nonpositive h²|Affine obstruction|Remaining root|
|---|---:|---:|---:|---:|---:|
|(2,4)|16,384|12|3|8|1|
|(2,5)|16,384|7|3|3|1|
|(3,5)|16,384|7|3|3|1|

The saved certificate includes every finite-eligible vector, its exact
height, and an explicit singular principal subset for each affine
exclusion. For the three survivors it verifies every nonempty principal
subset. The saved checker reconstructs the entire 49,152-vector finite
filter and its deterministic digest, recomputes all heights, checks each
affine witness using a different Schur elimination order, and exhausts
all 255 principal subsets of each surviving eight-root diagram.

For each pair the sole remaining contact vector is the same:

\[
c\cdot r_5=-2\cos(\pi/5),\qquad
c\cdot r_j=0\ (j\ne5),\qquad
h^2=\frac{1+3\sqrt5}{4}.
\]

The old [exact local witness packet](fresh_compact29_small_faces_20261015.md)
also contains each of these triangle recipes. The present necessary
classification does not require sufficiency or replay its complete local
polygon geometry.

## The coupling uses the common finite vertex

Suppose two different E7 factors each supply one of the three triangles.
Their new roots c_i and c_j have orthogonal projections p_i,p_j in T.
Both have the same positive height h because their halfspaces contain
the same finite vertex v strictly. Consequently

\[
c_i\cdot c_j=-h^2,
\qquad (c_i\cdot c_j)^2=
\frac{23+3\sqrt5}{8}.
\]

These are distinct roots: their nonzero contact supports lie in different
old factors. Their mutual product is not native. The exact bounds are

\[
\bigl(2\cos(\pi/11)\bigr)^2
<\frac{37}{10}
<\frac{23+3\sqrt5}{8}
<\frac{93}{25}
<\bigl(2\cos(\pi/12)\bigr)^2.
\]

The producer and checker use directed rational Machin-formula and cosine
Taylor enclosures, with different truncation lengths. Since squared
finite bonds increase strictly with their integer label and the displayed
value is less than4, no finite, parallel, or dotted Coxeter product can
have this value. The two triangular sides therefore cannot coexist in
the same chamber. At most one E7 factor contributes a designated triangle;
there are only three designated angles per factor. Simplicity makes the
nine marked pairs correspond to nine distinct incident two-faces, proving
the stated at-most-three / at-least-six conclusion.

## What remains for a compact-bound improvement

The [rank29 savings diagnostic](fresh_compact29_savings_frontier_20261018.md)
identifies E7³+E8 as one of the restrictive rank patterns. The present
result is a necessary condition on its actual two-face incidence. The
existing charge theorem includes quadrilaterals, so nontriangular does
not mean uncharged. A saving would still require an actual star charging
inequality, a restriction on those quadrilaterals, or a compatible change
to the global incidence argument. None is asserted by this packet.

## Verification and scope

```sh
python3 scripts/fresh_compact29_three_triangle_star_constraint_20261018.py --verify
python3 scripts/check_compact29_triangle_star_saved_20261018.py
```

Generation7924/ef0e85, separate complete producer replay57352/407f35,
and independent checker24256/ee6bb9 terminated exit0/PASS. The helper
uses only Python's standard library, imports no project module, and
independently verifies the exact finite domain, heights, principal
exclusions, and native gap. The written argument supplies the actual
corner reduction and common-vertex coupling; the helper does not prove
hyperbolic geometry from first principles. Finite Coxeter classification,
acute finite-vertex simplicity, and the frozen source identities remain
explicit mathematical inputs.

The root and order16 agents each completed full mathematical reads of
the producer, report and stand-alone helper. Both checked untouched
component exclusion, high-label completeness, the positive-height and
compact affine arguments, and the common-vertex product and native gap.
The order16 reader additionally read the immediate finite-tree classifier,
ordered-field/Gram/inverse/PSD interfaces and cosine enclosures. Neither
reviewer executed the packet, replayed its finite domain, inverses,
principal tests or Taylor bounds, or independently reproved the published
finite classification or compact simplicity. Neither review claims a new
charge saving. Both passed without correction.

The [standalone TeX exposition](../manuscripts/compact29_triangle_constraints.tex)
and [lemma/reproduction map](../docs/compact29_triangle_constraints_map.md)
were prepared from this proof; the three-page PDF build passed9d7889.
The two reviews above concern this packet, not a separate full review of
the manuscript. This report is final and frozen after the completed
generation, independent replay/checker and mathematical reviews.

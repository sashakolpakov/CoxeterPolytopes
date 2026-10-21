# Affine closure makes every rank4 leaf refund safe

18 October 2026. **New written argument, awaiting independent review.**
No proof program, enumeration, optimizer, Lean process or TeX build was
run. Earlier reports and manuscripts are unchanged.

The surviving-end refund extends from B4/F4/H4 to A4 and D4. The new
point is that an exceptional E8 continuation cannot make both remaining
corners of the actual triangle unpaid: reinstating the deleted old leaf
would produce an explicit parabolic diagram of actual walls. This avoids
any assertion about the positions of Vinberg's four omitted E8 pairs.

With the original-heavy selector, the resulting valid face rule removes
every rank4-containing forest from the heavy domain. Relative to the
inherited complete 2,126-type domain, only 98 abstract types remain.
The compact dimension ceiling remains 29. The 98 types are obligations,
not excluded or asserted realizable vertex types.

## Inputs and relation to the primary rule

We use compact simplicity, positive definiteness of every actual finite
corner, the complete classification of irreducible finite Coxeter
diagrams, and the acute Gram convention: distinct outward wall normals
have nonpositive products. Normalize each root to square 2, so a simple
Coxeter edge has product -1. Every wall is strictly negative against an
interior future timelike point.

The required features of the original Vinberg weights w0 are:

* adjacent marked roots have positive weight;
* all pairs in a whole connected component of rank at most 7 have
  positive weight;
* pairs in a classical end segment of order at most 6 have positive
  weight, with the stated whisker invariance in type D;
* every pair in a whole rank4 component has weight 2;
* every actual k-gon has total weight at least 5-k.

These are the inputs of the
[reviewed surviving-end argument](compact29_rank4_surviving_end_refund_20261018.md).
The primary reference is
[Vinberg, Section 6, printed pp.80–81 and the subsequent face proof](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).
Rule IV does have four omitted E8 pairs. This proof neither transcribes
their image-only positions nor enlarges its positive set. The primary
PDF text was read; attempts to retrieve a legible p80 diagram through
the web screenshot tool failed, and a bounded read-only download returned
HTTP403. The argument below applies even if every nonadjacent E8 pair
were regarded as potentially unpaid.

## 1. An explicit positive-null obstruction for compact acute chambers

Let R be a set of actual outward normals. Suppose its Gram matrix G
has a nonzero vector v with strictly positive coordinates and Gv=0.
Put z=sum(v_i r_i). If x is an interior timelike point, then x.z<0, so
z is nonzero. Also z.z=v^T Gv=0. With the time orientation chosen by x,
z is future null. For every other wall normal s, acuteness gives

    s.z = sum(v_i s.r_i) <= 0.

Thus the ideal ray of z satisfies every wall inequality. It belongs to
the ideal closure of the chamber, contradicting compactness. In
particular the chamber cannot contain any of the following explicitly
displayed positive-null wall diagrams.

For a simply laced tree with Gram diagonal 2 and edge entries -1,
the equation Gv=0 says that twice the coefficient at a vertex equals
the sum of its neighboring coefficients. All certificates below are
checked by that identity alone. No determinant computation or affine
classification is needed.

## 2. The only A3-to-E8 continuation, including all outside factors

Let C be a whole A4 or D4 component at a finite corner. Delete a leaf b,
and write K=C minus b=A3. Let U be the union of every other whole old
component. In an actual triangle marked by a,b, where a belongs to K,
write c for the third side wall. The new corner which omits b contains
all of K, U and c and is positive definite.

Assume c contacts K. A finite diagram is a tree, so c contacts exactly
one vertex of K. All K-to-U products are zero because C was a whole
old component. If the connected component containing K+c is E8, write
the E8 tree using its arm lengths 1,2,4. A contact at the middle of K
would make that old middle vertex a fork with two terminal arms of
length 1. E8 has no such fork. Therefore c contacts an endpoint p of K;
write the other two vertices in order p-q-r. These form a terminal
three-vertex chain in E8, since no old K vertex has a contact with U.

Only the length-four arm can contain this terminal chain. Consequently
the entire connected component has exactly these edges:

    z--o--y--x
       |
       c--p--q--r

Here the arms at o have lengths 1,2,4: z is the short leaf, y-x is the
medium arm, and c-p-q-r is the long arm. All edges have label 3. The
four old outside vertices {x,y,o,z} form one whole A4 factor, not a
selected subdiagram inside a larger old factor. Indeed the whole
connected component is exactly E8, and removing c leaves precisely K
and this A4. Any further outside factors remain orthogonal to it.

This proves the full placement statement. It does not assume c is
orthogonal to the old outside factors; its necessary contact with the
outside A4 is precisely c.o=-1.

## 3. Reinstating b creates a positive-null diagram

Suppose both new triangle corners had original weight zero. At the
corner which omits a, the two marked roots b,c cannot be adjacent.
Adjacency positivity therefore forces b.c=0. The remaining products
of b with K and U are the original ones. There are just three cases.

**A4, with b attached to the far endpoint r.** The complete nine-node
diagram is the tree with arms 1,2,5 at o. Assign these coefficients:

|Vertex|o|z|y|x|c|p|q|r|b|
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|Coefficient|6|3|4|2|5|4|3|2|1|

At every vertex twice its coefficient is the sum of its neighbors.
For example 12=3+4+5 at o. These coefficients are strictly positive,
so Section 1 gives a contradiction. This is the usual affine E8 tree.

**A4, with b attached to the near endpoint p.** Take the seven-node
subdiagram on {o,c,p,y,z,b,q}, omitting x and r. The two forks o,p
are distance 2 apart, joined by o-c-p. Their two leaves are y,z and
b,q respectively. Give each vertex of o-c-p coefficient 2, and each
of the four leaves coefficient 1. The same neighbor-sum identity holds.
This is a positive-null subdiagram (the usual affine D6).

**D4.** Now b is attached to the middle vertex q of K. Take the
eight-node subdiagram on {o,c,p,q,y,z,r,b}, omitting x. The forks o,q
are joined by o-c-p-q and have leaves y,z and r,b. Give every vertex
on that joining path coefficient 2, and the four leaves coefficient 1.
Again the neighbor-sum identity holds at every vertex. This is the
usual affine D7 positive-null subdiagram.

The listed subdiagrams contain no unrecorded edges. All old C-to-U
products are zero, the new contacts are those of the actual E8 corner,
and b.c=0 was forced by the other corner's zero payment. The diagrams
therefore concern actual wall normals, rather than a formal tree which
need not embed into the chamber. Section 1 excludes every case.

Thus an E8 continuation cannot make both remaining triangle corners
unpaid, regardless of the locations of its exceptional marked pairs.

## 4. The complete leaf-pair triangle lemma

**Lemma.** Let C be a whole A4 or D4 component at a finite vertex of
a compact hyperbolic Coxeter polytope. If b is any leaf of C and a is
any other vertex of C, every actual triangle at (a,b) has another
corner of positive original Vinberg weight.

**Proof.** Keep all of U, K and c as above. If c contacts K, its finite
connected component is either classical, one of the exceptional types
of rank at most 7, F4/H4, or E8. The small whole components give positive
weight directly.

In a classical component, the old A3 has no contacts outside itself
except the contact with c. If c attaches at an endpoint, K+c is a
terminal four-node chain. It is a classical end segment; in a D-series
component the chosen set is invariant under interchange of the two
whiskers. If c attaches at the middle, that middle becomes the unique
D-series fork, and K+c is the four-node fork-end segment including both
whiskers. Hence the marked pair a,c has positive weight in either
classical case. A second contact with K would create a cycle and is
incompatible with finite classification. These statements include
outside factors attached through c.

Only E8 could evade this positivity argument, and Sections 2–3 show
that both new corners cannot have weight zero there.

Finally suppose c is orthogonal to K. If the other new corner also
had weight zero, adjacency positivity would give c.b=0. The actual
positive-definite corner K+(U+c) shows that U+c is positive definite.
Since c is now orthogonal to C, the complete root set has positive-
definite block decomposition C+(U+c), of order n+1 in Lorentz signature
(n,1). This contradicts positive index n. Therefore another positively
weighted corner exists in this last case too. QED.

The proof uses compactness at the positive-null obstruction. It is not
a finite-volume theorem in which ideal rays are permitted.

## 5. A simultaneous face rule for all five rank4 types

At any chosen set of vertices make the following changes from 2 to 1
inside the indicated whole rank4 components:

|Component|Selected marked pairs|Exact saving|
|---|---|---:|
|A4 path 0-1-2-3|All five pairs meeting an end 0 or 3|5|
|D4|All six pairs|6|
|B4, non-simple bond 0-1|Three pairs meeting the far leaf 3|3|
|F4 path 0-1-2-3|All five pairs meeting an end 0 or 3|5|
|H4, non-simple bond 0-1|Three pairs meeting the far leaf 3|3|

The last three rows are the reviewed surviving-B3/H3 rule. For A4
each selected pair contains a leaf. Every D4 pair contains a leaf.
The new lemma therefore covers every newly selected pair. The selections
are invariant under all automorphisms of the finite component. They
are made once per original vertex, not independently for its consumers.

Here is the simultaneous face proof. Every originally positive weight
stays at least 1, and every originally zero weight stays zero. An
affected triangle has a changed corner of weight 1 and, by the applicable
leaf lemma, a different originally positive corner, whose new weight
is still at least 1. Even if several corners change, the triangle
therefore pays at least 2. An affected quadrilateral retains its changed
weight 1, sufficient by itself. Longer faces retain nonnegative weights.
Unaffected faces keep the old inequality. Thus all face consumers are
covered by one rule, without a quadrilateral contact classification.

## 6. Apply only at original-heavy rank29 vertices

Use the inherited complete rank29 domain and call a vertex heavy when
its original total sigma0 exceeds 84. Apply Section 5 at exactly these
vertices. Original totals are at most 87. A heavy forest containing
any rank4 factor loses at least 3 and hence ends at most 84. Multiple
components contribute disjoint savings. All vertices originally at most
84 remain unchanged. The zero set is the original zero set.

This removes all seven rank4-containing rows of the inherited nine-row
domain. The remaining rows are:

|Ranks|Type choices|Original charge|
|---|---:|---:|
|7+7+7+8|binom(6,3)*4=80|87|
|9+10+10|3*binom(4,2)=18|85|
|Total|98||

The rank7 and rank8 types are A,B,D,E; the rank9 and rank10 types are
A,B,D. Repeated equal-rank components are unordered. These elementary
counts use the already complete typed domain; they are not a new census.
Relative to the reviewed 476 frontier, 378 further abstract types are
controlled. Relative to the original 2,126, 2,028 are controlled.

## 7. Composition with the E7 facewise refund

Use the new integer rule as the background for the
[E7-triple facewise refund](compact29_e7_triple_facewise_refund_20261018.md).
Its eligible vertices T=E7+E7+E7+E8 have no rank4 factor and are
unchanged. At a face F with m eligible corners and other weight B,
replace every eligible unit by max(0,(5-k-B)/m). The background face
inequality makes this number lie in [0,1] and verifies the final
inequality on that same actual face.

The protected weight2 source corners used by the active-state proof
have whole original charges 80,84,78. Their changing parts are
A2+H4+A1, A4+H3, D5+I2(5), respectively, added to the unchanged
2E7+E8 contribution 66. None is original-heavy, so no selected source
angle is changed. In particular the boundary case A4+H3 at exactly
84 is retained. The six active-state refunds still give charge at most
81 at T.

For an inactive T vertex, a chosen quadrilateral with zero other
background weight has precisely the old zero-weight pattern. The
old essential-quad theorem therefore applies and forces its companion
triangle with the unchanged D5+I2(7) weight2 source, of whole charge 78.
Other selected quadrilaterals either have another positive integer
payment or share at least two eligible units, giving the old full or
half refund. Longer faces fully refund the eligible unit. All three
factor refunds remain distinct and yield charge at most 171/2.

Consequently the unresolved simultaneous allowances are 79 types at
87, the T type at 171/2 (at most 81 when active), and 18 types at 85.
If N79, NT and N18 count these remaining actual excess contributors,
the old strict balance implies the necessary inequality

    3*N79 + (3/2)*NT + N18 >= 1050.

Here NT may be restricted to inactive T vertices. Vertices with
nonpositive excess have been discarded in obtaining this necessary
inequality; no lower bound on their charges is asserted.

This is a new valid weighting refinement. It does not yet bound every
vertex by 84 and does not exclude dimension 29. The still-open cases
are precisely the surviving rank7/8 and rank9/10 forests, together
with their actual simultaneous face incidences.

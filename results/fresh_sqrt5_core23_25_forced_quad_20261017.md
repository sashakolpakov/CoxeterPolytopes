# A forced baseline-two quadrilateral links the dominant prism profiles

Two specific prism types in the current rank13 arithmetic witness cannot
occur with the saved face and vertex frequencies. Their coexistence at a
common maximal simplex forces a compact quadrilateral with a vertex of type
`A4 + 3 H3`. The witness assigns that vertex type mass zero.

The resulting necessary incidence inequality is

\[
 2a+b-2s\ \le\ 3v_{A_4+3H_3}.
\tag{1}
\]

Here `a` counts the selected actual raw23 prism type (face26410), `b` the
selected raw25 type (face26880), and `s` maximal simplex type580. These are
counts of whole face occurrences, with the actual provenance and spectator
types specified below. At the frozen
[A3/D5 repaired witness](fresh_sqrt5_a3_d5_cap_joint_ports_20261017.md), the
left side is exactly `3363761/1000000000`, while the right side is zero.
Thus (1) rejects that particular incidence distribution. It does not give
a new dimension bound or an optimum for a strengthened model.

The [source](../scripts/fresh_sqrt5_core23_25_forced_quad_20261017.py) and
[exact certificate](fresh_sqrt5_core23_25_forced_quad_20261017.json) reconstruct
the literal input profiles, their common full basis, all five possible
common-vertex Schur blocks, the complete forced quadrilateral, its four
corner forests, and the exact violated inequality. No optimizer is used.

## The common simplex and the two exact profiles

The maximal4-simplex type580 has base type55301. Its complete14-root star
is the orthogonal union of the following five-root Lorentz core and
`2 A1 + A4 + H3`:

```
4 --5-- 0 --3-- 1 --3-- 2 --3-- 3
```

Every core deletion is spherical, and the labelled core has trivial
automorphism group. The two A1 supporting roots are denoted q1,q2. All
roots are initially normalized to norm2. Write phi=(1+sqrt5)/2.

The selected raw23 profile A is actual provenance483, base occurrence1513,
face26410. It has product−1 with core root3 and product−sqrt(tauA) with one
q, where

\[
 \tau_A=(7+\sqrt5)/2.
\]

Every other product with the full S-star is zero. Its base contains all
five S vertices. The resulting canonical prism has dimension5. The selected
raw25 profile B is provenance535, base occurrence1617, face26880. It has
product−1 with core root2, product−sqrt(tauB) with one q, and all other
S-star products zero, where

\[
 \tau_B=6+2\sqrt5.
\]

Its base consists of S vertices0,1,2,4, and its canonical prism has
dimension4. Both actual prism types have additional spectator forest
`A1 + A4 + H3`. Each has exactly one base with containing maximal simplex
type580; the other base has a different Lorentz core. Consequently each
selected Q occurrence contributes once to the relevant count at S.

The source checks these assertions directly through the raw diagrams,
canonical prism partitions, deleted-root/core isomorphisms, full spectator
forests, current actual-Q features, and the unchanged provenance-to-face
map. It does not replace a merged row by a convenient representative.

The full S-star is a nondegenerate basis of the14-dimensional Lorentz
space. Therefore a cap's complete contact vector determines that root
uniquely; there is no free distance or cross-product parameter. If vA,vB
are their projections to the core, exact inversion gives

\[
 v_A^2=(1-\sqrt5)/4,\qquad
 v_B^2=-1-\sqrt5,\qquad
 v_Av_B=-\phi.
\tag{2}
\]

Two A caps on different q's have product vA² in(−1,0), which is impossible
for distinct norm2 Coxeter normals. On the same q, a fixed profile gives
the same root. Thus at each S there is at most one A cap. A and B on the
same q have product

\[
 -\phi+\tfrac12\sqrt{\tau_A\tau_B}>0,
\]

contradicting acuteness. There is at most one B on each of the two q's.
If A is present, B can occur only on the other q. Spectator permutations
introduce exactly these two q choices; the other spectator automorphisms
act trivially on both complete contact vectors.

## Every mixed pair forces the same pointed quadrilateral

Suppose A and B both occur at one S, necessarily on distinct q's. Name
their cofaces qA,qB. Their product is exactly−phi by(2). The supporting
white set of the forced quadrilateral consists of core roots0,1,3,4 and
the original A4+H3 spectators. Its rank is11. The black side roots, in
cyclic order, are qA,qB,A,B.

The source checks every possible S vertex. After eliminating the whites
at S vertex2, the A,B residual Gram is

\[
 \begin{pmatrix}3/2&-\phi\\-\phi&2\end{pmatrix},
 \qquad \det=(3-\sqrt5)/2>0.
\tag{3}
\]

The corresponding determinant is negative at the other common base
vertices0,1,4; vertex3 is not in the B base. Thus vertex2 supplies the
required fourth finite corner, with no profile or length search.

More explicitly, all four rank13 blocks W+{qA,qB}, W+{qB,A}, W+{A,B}, and
W+{B,qA} are positive definite. In an actual acute chamber, projecting a
strict timelike center onto the orthogonal complement of any such block
produces a future timelike point strict against every other wall: the
block's inverse is entrywise nonnegative, and all off-diagonal wall
products are nonpositive. Hence these are four actual finite vertices,
including the new A,B vertex. This also rules out an unseen old wall
cutting off the new vertex.

After splitting off the white spectators, the minimal quadrilateral core
has roots(qA,qB,A,B,w), with w=core root3. Its only nonright bonds are
qA–A dotted tauA, qB–B dotted tauB, A–B label5, and A–w label3. Its upper
label list is

`[2,0,2,2,2,0,2,5,3,2]`.

Rescale qA by sqrt(tauA) and qB by sqrt(tauB). This gives a Gram over
K=Q(sqrt5), with those diagonal entries2tauA,2tauB, contacts−tauA,−tauB,
and the remaining diagonal entries2. It has kernel

\[
 (2,-1,4,-2,2)
\]

and, in order(qA,qB,w,A,B), exact congruence pivots

\[
 2\tau_A,\ 2\tau_B,\ 2,\ -\phi/2,\ 0.
\]

Its physical inertia is(3,1,1); the other embedding gives(4,0,1).
Every Cartan coefficient is integral in O_K. The four black signs of the
kernel alternate. The
[singular-star quadrilateral lemma](fresh_singular_quad_local_realization_20261017.md)
therefore identifies the common-white cone as exactly the compact
four-ray quadrilateral. Its four vertices are already feasible for every
ambient wall, so convexity makes every other inequality redundant on that
quadrilateral. This establishes the complete face, not merely a spherical
candidate for one corner.

The orthogonal spectators of this minimal core are `A4 + 2 H3`. The four
full corner forests are, in cyclic order,

1. `3 A1 + A4 + 2 H3`;
2. `A1 + A2 + A4 + 2 H3`;
3. `A4 + 3 H3`;
4. `3 A1 + A4 + 2 H3`.

Only the third corner has a nonright black pair. Its H3 component is
`w --3-- A --5-- B`, so its native framing5 bond has the established
baseline payment2; the other three corners have payment0. This is a
baseline2 quadrilateral, outside the completed baseline1 dotted catalogue.
The label5 refers to the original ambient wall normals. The induced angle
inside the quadrilateral has cosine squared phi²/3, and is not claimed to
be pi/5.

The construction is conditional on simultaneous actual cap occurrences.
It does not assert that the combined full star embeds in a compact ambient
chamber. Such an additional obstruction would only strengthen the incidence
conclusion.

## Counting without duplicate originating simplices

For an actual S occurrence, let Ai,Bi count the selected profiles and let
Xi=AiBi count mixed pairs. The preceding exact compatibility gives
Ai in{0,1}, Bi≤2−Ai. Therefore

\[
 X_i\ge 2A_i+B_i-2.
\]

The source checks every allowed integer case. Summing over S gives
X≥2a+b−2s.

The map from a mixed pair to its forced quadrilateral is injective. To
recover a possible origin, point the new quadrilateral at its native5
corner. In the H3 component, A is the center and B its5-edge leaf; these
roles cannot be interchanged. The3-edge leaf is w. The quadrilateral
determines its opposite old corner qA,qB. At that old corner, the same
ambient wall w is a free port of every possible originating S.

Distinct maximal simplex faces through a simple vertex have
[disjoint port sets](fresh_simplex_disjoint_ports_20261017.md). Thus there
is at most one originating S containing this w port. For that S, the
complete contact vectors determine A and B uniquely. Each forced
quadrilateral is counted once, with the ordered roles A=center and
B=5-leaf, not twice by exchanging the caps.

At a vertex of type Y=`A4+3H3`, exactly three native5 bonds are available,
one in each H3 component. Each marked pair determines one actual2-face
of a simple polytope. Thus the number of pointed forced quadrilaterals
at all Y vertices is at most3vY. Combining this upper bound with the
lower bound for X proves(1).

This incidence argument does not require a global choice of root gauges
or a simultaneous realization of the saved LP rows. It supplies a
necessary consistency condition precisely where the partial local model
had allowed independent face frequencies.

## Exact separation and scope

The unchanged witness has

```
a = 8550239/625000000
b = 23409967/5000000000
s = 71697493/5000000000
vY = 0  (vertex index5572)
```

The resulting gap in(1) is `3363761/1000000000 > 0`. This is a new
geometric obstruction to the saved witness of ratio12.128700376798244,
not a repaired ratio. No masses, old columns, or frozen files are changed.

The package certifies one complete forced-quad profile for one common
S type and two selected Q provenances. It does not enumerate every
baseline2 quadrilateral, every possible cap pair, or every global face
incidence. A future model must implement the signed face coefficients
and the vertex capacity together; no dimension conclusion follows here.

Internal review: construction_audit independently checked the conceptual
injection and counting argument, conditional on the literal root roles;
the source verifies those roles. This is distinct from independent
arithmetic replay or external review. **FINAL/FROZEN:** generation and a fresh
saved replay both passed. Direct input hashes, local report links, and source
syntax were checked. The mathematical injection proof and these exact
computations have separate scopes.

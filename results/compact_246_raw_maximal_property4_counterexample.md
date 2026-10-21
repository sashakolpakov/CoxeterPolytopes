# A raw maximal arithmetic counterexample to the current three-branch Property-4 compiler

**Status (2026-08-11): proved by a self-contained exact certificate.**

## 1. Bottom line

Let

\[
                       P=\Delta(2,4,6)\subset \mathbb H^2
\]

be the compact Coxeter triangle with angles
\(\pi/2,\pi/4,\pi/6\).  Then:

1. `P` is an arithmetic Coxeter chamber over \(\mathbb Q\);
2. its reflection group is maximal among **all** discrete hyperbolic
   reflection groups, not merely relative to a chosen lattice; and
3. `P` admits none of the three **exact-return Property-4 packages** obtained
   by combining Chamber Fidelity with the carrier alternatives in Theorem
   4.1 of
   [`canonical_chamber_rigidity_absw_transfer_theorem.md`](canonical_chamber_rigidity_absw_transfer_theorem.md):

   - direct rational full Néron--Severi;
   - certified weak-log rational; or
   - minimalizable canonical-orthogonal full trace.

Consequently the universal outer assertion

\[
 \text{every raw maximal arithmetic }P
 \longmapsto
 \text{one of those three exact carriers}
\]

is **false**, rather than merely unproved.

This conclusion has a precise escape clause.  It does not exclude a future
category which permits a proper invariant/eigen-nef slice even when the
minimal Vinberg field is \(\mathbb Q\), or which weakens literal full-nef
return.  Either change would alter Property 2 of the current contract and
would require a new AG rank theorem.  It therefore cannot be used to claim
the present \(n\leq25\) transfer for all raw maximal chambers.  Nor does the
argument exclude a **bare carrier-relative** canonical-orthogonal trace
injection which is not required to return `P`; Theorem 4.1's third row alone
has that weaker scope.

## 2. Integral arithmetic realization

Use negative wall normals \(r_1,r_2,r_3\) with Gram matrix

\[
 G=\begin{pmatrix}
 -12&6&0\\
 6&-6&3\\
 0&3&-2
 \end{pmatrix}.                                             \tag{2.1}
\]

Its determinant is \(36\), and exact LDL decomposition has pivots
\((-12,-3,1)\), so its signature is \((1,2)\).  The normalized squared
pairings are

\[
 {6^2\over12\cdot6}={1\over2},\qquad
 {3^2\over6\cdot2}={3\over4},\qquad 0,                       \tag{2.2}
\]

giving labels \(4,6,2\).  All three rank-two principal Grams are negative
definite: their determinants are \(36,3,24\).  Thus every vertex is finite
and the chamber is compact.  Gauss--Bonnet gives

\[
 \operatorname {Area}(P)
 =\pi\left(1-\frac12-\frac14-\frac16\right)
 =\frac\pi{12}.                                             \tag{2.3}
\]

In the wall basis, the reflections are

\[
 s_1=\begin{pmatrix}-1&1&0\\0&1&0\\0&0&1\end{pmatrix},\quad
 s_2=\begin{pmatrix}1&0&0\\2&-1&1\\0&0&1\end{pmatrix},\quad
 s_3=\begin{pmatrix}1&0&0\\0&1&0\\0&3&-1\end{pmatrix}.    \tag{2.4}
\]

They are integral involutions and satisfy \(s_i^tGs_i=G\).  Hence
\(W(P)\subset O(G,\mathbb Z)\).  The left side has compact fundamental
domain, hence finite covolume.  Since \(O(G,\mathbb Z)\) is an arithmetic
lattice, the inclusion has finite index.  Thus `P` is arithmetic, with
minimal Vinberg field \(\mathbb Q\).

The rational form is anisotropic by a one-line descent.  For \(q=-G/2\),
the associated integral quadratic polynomial is

\[
 q(x,y,z)=6x^2-6xy+3y^2-3yz+z^2.                            \tag{2.5}
\]

If \(q=0\), reduction modulo three first gives \(3\mid z\).  Put
\(z=3z_1\), divide by three, and reduce again.  The residual equation is

\[
                       2x^2+xy+y^2=0\pmod3.                 \tag{2.6}
\]

Its only zero in \(\mathbb F_3^2\) is \((0,0)\), so \(3\mid x,y\) as well.
Infinite descent excludes a nonzero rational isotropic vector.  This also
identifies the earliest rational-surface compiler failure: a full-rank
rational host would be \(I_{1,2}\), which is isotropic, and therefore is not
rationally similar to (2.1).

## 3. Absolute reflection maximality

Suppose \(W(P)\) were properly contained in a discrete reflection group
\(W(Q)\), of index \(N\geq2\).  After choosing compatible chambers, the
compact triangle `P` is tiled by \(N\) copies of a compact Coxeter polygon
`Q`, and

\[
                         \operatorname {Area}(Q)
                         ={\pi\over12N}\leq{\pi\over24}.     \tag{3.1}
\]

A compact Coxeter polygon with at least four sides has area at least
\(\pi/6\).  For four sides, positivity excludes four right angles and the
next-largest angle sum is
\(\pi/2+\pi/2+\pi/2+\pi/3=11\pi/6\); five or more sides give an even larger
bound.  Hence `Q` must be a triangle.

Write its ordered labels as \(2\leq a\leq b\leq c\).  The conditions

\[
 0<1-{1\over a}-{1\over b}-{1\over c}\leq{1\over24},
 \qquad
 {1/12\over1-1/a-1/b-1/c}\in\mathbb Z_{\geq2}               \tag{3.2}
\]

leave exactly one candidate:

\[
                         Q=\Delta(2,3,8),\qquad N=2.         \tag{3.3}
\]

The index-two case is elementary.  Two adjacent `Q` chambers are glued
across one side.  Their union can be a triangle only when the glued side
meets the right-angle vertex; that vertex then disappears.  Gluing the
\((2,8)\)-side gives \(\Delta(3,3,4)\).  Gluing the \((2,3)\)-side gives a
non-Coxeter angle \(2\pi/3\), and gluing the \((3,8)\)-side gives a
quadrilateral.  None gives \(\Delta(2,4,6)\).

Therefore no proper discrete reflection overgroup exists.  In particular,
`P` is a raw maximal arithmetic chamber in every maximality sense admitted
by the current transfer theorem.

## 4. The two rational branches are empty

For \(k=\mathbb Q\), Chamber Fidelity in the corrected contract requires

\[
                         P=\mathbf P(\operatorname {Nef}S),  \tag{4.1}
\]

not a proper eigenslice.  Equality with a two-dimensional hyperbolic chamber
forces \(\rho(S)=3\).

There is already a rational HOST obstruction.  Property 1 would give a
full-rank rational similarity embedding of the anisotropic form (2.1) in
\(\operatorname {NS}(S)_\mathbb Q\).  But a smooth rational surface of rank
three has the odd unimodular, isotropic numerical lattice \(I_{1,2}\).  Two
full-rank rational quadratic spaces related by a positive similarity have
the same isotropy, so such an embedding is impossible.

Every smooth projective rational surface of Picard rank at least two has a
nef isotropic class.  Run a rational-surface MMP and pull back a ruling
fibre; if the endpoint is \(\mathbb P^2\), retain the ruling after the first
blow-up.  Pullback preserves nefness and square zero.  Its ray is an ideal
point of the projectivized nef cone.

Thus every positive-dimensional rational full-nef chamber is noncompact.
The target `P` is compact, so both rational branches are empty:

\[
 \operatorname {Car}_{\rm direct}(P)
 =\operatorname {Car}_{\rm weak\mbox{-}log}(P)
 =\varnothing.                                             \tag{4.2}
\]

This uses neither the direct \(n\leq11\) theorem nor the weak-log
\(n\leq25\) inequality.  It rules out the carrier itself.

## 5. The exact-return minimalizable full-trace branch is empty

The rational trace block has rank \(3\).  In an exact Property-4 package,
literal full-nef equality again gives \(\rho(S)=3\).  If a supplied
finite-orbit contraction were nontrivial,
its minimal endpoint would have Picard rank at most two and could not contain
the descended rank-three trace block.  Hence any carrier in the third branch
is already minimal.

That branch places the full trace space in \(K_S^\perp\).  Since the trace
space is all of \(\operatorname {NS}(S)_\mathbb Q\), this forces

\[
                              K_S\equiv0.                    \tag{5.1}
\]

Let \(C\) be an irreducible negative curve defining a nef facet.  Adjunction
then gives

\[
                  C^2=2p_a(C)-2<0,
\]

so \(p_a(C)=0\) and \(C^2=-2\).  Thus every physical facet is a smooth
rational \((-2)\)-curve.  For two distinct such curves, \(C.D\) is a
nonnegative integer.  If their mirrors meet at a finite Coxeter angle, then

\[
                         \cos(\pi/m)={C.D\over2}<1.          \tag{5.2}
\]

Consequently \(C.D\in\{0,1\}\) and \(m\in\{2,3\}\).  The target has edges
of labels \(4\) and \(6\), an immediate contradiction.  Hence

\[
                 \operatorname {Car}^{\rm exact}_{\rm min\mbox{-}trace}(P)
                 =\varnothing.                              \tag{5.3}
\]

The superscript is essential: (5.3) does not assert emptiness of the bare
carrier-relative third row of the dimension-transfer theorem after its
exact chamber-return requirement is removed.

## 6. Consequence and exact no-go frontier

The obstruction is not another failure of effectivity for one chosen
decoration.  It starts with an undecorated raw chamber, proves absolute
reflection maximality, and empties every branch of the exact complex
Property-4 disjunction which, if exhaustive, would activate the conditional
\(n\leq25\) theorem.

There is also a broader common-class observation.  If all facet reflections
fix one canonical class, nonsingularity of (2.1) forces that class to vanish,
and Section 5 applies.  If they fix one log class \(K_S+B/2\) with effective
boundary, the same spanning argument gives \(B\equiv-2K_S\).  A ruled
surface is excluded by the nef isotropic class.  For nonnegative Kodaira
dimension, pseudo-effectivity of both \(K_S\) and \(-K_S\) forces
\(K_S\equiv0\); a blow-up is then impossible, and Section 5 again applies.
The finite-orbit descent alternative is exactly the rank argument already
used above.

Accordingly, the usual common-class realizations of the three
canonical/log/MMP alternatives do not evade this example.  This paragraph is
not promoted to a no-go theorem for an unstated broader log/MMP alternative
outside the frozen compiler.  The conservative theorem needed by the current
project is (4.2)--(5.3): **the presently advertised exact-return outer
three-branch compiler cannot be exhaustive.**

The correct next target is therefore not to prove this form of Property 4.
It is to alter the carrier category and prove a different exhaustive
disjunction.  The least invasive possible repair is a proper invariant-nef
return for rational-field inputs, but that loses the current rational MMP
rank bounds and needs a replacement AG ceiling.

## 7. Replay

```bash
python3 -m py_compile \
  scripts/certify_compact_246_property4_counterexample.py

PYTHONPATH=src:. python3 \
  scripts/certify_compact_246_property4_counterexample.py --check-archive

PYTHONPATH=src:. pytest -q \
  tests/test_compact_246_property4_counterexample.py
```

Frozen certificate:
[`compact_246_property4_counterexample.json`](../data/enumerations/compact_246_property4_counterexample.json).

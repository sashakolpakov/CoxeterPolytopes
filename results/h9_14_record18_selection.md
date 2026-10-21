# An audited abstract candidate nerve for H9/14 record 18

## Result and exact scope

**TOOT.** Second-layer survivor 18, with graph6 word

```text
M???????????o?K?_
```

has an explicit 110-facet **abstract audited nerve candidate**.  This is a
positive result at the incidence level.  It is not yet a PL-sphere theorem, a
convex-polytope construction, or a hyperbolic Coxeter example.

The missing edges are

\[
\{0,12\},\quad\{1,12\},\quad\{2,13\},\quad
\{3,13\},\quad\{12,13\}.
\]

Adjoin the eight missing five-faces

\[
\begin{aligned}
&02456,\quad 05678,\quad 03789,\quad 3789{10},\\
&2456{11},\quad 124{10}{11},\quad 139{10}{11},\quad 149{10}{11}.
\end{aligned}
\]

The complete blocker of these thirteen minimal nonfaces has exactly 110
members, all of order five.  Their complements generate a pure
eight-dimensional complex with

\[
f=(14,86,310,735,1204,1372,1052,495,110).
\]

## Blocker-duality and ridge reduction

Put

\[
L=\{0,1\},\quad R=\{2,3\},\quad I=\{4,\ldots,11\},
\quad P=\{12,13\},\quad W=L\sqcup R\sqcup I.
\]

The missing-edge graph has exactly three minimal covers:

\[
P,\qquad L\cup\{13\},\qquad R\cup\{12\}.
\]

Consequently every possible five-element facet complement has exactly one of
the following forms:

\[
Y_T=P\cup T,\quad |T|=3,\ T\subset W;
\]
\[
Z_S=L\cup\{13\}\cup S,\quad |S|=2,\ S\subset R\cup I;
\]
\[
X_S=R\cup\{12\}\cup S,\quad |S|=2,\ S\subset L\cup I.
\]

There are respectively 220, 45, and 45 possibilities.  The complete symbolic
ridge system has 735 nonempty six-set groups: 495 containing both special
vertices and 120 of each one-special kind.  Equivalently, for every
four-set \(Q\subset W\),

\[
\sum_{T\in {Q\choose3}}y_T
+\mathbf1_{L\subset Q}z_{Q\setminus L}
+\mathbf1_{R\subset Q}x_{Q\setminus R}\in\{0,2\},
\]

while on every triple in \(R\cup I\) (respectively \(L\cup I\)) the sum of
its three \(z\)-edges (respectively \(x\)-edges) is in \(\{0,2\}\).  The
certifier mechanically reconstructs all 735 equations, with total incidences

\[
(Y,Z,X)=(1980,405,405).
\]

In the selected completion the counts are

\[
(|Y|,|Z|,|X|)=(60,25,25).
\]

The two 25-edge graphs are the balanced complete bipartite cuts

\[
Z=K_{\{2,4,5,6,11\},\{3,7,8,9,10\}},
\]

\[
X=K_{\{0,5,6,7,8\},\{1,4,9,10,11\}}.
\]

This is the rigid algebro-combinatorial feature of the candidate: it realizes
the two balanced cuts forced by blocker involution and the Lanner-order bound,
rather than arising from an unconstrained facet search.

## Complete exact audits

The generated face table contains 5,379 faces including the empty face.  The
solver-independent audit verifies:

- the 110 facet complements are in the certified ridge-support core;
- every one of the 495 ridges has degree two;
- the dual facet graph is connected;
- every independent triple is a face;
- the only minimal nonfaces are the five edges and eight displayed five-sets;
- the face vector is exactly schedule \(a=0\);
- every strict local \(g\)-condition and applicable Macaulay inequality holds;
- every face link has spherical Euler characteristic; and
- the complex and every face link have the mod-two homology of the appropriate
  sphere.

The blocker and packed face-table hashes are respectively

```text
969da679ad87ca730e822ceb7c5a4eee7f14b461062cf71543293c50c9872afc
620c83d073c4a141ef2ec5cf9333326cdc8d7c7e94d5e48e647b19b030943eaf
```

A uniform rank-four Gale chirotope was also found.  Its 1,001 archived signs
pass all 45,045 Grassmann--Pluecker relations, and its positive five-circuits
are exactly the 110 facet complements.  This is a necessary oriented-matroid
gate only: the oriented matroid has not been proved realizable over
\(\mathbb R\).

## Vinberg finite-volume warning

The degree-two ridge test above is the closed-pseudomanifold incidence
condition for the proposed abstract compact nerve.  It is **not** Vinberg's
finite-volume criterion.  That criterion must be applied to an actual Coxeter
diagram/Gram matrix: elliptic subdiagrams must have the required elliptic or
parabolic extensions, and parabolic subdiagrams must have the required ranks
and maximality properties.  Record 18 currently has neither Coxeter labels nor
a Lorentzian Gram matrix, so this certificate makes no compactness or
finite-volume claim.

The compact branch has since been excluded for this selected nerve by an
independently replayed unrestricted Lanner/ellipticity and exact-inertia
gate; see `results/h9_14_record18_coxeter_label_gate.md`.  That later theorem
does not decide finite-volume noncompact realizability or any other
completion of the record-18 graph shadow.

For a positive finite-volume noncompact example, the remaining gates are:

1. prove PL sphericity or construct an ordinary convex realization;
2. find Coxeter labels compatible with finite and ideal vertex diagrams;
3. solve the exact Lorentzian Gram rank/signature equations; and
4. only then apply Vinberg's elliptic/parabolic finite-volume criterion.

## Reproduction

The frozen artifacts are

- `src/coxeter_polytopes/combinatorics/h9_14_record18.py`;
- `scripts/certify_h9_14_record18.py`;
- `scripts/independent_h9_14_record18_audit.py`;
- `data/enumerations/h9_14_record18_selection.json`; and
- `tests/test_h9_14_record18.py`.

The certificate digest is

```text
27193c530804237c8505b7ce7400e7551a6b152d5a4391f2dc624371d8f54406
```

Replay it with

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_record18.py
PYTHONPATH=.:src python3 scripts/independent_h9_14_record18_audit.py
PYTHONPATH=.:src pytest -q tests/test_h9_14_record18.py
```

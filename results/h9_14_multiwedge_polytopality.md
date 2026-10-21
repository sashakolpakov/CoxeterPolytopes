# Three multiwedge polytopality upgrades in H9/14

## Result

H9/14 records **43, 44, and 47** are not merely abstract
homology-sphere/oriented-matroid candidates. Each is the nerve of an ordinary
convex simple 9-polytope with fourteen facets.

Subsequent work proves that records 44 and 47 are each the unique residual
incidence orbit of their graph shadows under the compact missing-face
necessities and transfers their unrestricted label obstructions to the whole
shadows; see
[`h9_14_record44_whole_shadow_compact_no_go.md`](h9_14_record44_whole_shadow_compact_no_go.md)
and
[`h9_14_record47_whole_shadow_compact_no_go.md`](h9_14_record47_whole_shadow_compact_no_go.md).
This does not change the status of the construction as ordinary-polytopal,
nor does it exclude finite-volume noncompact chambers.

The common construction is

\[
K_R=\partial\Delta^4 * K_R^0(J_R),
\tag{1}
\]

where `*` is simplicial join and `K_R^0(J_R)` is a simplicial multiwedge of a
small explicitly polytopal sphere. Dually, the desired simple polytope is

\[
Q_R=\Delta^4\times W_R,
\tag{2}
\]

where `W_R` is obtained from the dual base polytope by repeated facet wedges.

With the subsequent record-37 and record-34 exclusions and record-35 rational
multiwedge included, the current theorem-facing partition of the 142 graph
shadows is

\[
121\text{ independently impossible}
\; +\;6\text{ explicit ordinary-polytopal}
\; +\;15\text{ theorem-undecided}.
\tag{3}
\]

The sixth realization is documented separately in
[`h9_14_record35_polytopality.md`](h9_14_record35_polytopality.md); it is not
one of the three Delzant constructions in this report.

Ordinary polytopality is not Coxeter realizability: no edge labels or
Lorentzian Gram matrices are supplied.  A subsequent
[exact metric-gate audit](h9_14_multiwedge_metric_gate.md) obstructs several
natural right-angle ansatzes for the displayed Delzant/projective
realizations only; it does **not** exclude a different realization of any of
the three combinatorial types.

## The wedge lemma

Let a simple polytope be written using facet slacks

\[
P=\{x:\ell_i(x)\geq0,\ 1\leq i\leq m\}.
\]

The facet wedge at `i` can be realized as

\[
\operatorname{Wedge}_i(P)
=\{(x,t):\ell_j(x)\geq0\ (j\ne i),\quad
0\leq t\leq\ell_i(x)\}.
\tag{4}
\]

It is a bounded simple polytope of dimension one higher with one additional
facet. In the dual boundary complex, the vertex `i` is replaced by two clone
vertices: every minimal nonface not containing `i` is unchanged, while every
minimal nonface containing `i` contains both clones. Iteration replaces `i`
by a clone block of any prescribed positive size. This is exactly the
simplicial multiwedge or `J`-construction.

This preservation of polytopality and the duality between simplicial and
facet wedges are stated in [Choi--Park, *Wedge Operations and Torus
Symmetries*, Proposition 2.2(3) and Lemma
2.3](https://www.jstage.jst.go.jp/article/tmj/68/1/68_1458248864/_pdf).

Finally, joins of polytopal boundary complexes are polytopal; dually, they
correspond to products of simple polytopes. Thus (1) and (2) are explicit
polytopality constructions.

## Record 44: a multiwedged hexagon

Put

\[
A=\{2,3,4,5,6\},\qquad B=\{7,8,9,10\}.
\]

The set `A` occurs as one minimal nonface and nowhere else, giving the join
factor \(\partial\Delta^4\) in (1). Compress `B` to one base vertex `b`. The
remaining base is the boundary cycle of the hexagon

\[
b-0-1-11-12-13-b.
\tag{5}
\]

Its six nonedges not containing `b` are exactly

\[
(0,11),(0,12),(0,13),(1,12),(1,13),(11,13),
\]

the record-44 missing edges. Its three nonedges containing `b` are

\[
(b,1),(b,11),(b,12).
\]

Replacing `b` by the four-clone block `B` turns these into the three A5s

\[
B\cup\{1\},\quad B\cup\{11\},\quad B\cup\{12\}.
\]

Together with `A`, these are exactly the four recorded A5s. Hence

\[
K_{44}=\partial\Delta(A)*C_6(1,1,1,1,1,4).
\]

The hexagon has dimension two; three facet wedges give a simple 5-polytope
with nine facets, and its product with \(\Delta^4\) is the required simple
9-polytope with fourteen facets.

## Records 43 and 47: multiwedged stacked spheres

For both records put

\[
A=\{3,4,5,6,7\},\qquad B=\{8,9,10\},
\]

and compress `B` to `b`. Again `A` is the separate \(\partial\Delta^4\) factor.
The remaining bases are seven-vertex triangulated 2-spheres, explicitly
constructed by stacking tetrahedral facets.

For record 43, start with the tetrahedron on

\[
\{1,2,11,b\}
\]

and stack, in order,

\[
12\text{ on }\{2,11,b\},\qquad
0\text{ on }\{1,2,b\},\qquad
13\text{ on }\{11,12,b\}.
\tag{6}
\]

The resulting sphere has precisely the six record-43 missing edges and three
additional minimal triples

\[
\{1,2,b\},\quad\{2,11,b\},\quad\{11,12,b\}.
\]

Replacing `b` by `B` produces exactly the three non-join A5s of record 43.

For record 47, start with the same tetrahedron and stack

\[
0\text{ on }\{1,2,b\},\qquad
12\text{ on }\{2,11,b\},\qquad
13\text{ on }\{1,11,b\}.
\tag{7}
\]

Its extra minimal triples are

\[
\{1,2,b\},\quad\{1,11,b\},\quad\{2,11,b\},
\]

which expand to exactly the three non-join A5s of record 47.

Each base is the boundary of an explicit stacked simplicial 3-polytope. Two
facet wedges at the dual facet `b` give a simple 5-polytope with nine facets;
product with \(\Delta^4\) again gives dimension nine and fourteen facets.

## Algebraic-geometric realization

The construction can be made **Delzant**, so each record also comes from a
smooth projective toric 9-fold.

- For record 44 choose an integral Delzant hexagon, for example the moment
  polygon of the blow-up of \(\mathbb P^2\) at its three torus fixed points.
- For records 43 and 47, dualize (6) and (7) to toric corner cuts of a large
  integral Delzant tetrahedron. At a corner with primitive inward normals
  \(a_1,\ldots,a_d\), the new primitive normal is
  \(a_1+\cdots+a_d\). The last record-43 cut is at the newly created corner
  \(\{b,11,12\}\), so it is an iterated torus-fixed-point blow-up.
- In (4), use a *primitively normalized* integral slack \(\ell_i\). The new
  inward normals are

  \[
  (a_j,0),\qquad e_{d+1},\qquad(a_i,-1).
  \]

  At a vertex away from the wedged facet, the old lattice basis together with
  either new normal is unimodular. At a vertex on the facet, adding
  \(e_{d+1}\) to \((a_i,-1)\) recovers \((a_i,0)\), again giving the old basis
  plus \(e_{d+1}\). Thus the wedge is Delzant. The coefficient-one primitive
  normalization is essential: replacing \(\ell_i\) by \(k\ell_i\) can create
  lattice index \(k\).
- Products with a simplex are Delzant.

Let \(X_R\) be the resulting smooth projective toric variety, and let

\[
D_R=\sum_{i=0}^{13}D_i
\]

be its reduced toric boundary. It is a simple-normal-crossings divisor and

\[
K_{X_R}+D_R=0.
\tag{8}
\]

The fourteen components correspond to the fourteen facets of \(Q_R\), and

\[
D_{i_1}\cap\cdots\cap D_{i_k}\ne\varnothing
\quad\Longleftrightarrow\quad
\{i_1,\ldots,i_k\}\in K_R.
\tag{9}
\]

Thus the dual intersection complex \(\operatorname{Dual}(D_R)\) is exactly
the record nerve. In particular, \((X_R,D_R)\) is a smooth projective toric
log-Calabi--Yau pair carrying the desired combinatorics.

Equivalently, its Chow ring has the presentation

\[
A^*(X_R)=
\mathbb Z[D_0,\ldots,D_{13}]/(I_{\rm SR}+I_{\rm lin}),
\tag{10}
\]

where the Stanley--Reisner ideal `I_SR` is generated by the six quadratic and
four quintic monomials listed in the record certificate. The Picard rank is

\[
14-9=5.
\]

This gives a precise algebro-geometric answer to the divisor question: the
distinguished toric boundary pair, or equivalently its labelled
Stanley--Reisner presentation, recovers the candidate nerve automatically.
The abstract variety or abstract graded Chow ring alone need not canonically
recover the boundary labels. The bare Neron--Severi group has rank only five,
so it cannot itself be the rank-ten Lorentzian Coxeter lattice. Neither the
boundary complex nor its intersection algebra supplies Coxeter angle data;
that remains an additional Gram-matrix or reflection-lattice constraint.

There is nevertheless a clean algebro-geometric meaning for the *ordinary*
compactness and invariant measure. If \(H\) is the ample toric divisor whose
moment polytope is \(Q_R\), projectivity makes \(Q_R\) compact and

\[
H^9=9!\,\operatorname{Vol}_{\rm lattice}(Q_R).
\tag{11}
\]

Likewise, \(H^8D_i\) is the normalized volume of the corresponding facet,
up to the standard factor \(8!\). Symplectic volume pushes forward under the
moment map to the Duistermaat--Heckman/Lebesgue measure on \(Q_R\). Thus the
ordinary polytope volume is encoded by classical intersection numbers of the
polarized pair \((X_R,H)\). This should not be confused with the hyperbolic
volume of a Coxeter chamber, which still depends on the missing Lorentzian
metric and angle labels.

Numerically, the divisor relation space also explains the Gale rank. There
are fourteen invariant divisors and nine character relations, giving Picard
rank five; imposing the affine support/ample class leaves a four-dimensional
relation space, the expected Gale rank \(14-9-1=4\). This is a useful bridge
between the toric divisor presentation and the realizable Gale transform,
but it is not an identification with the Coxeter Gram space.

## Certificate

The machine audit constructs every base facet, derives every base minimal
nonface, performs the clone substitution, adds the simplex-boundary factor,
and checks exact equality with each archived ten-element missing-face clutter.
It also checks all dimensions and facet counts. A second implementation uses
explicit integer coordinates for all three small bases, reconstructs their
supporting facets by exact determinant signs, and independently recovers the
common residual face polynomial

\[
1+9t+30t^2+50t^3+45t^4+18t^5.
\]

Multiplication by the face polynomial of \(\partial\Delta^4\),
\((1+t)^5-t^5\), gives the full face vector in every case.

The frozen certificate is
`data/enumerations/h9_14_multiwedge_polytopality.json`, with digest

```text
d5efd409cf5d5a5a62e9784f8c9e71eebb5a08ee20e8ea7f1d05e894d305a183
```

Reproduce it with

```sh
PYTHONPATH=src:. python3 scripts/certify_h9_14_multiwedge.py
PYTHONPATH=src:. python3 scripts/independent_h9_14_multiwedge_audit.py
PYTHONPATH=src:. python3 -m pytest -q tests/test_h9_14_multiwedge.py
```

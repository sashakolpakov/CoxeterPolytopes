# Golden Lanner blocks and the global 9/14 Gram constraint

## Outcome

For the norm-two \(\mathbb Z[\phi]\) K3 construction, the statement

> every nine-facet vertex becomes ADE after replacing both 3- and 5-edges
> by simple curve edges

is necessary but gives no new vertex-level pruning.  It is already a
consequence of ellipticity and the classification of finite Coxeter diagrams.

There are, however, two exact global reductions.

1. Every hypothetical compact \(\mathbb H^9\) polytope with 14 facets on this
   golden track contains one of **eight** explicitly classified Lanner
   missing-face blocks.
2. At every finite vertex, the five complementary facets have a rank-one
   Schur complement.  This is the algebraic condition that fuses all local
   hyperbolic missing-face directions into the single negative direction of
   a rank-ten Lorentzian Gram matrix.

The local conditions alone do not force a contradiction.  An exact
14-facet countercertificate based on

\[
\Delta_3\times\Delta_3\times\Delta_1^3
\]

passes every local elliptic, trace-ADE, pairwise curve-intersection-sign,
arithmetic-conjugate, and missing-face test, but its global field Gram inertia
is \((9,5,0)\), not
\((9,1,4)\).  Thus the rank-one fusion equations, rather than the vertex-ADE
condition, are the useful next constraint.

## 1. Why vertex-ADE is redundant

Let \(I\) be the nine facets through a finite vertex.  The Coxeter diagram on
\(I\) is elliptic.  With finite labels restricted to \(2,3,5\), every
irreducible component is one of

\[
A_n,D_n,E_6,E_7,E_8,H_3,H_4,I_2(5),A_1.
\]

After forgetting whether a non-right edge was labelled 3 or 5, the last
three exceptional types become respectively \(A_3,A_4,A_2\).  Hence the
undecorated graph is automatically a disjoint union of ADE diagrams.

Equivalently, for the K3 curve matrix

\[
Q_I=-2I+\operatorname{Adj}(\Gamma_I),
\]

ellipticity of the golden Coxeter matrix already implies that \(Q_I\) is
negative definite.  The K3 statement supplies the geometric interpretation
as a contractible rational-curve configuration, but it does not remove any
elliptic \(\{2,3,5\}\)-labelled vertex diagram.

The same observation explains the first global graph consequence.  Burcroff's
absence of missing triples implies that a triangle in the global 3/5 graph
would extend to a vertex, where it is not ADE.  Thus the global 3/5 graph is
triangle-free.  But this is again already forced by ellipticity: the all-3
triangle is affine and inserting a 5 makes it indefinite.

## 2. A missing face of order four or five is forced

Burcroff's Lemma 10.10 says that a hypothetical compact 9-polytope with 14
facets has no missing face of order three and no finite label \(m\geq6\).
If all its missing faces had order two, it would be 3-free in Esselmann's
sense.  Esselmann's facet bound would then give at least \(2d=18\) facets
(with equality only for a cube), contradicting 14 facets.

Missing faces of a compact Coxeter polytope are Lanner diagrams, and a Lanner
diagram has order at most five.  Therefore a 9/14 candidate has a missing
face of order four or five.  Equal norm two excludes label four over
\(\mathbb Q(\sqrt5)\), so all of its ordinary edges have labels 3 or 5.

## 3. Exact eight-block theorem

Let \(M\) be a missing face of order four or five.  Every proper subset of
\(M\) is a face and extends to a finite vertex.  Consequently every proper
induced subgraph of the undecorated 3/5 graph is ADE.

For four nodes, a connected triangle-free graph with this property is
\(A_4\), \(D_4\), or the four-cycle \(C_4=\widetilde A_3\).  For five nodes,
the only additional minimally non-ADE possibilities are the five-cycle
\(C_5=\widetilde A_4\) and the four-leaf star
\(K_{1,4}=\widetilde D_4\).

The decorations are then forced by rank-three and rank-four ellipticity:

- on a path, adjacent 5-edges are impossible;
- on \(D_4\), at most one arm can carry a 5;
- a \(C_4\) must have one 5 or two opposite 5s, since all-3 is affine;
- a \(K_{1,4}\) with a 5 contains a proper \(D_4(5,3,3)\) Lanner
  subdiagram, while all-3 is affine;
- a \(C_5\) with a 5 has a four-node deletion in which that 5 is internal,
  producing the Lanner path \(3-5-3\), while all-3 is affine.

Thus \(C_5\) and \(K_{1,4}\) cannot be order-five Lanner diagrams.  Exactly
the following eight diagrams remain.  Here the determinant and Smith group
refer to the full restriction-of-scalars lattice on
\(e_i,\phi e_i\), with scalar \(-\phi/\sqrt5\).

| order | labelled Lanner graph | curve graph | field determinant | trace determinant | trace discriminant group |
|---:|---|---|---:|---:|---|
| 4 | \(A_4(3,5,3)\) | \(A_4\) | \(5-4\phi\) | \(-11\) | \(\mathbb Z/11\) |
| 4 | \(A_4(5,3,5)\) | \(A_4\) | \(6-5\phi\) | \(-19\) | \(\mathbb Z/19\) |
| 4 | \(D_4(5,3,3)\) | \(D_4\) | \(4-4\phi\) | \(-16\) | \((\mathbb Z/2)^4\) |
| 4 | \(C_4\), one 5 | \(\widetilde A_3\) | \(2-5\phi\) | \(-31\) | \(\mathbb Z/31\) |
| 4 | \(C_4\), opposite 5s | \(\widetilde A_3\) | \(1-7\phi\) | \(-55\) | \(\mathbb Z/55\) |
| 5 | \(A_5(5,3,3,3)\) | \(A_5\) | \(6-4\phi\) | \(-4\) | \((\mathbb Z/2)^2\) |
| 5 | \(A_5(5,3,3,5)\) | \(A_5\) | \(8-6\phi\) | \(-20\) | \(\mathbb Z/2\oplus\mathbb Z/10\) |
| 5 | \(D_5\), terminal 5 | \(D_5\) | \(4-4\phi\) | \(-16\) | \((\mathbb Z/4)^2\) |

Every identity-place matrix in the table has inertia \((r-1,1)\), and every
conjugate matrix is positive definite.  The full trace lattice has inertia
\((1,2r-1)\).  The table is produced by an exact enumeration of every graph
on four or five nodes and every 3/5 decoration; no floating-point spectral
comparison enters the certificate.

This is a constructive result as much as an obstruction.  In particular,
the two path blocks of trace discriminants 11 and 19 line up with the current
discriminant-11 and discriminant-19 K3 candidates.  A missing face can be used
as the initial Lorentzian divisor sublattice from which the rest of the
surface is glued.

### The affine four-cycle branch

For either four-cycle, the sum

\[
F=\delta_1+\delta_2+\delta_3+\delta_4
\]

is the primitive null vector of the curve matrix.  If the four classes are
effective irreducible curves, \(F\) is a nef isotropic fiber class.  The
second pairing nevertheless gives

\[
B(F,F)=2\quad\text{or}\quad4
\]

for one or two opposite 5-edges.  Hence \(F\) is not fixed by the golden
endomorphism: \((\tau F)^2=F\cdot\tau F=B(F,F)>0\).  The affine curve fiber
therefore does not by itself put an ideal point in the golden eigenslice.
It is a geometric feature, not a compactness contradiction.

## 4. Rank-one vertex-complement lemma

Let \(h\) be a 14-wall norm-two Gram matrix of rank 10 with identity inertia
\((9,1,4)\) and conjugate inertia \((10,0,4)\).  Let \(I\) be the nine facets
through a finite vertex and let \(J\) be the five remaining facets.  Write

\[
h=\begin{pmatrix}H&C\\ C^t&D\end{pmatrix},
\qquad H=h_I,
\]

and form the Schur complement

\[
S_I=D-C^tH^{-1}C.
\]

Since \(H\) is positive definite at both real places, block Gaussian
congruence gives

\[
h\sim H\perp S_I.
\]

It follows immediately that

\[
\operatorname{rank}S_I=1,
\]

with inertias

\[
\operatorname{inertia}(S_I)=(0,1,4),\qquad
\operatorname{inertia}(S_I^\sigma)=(1,0,4).
\]

In particular every two-by-two minor of \(S_I\) vanishes.  At the identity
place it has the form \(-uu^t\).  The inward-normal orientation supplies an
additional, unsquared datum: the five entries of \(u\) have one common sign,
because all five nonincident facets are strict on the same side of the
vertex.  Hence every off-diagonal entry of \(S_I\) is strictly negative.
Written entrywise, every pair of outside facets \(j,k\in J\) satisfies

\[
\left(h_{jk}-h_{jI}H^{-1}h_{Ik}\right)^2
=
\left(2-h_{jI}H^{-1}h_{Ij}\right)
\left(2-h_{kI}H^{-1}h_{Ik}\right).
\]

The squared equations alone are not the complete realization gate: one must
also impose

\[
h_{jk}-h_{jI}H^{-1}h_{Ik}<0
\]

at the identity place.  If \(j_0\) is used as the outside anchor, then the
coordinate of \(r_k\) along \(r_{j_0}\) is the ratio of this off-diagonal
Schur entry to the negative diagonal Schur entry for \(j_0\).  Thus the
inequality is exactly the positive vertex-anchor same-side test.

These are exact coupled equations among all five outside walls.  They are
strictly stronger than checking each ten-wall extension \(I\cup\{j\}\)
separately.

### Spectral formulation and the limit of interlacing

At the identity place write \(h=2I-A\), where \(A\) is the nonnegative
weighted adjacency matrix: a 3-edge has weight 1, a 5-edge weight \(\phi\),
and a dotted edge weight greater than 2.  The target inertia says that \(A\)
has

\[
1\text{ eigenvalue above }2,qquad
4\text{ eigenvalues equal to }2,qquad
9\text{ eigenvalues below }2.
\]

For a vertex, \(\rho(A_I)<2\).  Cauchy interlacing is consistent with this:
deleting the five complementary nodes can remove the Perron eigenvalue and
the four eigenvalues equal to 2.  Thus interlacing alone gives no
contradiction.  The fact that the deletion size is exactly five explains why
the rank-one Schur condition is the sharp equality-case datum to retain.

The known compact Bugaenko chamber in \(\mathbb H^8\) provides an exact
calibration: its 11-wall Gram matrix has inertia \((8,1,2)\), and at each
eight-wall vertex the three-by-three complement Schur matrix has rank one and
inertias \((0,1,2)\), \((1,0,2)\) at the two places.

## 5. Exact 14-node countercertificate to local sufficiency

Consider the simple combinatorial 9-polytope

\[
P=\Delta_3\times\Delta_3\times\Delta_1\times\Delta_1\times\Delta_1.
\]

It has 14 facets.  Partition them into factor blocks of sizes

\[
4,4,2,2,2.
\]

Its minimal missing faces are exactly those five blocks, so their orders are
\(4,4,2,2,2\); in particular there is no missing triple.  A vertex omits one
facet from each four-block and chooses one facet from each two-block, giving
\(4\cdot4\cdot2^3=128\) simple nine-facet vertices.

Place the Lanner path \(A_4(3,5,3)\) on each four-block, put dotted product

\[
h_{ij}=-\phi^2=-1-\phi

\]

on each two-block, and make different blocks mutually orthogonal.  The
conjugate of a dotted two-block is positive definite, while its identity
matrix has inertia \((1,1)\).  Its curve intersection is

\[
\operatorname{Tr}\left(-\frac{\phi}{\sqrt5}(-\phi^2)\right)=2.
\]

The exact audit verifies:

- both four-blocks and all three dotted pairs are Lanner;
- every one of the 128 vertex diagrams is elliptic at both real places;
- every vertex curve graph is ADE;
- every off-diagonal formal curve intersection is a nonnegative integer; and
- the conjugate 14-wall Gram matrix is positive definite.

Nevertheless the block diagonal identity Gram matrix has

\[
\operatorname{inertia}(h)=(9,5,0),
\]

and the rank-28 restriction-of-scalars lattice has

\[
\operatorname{inertia}(q)=(5,23,0).
\]

At every vertex, the five-by-five Schur complement is negative definite of
rank five, rather than rank one.  The five independent negative directions
come from the five Lanner factor blocks.  A genuine \(\mathbb H^9\) design
must introduce cross-block couplings that fuse those five directions into
one while simultaneously creating four field-linear relations.

This countercertificate proves that no argument using only no-missing-three,
vertex-ADE, local Lanner blocks, conjugate positivity, and the pairwise
effectivity-sign test can establish nonexistence.  It is not an actual K3
curve configuration: its failed global signature is precisely what prevents
that interpretation.

## 6. Two bad neighbours are also locally compatible

Burcroff's Theorem 10.11 forces an \(H_4=[5,3,3]\) inside every order-five
missing face to have at least two bad neighbours (the \(F_4\) alternative is
unavailable on the norm-two golden track).  This does not contradict the K3
trace condition.

Attach two new 3-edges to the terminal node of \([5,3,3]\).  Each new node
separately completes the \(H_4\) to the order-five Lanner path
\([5,3,3,3]\).  On all six nodes the undecorated curve graph is \(D_6\), so
its curve matrix is negative definite.  Exact leading minors of the
identity-place field matrix are

\[
2, 3-\phi, 4-2\phi, 5-3\phi, 6-4\phi, 4-4\phi.
\]

Thus its inertia is \((5,1)\); the conjugate matrix is positive definite and
has determinant \(4\phi\).  Its rank-12 trace lattice has determinant
\(-16\) and discriminant group \((\mathbb Z/2)^4\).  Hence even the
published two-bad-neighbour condition and a full local K3 trace lift coexist
without a contradiction.

## 7. Consequence for the inverse search

The useful order is now:

1. enumerate a 9/14 nerve and its order-four/five missing faces;
2. assign one of the eight exact golden Lanner blocks to each such face;
3. enforce ellipticity on all nine-facet vertices;
4. for every vertex, impose the rank-one five-facet Schur equations;
5. only then construct the doubled integral \((Q,B)\) presentation and test
   the rank-20 K3 discriminant form.

The vertex-ADE check remains valuable as a semantic and implementation
cross-check, but it should not be counted as an independent pruning theorem.
The discriminants in the eight-block table are better design data: they turn
the missing faces into small Lorentzian divisor seeds before the rest of the
surface is constructed.

## Reproduction

The exact census, audits, and countercertificate are implemented in
`coxeter_polytopes.geometric.golden_lanner` and tested by
`tests/test_golden_lanner.py`.

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_golden_lanner.py
```

The combinatorial inputs are Burcroff, Lemma 10.10 and Theorem 10.11, and
Esselmann's 3-free facet bound as quoted in Burcroff, Lemma 10.6.  The eight
Lanner classification and all lattice calculations above are independently
reproduced by the exact implementation.

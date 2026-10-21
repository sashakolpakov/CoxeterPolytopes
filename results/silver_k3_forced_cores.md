# Silver K3 forced cores for the H9/14 inverse problem

## Outcome

The divisor-lattice construction has a second exact quadratic branch.  Let

\[
F=\mathbb Q(\sqrt2),\qquad \mathcal O_F=\mathbb Z[\sqrt2],
\]

and let a rank-ten norm-two \(\mathcal O_F\)-lattice have signatures
\((9,1)\) and \((10,0)\) at the two real embeddings.  With

\[
\alpha=-\frac{2+\sqrt2}{4},
\qquad
N=\left(M,\operatorname{Tr}_{F/\mathbb Q}(\alpha h)\right),
\]

the trace lattice has rank twenty and signature \((1,19)\), exactly the
signature and maximal rank of a complex K3 Néron--Severi lattice.

For a hypothetical compact \(\mathbb H^9\) Coxeter polytope with fourteen
facets on this norm-two **silver** branch, every missing face of order four or
five is forced into one of only two arithmetic Lannér types:

\[
\boxed{C_4(4,3,3,3)},\qquad
\boxed{C_5(4,3,3,3,3)}.
\]

Restriction of scalars identifies their complete trace lattices by explicit
unimodular congruences:

\[
C_4(4,3,3,3)\longmapsto U\oplus A_6(-1),
\]

\[
C_5(4,3,3,3,3)\longmapsto E_8(-1)\oplus U(2).
\]

Moreover, in either order the visible curve cycle has null vector
\(F_0=C_1+\cdots+C_r\), and

\[
F_0^2=(\sqrt2F_0)^2=0,
\qquad
F_0\cdot\sqrt2F_0=2.
\]

Thus the missing face itself supplies a primitive \(U(2)\) pair inside its
free trace root lattice.  If the curve classes are effective irreducible
curves on a K3 surface, \(F_0\) is the cycle class of an \(I_4\) or \(I_5\)
fiber configuration.  This is a rigid algebro-geometric input for the next
surface design, not a compact-polytope existence claim.

## 1. The integral trace dictionary

Write a field product as \(h_{ij}=a_{ij}+c_{ij}\sqrt2\), and define the two
integral divisor pairings

\[
Q_{ij}=\operatorname{Tr}(\alpha h_{ij}),
\qquad
B_{ij}=\operatorname{Tr}(\alpha\sqrt2 h_{ij}).
\]

Direct calculation gives

\[
Q_{ij}=-(a_{ij}+c_{ij}),\qquad
B_{ij}=-(a_{ij}+2c_{ij}),
\]

and hence

\[
h=(B-2Q)+(Q-B)\sqrt2.                                      \tag{1}
\]

The scalar \(\alpha\) is totally negative, has trace \(-1\), and lies in the
codifferent because

\[
2\sqrt2\,\alpha=-1-\sqrt2\in\mathcal O_F.
\]

Consequently a norm-two field root has divisor square \(-2\), and the finite
labels are

\[
\begin{array}{c|ccc}
m&2&3&4\\ \hline
(Q_{ij},B_{ij})&(0,0)&(1,1)&(1,2).
\end{array}                                                  \tag{2}
\]

A label five is unavailable because its normalized product
\(-2\cos(\pi/5)=-\phi\) does not belong to \(\mathbb Q(\sqrt2)\).  This
exclusion is specific to the norm-two silver factory; it is not a restriction
on arbitrary Coxeter polytopes.

## 2. The finite cores are exactly the golden trace cores

The standard missing-face reduction says that an order-four Lannér diagram
contains \(H_3\) or \(B_3\), while an order-five Lannér diagram contains
\(H_4\) or \(F_4\).  Since the silver field excludes label five, it forces
\(B_3\) and \(F_4\), respectively.

The exact trace lifts are

| silver finite type | field determinant | trace rank | trace determinant | integral trace type |
|---|---:|---:|---:|---|
| \(I_2(4)\) | \(2\) | 4 | 4 | \(D_4(-1)\) |
| \(B_3\) | \(2\) | 6 | 4 | \(D_6(-1)\) |
| \(F_4\) | \(1\) | 8 | 1 | \(E_8(-1)\) |

The implementation stores a determinant-\(\pm1\) integral change of basis
for every row.  In particular, an order-five missing face forces a primitive
unimodular \(E_8(-1)\) divisor summand, which splits orthogonally in any
ambient integral Néron--Severi lattice.

The striking coincidence is

\[
B_3\ \text{(silver)}\quad\text{and}\quad H_3\ \text{(golden)}
\quad\longmapsto\quad D_6(-1),
\]

\[
F_4\ \text{(silver)}\quad\text{and}\quad H_4\ \text{(golden)}
\quad\longmapsto\quad E_8(-1).
\]

Thus the two real-quadratic branches have different Coxeter decorations but
the same forced crystallographic K3 cores.

## 3. Exhaustive silver Lannér census

The census considers every graph on four and five labelled vertices and
every independent decoration of each active edge by 3 or 4.  A branch is
retained exactly when

1. the identity-place Gram matrix has inertia \((r-1,1,0)\);
2. every one-node deletion is positive definite; and
3. the conjugate-place Gram matrix is positive definite.

Isomorphic decorations are identified by the least complete pair-label code.
Exactly one type remains in each order:

| order | diagram | field determinant | trace determinant | trace lattice |
|---:|---|---:|---:|---|
| 4 | one-4 cycle \(C_4(4,3,3,3)\) | \(-1-2\sqrt2\) | \(-7\) | \(U\oplus A_6(-1)\) |
| 5 | one-4 cycle \(C_5(4,3,3,3,3)\) | \(-2-2\sqrt2\) | \(-4\) | \(E_8(-1)\oplus U(2)\) |

The literal enumeration counts are

| order | graph masks | ADE-prefilter survivors | decorated candidates | exact labelled survivors | isomorphism types |
|---:|---:|---:|---:|---:|---:|
| 4 | 64 | 19 | 176 | 12 | 1 |
| 5 | 1,024 | 137 | 2,384 | 60 | 1 |

The integral trace identifications again come with explicit unimodular
congruence matrices.  No floating-point eigenvalue comparison or lattice
identification from determinant alone enters the certificate.

## 4. The forced genus-one pair

By (2), both a 3-edge and a 4-edge lift, in the \(Q\)-pairing, to two
\((-2)\)-classes meeting once.  Therefore the visible curve Gram of either
one-4 cycle is the affine matrix \(-\widetilde A_{r-1}\).  For

\[
F_0=C_1+\cdots+C_r
\]

one has \(F_0^2=0\).  There is exactly one 4-edge, so the \(B\)-pairing adds
one extra intersection relative to \(Q\), giving

\[
B(F_0,F_0)=F_0\cdot\sqrt2F_0=2.
\]

Self-adjointness and \((\sqrt2)^2=2\) give

\[
(\sqrt2F_0)^2=F_0\cdot2F_0=0.
\]

In the interleaved basis \(C_i,\sqrt2C_i\), the two coordinate columns of
\(F_0,\sqrt2F_0\) have a unit minor.  Their span is therefore primitive in
the free trace lattice generated by the missing-face roots and has Gram

\[
\begin{pmatrix}0&2\\2&0\end{pmatrix}=U(2).                    \tag{3}
\]

Primitivity in a later saturated ambient K3 overlattice must still be checked
after gluing; equation (3) itself survives unchanged.

## 5. How this routes the inverse search

The exact H9/14 label gate should now divide candidates before any surface is
chosen.

1. A norm-two diagram with 5-edges but no 4-edges goes to the golden
   \(\mathbb Z[\phi]\) factory.
2. A norm-two diagram with 4-edges but no 5-edges goes to the silver
   \(\mathbb Z[\sqrt2]\) factory and must contain the cycle and \(U(2)\)
   configuration above.
3. A diagram using both 4- and 5-edges belongs to neither quadratic factory.
   A uniform rank-ten restriction-of-scalars construction over
   \(\mathbb Q(\sqrt2,\sqrt5)\) would have integral rank forty, beyond the
   K3 Picard bound.  This is a no-go for that uniform K3 factory, not for an
   arbitrary nonuniform divisor realization or an arbitrary Coxeter Gram
   matrix.

The most constrained new surface target is therefore a singular K3 carrying
a self-adjoint square-root-of-two action, an \(I_4\) or \(I_5\) fiber cycle,
and the forced \(D_6(-1)\) or split \(E_8(-1)\) core.  The remaining task is
to enumerate the negative-definite trace complement and its discriminant-form
gluing, then run wall closure and Vinberg compactness.  The Hodge index theorem
supplies hyperbolicity; the cycle, trace core, and global Gram equations supply
the nontrivial design restrictions.

## Reproduction and scope

The exact arithmetic and certificates are in

- `coxeter_polytopes.exact.silver`,
- `coxeter_polytopes.geometric.silver_finite`, and
- `tests/test_silver_finite_trace_lifts.py`.

Run

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_silver_finite_trace_lifts.py
```

This is a new carrier-level theorem and a finite seed classification.  It
does not assert that the trace lattices embed primitively in the K3 lattice,
that all displayed roots are effective walls of one chamber, that the wall
list is complete, or that the resulting eigen-nef polytope is compact.

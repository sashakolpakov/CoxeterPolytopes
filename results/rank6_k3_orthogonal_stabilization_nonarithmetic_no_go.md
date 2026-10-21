# The first rank-six K3 stabilization is compact but nonarithmetic

**Status (2026-08-11): proved \(n=2\) no-go and all-orthogonal-stabilization
corollary.**  Orthogonally
duplicating the negative field line in the rank-four construction produces a
primitive rank-six K3 trace carrier whose distinguished eigen-nef component
is compact.  However, the mixed K3 root required to close the new direction
has a nonintegral cyclic Cartan product.  Its reflection cannot preserve any
order-lattice or commensurable overlattice.  Thus the most direct
stabilization from \(\mathbf H^1\) to \(\mathbf H^2\) fails the arithmetic
reflection-saturation gate.

The same two-coordinate obstruction occurs in every higher orthogonal
stabilization.  This does not rule out nonorthogonal couplings.  It
identifies the exact invariant such a coupling must repair.

## 1. The orthogonal duplicate

For every integer \(a\ge2\), let

\[
 d=a-1,\qquad
 t=a+\sqrt{a^2-1},\qquad
 A=\mathbf Z[t],
\tag{1.1}
\]

and retain

\[
 B_0=\begin{pmatrix}-1&0\\0&1\end{pmatrix},
\qquad
 B_a=\begin{pmatrix}1-2a&-1\\-1&-1\end{pmatrix},
\qquad
 C_a=\begin{pmatrix}0&-1\\1&2a\end{pmatrix}.
\tag{1.2}
\]

The first stabilization is

\[
 L_a^{(2)}=2B_0\oplus2B_a\oplus2B_a,
\qquad
 \tau_a^{(2)}=C_a\oplus C_a\oplus C_a.
\tag{1.3}
\]

After the same unimodular change on both negative blocks,

\[
 G_a^{(2)}=-2q_d^{(2)},\qquad
 q_d^{(2)}
 =-Y^2+X^2+U^2+2dZ^2+V^2+2dW^2.
\tag{1.4}
\]

Thus \(L_a^{(2)}\) has signature \((1,5)\), is free of rank three over
\(A\), and its \(t\)-eigenspace has form

\[
 (t+1)\bigl(-(t-1)A_0^2+2dB^2+2dC^2\bigr).
\tag{1.5}
\]

The conjugate form is definite.

### Explicit primitive K3 embedding

Choose three pairwise nonadjacent simple roots
\(v_1,v_2,v_3\) in one \(E_8(-1)\) summand.  In the diagonal basis
\((Y,X,U,Z,V,W)\), send

\[
\begin{aligned}
Y&\longmapsto e_1+f_1,&
X&\longmapsto v_1,&
U&\longmapsto v_2,\\
Z&\longmapsto e_2-2df_2,&
V&\longmapsto v_3,&
W&\longmapsto e_3-2df_3.
\end{aligned}
\tag{1.6}
\]

The images are orthogonal and have squares

\[
 2,\ -2,\ -2,\ -4d,\ -2,\ -4d.
\tag{1.7}
\]

Because \(v_1,v_2,v_3\) are a subset of a simple-root basis, the six rows
\(e_1,v_1,v_2,e_2,v_3,e_3\) give a determinant-one minor.  The embedding in
\(\Lambda_{K3}\) is primitive.  A generic period in its signature-\((2,14)\)
complement gives a minimal projective K3 surface with exact Picard lattice
\(L_a^{(2)}\), and every returned NS operator is induced by a rational
divisor-product correspondence.

## 2. Five roots close the visible eigenslice

Use slope coordinates

\[
 s={B\over A_0},\qquad r={C\over A_0}.
\tag{2.1}
\]

The Klein disk is

\[
 s^2+r^2<R^2,\qquad
 R^2={t\over t-1}.
\tag{2.2}
\]

Since \(1<R^2<2\), the five halfspaces

\[
 s\ge0,\qquad r\ge0,\qquad
 s\le1,\qquad r\le1,\qquad
 s+r\le R^2
\tag{2.3}
\]

cut out a compact right-angled pentagon.  They come from the K3 roots

\[
\begin{aligned}
\alpha_s&=(0,0,-1,0,0,0),&
\alpha_r&=(0,0,0,0,-1,0),\\
\beta_s&=(1,-1,1,0,0,0),&
\beta_r&=(1,-1,0,0,1,0),\\
\gamma&=(1,0,1,0,1,0).
\end{aligned}
\tag{2.4}
\]

Each has \(q_d^{(2)}\)-norm one.  The last root is the mixed closing wall:

\[
                         s+r=R^2.
\tag{2.5}
\]

The complete restricted K3 chamber can have further walls, but it is
contained in the compact pentagon (2.3), hence is compact.  More importantly,
\(\gamma\) is a genuine exposed facet of the same chamber as the coordinate
wall \(\alpha_s\), not a redundant inequality.

### Proposition 2.1 (the mixed wall is exposed)

The midpoint

\[
                         s=r={R^2\over2}
\tag{2.6}
\]

lies on \(\gamma^\perp\) and strictly inside every other restricted K3 root
halfspace of the selected chamber.

### Proof

For a positive-height root, changing the signs of the \(X,U,V\)
coordinates independently can only increase its functional on the selected
sector when they are replaced by \(-|X|,|U|,|V|\).  It therefore suffices to
test the dominating normalization

\[
 \rho=(h,-m,u,z,v,w),\qquad h>0,\quad m,u,v\ge0,
\tag{2.7}
\]

with

\[
 -h^2+m^2+u^2+v^2+2d(z^2+w^2)=1.
\tag{2.8}
\]

At (2.6), twice the root functional is

\[
                         2L_\rho=B-tA,
\tag{2.9}
\]

where

\[
 A=2h-u-v+z+w,\qquad B=2m+z+w.
\tag{2.10}
\]

First \(A\ge0\).  Indeed, weighted Cauchy--Schwarz gives

\[
 (u+v-z-w)^2
 \le\left(2+{1\over d}\right)
       \bigl(u^2+v^2+2d(z^2+w^2)\bigr).
\tag{2.11}
\]

For \(h\ge2\), its right side is less than \(4h^2\).  At \(h=1\), (2.8)
is completely elementary.  If \(d\ge2\), then \(z=w=0\) and exactly two
of \(m,u,v\) equal one; if \(d=1\), the only additional cases have exactly
one of \(z,w\) equal to \(\pm1\) and all of \(m,u,v\) zero.  These cases
give \(A\ge0\), with \(A=0\) only for \(\rho=\gamma\).

Suppose \(A\ge1\) and \(L_\rho>0\).  Since
\(t>2d+1\) and \(A,B\) are integers,

\[
                         B\ge(2d+1)A+1.
\tag{2.12}
\]

After substituting (2.10), this says

\[
 F:=2m+(2d+1)(u+v)-2d(z+w)
       \ge2(2d+1)h+1.
\tag{2.13}
\]

Another weighted Cauchy estimate gives

\[
 F^2\le
 (8d^2+12d+6)(h^2+1).
\tag{2.14}
\]

For \(h\ge2\),

\[
 (8d^2+12d+6)(h^2+1)
 <4(2d+1)^2h^2,
\tag{2.15}
\]

contradicting (2.13).  In the \(h=1\) list above, the two cases with
\((m,u,v)\) supported in two coordinates give either \(\gamma\), or
\(B-tA=2-t<0\); the extra \(d=1\) cases give
\(B-tA=1-3t<0\) or \(-1-t<0\).  Hence equality occurs only for
\(\gamma\).  Thus
\(L_\rho\le0\), as claimed. \(\square\)

At height zero, the norm equation leaves only the \(X,U,V\) coordinate
roots.  The \(X\)-root has no projective eigenslice wall, and the other two
are \(s=0\) and \(r=0\).  The origin is strict for every positive-height
root; local finiteness of the integral root arrangement therefore gives a
short segment on \(\alpha_s^\perp\) exposing the coordinate facet.
Convexity joins that segment to (2.6) inside the same complete chamber.

### Corollary 2.2 (every orthogonal stabilization fails)

Replace (1.3) by

\[
 2B_0\oplus(2B_a)^{\oplus r},\qquad r\ge2.
\tag{2.16}
\]

For every pair of slope coordinates \(s_i,s_j\), the vector with
\(Y=U_i=U_j=1\) and every other coordinate zero is a norm-one lattice
root, and its restricted wall is

\[
                         s_i+s_j=R^2.
\tag{2.17}
\]

Set every other slope to zero.  The proof of Proposition 2.1 applies
verbatim: the unused squares only decrease the active energy in both
Cauchy estimates.  At the midpoint the only additional equalities are the
height-zero coordinate walls for the unused slopes.  Root-arrangement local
finiteness permits a small positive perturbation in those slopes, within
the mixed wall, so (2.17) is an exposed facet of the full selected chamber.
Its cyclic product with the \(s_i=0\) wall is still (3.4).  Hence every
orthogonal stabilization with \(r\ge2\) fails integral reflection
saturation.

This corollary is a statement about the orthogonal trace lattice and its
restricted root packet.  The explicit primitive K3 embedding and compact
pentagon asserted above remain rank-six statements; no uniform K3 embedding
or full-chamber compactness is asserted for arbitrary \(r\).

## 3. The unavoidable nonintegral cyclic product

Drop the harmless scalar \(t+1\) in (1.5), and write

\[
 c=t-1,\qquad b=2d,\qquad
 h=\operatorname {diag}(-c,b,b).
\tag{3.1}
\]

A normal for \(\alpha_s\) is

\[
 n_s=(0,1,0),\qquad n_s^2=b.
\tag{3.2}
\]

A normal for the mixed wall is

\[
 n_\gamma=(1,1,1),\qquad
 n_\gamma^2=2b-c=4d+1-t,
\tag{3.3}
\]

and \(n_s.n_\gamma=b\).  Their scale-invariant cyclic Cartan product is

\[
\begin{aligned}
 \chi_d
 &= {4(n_s.n_\gamma)^2\over n_s^2n_\gamma^2}\\
 &= {8d\over4d+1-t}
  ={4(2d-1+t)\over4d-1}.
\end{aligned}
\tag{3.4}
\]

Its field norm is

\[
 N_{k/\mathbf Q}(\chi_d)
 ={32d\over4d-1}
 =8+{8\over4d-1}.
\tag{3.5}
\]

This is never an integer: \(4d-1\equiv3\pmod4\) cannot divide \(8\).
Therefore \(\chi_d\) is not an algebraic integer.

This failure cannot be repaired by rescaling the normals.  Nor can it be
repaired by changing the order or passing to any commensurable full
restriction-of-scalars lattice.  Indeed, the product of the two field
reflections has eigenvalues \(1,\lambda,\lambda^{-1}\), and its trace on the
two-normal plane is \(\chi_d-2\).  If the restriction-of-scalars operator
preserved a \(\mathbf Z\)-lattice, both \(\lambda\) and \(\lambda^{-1}\)
would be eigenvalues of an integral matrix and hence algebraic integers.
Consequently
\(\chi_d=\lambda+\lambda^{-1}+2\) would be an algebraic integer,
contradicting (3.5).  This eigenvalue argument, unlike a rational-trace
argument, applies even when the commensurable \(\mathbf Z\)-lattice is not
itself presented as an \(A\)-module.

Equivalently, the mixed reflection matrix contains the unavoidable
denominator \(4d-1\).  For example,

\[
 {2b\over2b-c}
 ={4d-2+2t\over4d-1}.
\tag{3.6}
\]

## 4. Exact conclusion and remaining coupling problem

The orthogonal duplicate achieves:

- a primitive exact-Picard rank-six K3 carrier;
- a full rank-three \(A\)-trace module;
- a distinguished hyperbolic plane;
- a compact complete restricted eigen-nef component; and
- actual mixed K3 roots closing the extra direction.

It fails precisely where the strong non-rational return needs arithmetic
reflection integrality.  The exposed mixed facet cannot join the coordinate
facet in any integral arithmetic reflection host.

Thus a successful \(n\ge2\) construction cannot be a product or any
orthogonal stabilization of the rank-four interval.  It must alter the
coupling so that
every mixed closing wall has algebraic-integral cyclic products while
retaining:

1. a primitive K3 trace embedding;
2. compact or finite-volume complete restricted return;
3. one complete mirror orbit; and
4. reflection maximality.

No impossibility theorem for all such nonorthogonal couplings is claimed.
That is the exact surviving \(n\ge2\) field-height problem: construct such a
coupling, or bound the discriminant/height of every integral mixed-wall
packet from finite volume and maximality alone.
The exact replay is supplied by the
[certifier](../scripts/certify_rank6_k3_orthogonal_stabilization_nonarithmetic_no_go.py),
[archive](../data/enumerations/rank6_k3_orthogonal_stabilization_nonarithmetic_no_go.json),
[tests](../tests/test_rank6_k3_orthogonal_stabilization_nonarithmetic_no_go.py),
and
[hostile scope audit](rank6_k3_orthogonal_stabilization_nonarithmetic_no_go_hostile_scope_audit.md).

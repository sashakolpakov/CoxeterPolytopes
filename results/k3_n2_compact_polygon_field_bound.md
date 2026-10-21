# Every compact K3 Coxeter polygon has a finite field set

**Status (2026-08-11): proved for arbitrary compact Coxeter polygons in
\(\mathbf H^2\) with at least four facets.**  Let a rank-three module over a
totally real field return an exact compact Coxeter polygon.  The primitive
full trace has rank \(3d\), hence

\[
                              d=[k:\mathbf Q]\le6.
\tag{0.1}
\]

Integral adjacent reflections force every Coxeter angle to be at least
\(1/64\).  A direct hyperbolic sector-packing argument then finds two
nonadjacent facets whose dotted pair cycle satisfies

\[
                           4<p<69,798,469,636.
\tag{0.2}
\]

At every definite embedding the same pair cycle lies in \((0,4)\).
Therefore \(p\) has a unique large conjugate and generates the full field
\(k\).  It follows that

 \[
 |\operatorname {Disc}(k)|
 \le
 \bigl(69,798,469,636\bigr)^{30}.
\tag{0.3}
\]

Hermite--Minkowski gives a finite field set.  This closes the compact
\(n=2\) field-height gate for every facet count, without a Coxeter polygon
classification, reflective-lattice classification, imposed correspondence
height, saturation, or reflection maximality.

Compact triangles are covered separately by the finite-volume simplex
theorem; they have no nonadjacent facet pair.  Everything below assumes at
least four facets.

## 1. Exact carrier package

Let \(k\) be totally real of degree \(d\), let \(\mathcal O\subset k\) be a
returned order, and let \(M\) be a finite projective rank-three
\(\mathcal O\)-module with an admissible symmetric \(k\)-form \(h\).  At one
distinguished embedding \(\sigma_0\), the form has signature \((2,1)\); at
every other embedding it is positive definite.  Assume:

1. a minimal complex projective K3 surface \(X\) has Picard lattice equal
   to the primitive full trace of \((M,h)\), and the returned order acts
   self-adjointly by divisor-product correspondences;
2. the prescribed projectivized eigen-nef chamber is exactly a compact
   Coxeter polygon \(P\subset\mathbf H^2\) with \(N\ge4\) facets;
3. every facet reflection preserves \(M\); and
4. \(k\) is the minimal cyclic Cartan ground field of the returned packet.

Completeness, saturation, and declared-host reflection maximality may be
retained but are not used.  Exact full trace gives
\(3d=\rho(X)\le20\), proving (0.1).

Write the supporting geodesics cyclically as
\(H_0,\ldots,H_{N-1}\), their compact side facets as \(F_i=P\cap H_i\),
and the vertices as

\[
                              v_i=F_{i-1}\cap F_i
                              \qquad(i\bmod N).
\tag{1.1}
\]

Every interior angle \(\alpha_i\) is a Coxeter angle in \((0,\pi/2]\), so
the polygon is non-obtuse.

## 2. Acute support and projection lemmas

We record the elementary non-obtuse geometry used by the packing argument.

### Lemma 2.1 (facet projection)

Let \(F=P\cap H\) be a side of a convex non-obtuse hyperbolic polygon.

1. The orthogonal projection of every point of \(P\) onto \(H\) belongs
   to \(F\).
2. The metric projection onto \(P\) of every point of \(H\) belongs to
   \(F\).

For the first assertion, let \(x\in P\) and let \(z\in H\) be its
orthogonal foot.  If \(z\) lay beyond an endpoint \(q\) of the segment
\(F\), let \(K\) be the other side through \(q\).  In the hyperbolic
two-plane, the intersection of the inward half-planes bounded by \(H\)
and \(K\) is a wedge of angle at most \(\pi/2\).  Fermi coordinates along
\(H\) (equivalently, the sign rule in the resulting right triangle) show
that every point of this wedge projects orthogonally to the \(F\)-ray of
\(H\) issuing from \(q\).  Since \(x\) lies in both inward half-planes, its
foot cannot lie beyond \(q\).  Applying this at both endpoints proves the
first assertion.

For the second, let \(F^\perp\) be the union of all geodesics perpendicular
to \(H\) through points of \(F\).  The first assertion gives
\(P\subset F^\perp\).  This set is convex and invariant under reflection in
\(H\), with \(F^\perp\cap H=F\).  If \(x\in H\), uniqueness of metric
projection makes its nearest point in \(F^\perp\) reflection-fixed, hence
in \(F\subset P\).  Since \(P\subset F^\perp\), that same point is also the
nearest point of \(P\).  This proves the second assertion.

### Corollary 2.2 (nonadjacent supports are ultraparallel)

If two facets \(F_i,F_j\) are nonadjacent, their supports \(H_i,H_j\) are
disjoint.  Otherwise an intersection point \(x\in H_i\cap H_j\) would have
a unique metric projection to \(P\); Lemma 2.1(2) would put it in both
\(F_i\) and \(F_j\), a contradiction.

The compact facets have a positive attained distance.  If
\(a\in F_i,b\in F_j\) realize it, project \(a\) orthogonally to \(H_j\).
Lemma 2.1(1) puts its foot in \(F_j\), so uniqueness and minimality force
that foot to be \(b\).  Reversing the roles forces \(a\) to be the foot of
\(b\) on \(H_i\).  Thus \(ab\) is the common perpendicular and

\[
                        d(H_i,H_j)=d(F_i,F_j)>0.
\tag{2.1}
\]

Thus every nonadjacent pair is genuinely dotted, not merely asymptotic.
No structural classification is used in Lemma 2.1 or Corollary 2.2.

## 3. Vertex separation and contained sectors

Let

\[
 \delta=\min\{d(H_i,H_j):F_i,F_j\text{ nonadjacent}\}>0.
\tag{3.1}
\]

The following two distance facts include all cyclic edge cases.

### Lemma 3.1

1. Every vertex \(v_i\) is at distance at least \(\delta\) from every
   nonincident support \(H_j\).
2. Every two distinct vertices are at distance at least \(\delta\).

For (1), among the two incident facets \(F_{i-1},F_i\), at least one is
nonadjacent to the prescribed nonincident \(F_j\).  Since \(v_i\) lies on
its support,

\[
 d(v_i,H_j)\ge d(H_a,H_j)\ge\delta.
\]

For (2), first suppose \(v_i,v_{i+1}\) are adjacent vertices.  The outer
facets \(F_{i-1}\) and \(F_{i+1}\) are nonadjacent, so

\[
 d(v_i,v_{i+1})\ge d(H_{i-1},H_{i+1})\ge\delta.
\]

For nonadjacent vertices, orient the cyclic indexing so that
\(j=i+\ell\), where \(2\le\ell\le\lfloor N/2\rfloor\) is the shorter
cyclic separation.  If \(\ell=2\), choose the incident facets
\(F_{i-1}\ni v_i\) and \(F_{i+1}=F_{j-1}\ni v_j\).  Their indices differ
by two in both cyclic directions when \(N=4\), and by at least two in the
shorter direction when \(N>4\), so the facets are nonadjacent.  If
\(\ell\ge3\), choose \(F_i\ni v_i\) and \(F_{j-1}\ni v_j\).  Their two
cyclic index gaps are \(\ell-1\ge2\) and
\(N-(\ell-1)\ge N-\lfloor N/2\rfloor+1\ge3\), so again they are
nonadjacent.  The chosen supports contain the two vertices, hence their
support distance is at most the vertex distance and (2) follows.

Put \(r=\delta/2\).  By Lemma 3.1(2), the open balls
\(B(v_i,r)\) are pairwise disjoint.  By Lemma 3.1(1), no nonincident support
meets \(B(v_i,r)\).  Therefore

\[
                         P\cap B(v_i,r)
\tag{3.2}
\]

is the full hyperbolic sector of radius \(r\) and angle \(\alpha_i\).
These \(N\) sectors have pairwise disjoint interiors and are contained in
\(P\).  A sector of angle \(\alpha\) and radius \(r\) has area

\[
                         \alpha(\cosh r-1).
\tag{3.3}
\]

## 4. Reflection integrality supplies a uniform angle floor

Let \(s_i,s_{i+1}\) be the reflections in adjacent facets and put
\(P_i=1-s_i\).  On the finite projective \(\mathcal O\)-module,
\(P_iP_{i+1}\) is an \(\mathcal O\)-linear endomorphism and

\[
 q_i:=\operatorname {Tr}_{\mathcal O}(P_iP_{i+1})
     =4\cos^2\alpha_i
     \in\mathcal O.
\tag{4.1}
\]

Set

\[
                         z_i=4-q_i=4\sin^2\alpha_i\in\mathcal O.
\tag{4.2}
\]

At the distinguished place \((z_i)_0>0\).  At every definite embedding,
positive-definite Cauchy--Schwarz gives \(0<\sigma(z_i)\le4\); the endpoint
\(4\) occurs only when the adjacent pairing is zero.  Thus the nonzero
integer norm satisfies

\[
 1\le N_{k/\mathbf Q}(z_i)
   =(z_i)_0\prod_{\sigma\ne\sigma_0}\sigma(z_i)
 \le (z_i)_0\,4^{d-1},
\]

and hence

\[
 (z_i)_0\ge4^{-(d-1)},
 \qquad
 \sin^2\alpha_i\ge4^{-d}\ge4^{-6},
 \qquad
 \alpha_i\ge\sin\alpha_i\ge{1\over64}.
\tag{4.3}
\]

## 5. Sector area bounds the shortest dotted pair

Summing the disjoint sector areas from Section 3 and using (4.3),

\[
 \operatorname {Area}(P)
 \ge {N\over64}\bigl(\cosh(\delta/2)-1\bigr).
\tag{5.1}
\]

Gauss--Bonnet gives

\[
 \operatorname {Area}(P)
 =(N-2)\pi-\sum_i\alpha_i
 <N\pi.
\tag{5.2}
\]

Combining (5.1)--(5.2), and using \(\pi<4\), yields

\[
 \cosh(\delta/2)<1+64\pi<257.
\tag{5.3}
\]

Therefore

\[
\begin{aligned}
 \cosh\delta
 &=2\cosh^2(\delta/2)-1<2\cdot257^2-1=132,097,\\
 4\cosh^2\delta
 &<4\cdot132,097^2=69,798,469,636.
\end{aligned}
\tag{5.4}
\]

The pair realizing \(\delta\) proves (0.2).

## 6. The shortest dotted cycle generates the full field

Let \(e_i,e_j\) be normals for the minimizing nonadjacent facets.  Their
pair cyclic product is exactly

\[
 p=\operatorname {Tr}_{\mathcal O}(P_iP_j)
  ={4h(e_i,e_j)^2\over h(e_i,e_i)h(e_j,e_j)}
  =4\cosh^2\delta
  \in\mathcal O.
\tag{6.1}
\]

At \(\sigma_0\), (0.2) holds.  At every other embedding, the form is
positive definite, so

\[
                              0<\sigma(p)<4.
\tag{6.2}
\]

Put \(F=\mathbf Q(p)\).  If \([k:F]>1\), the restriction of \(\sigma_0\)
to the totally real field \(F\) has \([k:F]\) distinct real extensions to
\(k\).  Every extension takes the same value \(p_0>4\), whereas only
\(\sigma_0\) is Lorentzian and every other embedding must satisfy (6.2).
This contradiction proves

\[
                                  \mathbf Q(p)=k.
\tag{6.3}
\]

All conjugates of this primitive algebraic integer lie in the interval
\((0,69,798,469,636)\).  Its Vandermonde discriminant and (0.1) give

\[
\begin{aligned}
 |\operatorname {Disc}(k)|
 &\le|\operatorname {Disc}(\mathbf Z[p])|\\
 &<\bigl(69,798,469,636\bigr)^{d(d-1)}\\
 &\le\bigl(69,798,469,636\bigr)^{30},
\end{aligned}
\tag{6.4}
\]

which proves (0.3).

### Theorem 6.1 (general compact-polygon field theorem)

Every exact full-trace complex-K3 carrier satisfying Section 1 and having
an arbitrary compact Coxeter polygon with at least four facets as its
prescribed eigen-nef chamber belongs to the finite set of totally real
fields bounded by (0.3).

## 7. Scope

The proof uses compactness twice: facets are compact and their disjoint
distances are positive, and Gauss--Bonnet has no ideal-vertex sectors.
It does not cover finite-area noncompact polygons, ultraideal/truncated
polygons, or higher-dimensional chambers.  It does not construct the K3
carrier, establish essential surjectivity, or prove the separate
Property-4 correspondence return.  Saturation and maximality are not used.

This is a dimension-two automatic-height theorem in intrinsic chamber
geometry: the K3 trace bound controls the angles, and area forces one short
dotted pair which already carries the whole field.

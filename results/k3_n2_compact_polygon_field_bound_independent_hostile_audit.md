# Independent hostile audit: the compact-polygon field theorem

**Verdict (2026-08-11): PASS in the exact compact, rank-three, primitive
full-trace K3 scope with \(N\ge 4\).** The argument really produces a
uniformly bounded primitive algebraic integer from the polygon itself. No
classification, imposed correspondence height, saturation, or reflection
maximality enters the proof.

This audit is independent of the owner and local-hostile proof code. It
checks the frozen owner report, certifier, archive, and test bytes, then
re-derives the geometric and arithmetic pressure points below.

## 1. The non-obtuse projection lemma

The projection statement needs both half-planes at an endpoint, not the
other-side half-plane alone. The frozen proof does use their intersection,
and the following Fermi-coordinate calculation verifies the precise claim.

Place the support \(H\) at \(r=0\), an endpoint \(q\) at \(t=0\), and the
facet ray at \(t\ge 0\). In hyperboloid Fermi coordinates,

\[
 x(t,r)=(\cosh r\cosh t,\ \cosh r\sinh t,\ \sinh r).
\]

If the other support makes interior angle
\(0<\alpha\le\pi/2\), the inward non-obtuse wedge is described, after
choosing the inward side, by

\[
 r\ge 0,\qquad
 -\sin\alpha\,\cosh r\,\sinh t+\cos\alpha\,\sinh r\le 0.
\tag{1.1}
\]

For \(t<0\), the first summand on the left of (1.1) is strictly positive
and the second is nonnegative. Thus (1.1) forces \(t\ge 0\): orthogonal
projection to \(H\) lands on the facet ray. Applying this at both endpoints
proves \(\operatorname{proj}_H(P)\subset F\).

There is no hidden convexity gap in the reverse projection. If
\(F=\{(t,0):a\le t\le b\}\), then the full perpendicular strip is

\[
\begin{aligned}
 F^\perp
 &=\{a\le t\le b\}\\
 &=\{x_1\cosh a-x_0\sinh a\ge 0\}
   \cap\{x_1\cosh b-x_0\sinh b\le 0\}.
\end{aligned}
\tag{1.2}
\]

It is therefore an intersection of two geodesic half-planes, hence convex,
and it is invariant under reflection in \(H\). Since
\(P\subset F^\perp\), uniqueness and reflection invariance put the
projection of any \(x\in H\) to \(F^\perp\) in
\(F^\perp\cap H=F\subset P\). The inclusions
\(P\subset F^\perp\) and \(F\subset P\) show that this same point is the
metric projection to \(P\).

Consequently, intersecting nonadjacent supports would project their
intersection into two disjoint facets. They cannot intersect. Minimizing
the positive distance between the two compact facets and applying the
first projection in both directions puts both feet on the facets; the
joining segment is the common perpendicular. Hence every nonadjacent
support pair is ultraparallel and

\[
 d(H_i,H_j)=d(F_i,F_j)>0.
\tag{1.3}
\]

## 2. Vertex separation and sector packing for every \(N\ge 4\)

Let \(\delta\) be the least distance in (1.3). For a vertex and a
nonincident wall, at least one of the two incident facets is nonadjacent
to that wall. For adjacent vertices, their two outer facets are
nonadjacent. If the shorter cyclic vertex separation is \(2\), use
\(F_{i-1}\) and \(F_{i+1}\); this also covers the quadrilateral
wraparound. If it is \(\ell\ge 3\), use \(F_i\) and \(F_{j-1}\); their
cyclic gaps are \(\ell-1\ge 2\) and
\(N-(\ell-1)\ge 3\). These witnesses prove

\[
 d(v_i,H_j)\ge\delta
 \quad\text{for every nonincident }H_j,\qquad
 d(v_i,v_j)\ge\delta
 \quad(i\ne j).
\tag{2.1}
\]

Thus the open balls \(B(v_i,\delta/2)\) are pairwise disjoint. Each such
ball stays in every nonincident inward half-plane: its connected interior
contains \(v_i\) and meets none of their supports. Its intersection with
\(P\) is therefore the untruncated sector determined by the two incident
supports. A radius-\(r\), angle-\(\alpha\) sector has area

\[
 \int_0^\alpha\int_0^r\sinh t\,dt\,d\theta
 =\alpha(\cosh r-1).
\tag{2.2}
\]

The sectors are contained in \(P\) and have disjoint interiors, so their
areas may be summed.

## 3. The norm gap gives the angle floor

Exact full trace gives \(3d=\rho(X)\le 20\), hence \(d\le 6\). For adjacent
reflections set \(P_i=1-s_i\). The canonical trace on the finite
projective \(\mathcal O\)-module gives

\[
 q_i=\operatorname{Tr}_{\mathcal O}(P_iP_{i+1})
     =4\cos^2\alpha_i\in\mathcal O,\qquad
 z_i=4-q_i=4\sin^2\alpha_i\in\mathcal O.
\tag{3.1}
\]

At the distinguished place \(z_i>0\). At every definite place strict
Cauchy--Schwarz for the two linearly independent facet normals gives
\(0<\sigma(z_i)\le 4\). Thus \(z_i\) is a nonzero algebraic integer and

\[
 1\le N_{k/\mathbf Q}(z_i)
 \le (z_i)_0\,4^{d-1}.
\]

It follows that

\[
 \sin^2\alpha_i\ge 4^{-d},\qquad
 \sin\alpha_i\ge 2^{-d}\ge\frac1{64},\qquad
 \alpha_i\ge\frac1{64}.
\tag{3.2}
\]

The exponent is \(d\), not \(d-1\), after dividing \(z_i\) by \(4\).

## 4. Exact constants and the primitive generator

Equations (2.2) and (3.2), followed by Gauss--Bonnet, give

\[
 \frac N{64}\bigl(\cosh(\delta/2)-1\bigr)
 \le\operatorname{Area}(P)
 =(N-2)\pi-\sum_i\alpha_i<N\pi.
\]

Using \(\pi<4\),

\[
\begin{aligned}
 \cosh(\delta/2)&<257,\\
 \cosh\delta&<2\cdot257^2-1=132\,097,\\
 p=4\cosh^2\delta&<4\cdot132\,097^2
                  =69\,798\,469\,636.
\end{aligned}
\tag{4.1}
\]

The minimizing supports are ultraparallel, so \(p>4\). Projective-module
trace integrality gives \(p\in\mathcal O\). At every definite embedding,
strict Cauchy--Schwarz gives \(0<\sigma(p)<4\); strict positivity follows
because the normal pairing is a nonzero element of \(k\) at the
distinguished place. Hence \(p\) has exactly one conjugate larger than
\(4\).

If \(F=\mathbf Q(p)\) were proper in \(k\), the distinguished embedding
of \(F\) would have \([k:F]\ge 2\) real extensions to the totally real
field \(k\). Every extension would give \(p\) the same value \(p_0>4\),
contradicting the unique Lorentzian place. Therefore
\(\mathbf Q(p)=k\). All conjugates of this primitive algebraic integer lie
in \((0,B)\), where \(B=69\,798\,469\,636\), and Vandermonde gives

\[
 \lvert\operatorname{Disc}k\rvert
 \le\lvert\operatorname{Disc}\mathbf Z[p]\rvert
 <B^{d(d-1)}
 \le B^{30}.
\tag{4.2}
\]

This is a discriminant bound and therefore a literal finite field set,
not merely a degree bound.

## 5. Exact scope

The proof assumes a finite projective rank-three order-module, one
Lorentzian and all other positive-definite places, a primitive exact
full-trace Picard lattice on a complex projective K3 surface, integral
facet reflections, the returned self-adjoint correspondence action, and
the exact same compact Coxeter eigen-nef chamber with at least four facets.
The owner also retains the minimal cyclic Cartan ground-field hypothesis;
the unique-large-place argument actually forces \(\mathbf Q(p)=k\) without
using that minimality. The theorem does not construct the carrier or prove
essential surjectivity or Property 4. Compactness is essential to the
positive attained facet distances and the vertex-sector argument.
Finite-area noncompact, ultraideal/truncated, and higher-dimensional
chambers remain outside this theorem. Triangles require the separate
simplex theorem.

The independent executable audit is
[audit_k3_n2_compact_polygon_field_bound_independent_hostile.py](../scripts/audit_k3_n2_compact_polygon_field_bound_independent_hostile.py),
with separate tests in
[test_k3_n2_compact_polygon_field_bound_independent_hostile.py](../tests/test_k3_n2_compact_polygon_field_bound_independent_hostile.py).

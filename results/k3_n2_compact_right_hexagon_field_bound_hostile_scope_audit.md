# Hostile audit: the compact right-angled K3 hexagon field bound

**Verdict (2026-08-11): PASS in exact compact right-angled scope.**  The
rank-three closure does not bound every hexagon modulus: one alternating
side can escape while another derived side also grows.  The theorem does
not need that stronger statement.  It proves that at least one dotted pair
is uniformly short, and the admissible signed-place pattern makes that one
pair a primitive field generator.

## 1. Cofactor closure audit

For alternating side lengths put

\[
 x_0=\sinh^2\ell_0,
 \quad x_2=\sinh^2\ell_2,
 \quad x_4=\sinh^2\ell_4,
\]

and \(A=x_0+1\), \(C=x_2+1\), \(E=x_4+1\).  Inverting the Gram matrix of
the three alternating normals gives

\[
\begin{aligned}
 S&=A+C+E-1+2\sqrt{ACE},\\
 x_1&={S\over x_0x_2},
 &x_3&={S\over x_2x_4},
 &x_5&={S\over x_4x_0}.
\end{aligned}
\]

The replay checks these cofactor identities and all six compatible
four-consecutive Gram-minor identities on exact rational grids.  Positive
square roots occur only at the distinguished real place; no square root is
silently asserted to lie in \(k\).

## 2. Short-pair inequality audit

If an alternating \(x_i<16\), the conclusion is immediate.  Otherwise
order the alternating values \(x\le y\le z\), all at least \(16\), and use
the odd variable with denominator \(yz\).  Monotonicity gives

\[
 x_{\mathrm{odd}}
 <{51+34\sqrt{17}\over256}
 <{221\over256}<1.
\]

Thus some second-neighbour shift is always below \(16\), even though other
coordinates can escape.  Its pair cycle is below \(68\).

## 3. Primitive-field audit

Integral facet reflections put the short pair cycle \(p\) in the returned
order.  At the distinguished place,

\[
                              4<p_0<68,
\]

whereas positive-definite Cauchy--Schwarz puts every other conjugate in
\((0,4)\).  If \(p\) lay in a proper subfield \(F\subsetneq k\), the
distinguished embedding of \(F\) would have \([k:F]>1\) real extensions to
the totally real field \(k\).  Each extension would give the same value
\(p_0>4\), contradicting uniqueness of the Lorentzian place.  Hence
\(\mathbf Q(p)=k\).

This subfield argument is the decisive point.  Bounded degree alone is not
being mistaken for a finite field set, and no pair-generation assumption is
silently imported.

## 4. Discriminant and scope

The full trace gives \(d\le6\).  Every conjugate of the primitive integral
generator \(p\) lies in \((0,68)\), so the Vandermonde bound is

\[
                         |\operatorname {Disc}(k)|\le68^{30}.
\]

Hermite--Minkowski proves the finite field set.  The theorem neither uses
maximality nor proves carrier existence.  It does not cover non-right or
noncompact hexagons, arbitrary six-facet incidence, or higher-facet
polygons.  Its general output is the reduction of the latter problem to a
uniform short-dotted-pair theorem.  No reflective-lattice classification
enters.

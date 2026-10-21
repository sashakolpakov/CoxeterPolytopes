# Hostile audit: the general compact K3 Coxeter-polygon field bound

**Verdict (2026-08-11): PASS in the literal compact, rank-three,
full-trace scope with at least four facets.**  The proof supplies its own
height: an integral norm gap prevents arbitrarily small Coxeter angles,
and disjoint vertex sectors force one uniformly short dotted pair.  The
signed-place pattern then makes that pair a primitive generator of the
whole field.

## 1. Acute geometry audit

For a side \(F\subset H\), the two inward half-planes at either endpoint
form a wedge of angle at most \(\pi/2\).  Orthogonal projection of that
wedge to \(H\) lands on the ray occupied by \(F\).  Applying this at both
endpoints proves that the projection of \(P\) to \(H\) lands in \(F\).
The perpendicular strip over \(F\) is convex and reflection-invariant, so
the reverse projection, from \(H\) to \(P\), also lands in \(F\).

If nonadjacent supports met, the projection of their intersection to
\(P\) would lie in two disjoint facets.  Thus the supports are disjoint.
Projecting a minimizing point of either compact facet to the other support
shows that the common perpendicular has both feet on the facets.  Hence

\[
                 d(H_i,H_j)=d(F_i,F_j)>0
\]

for every nonadjacent pair.  This proves ultraparallelity and positivity;
neither is silently inferred merely from disjoint compact side segments.

## 2. Cyclic separation and sector audit

Let \(\delta\) be the least nonadjacent-support distance.  A vertex and a
nonincident wall admit a nonadjacent incident-side witness.  Two adjacent
vertices use their two outer sides.  For vertices at shorter cyclic
separation \(\ell\ge2\), separation two uses \(F_{i-1},F_{i+1}\), while
separation at least three uses \(F_i,F_{j-1}\).  Both cyclic gaps are then
at least two.  The exact replay checks these selectors for all polygon
sizes from four through 128, including every quadrilateral wraparound.

It follows that distinct vertices are at least \(\delta\) apart and that
each vertex is at least \(\delta\) from every nonincident support.  The
open radius-\(\delta/2\) vertex balls are therefore pairwise disjoint, and
all remote half-planes contain each such ball.  Their intersections with
\(P\) are full sectors, not truncated pieces.  Their total area is at
least

\[
 {N\over64}\bigl(\cosh(\delta/2)-1\bigr).
\]

## 3. Angle norm audit

On the finite projective \(\mathcal O\)-module, the trace of the
\(\mathcal O\)-linear endomorphism \((1-s_i)(1-s_{i+1})\) is an algebraic
integer

\[
 q_i=4\cos^2\alpha_i.
\]

Thus \(z_i=4-q_i=4\sin^2\alpha_i\) is integral.  It is positive at every
real place and at most four at every definite place.  The norm inequality
is exactly

\[
 1\le N_{k/\mathbf Q}(z_i)
   =(z_i)_0\prod_{\sigma\ne\sigma_0}\sigma(z_i)
 \le (z_i)_0\,4^{d-1}.
\]

Consequently \((z_i)_0\ge4^{-(d-1)}\), so
\(\sin\alpha_i\ge2^{-d}\ge1/64\) because the K3 full trace gives
\(3d\le20\), hence \(d\le6\).  This is the correct exponent and direction;
bounded degree alone was not substituted for the norm gap.

## 4. Area, height, and primitive-field audit

Gauss--Bonnet gives

\[
 \operatorname {Area}(P)=(N-2)\pi-\sum_i\alpha_i<N\pi.
\]

Comparison with the sectors and \(\pi<4\) yields

\[
 \cosh(\delta/2)<257,
 \qquad
 \cosh\delta<132,097,
\]

so the minimizing dotted pair cycle satisfies

\[
 4<p=4\cosh^2\delta<69,798,469,636.
\]

The same finite-projective-module trace makes \(p\in\mathcal O\).  At
every definite place, positive-definite Cauchy--Schwarz gives
\(0<\sigma(p)<4\).  If \(p\) belonged to a proper subfield \(F\), the
distinguished embedding of \(F\) would have \([k:F]>1\) real extensions,
all taking the repeated value \(p_0>4\).  This contradicts the unique
Lorentzian place.  Therefore \(\mathbf Q(p)=k\), and the Vandermonde bound
is

\[
 |\operatorname {Disc}(k)|
 <69,798,469,636^{d(d-1)}
 \le69,798,469,636^{30}.
\]

Hermite--Minkowski now gives a finite field set.

## 5. Scope verdict

The theorem is classification-free and does not assume an imposed mixed
height, saturation, or reflection maximality.  It is not a theorem for
finite-area noncompact polygons, ultraideal or truncated polygons,
higher-dimensional chambers, or arbitrary abstract correspondence
packets.  It assumes the exact primitive full-trace K3 carrier and proves
neither carrier existence, essential surjectivity, nor Property 4.
Triangles are intentionally separate because they have no nonadjacent
facet pair.


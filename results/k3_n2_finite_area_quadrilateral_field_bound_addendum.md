# Addendum: finite-area K3 Coxeter quadrilaterals

**Status (2026-08-11): proved.**  The frozen compact-quadrilateral field
theorem extends unchanged when one or more vertices are ideal.  The field
degree remains at most six and the same bounds hold:

\[
 |\operatorname{Disc}(k)|
 \le(590,584,093,515,710,595,072)^{30},
\]

with the sharper quadratic bound

\[
                         |\operatorname{Disc}(k)|\le4,857,950,597.
\]

## The endpoint extension

Retain every hypothesis of the frozen theorem—primitive exact full trace,
integral facet reflections, exact eigen-nef return, and the requirement that
the cyclic Cartan packet generate the returned field—but let the Coxeter
quadrilateral have finite or ideal vertices.  An acute-angled polygon has
disjoint supporting lines for nonadjacent facets, so the two opposite wall
pairs remain ultraparallel: otherwise an intervening side and the two
meeting support lines would cut out a hyperbolic triangle with two angles
supplementary to polygon angles, hence both at least \(\pi/2\), contradicting
the strict hyperbolic triangle angle sum. Thus their normalized magnitudes
satisfy \(U,V>1\), while the four adjacent magnitudes now satisfy

\[
                         0\le a,b,c,d_0\le1.
\]

The opposite pair cycles are still \(4A=4U^2\) and \(4C=4V^2\). At every
definite place their squared ratios lie strictly between zero and one.
The same nonzero-norm argument therefore gives the strict floor

\[
                         A-1,C-1>4^{-d}\ge4^{-6}.           \tag{1}
\]

Only one line of the determinant proof changes.  The coefficient estimates
become weak:

\[
 0\le P,Q,R\le2,\qquad |K|\le8,
\]

and hence

\[
 (U^2-1)(V^2-1)\le4U+4V+4UV+8.                            \tag{2}
\]

This does not weaken the final diagonal bound. From (1),
\(V>\sqrt{4097}/64\) strictly. Both coefficient functions obtained after
dividing (2) by \(V^2-1\), namely

\[
                  {4\over V-1},\qquad {4V+8\over V^2-1},
\]

are strictly decreasing for \(V>1\). Evaluating them at the excluded lower
endpoint therefore restores the strict inequality

\[
 U^2<(16384+256\sqrt{4097})U
        +(32769+256\sqrt{4097}).                           \tag{3}
\]

The frozen exact trial remains valid:

\[
 537034767-8389888\sqrt{4097}>0,
\]

because its square difference is \(17{,}606{,}953{,}271{,}521\).
Consequently \(U,V<32772\), exactly as before.

At an ideal vertex an adjacent normalized factor may equal one, so the
simple-cycle height is now a weak rather than strict bound.  Its numerical
value is unchanged:

\[
 H=16\cdot32772^4=18,455,752,922,365,956,096.
\]

The complete cyclic packet still consists of cycles of length at most four.
The same two-generator subpacket, \(16\)-point primitive grid, and
Vandermonde argument give the displayed discriminant bounds and hence a
finite field set.

This addendum covers finite-area quadrilaterals whose vertices are finite or
ideal.  It does not cover ultraideal/truncated or infinite-area cases,
arbitrary higher-sided polygons, gratuitous operator-field extensions, or
carrier construction.  Completeness, saturation, and declared-host
maximality are retained but unused in the numerical estimate.

The exact replay is
[certify_k3_n2_finite_area_quadrilateral_field_bound_addendum.py](../scripts/certify_k3_n2_finite_area_quadrilateral_field_bound_addendum.py),
with tests in
[test_k3_n2_finite_area_quadrilateral_field_bound_addendum.py](../tests/test_k3_n2_finite_area_quadrilateral_field_bound_addendum.py).

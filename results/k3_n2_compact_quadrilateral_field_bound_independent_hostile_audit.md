# Independent hostile audit: compact K3 Coxeter quadrilaterals

**Verdict (2026-08-11): PASS under the exact frozen carrier hypotheses.**
The conclusion is a finite set of totally real fields, not merely bounded
degree.  For degrees above two the proof correctly uses the complete bounded
cyclic Cartan packet; it does not pretend that one opposite-pair coordinate
generates the field.

## 1. Full trace and the norm gap

A rank-three order module has full trace rank (3d).  Equality with the
Picard lattice of a complex projective K3 gives

\[
                         3d=\rho(X)\le20,
\]

and hence (d\le6).  For an opposite pair let
(A=h_{ij}^2/(h_{ii}h_{jj})).  At the Lorentzian place (A>1), whereas at
every definite place strict Cauchy--Schwarz gives (0<\sigma(A)<1).
Nonzero pairings exclude zero, and a nonzero wedge excludes equality.

The product of the two Cartan entries is (4A\in\mathcal O): equivalently,
it differs by an integer from the module trace of the product of the two
integral reflections.  Therefore (4(A-1)) is a nonzero algebraic integer.
Its distinguished value is positive and its other (d-1) conjugates have
absolute value strictly below four, so

\[
 1\le |N(4(A-1))|<4^d(A_0-1),\qquad A_0-1>4^{-d}.
\]

The same holds for the other opposite ratio (C).  This rederives the
power (4^{-d}), including the strict endpoint.

## 2. Determinant and diagonal bound

I recomputed the (4\times4) determinant as a formal sparse polynomial,
not by importing or numerically sampling the owner formula.  With adjacent
magnitudes (a,b,c,d_0<1), opposite magnitudes (U,V>1), and
(r=U^2-1, s=V^2-1), it is exactly

\[
 0=rs+K-2PU-2QV-2RUV,
\]

where

\[
\begin{aligned}
P&=ab+cd_0,&Q&=ad_0+bc,&R&=ac+bd_0,\\
K&=-a^2-b^2-c^2-d_0^2+a^2c^2+b^2d_0^2-2abcd_0.
\end{aligned}
\]

Thus (P,Q,R<2, |K|<8) imply

\[
 (U^2-1)(V^2-1)<4U+4V+4UV+8.                 \tag{2.1}
\]

At the weakest degree-six norm gap,
(V>\sqrt{4097}/64).  Dividing (2.1) by (V^2-1), and using the fact that
both resulting coefficient functions decrease on (V>1), gives

\[
U^2<(16384+256\sqrt{4097})U+(32769+256\sqrt{4097}).
\]

At (U=32772) the remainder is

\[
537034767-8389888\sqrt{4097}>0,
\]

with exact square difference (17,606,953,271,521).  The constant term is
negative in the associated quadratic, so this trial lies past its unique
positive root.  Cyclic rotation yields

\[
                     U,V<32772,qquad A,C<1,074,003,984.
\]

## 3. The complete field-generating packet

For an integral reflection (s_i), the endomorphism (1-s_i) is integral.
The trace of a cyclic product of such rank-one endomorphisms is precisely
the corresponding cyclic Cartan product, so every cycle coordinate is an
algebraic integer.  A closed walk factors at repeated vertices into pair
cycles and simple cycles.  Four facets therefore require only lengths two,
three, and four.

At definite places every normalized factor is below one; at the distinguished
place every factor is below (32772).  Hence every generator in the simple
cycle packet has all conjugates bounded by

\[
 H=16\cdot32772^4=18,455,752,922,365,956,096.
\]

The frozen theorem explicitly assumes that this whole packet generates
(k).  A minimal generating subpacket has at most two elements, because
three strict adjunctions would force degree at least (2^3>6).  For two
generators (\theta_1,\theta_2), each unordered pair of the at most six
embeddings forbids at most one coefficient in
(\theta_1+c\theta_2).  There are at most ({6\choose2}=15) bad values,
so some (c\in\{0,\ldots,15\}) gives a primitive integral element with
all conjugates bounded by (16H).

The Vandermonde estimate is consequently

\[
 |\operatorname{Disc}(k)|
 \le(32H)^{d(d-1)}
 \le(590,584,093,515,710,595,072)^{30}.
\]

Hermite--Minkowski now gives an actual finite field set.  A single-coordinate
generation assertion would be unjustified for (d>2), but none is used.

## 4. Degree-two specialization and scope

For (d=2), the gap is (1/16).  Repeating the exact quadratic comparison
at (U=132) gives the positive remainder

\[
8847-2128\sqrt{17},qquad
8847^2-17\cdot2128^2=1,286,881.
\]

Thus (U,V<132).  Now (alpha=4A) has its two conjugates in the disjoint
intervals ((4,69696)) and ((0,4)), so it is nonrational and generates
the quadratic field.  Its integral trace is at most (69,699), its norm is
at least one, and

\[
                    |\operatorname{Disc}(k)|\le4,857,950,597.
\]

The theorem is conditional on a primitive exact full-trace K3 Picard
carrier, a compact four-facet eigen-nef chamber, integral facet reflections,
and the literal cyclic-ground-field hypothesis.  Completeness, saturation,
and declared-host reflection maximality are retained but are not used in
the numerical bound.  The result does not cover noncompact quadrilaterals,
gratuitous operator-field extensions, or arbitrary compact polygons.

The independent executable audit is
[audit_k3_n2_compact_quadrilateral_field_bound_independent_hostile.py](../scripts/audit_k3_n2_compact_quadrilateral_field_bound_independent_hostile.py),
with separate tests in
[test_k3_n2_compact_quadrilateral_field_bound_independent_hostile.py](../tests/test_k3_n2_compact_quadrilateral_field_bound_independent_hostile.py).

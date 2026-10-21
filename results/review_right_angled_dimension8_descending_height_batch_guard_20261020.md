# Independent adversarial review of the descending-height proof

20 October2026. **Full written mathematical review PASS**, relative to
the explicitly published six-dimensional27-facet theorem. No computation,
optimizer, scientific check, Lean compilation or build was run.

Reviewed in full:
`right_angled_dimension8_descending_height_proof_20261020.md`, SHA256
`9470fa1fa3a973a12c412fbc1257537655536dba4ffc4e9ac23c320615babdbe`.
The reviewed source remains unchanged. This review independently derives
the real projection, the tensor identities and the strict endpoint.

## 1. The external facet minimum has the required universal scope

The needed input is that EVERY finite-volume right-angled hyperbolic
six-polytope has at least27 facets, allowing finite and ideal vertices.
It is not the construction of a particular27-facet P6.

I directly opened Alexandrov's original2023 mathematical text,
[Section4](https://arxiv.org/html/2303.09533v1#S4), which states that
universal six-dimensional assertion and attributes it to Dufour,
Lemma1 and Proposition4. Section6 separately records the two steps:
Proposition6.4 gives nu5>=26 and Proposition6.5 gives
F_n>=1+nu_(n-1). The PDF has the same explicit assertion on its
fourth page. Egorov's original2026 Section2 and Section3.2 also
explicitly restate the resulting27-facet bound.

The original Dufour full text was not retrieved in this review;
the publisher abstract alone would not have verified this particular
input. The boundary is therefore the named published Dufour result,
with its precise statement independently checked in the primary
mathematical texts of Alexandrov and Egorov. Potyagailo--Vinberg2005,
Proposition3.1, constructs P6 with27 facets but its final question
does NOT by itself establish the universal minimum. The author note
correctly avoids that misattribution.

## 2. Three actual projections and the strict real bound

Let q be a future unit point in the relative interior of an actual
right-angled polytope, and let n maximize t=-q.n>0 over its finitely
many intrinsic unit facet normals. The proposed numerator q+t*n
has square-1-t^2 and product0 with n. Its product with q is
-1-t^2<0, so its normalized point q' is future.

For EVERY other facet m, q.m is strictly negative and n.m<=0 by
all-wall acuteness. Consequently q'.m<0. Thus q' lies in the relative
interior of the ACTUAL facet, not just its support plane. This strict
feasibility is preserved at every stage. Each intrinsic facet normal
is an actual orthogonal ridge neighbor, so it remains unit after
restriction and has new height

    (-q.m)/sqrt(1+t^2)<=t/sqrt(1+t^2).

The positive function x/sqrt(1+x^2) is strictly increasing, with
limit1. Starting with a finite maximum, the first bound is strictly
less than1, the second strictly less than1/sqrt2, and the third
strictly less than1/sqrt3. Ties for a maximum do not remove strictness.
Equivalently its k-fold iterate at a finite a>0 equals
a/sqrt(1+k*a^2), which is strictly less than1/sqrt(k).

The selected face dimensions are8,7,6 in a hypothetical nine-polytope.
Right-angle and finite-volume inheritance apply throughout. Restriction
does not introduce a new normalization factor: all selected intrinsic
facet normals are already perpendicular to the earlier supporting
normals. There is no descent to an ideal polygon and no substitution
of an average facet for an actual facet.

In the final Lorentz seven-space, write each actual normal as
(t_i,sqrt(1+t_i^2)*u_i). Here u_i is a unit vector in R6, and
0<t_i<1/sqrt3. Pairwise all-wall signs give

    u_i.u_j<=t_i*t_j/sqrt((1+t_i^2)(1+t_j^2))<1/4.

Repeated directions would give a positive Lorentz wall product and
are excluded. The resulting cardinality is exactly the actual
six-face facet count.

## 3. Independent tensor-square and cardinality calculation

For unit u,v in R6 set x=u.v, A(u)=u tensor u-I/6, and

    S(u)_abc=u_a*delta_bc+u_b*delta_ac+u_c*delta_ab,
    T(u)=u tensor u tensor u-S(u)/8.

Direct contraction gives A(u).A(v)=x^2-1/6. Each contraction of
u tensor u tensor u with S(v), or S(u) with v tensor v tensor v,
is3x. The nine contractions in S(u).S(v) comprise three6x terms
and six x terms, totaling24x. Thus

    T(u).T(v)=x^3-(6/8)x+(24/64)x=x^3-(3/8)x.

Substituting independently verifies the exact cubic identity

    (x+1/2)^2(x-1/4)
      =1/16+(3/8)x+(3/4)(x^2-1/6)+(x^3-(3/8)x).

Summing all ordered pairs indeed gives N^2/16 plus the three
nonnegative squared tensor sums with coefficients3/8,3/4,1.
There is no negative coefficient, dimension mismatch, uncontrolled
diagonal term or numerical approximation. The equivalent normalized
Gegenbauer coefficients are1/16,3/8,5/8,5/8 in degrees0,1,2,3.

The diagonal polynomial value is27/16. Every off-diagonal value is
nonpositive for x<1/4, so N^2<=27N and N<=27 when N>0. For N=27,
both finite sums must be equal. Every off-diagonal value must then
vanish separately. The only remaining polynomial root under the
STRICT constraint is x=-1/2. Four distinct code vectors with that
common pair product would have sum square4+12*(-1/2)=-2, impossible.
Hence N<=26. This endpoint argument needs neither uniqueness of
the27-point Schlaefli configuration nor any stronger code theorem.

## 4. Conclusion and verification boundary

The actual six-face has at most26 facets and at least27 by the
published universal theorem. This excludes a finite-volume right-angled
nine-polytope, with no arithmetic or vertex-type restriction. An actual
nine-face excludes every higher dimension. The known P8 construction
is a separate published sharpness example, not needed for the exclusion.

No hidden compactness assumption appears: a finite interior point
exists for every full-dimensional finite-volume polytope, including
entirely ideal ones. Finite volume is used for face inheritance and
the universal six-dimensional minimum; it is not inferred from a
collection of feasible local stars. I found no mathematical gap.

This PASS is a written proof review, not a formalization of real
Lorentz geometry or the universal tensor algebra. Those interfaces
must remain explicit in any conditional Lean module. It is also not
a novelty or published-priority audit. No earlier proof bytes or
current status table were changed by this review.

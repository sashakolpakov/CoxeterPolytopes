# Descending maximum facet height and a dimension-eight obstruction

20 October2026. Complete paper argument for independent review; no
scientific execution, optimizer, numerical certificate, build or Lean
run. No dimension table or existing manuscript was changed. The proof
is separate from the ideal-cusp argument and does not assume its result.

## Statement and the one numerical geometric input

**Proposed theorem.** Every finite-volume right-angled hyperbolic
polytope has dimension at most eight, with finite and ideal vertices
allowed in any combination.

The numerical geometric input is the published statement that a
six-dimensional finite-volume right-angled hyperbolic polytope has
at least27 facets. It is attributed to Lemma1 and Proposition4 of
G. Dufour, *Notes on right-angled Coxeter polyhedra in hyperbolic
spaces*, Geometriae Dedicata147 (2010),277–282,
[doi:10.1007/s10711-009-9454-2](https://doi.org/10.1007/s10711-009-9454-2).
The exact six-dimensional statement was independently checked in
the actual mathematical text of S. Alexandrov,
[*On ideal vertices of right-angled hyperbolic polyhedra*,
Section4](https://arxiv.org/html/2303.09533v1#S4), which cites those
two Dufour results. The original Dufour full text was not retrieved
in this pass; that bibliographic limitation does not change the
precise input being used.

As foundations we use all-wall acuteness and right-angle/finite-volume
inheritance for faces of dimension at least three. All-wall acuteness
is the nonpositive off-diagonal Gram condition in Vinberg,
[*Hyperbolic reflection groups*, printed page42, §2.1](https://doi.org/10.1070/RM1985v040n01ABEH003527).
It applies to actual right-angled polytopes, including nonadjacent
facet walls. The proof below checks actual points against every wall.
It never descends to a polygon, where ideal vertex angles would not
be right angles.

## 1. Farthest-facet projection from a finite point

Use Lorentz signature `(-,+,...,+)` and outward UNIT spacelike
normals. Let `P` be a right-angled polytope of dimension at least
three and `q` a finite unit future point in its relative interior.
For every intrinsic facet normal `m`, define

    t(m) = -<q,m> > 0.

Choose `n` with the largest value `t>0`. Set

    q' = (q+t*n)/sqrt(1+t^2).

The numerator has square `-1-t^2`, is orthogonal to `n`, and has
negative product with the future timelike `q`; hence `q'` is unit
future timelike. For ANY other actual normal `m`,

    <q',m> = (<q,m>+t*<n,m>)/sqrt(1+t^2) < 0,

by strictness at `q` and all-wall acuteness. Thus `q'` is in the
relative interior of the ACTUAL facet `Q=P intersect n-perp`.

Every intrinsic facet of `Q` comes from an actual ridge neighbor
`m_i` of `n`. Right angles give `<m_i,n>=0`, so the restriction
remains unit, and its height at `q'` is

    t'_i = t(m_i)/sqrt(1+t^2) <= t/sqrt(1+t^2).       (1)

If all old heights are at most `A`, monotonicity of
`x/sqrt(1+x^2)` gives a new uniform bound

    A' <= A/sqrt(1+A^2).

For the first projection `A=t` is merely finite, so every new height
is STRICTLY less than1. Repeat intrinsically in the selected facet.
After two projections every height is less than `1/sqrt2`, and
after THREE it is less than `1/sqrt3`. These are strict inequalities
even when several facets tie for the maximum at any step.

Assume a nine-dimensional finite-volume right-angled polytope exists.
Starting from an arbitrary finite interior point and making those
three projections produces an actual finite-volume right-angled
SIX-face `R`, and a relative-interior future unit point `z`, such that
all its intrinsic unit facet normals have heights

    0 < t_i = -<z,m_i> < 1/sqrt3.                    (2)

All intermediate faces have dimension at least six, so there is no
polygon-angle inheritance issue.

## 2. The strict spherical code in Euclidean six-space

Use `z` as the time axis in the Lorentz seven-space of `R` and write

    m_i = (t_i,sqrt(1+t_i^2)*u_i),    u_i in S^5.

For two distinct actual facet normals, acuteness gives

    u_i.u_j <= t_i*t_j/sqrt((1+t_i^2)*(1+t_j^2)) < 1/4.    (3)

The strict last inequality is exactly (2). Equal directions would
instead give positive Lorentz product for two unit normals, so all
directions are distinct. Their number `N` is the actual number of
facets of `R`, not a set of formal candidate walls.

## 3. A self-contained cubic spherical-code bound

We prove that any finite set of unit vectors in Euclidean SIX-space
with distinct-pair products strictly less than `1/4` has at most26
members. No external linear-programming bound or uniqueness of a
spherical configuration is used.

Put

    F(x) = (x+1/2)^2*(x-1/4) = x^3+(3/4)x^2-1/16.

For a unit vector `u`, define ordinary Euclidean tensors

    A(u) = u tensor u - I/6,
    S(u)_(abc) = u_a*delta_bc+u_b*delta_ac+u_c*delta_ab,
    T(u) = u tensor u tensor u - S(u)/8.

Their tensor-product inner products satisfy, for `x=u.v`,

    A(u).A(v) = x^2-1/6,
    T(u).T(v) = x^3-(3/8)x.                        (4)

For the second identity, the two cross contractions with the cubic
tensors are each `3x`; in dimension six `S(u).S(v)=24x`.
Consequently

    F(u.v) = 1/16 + (3/8)(u.v)
             +(3/4) A(u).A(v) + T(u).T(v).

For `u_1,...,u_N`, summing ordered pairs gives the exact identity

    sum_(i,j) F(u_i.u_j)
      = N^2/16 + (3/8)|sum_i u_i|^2
        +(3/4)|sum_i A(u_i)|^2 + |sum_i T(u_i)|^2
      >= N^2/16.                                   (5)

But `F(1)=27/16`, while `F(x)<=0` for every `x<1/4`. Thus

    N^2/16 <= sum_(i,j) F(u_i.u_j) <= 27N/16,

and `N<=27`. If `N=27`, equality holds throughout, so every
off-diagonal polynomial value is zero. The STRICT upper bound in
(3) then makes every off-diagonal product `-1/2`, the only remaining
root. This is impossible already for four such vectors, since

    |sum_(i=1)^4 u_i|^2 = 4+12*(-1/2) = -2.

Hence `N<=26`, as claimed. There is no numerical tolerance or
boundary-case appeal here.

## 4. Contradiction and higher dimensions

The actual six-face in Section1 has at most26 facets by Sections2–3.
The published six-dimensional facet minimum is27. This contradicts
the existence of the nine-polytope.

For any higher-dimensional right-angled finite-volume polytope,
an actual nine-face inherits finite volume and right angles. It is
excluded by the same argument. The proposed ceiling is therefore
eight, without any vertex-type restriction or arithmetic assumption.
Known finite-volume right-angled examples in dimension eight would
make the ceiling sharp; their existence is a separate standard
construction input, not required for the exclusion.

## 5. Scope of the proof packet

The only numerical geometric input is the lower bound27 for facets
in dimension six. The spherical polynomial inequality is proved by
the displayed tensor-square identity, and the strict endpoint is
handled explicitly. No spherical-code classification, Musin theorem,
27-point uniqueness statement, local search, optimizer or simple
polytope count relaxation enters this proof.

This source is submitted for full independent review. It does not
claim that a manuscript has been built, that a checker has run, or
that Lorentz geometry is formalized in Lean. Earlier opposition and
incidence research sources remain unchanged.

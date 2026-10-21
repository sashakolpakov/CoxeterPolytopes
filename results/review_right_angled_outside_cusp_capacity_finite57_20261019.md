# Independent written review of outside-cusp capacity and ideal counts

**PASS, relative to the explicitly inherited inputs.** No correction is
requested to
[the outside-cusp report](right_angled_outside_cusp_capacity_20261019.md),
SHA256 `ffbae765aa2594283e347c75e2a2bf41bfbe5ded217c66c3e53746c641c7ccb0`.
This review covers the complete written argument, including the optional
mixed-face fractions and mixed-vertex inequalities. It ran no scientific
calculation, checker, optimizer, Lean process or build. The arithmetic
below was checked on paper. The stated bounds remain necessary conditions,
with no existence or dimension-exclusion conclusion.

## Projection is against every wall

For a positive definite selected Gram M with nonnegative inverse, the
strict interior vector b has negative coordinates. Thus u=M^(-1)b is
strictly negative; every row has a positive diagonal entry. Subtracting
the selected projection gives a timelike point p orthogonal to every
selected normal. Its product with any other actual facet normal remains
strictly negative, because the correction is a nonnegative combination
of normals and all off-diagonal wall products are nonpositive. Also
(q,p)<0 fixes the future component. This proves actual face feasibility,
rather than only feasibility for a chosen cusp's incident walls.

The two consequences used are valid: n mutually orthogonal actual walls
force an actual finite vertex; and a unit product in(-1,0) would force an
actual nonright codimension-two face. The right-angled product gap0 or<=-1
therefore has the required global meaning. The ordinary acute-wall fact
is explicitly inherited and not inferred from the new model counts.

## The capacity-two lemma

At an outside vertex w sent to infinity, a facet F not containing w must
be a hemisphere. A vertical supporting hyperplane would contain w and
therefore make w a vertex of F as a face of P. Its orientation is the
exterior-of-ball side because w belongs to P. The products with both
vertical walls in direction i give precisely
c_i>=0 and L_i-c_i>=0; each positive quantity is at least R.

Every edge from w lies over one cusp-box vertex xi. Entire ideality is
essential here: its other endpoint is the ideal boundary point xi, rather
than a finite point at positive height. If that endpoint lies in F, the
sphere equation is sum_i(xi_i-c_i)^2=R^2. Every nonzero summand is at least
R^2, so exactly one coordinate differs and has absolute difference R.
This places c on at least n-2 box boundary coordinate planes.

The possible box-corner center is correctly excluded using the full
projection lemma: F and the n-1 vertical facets through that corner are
n mutually orthogonal actual walls. Their actual finite intersection is
strict for every remaining chamber wall. Thus an untested wall cannot
remove the forbidden vertex. The remaining center lies in one box edge's
relative interior, and only its two endpoints can satisfy the sphere
equation. The report explicitly makes this placement conditional on
having at least one neighbor, so it does not make an assertion about an
arbitrary sphere with no such endpoint.

Consequently an outside ideal vertex has at most two graph-neighbors in
F. The outgoing-edge count is2^(n-2) at each vertex of F: it is the
difference between the numbers of vertices of the(n-1)- and(n-2)-box
links. An edge with both endpoints in F is contained in F by convexity.
Thus all these outgoing edges have outside endpoints, giving

    V(P)>=(1+2^(n-3))*V(F).

No finite-vertex version or general mixed-polytope version is asserted.

## The two carrier estimates are independent

The reviewed Figure7/small-spectrum theorem and the positive-kappa
carrier identity give a three-facet R of Q with m=V(R)>=11. Applying the
capacity lemma in dimension4 gives V(Q)>=3m>=33.

Separately, ideal right-angled R has degree4, hence F(R)=m+2 by Euler.
The actual-facet injection gives F(Q)>=1+F(R)+V(R)=2m+3>=25. Ridge
neighbors are distinct because two convex facets have a single face as
their intersection. Cusp-opposite facets are distinct because their
supporting hyperplanes are asymptotic to R at one unique ideal point.
They are disjoint from the ridge-neighbor list, whose hyperplanes are
orthogonal to R. The argument therefore counts actual facets rather
than potential extensions. Applying this injection again gives F(P5)>=59.
The V(Q) estimate alone was not used to infer F(Q).

For ideal P5 the box link gives f1=8v and total three-face vertex incidence
24v. Each polygon has three three-cofaces, so3f2=24v+2g. Euler then yields
g=3v+3f-6. The published minimum24 for all Q and the one actual Q with at
least25 facets give2g>24f. Hence v>3f+2, and integrality gives

    V(P5)>=3F(P5)+3>=180.

The strict step is valid and does not require uniqueness of the24-facet
minimum or the second-volume equality classification. In fact the total
facet excess is the positive multiple6(v-3f-2), as recorded.

## Six-dimensional transport and the exact constants

Applying the facet injection to an ideal five-facet gives F(P6)>=240.
The inherited transport relation a4-F<=136V/99 is universal before its
old final facet-minimum substitution. It may therefore be reused with
2a4>=59F, giving

    (57/2)F<=136V/99,
    V>=(5643/272)F.

The coefficient is57*99/(2*136), with neither a lost coface factor nor
a change to the triangular-corner convention. At F240 the numerator is
1,354,320. The adjacent multiples are272*4979=1,354,288 and
272*4980=1,354,560, so the required integer bound is4980. Equivalently
the unrounded quotient is4979+2/17.

This uses the already reviewed sharp reciprocal maximum103/132 and its
actual-cusp transport theorem; it is not a new maximizer claim or an
optimization. The old proof's intermediate136/99 is the input actually
needed, not its final old4462 numerical conclusion.

## Mixed consequences and scope

I read the inherited density report's NK ranges, strict average bounds,
general facet minima and edge incidences. Replacing the ideal-five minimum
56 by59 gives I5/N5<4/43 in dimensions11,12 and<14/43 in dimensions9,10.
Replacing the ideal-six minimum227 by240 gives I6/N6<15/213=5/71 in
dimensions11,12. All inequalities retain NK strictness.

For the optional vertex inequalities, substitution of c=39/43 into
binom(n,5)*v_f>c*N5, followed by the unchanged strict edge/NK bounds, gives

    (39459/86)*v_f>(139776/473)*v_i, n11,
    (33822/43)*v_f>(39936/43)*v_i, n12.

Reducing gives exactly20669v_f>13312v_i and5637v_f>6656v_i. This is valid
without dividing by v_i, including the v_i=0 case. No contradiction follows.

The geometric dependencies remain the foundational acute-normal/cusp-box
facts, the new Figure7 theorem and external small-three-polyhedron list,
the positive-kappa carrier identity, the published ideal4 facet minimum,
the earlier actual-facet injection, and the sharp transport theorem.
The optional mixed statements also retain the documented NK/Dufour
literature and source-access limits. This review did not independently
reprove those inputs or run their certificates. No full kernel proof or
external human review is claimed, and the ideal/mixed dimension ceilings
remain6 and12.

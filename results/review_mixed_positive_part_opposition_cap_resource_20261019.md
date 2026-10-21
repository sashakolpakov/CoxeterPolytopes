# Independent review of the positive-part opposition obstruction

19 October 2026. **PASS under the stated standard geometric inputs and
the inherited published eight-dimensional facet lower bound. No flaw
found.** This is an independent full written mathematical review, not
an inference from other reviews. I did not read the root's new review
before completing the derivations below. No scientific calculation,
checker import, enumeration, optimization, Lean run or TeX build was
performed. Streamed source hashing was the only executed arithmetic
operation associated with this review.

The main source is
`results/right_angled_mixed_positive_part_opposition_cap_research_20261019.md`,
SHA256 `19ec3d814759cd2d9b9b1df8ef128f1df7d7afee68590aa3310479c954a65c9f`.
I read its entire proof, the full opposition/moment and polygon-coupling
packets, the full weighted-triangle packet, the earlier reciprocal-square
identity packet, and the attachment/three-dimensional inputs needed for
the Q8 and exceptional-face arguments. None was changed.

|Input source under `results/`|SHA256|
|---|---|
|`right_angled_mixed_opposition_moments_research_20261018.md`|`a95079d23556e82bcfbb369abccfb751fe96655ff03d95682b0a7e0c8f6bdd7e`|
|`right_angled_mixed_polygon_three_face_coupling_research_20261018.md`|`f2949232322ec28ba5c2d10a97a30a0a47a519c3dc6061e88dd95c4cf0b31a41`|
|`right_angled_mixed_weighted_triangles_research_20261018.md`|`2e9fac6d7c82493e8e7e1e0172235cc3b910f5f3b05d7ded0d7056bdae217a60`|
|`right_angled_mixed_low_facet_rigidity_research_20261018.md`|`a57761cb4e432a2776018509ac63ce9ce6bab162343f5c759da9b75a07a0b899`|
|`right_angled_mixed_opposition_power_two_research_20261018.md`|`ab7ed2a9e9c04947538589aeee11accb238551d20fca017e58aa9bd892701a2f`|

## 1. Actual oppositions, invariance and the componentwise choice

At a finite vertex q the incident twelve unit outward normals form an
orthonormal basis of q-perpendicular space. For a bad-graph triangle
S={i,j,k}, the three actual opposite walls are
m_ij=q-n_i-n_j. They are orthonormal. The candidate
q_S=2q-n_i-n_j-n_k has square -1 and q.q_S=-2, hence is future.

Its feasibility is an all-wall statement, not just a Gram calculation.
For every other wall r=a*q-sum c_l*n_l, strictness at q gives a>0
and acuteness gives c_l>=0. Its products with the three m walls give
c_i+c_j<=a, c_i+c_k<=a, c_j+c_k<=a, so

    r.q_S=-2a+c_i+c_j+c_k<=-a/2<0.

The selected m walls and unchanged n_l vanish there; the three replaced
n_i are strictly negative. This is an actual simple finite vertex.
The three known ideal endpoints are actual, and the six known walls
enclose their finite-volume bipyramid with q,q_S. Convexity identifies
the whole actual three-face P3. Each P3 has exactly two finite vertices.
Linear independence of the n_i makes q_S distinct for distinct triples.
Consequently the finite opposition graph is undirected, has B edges,
and has degree sum exactly 2B, with no P3 multiplicity lost.

The new incident-normal sum is 3q+sum_all n_i-3sum_S n_i. Hence
3q_S minus that sum is exactly h=3q-sum_all n_i. This proves h is
constant throughout each actual opposition component, without any
complete-corner hypothesis. Also h²=3 and h.q=-3 everywhere.

A spacelike h has a Lorentzian orthogonal complement, containing a
nonzero future null Z. A single such Z is fixed on each entire finite
component, including isolated vertices. For every future unit q,
E(q)=-Z.q>0. There is no normalization obstruction: rescaling Z by
a positive componentwise constant changes g=E^-2 by a constant and
leaves every actual ratio unchanged. Decomposition in the full q,n_i
frame gives sum p_i=3 and sum p_i²=1 for
p_i=-Z.n_i/E(q). No sign of p_i is required.

## 2. Full formal rows and pointwise moments

For x_i=p_i-1/4, sum x_i=0 and sum x_i²=1/4. Every one of the220
formal triples gives a future unit timelike vector 2q-sum_S n_i,
even if it is not feasible in the chamber. Thus all its energy ratios
r=5/4+d are positive. The centered three-subset indicator has squared
norm 9/4, so |d|<=3/4 and 1/2<=r<=2 on the whole formal row.
No absent formal triple is silently declared an actual opposition.

I reexpanded the subset moments using inclusion probabilities
1/4,1/22,1/220. Their coefficients give

    E d²=9/176,
    E d³=-(27/220)sum x_i³,
    E d⁴=3/440-(3/220)sum x_i⁴<=3/440.

In particular the fourth power-sum coefficient is negative, whereas
the remaining square term is (6/55)(sum x_i²)²=3/440. The norm
bound |sum x_i³|<=1/8 gives -E d³<=27/1760. The new polynomial's
cubic coefficient is negative and its quadratic and quartic
coefficients positive, so every moment bound is used in the correct
direction. These are pointwise statements at every finite vertex,
not expectations with an unproved vertex distribution.

## 3. Full-domain positive-part majorant

I checked both required inequalities P>=0 and P>=r^-2-1/2;
neither alone suffices for the later row deletion.

For d<=0 all five terms of P have nonnegative signs. For d>=0,
the stated lower bound P''>=12/5 follows from
36/625-2d+60d²>0. The exact values at1/4 are
P=1627/1200000 and P'=-10243/150000. Strong convexity yields
the uniform positive lower bound67/240000. This remains valid
for the whole nonnegative half-interval, not just near1/4.

The reciprocal-square identity has the displayed exact remainder;
4/r+(5/4)/r² decreases for r>0 and is at most13 on r>=1/2.
The difference P-P0 is exactly
3/250+(9/200)d+d³+(106/75)d⁴, nonnegative for d>=0.

For d=-t I independently multiplied by (5/4-t)². The resulting
degree-six coefficients are exactly those of N in the source,
including the constant3/160 and the negative linear coefficient
-321/3200. I also reexpanded all17 Bernstein coefficients in the
three displayed derivative rows using the supplied binomial formula.
Their coarse strict rational bounds have the stated signs. Thus
the argument controls complete intervals; it is not a sample test.

The central interval deserves a sign check. N'''' is convex and
negative at both endpoints1/5,3/10, so it is negative throughout.
Therefore N''' decreases and is bounded below by its positive value
at3/10. Thus N'' increases and is bounded below by its value at1/5,
which exceeds4/5. Using the weaker bound N''>=1/2 with the exact
N(1/4)=679/1200000 and N'(1/4)=3061/600000 gives the asserted
strict positive lower bound. The first Bernstein row propagates
positivity backward from N(1/5)>1/1000; the second propagates it
forward to1/2. The third proves concavity on[1/2,3/4], where
positive endpoints, including N(3/4)=57/64000, suffice. All four
intervals cover the entire required negative-d domain.

No root-location or approximate fit is a premise. The root separately
reports a fresh exact coefficient generation/replay, but this review
does not claim to have executed it or substitute its status for the
written polynomial proof.

## 4. Actual-subset deletion and the unweighted global sum

For each actual opposition the positive-part inequality gives

    g(q_S)/g(q)-1 <= P(d_S)-1/2.

Summing over actual triples and then inserting absent formal triples
is valid because P is nonnegative. It gives

    sum_actual(g(q_S)/g(q)-1)<=sum_all220 P-t(q)/2.

Summing this over ordinary finite vertices requires no selected
maximum-energy vertex and no weighted-to-unweighted conversion.
Each actual edge contributes r^-2+r²-2>=0, since both endpoints
use the same fixed componentwise Z. Different components cause no
problem. The moment estimate is exactly

    E_all P <=1420561/5500000.

The coefficient sum checks directly after clearing that denominator:
836000+345600+26217+212744=1420561. With sum t(q)=2B this gives

    2B<=(1420561/12500)f.

Here f counts finite vertices only and B the same actual canonical
P3 three-faces used below. No polygon-triangle count is substituted.

## 5. Independent reconstruction of the opposing lower bound

I checked the weighted-polygon/three-face chain rather than importing
its final number alone. Its external inputs are the standard finite
simplex and ideal box links, finite-volume/right-angle inheritance
for faces of dimension at least three, simple truncation and
Dehn--Sommerville, the published facet lower bound82 in dimension8,
the compact right-angled ceiling4, and the elementary at-least11
facets for a ten-polytope. Two-faces instead have finite angles pi/2
and ideal angles zero; they are not claimed to inherit pi/2 at cusps.

The truncation convention b0=f+2048c is essential. The displayed
DS coefficient identities give

    66f+11264*(2123/1280)c-(22/5)a2
      =(-22a7+88a8-198a9+308a10-308a11)/5<0.

Indeed5a7>=82a8 and3a9>=11a10 make the two grouped terms
strictly negative, as does the remaining facet term. The finite and
cusp coefficients cancel correctly. Thus a2>15f+4246c.

The three-dimensional Euler/area inequalities give F>=6,
F+2s>=12 and nu=F+s>=9, with equality nu=9 only for P3.
The six-facet equality proof supplies the actual P3 geometry rather
than an unattested classification. Summing nu over three-faces gives
B>=10a3-10a2-28160c. The DS relation
a3=5a2-55f-14080c and the weighted row yield

    B>50f+880c.

For the stronger coupling, a Q8 face has F=8,s=2 and eight finite
vertices. Its polygon excess vanishes. A bad triangle would share
its ideal/ideal edge with another bad triangle; the explicit normal
attachment lemma forces P3, a contradiction. Hence all its eight
polygons are(3,1) quadrilaterals. A finite Q8 vertex can have at
most one incident ideal-ending edge. No cusp-cusp edge is possible,
and the two degree-four cusps supply eight finite/ideal edges, so
every finite Q8 vertex has exactly one such edge.

At an ambient finite vertex with k ideal-ending directions, Q8
three-faces therefore occupy at most k*C(12-k,2) triples. P3 faces
occupy t(q)<=C(k,3). The elementary finite inequality

    14*(220-C(k,3))-27*k*C(12-k,2)>=0

holds for k=0,...,12. Independently, its successive integer values
are3080,1595,650,150,0,105,370,700,1000,1175,1130,770,0.
It implies X(q)>=(13/27)(220-t(q)) for exceptional three-face
incidences X(q), where exceptional means neither P3 nor Q8.

For each exceptional three-face, p=4F-2nu-4 satisfies

    p<=22(nu-10)+4(8-F).

The difference is24nu-184-8F. At nu=10 only F=7 remains after
excluding P3/Q8 and it gives equality. At nu=11 one has F<=10;
at nu>=12 use F<=nu. This proves the full budget. Q8 contributes
zero to both right-hand excesses and is not erroneously counted
among exceptional incidences.

Define E=sum_nonP3(nu-10) and R=sum_nonP3(8-F), the latter signed.
Actual coface counts give

    E=B-10a3+10a2+28160c,
    R=8a3-10a2-2B.

Summing the budget and the local lower bound gives
22E+4R>=(13/27)(220f-2B). Substitution of the DS identity yields
exactly

    101B-5130a2+69080f+22049280c>=0.

No low-facet rigidity in dimensions4/5, graph symmetrization,
ideal-only carrier spectrum or higher-dimensional Garland claim is
needed in this chain. The old triangle attachment lemma and its
all-wall normal-product gap are the only Q8 rigidity inputs.

Combining the strict weighted row gives
101B>7870f-267300c. Its positive multiple880 plus the positive
multiple267300 of B>50f+880c cancel the cusp terms and give

    356180B>20290600f,
    2B>(184460/1619)f.

The inequality is strict, and c still counts ambient cusps while f
counts only finite vertices. Thus the two sides compare exactly the
same B and f as the opposition upper bound.

## 6. Verdict and exact conclusion

The cross-product difference is the positive integer5861741:

    184460*12500-1420561*1619=5861741.

The new upper bound and the reconstructed strict lower bound are
therefore incompatible whenever f>0. This proves nonexistence of
a mixed right-angled finite-volume hyperbolic12-polytope under the
explicit standard/published inputs above. Combined with the prior
general ceiling12 it yields mixed ceiling11. Dimension11 and any
new example are outside the proof.

The result remains a written geometric proof with small separately
replayable coefficient obligations. This review finds no flaw; it
does not assert external refereeing, worldwide priority, end-to-end
Lean verification or a new lower-dimensional construction. Source
hashes identify the exact proof reviewed and do not themselves
prove its geometric interfaces.

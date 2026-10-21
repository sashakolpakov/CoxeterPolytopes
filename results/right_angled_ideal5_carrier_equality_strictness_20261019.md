# Extra triangle credit at the smallest ideal-five carrier

19 October 2026. **New written proof, awaiting independent review.** No
scientific script, enumeration, optimizer, Lean process or build ran.
Earlier source-bound proofs are unchanged. The new local assertion below
implies an entirely ideal right-angled five-polytope has at least208
vertices. No dimension exclusion or existence result is claimed.

The inputs are the reviewed exact cusp-defect identities in
`right_angled_ideal4_cusp_defect_research_20261019.md` and the reviewed
five-dimensional propagation proof in
`right_angled_ideal5_facet_cusp_propagation_20261019.md`. Their actual-wall,
fixed-null-scale, small-polygon-corner, and incidence arguments are used
with their full hypotheses. This proof does not assume a garland
classification for an arbitrary ideal four-polytope.

## 1. The equality carrier and at least five distinguished triangles

Let Q be a four-facet of an entirely ideal right-angled five-polytope P.
Assume Q has36 vertices and contains a three-facet R with11 vertices.
Write q_R for R-edges with nontriangular outgoing two-face in Q, N0_R
for Q-vertices without an R-neighbor, and G_R for good cusps with the
INTERNAL R normalization tau_u=r+g_u. Here r is R's unit normal in Q
and g_u its actual opposite facet in Q at u.

The reviewed identities and inequalities are

    36=3*11+N0_R+q_R,
    G_R>=2+max(0,4-q_R),    2N0_R>=G_R.

Their integer consequences force q_R=0, N0_R=3, G_R=6. For q_R=1,2,
the defect is at least4; for q_R>=3 it is at least5. At q_R=0 the
lower bound6 and upper bound2N0_R=6 give the exact good-cusp count.

Every edge of R consequently has -(tau_u,tau_v)=2. Every triangular
face of R has all three internal corner coefficients equal1. There
are at least8 such triangles. A nongood R-cusp has at most two small
polygon corners, so at most two triangle incidences. The five nongood
vertices supply at most10 triangle incidences, whereas eight triangles
supply24. Thus at least14 triangle incidences occur at good R-cusps.
At least

    K>=ceil(14/3)=5

DISTINCT triangles of R contain an R-good cusp. Call them distinguished.
This bound does not require a classification of11-vertex polyhedra.

## 2. Two normalizations and the small-scale case

Write n for Q's ambient normal in P and h_u for its opposite ambient
facet at u. The Q-in-P normalization is t_u=n+h_u. The six facet normals
of Q through u form three opposite pairs whose sums are

    alpha0_u*t_u, alpha1_u*t_u, alpha2_u*t_u.

At u in R one pair is r,g_u, so

    tau_u=r+g_u=alpha0_u*t_u,
    t_u=lambda_u*tau_u,    lambda_u=1/alpha0_u>0.

The two INTERNAL R coefficients are alpha1/alpha0 and alpha2/alpha0.
All supporting normals remain unit after restriction because the
relevant actual walls are orthogonal. There is no extra normalization
factor. For an R-edge uv, q_R=0 and the ambient edge inequality give

    -(t_u,t_v)=2*lambda_u*lambda_v>=2,
    lambda_u*lambda_v>=1.

For every triangle of R its internal coefficient is1 at every corner.
The corresponding ambient polygon coefficient is therefore
alpha0=1/lambda at that corner. If it is at least2, TWO of Q's three
coefficients coincide at a value at least2. In particular the possibly
small lambda<=1/2 case is not assumed away; it will be charged explicitly.

Let G count Q-cusps good in the ambient normalization (all three
coefficients<2). Let B count Q-cusps with at least TWO coefficients>=2.
A vertex in B lies in at most four triangles of R, since R has degree4.
Consequently at most4B distinguished triangles fail to be all-small.

Call an edge of Q exceptional when its outgoing two-face in P is not a
triangle; let q be their number. A triangle of Q with no exceptional
edge has all three ambient corner coefficients equal1. If a distinguished
triangle is of this kind, at an R-good corner it has alpha0=1 and both
other ambient coefficients<2. That corner is Q-good. Each Q-good vertex
lies in at most four R triangles, so at most4G distinguished triangles
are unexceptional.

Therefore at least

    J>=K-4B-4G

distinguished triangles are BOTH all-small and exceptional. A negative
right-hand side is allowed: the resulting lower bound remains valid.
These triangles provide extra credit distinct from the usual credit
for unexceptional triangles. No assumption lambda>1/2 is used.

## 3. Exact credit accounting

Put m=V(Q)=36 and f=F(Q), let r2 be Q's number of polygonal two-faces,
T its number of triangles, and S its total small polygon corners.
The reviewed identities give r2=3m+f, T>=4f. At most3q triangles are
exceptional. Every polygon has at least two small corners, every
unexceptional triangle has three, and the J exceptional all-small
triangles supply separate further units. Hence

    S>=2r2+(T-3q)+J
      >=6m+6f-3q+K-4B-4G.

At a Q-cusp each of the three pair coefficients occurs four times.
A good cusp contributes12 small corners; an ordinary nongood one at
most8; a cusp in B at most4. Therefore

    S<=8m+4G-4B.

Combining and canceling the B losses proves the new local inequality

    8G+3q>=6f-2m+K>=6f-2m+5.                  (A)

The exact ambient partition and reviewed good-opposite-facet capacity
are

    V(P)=5m+N0+q,      9G<=2N0.

Since8G<=16N0/9<=3N0, (A) gives

    3V(P)>=13m+6f+5                           (B)

for this equality carrier. This is actual extra triangle credit, not
a substitution into a previously relaxed scalar bound.

## 4. Universal numerical consequences and scope

The reviewed positive carrier theorem gives Q containing R with at
least11 vertices, with m>=36 and f>=25. If m>=37, the old universal
propagation already gives3V(P)>=13*37+6*25=631, hence V(P)>=211.
If m=36, the internal +3 inequality forces V(R)=11, so (B) applies
and gives3V(P)>=623. Thus every ideal right-angled five-polytope has

    V(P)>=208.

The existing facet minimum62 is unchanged by this argument. Ambient
facet injection in dimension6 then gives F6>=1+62+208=271. The reviewed
unchanged transport relation68V6>=1485F6 yields

    68*5918=402424<1485*271=402435<=402492=68*5919,
    V6>=5919.

These consequences await independent review. The geometric foundations
and the11-vertex carrier bridge remain the existing explicitly stated
inputs; no new small-polyhedron classification is used. The ideal and
mixed dimension ceilings remain6 and12, and no Lean or executable
verification is asserted for this new credit argument.

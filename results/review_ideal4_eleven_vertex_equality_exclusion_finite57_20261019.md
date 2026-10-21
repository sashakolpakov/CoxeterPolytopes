# Independent review of the eleven-vertex carrier equality exclusion

19 October 2026. Complete mathematical read of
`right_angled_ideal4_eleven_vertex_equality_exclusion_20261019.md`, bound
to SHA256
`b9af9b497c69da69a25faced7067b854f0e3b8332f89da57daa12a5a968c2121`,
and of its full262-line inherited cusp-defect proof. No scientific job,
matrix calculation program, census, numerical test, Lean or build ran.
No pinned source was edited.

**Verdict: PASS**, conditional on the already reviewed actual-wall,
capacity-two and cusp-defect inputs stated in the report. The new written
argument proves that an ideal RA4 polytope containing an eleven-vertex
three-facet has at least37 vertices. No omitted equality branch, cusp
normalization change or unjustified adjacency was found.

## Equality saturation and the six octahedra

At hypothetical equality V(Q)=36 with V(R)=11, the exact identity
V(Q)=3V(R)+N0+q gives N0+q=3. The good-cusp inequalities force q=0,
N0=3 and G=6. Counting small polygon corners then gives T<=8, while
the degree-four Euler identity gives T>=8. Equality forces precisely
eight triangular and five quadrilateral faces and exactly two small
corners on every quadrilateral. Every nongood cusp has exactly one
small opposite-pair coefficient, repeated twice; both-large and both-small
nongood alternatives are excluded by saturation.

Each good h_u has at least one extra N0 vertex. Six such required
incidences saturate the upper2N0=6, so every good h_u has exactly six
vertices. The octahedral graph follows from its simple degree-four graph
on six vertices; Euler forces eight triangles. The report's hemisphere
calculation is a direct proof of square cusp shape: sending the unique
opposite vertex to0 inside the four vertical walls gives x_+²=−y_-y_+
and x_-²=−y_-y_+, then the analogous horizontal equations. They force
the centered rectangle to be a square.

The four relevant unit normals are the SAME ambient normals for R and
h_u: at their common cusp they are orthogonal to both n and h_u. Thus
restriction changes neither their norms nor their sums. Both sums use
the fixed actual null vector t_u=n+h_u. Square shape really gives
alpha_u=beta_u, without an independent cusp rescaling. This is the
normalization point that a combinatorial octahedron alone would not prove.

## Quadrilateral normalization and nongood orthogonality

Since q=0, every consecutive null product is−2. Their spinor lifts have
unit consecutive determinants, so an SL2 transformation may put the
first two at(1,0),(0,1) without changing any actual t scale. The other
two are(−1,t),(−2/t,1), t>0. Direct determinants give corner
coefficients2/t,t,2/t,t and the opposite product8/B². This verifies
the fixed-scale formula used in Section4.

If h_u is orthogonal to h_v, then t_u·t_v=−1. The vertices cannot be
adjacent. At a nongood u write A<2<=B. If both facets from the A pair
missed v, their products with t_v would each be at most−1, forcing
A>=2. Thus one such facet contains v. At its polygonal face the corner
coefficient at u is the OTHER coefficient B. It cannot be a triangle,
whose coefficient is1. It must be a quadrilateral with v opposite u.
The normalized product gives B=2sqrt2 and adjacent coefficients1/sqrt2.

At each of the two other vertices, goodness makes the second coefficient
equal1/sqrt2, while nongoodness makes it at least2. Neither coefficient
is1, so no triangular face is incident there. Each vertex consequently
belongs to four of the five quadrilateral faces. Those two four-element
sets share at least three facets. Any three facets through an ideal
rectangular cusp include an opposite pair with that cusp as their unique
common ideal point. The second distinct vertex cannot be on all three.
This excludes orthogonality for nongood u as well as good u; it does
not silently apply the earlier good-cusp lemma outside its hypothesis.

## Actual-edge step for all eleven opposite facets

Every h_u contains a vertex w beyond u and its four neighbors. It lies
outside R, since R and h_u are tangent only at u. If w had an R-neighbor
v, the outgoing edge wv would lie on h_v. The two distinct incident
h_u,h_v are nonorthogonal and hence form an opposite pair at w's box
link. Sending w to infinity makes them parallel vertical planes. R's
hemisphere is tangent to both at the actual points u,v. A box-corner
tangency v forces the other tangency u to have the same other coordinates,
so u is also a box corner.

This last observation really gives an actual ambient edge uw: the three
common vertical walls cut out its vertical geodesic. Every other ambient
wall through w is a containing vertical halfspace; every wall missing w
is a hemisphere whose outside halfspace contains u, and moving vertically
up from u only increases the squared distance from its center. Thus no
untracked wall cuts that geodesic. This contradicts w's selection as a
nonneighbor of u. The report's extension from good to all u uses exactly
the new pairwise-nonorthogonality statement, with no residual goodness
assumption.

All11 facets h_u therefore require an N0 incidence. At an N0 cusp, any
three distinct incident box facets include an orthogonal pair, so at
most two of these pairwise nonorthogonal h_u occur there. The contradiction
11<=2N0=6 is valid. The same counting argument applies over7<=V(R)<=11
at hypothetical +3 equality; no arbitrary RA4 decomposition is asserted.

## Numerical and dependency scope

For the reviewed positive-kappa carrier, V(R)>=11. The equality exclusion
gives V(Q)>=37 at11; at least12 gives V(Q)>=39 by the older +3 theorem.
The independent facet injection still gives F(Q)>=25. Consequently
F5>=63 and the reviewed propagation3V5>=13V(Q)+6F(Q) gives3V5>=631,
so V5>=211. The next injection gives F6>=275. Combining
2f4>=63F6 with99f4<=99F6+136V6 gives272V6>=6039F6 and V6>=6106.
The integer products in the report were checked by hand.

These numerical corollaries additionally inherit the large-facet carrier,
facet injections, five-dimensional propagation and sharp six-dimensional
transport theorem. The new equality proof uses no eleven-vertex
classification, new enumeration, universal garland claim or formalized
geometry. This independent review adds no execution and proves no
ideal-dimension6 or mixed-dimension12 exclusion.

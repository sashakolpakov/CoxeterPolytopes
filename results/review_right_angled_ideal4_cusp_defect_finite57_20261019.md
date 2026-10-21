# Independent review of the ideal-four cusp defect

19 October 2026. **Written mathematical review: PASS, relative to the stated
foundational inputs.** No proof script, matrix calculation, optimizer, Lean
process or build was executed. The reviewed report is
`results/right_angled_ideal4_cusp_defect_research_20261019.md`, SHA256
`a50bb5ecc0947870c969098445618de8504c9ade119bd0e4706ccd010769f336`.
Its source was preserved.

The new theorem is universal within its entirely ideal finite-volume
right-angled scope:

    V(Q4) >= 3 V(R3)+3

for every three-facet R of Q. I found no omitted geometry or counting case.
The withdrawn +6 claim is neither needed nor reinstated. The resulting
62/189 and252/5504 necessary counts follow from the explicit older carrier,
facet-injection, strict Euler and136/99 transport interfaces; they do not
give new dimension ceilings.

## 1. The null normalization and exact outgoing count

At u in R, n and h_u are the actual opposite pair of unit normals in its
three-box cusp. Their product is -1, so t_u=n+h_u is null and orthogonal
to n. It is FUTURE null: its product with any strict interior future point
is negative. Thus these t_u have a common, geometrically fixed scale.
They are not independently rescaled vertex representatives.

Distinct u give distinct h_u, since the two hyperplanes of n,h_u are
asymptotic at that unique point. For an edge uv of R the other two ambient
walls a,b are mutually orthogonal and orthogonal to both h_u,h_v. The
alternative h_u.h_v=0 would give four orthogonal actual walls, hence an
actual finite vertex by the all-wall projection lemma. Entire ideality
excludes this. Therefore -t_u.t_v=1-h_u.h_v>=2.

The d_uv=2 direction really gives a triangle: h_u+h_v is future null and
feasible against every actual wall by the positive-sum inequality. It is
not in R, since its product with n is -2. The three halfplanes n,h_u,h_v
in the a,b perpendicular plane bound that ideal triangle; all other walls
contain its three ideal vertices and therefore their convex hull.

Conversely, the reviewed outside-cusp box argument puts the two R-neighbors
of an outside vertex at the endpoints of one box edge. The two common
vertical walls cut out the same actual outgoing two-face at their R edge.
The three ideal vertices already fill the three-halfplane triangle, so a
fourth side cannot intervene. This is a bijection, not just an injection.
Consequently N2=2m-q, N1=2q, and

    V(Q)=3m+N0+q.

The four outgoing edges at each R vertex come from the four vertices of
the cusp box's opposite face. No outgoing edge can end again in R: an edge
with both endpoints in a face lies in that face. Thus no edge was lost in
N1+2N2=4m, and no outside vertex has more than two R-neighbors by the
inherited capacity theorem.

## 2. The polygon lemma keeps the fixed scale

The opposite pairs among the four actual normals of R through u have
positive sums proportional to the fixed t_u. In each polygonal face the
intrinsic side normals are the same unit normals: right angles make them
orthogonal to the polygon's own normal. Hence their sums have exactly one
of the two coefficients alpha_u,beta_u, occurring on opposite polygons.
The paired-corner count uses this fixed normalization correctly.

For T(a,b)=(a^2+b^2,a^2-b^2,2ab), direct expansion gives

    T(v).T(w)=-2 det(v,w)^2.

The signs of cyclic spinor representatives can be chosen on one lift of
RP1, with v0=-vk and v_(k+1)=-v1. Successive determinants, including the
last lifted gap, are positive. Since adjacent t_i products are at most-2,
all ell_i are at least1. The note's local side-normal calculation gives
the stated gamma_i at the SAME t_i scale and the vector recurrence.

The recurrence matrix is symmetric, with the two POSITIVE wrap entries;
the endpoint minus signs are what produce them. Its two coordinate
columns are independent kernel vectors. If all but at most one gamma
were at least2, deleting that position would leave a path with diagonal
at least2 and off-diagonal magnitude at most1. Its quadratic form dominates
the positive path expression in the report, regardless of the sign of
the one remaining wrap bond. The resulting positive-definite principal
matrix contradicts nullity at least two. This also covers k=3.

If all three outgoing edge-faces of a triangular polygon are triangles,
all ell_i=1. Pairing its three positive side-normal sums gives
gamma_i gamma_j=1. Positivity then gives gamma_i=1 at all three corners.
Thus the extra small corner used in the triangle count is justified;
it is not an independent octahedron-normalization assumption.

The degree-four Euler identity gives at least eight triangular polygons
of R. An exceptional edge can spoil at most two. Combining their extra
small corners with the polygon lemma and the paired coefficients gives

    S>=2m+4+max(0,8-2q),
    S<=2m+2G,
    G>=2+max(0,4-q).

These counts do not assume a bounded polygon size or use a census.

## 3. Good opposite walls force N0 vertices: the delicate adjacency step

Suppose a good h_u were orthogonal to h_v for v in R. Then t_u.t_v=-1,
so u,v are nonadjacent by the edge inequality. Distinct nonadjacent
vertices of a convex three-polyhedron share at most one facet. One of
the two opposite pairs through u is therefore wholly nonincident to v.
For each normal f in that pair,

    f.t_v=f.h_v<=-1.

Indeed a zero product with t_v would make f contain the actual vertex v;
the nonzero product then uses the actual zero-or-at-most-minus-one wall
gap. Pairing its positive sum with t_v forces alpha_u>=2 or beta_u>=2,
contrary to goodness. This proves the global nonorthogonality assertion.

The ideal right-angled facet h_u has at least six vertices. Four of its
edges at u are precisely the four edges leaving R. It therefore has
another vertex w, outside R and not adjacent to u. If w had an R-neighbor
v, its outgoing edge wv would lie in h_v. At w the two incident walls
h_u,h_v have product0 or-1. The former was just excluded, so they are an
opposite pair in its rectangular three-box link.

Send w to infinity. These walls become parallel coordinate planes
x1=0 and x1=L. The hemisphere of n is tangent to each, at the KNOWN
actual points u,v. Tangency to both forces radius L/2 and center first
coordinate L/2. The tangency points share the other two coordinates.
Since v is joined to w, it is a box corner; the other tangency point u is
therefore also a box corner.

Here the asserted u-w adjacency does include every other actual wall.
The three orthogonal vertical walls at that box corner contain both
actual ideal vertices u and w. Every chamber halfspace contains their
geodesic by convexity. The intersection of those three wall hyperplanes
is exactly that geodesic, so it gives a one-dimensional face and hence an
actual edge u-w. No additional hemisphere can remove an interior point
while preserving both endpoints. This contradicts w's choice. Thus w
has NO R-neighbor and belongs to N0.

Finally any three distinct faces of a rectangular three-box contain an
orthogonal pair. Since no good h_u is orthogonal to another h_v, at most
two good opposite facets can contain a fixed N0 vertex. Choosing one
witness for each good h_u and double counting proves 2N0>=G. The witnesses
need not be distinct; precisely this capacity handles coincidences.

## 4. Defect, applications and formal scope

Together the exact partition and the two inequalities give

    V(Q)-3m=N0+q
      >=q+ceil((2+max(0,4-q))/2)>=3.

At q=0,1,2,3 the displayed bounds are3,4,4,5; at q>=4 they are at least5.
The theorem therefore follows without the withdrawn stronger claim.

For the independently reviewed m>=11 positive carrier, the TWO separate
conclusions are V(Q)>=36 and F(Q)>=25. Actual facet injection gives
F(P5)>=62. The old strict Euler argument still uses only the published
24-facet minimum plus one carrier with25 facets, and gives
V(P5)>=3F(P5)+3>=189. The next-dimensional injection gives F(P6)>=252.
Since2a4>=62F and a4-F<=136V/99, one obtains

    V>=(1485/68)F,
    1485*252=374220,
    68*5503=374204<374220<=68*5504.

Hence V(P6)>=5504. The stated strict fractions2/23 and1/15 are consistent
with the unchanged NK substitutions at the new facet minima; they remain
conditional on those previously stated external incidence/average inputs.

I also read `formal/CoxeterProofs/RightAngledCuspCapacity.lean`. At this
review boundary it proves the EARLIER conditional integer consequences
59/180 and240/4980, taking the outside-capacity, carrier, injections, Euler
and transport relations as explicit premises. It does not yet formalize
this new +3 geometry or its stronger counts. The coordinator's reported
kernel/audit PASS for that module must not be described as formal
verification of the new theorem reviewed here.

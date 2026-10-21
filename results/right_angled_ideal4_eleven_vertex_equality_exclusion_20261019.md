# Excluding the eleven-vertex equality carrier in an ideal four-polytope

19 October 2026. **New written theorem candidate, awaiting independent
review.** No scientific script, enumeration, optimizer, Lean process or
build ran. All earlier pinned sources are unchanged. This proof uses
actual facets and a small equality argument, not a classification of an
arbitrary ideal four-polytope or of eleven-vertex three-polyhedra.

**Local assertion.** If Q is an entirely ideal finite-volume right-angled
four-polytope and a three-facet R has eleven vertices, then

    V(Q)>=37.                                        (A)

The reviewed universal bound was V(Q)>=3V(R)+3=36. The proof below excludes
its equality case. The same argument works for7<=V(R)<=11, but only the
eleven-vertex instance is needed for the numerical application.

## 1. Explicit inherited interfaces and all equality consequences

Use the reviewed proof
`right_angled_ideal4_cusp_defect_research_20261019.md`, SHA256
`a50bb5ecc0947870c969098445618de8504c9ade119bd0e4706ccd010769f336`.
Let n be R's unit normal, h_u the actual opposite facet of Q at u in R,
and t_u=n+h_u. Let q count R-edges with nontriangular outgoing two-face,
N0 vertices outside R without an R-neighbor, G good R-cusps, and S small
polygon corners of R. Good means both fixed-scale coefficients alpha,beta
in the two opposite-pair sums are<2.

Suppose V(Q)=36. The exact identities and inequalities in that proof give

    N0+q=3,
    G>=2+max(0,4-q),    2N0>=G.

They force q=0, N0=3, G=6. Thus every R-edge has -(t_u,t_v)=2, and all
triangular faces of R have all three coefficients equal1. If T is the
number of triangular faces, the polygon lemma and paired-corner bound give

    2(11+2)+T<=S<=2*11+2*6=34.

Hence T<=8. The degree-four Euler identity forces T>=8, with equality
only when all other faces are quadrilaterals. Therefore R has precisely
eight triangles and five quadrilaterals, S=34, and every nongood cusp
has exactly ONE coefficient<2 and the other>=2. Every quadrilateral has
exactly two small corners.

At each good u, all vertices of h_u other than u and its four neighbors
belong to N0, by the reviewed actual tangency-to-adjacency argument. Each
h_u is an ideal right-angled three-polyhedron with at least six vertices.
The total number of incidences of GOOD h_u with N0 is therefore at least
G=6 and at most2N0=6. Each good h_u has exactly one such additional
vertex and consequently EXACTLY six vertices. No statement about additional
vertices of a nongood h_u has been inferred from this tight count.

## 2. Six vertices force a square cusp: elementary geometry

An ideal right-angled three-polyhedron with six vertices has degree4.
Its simple graph is the complete graph on six vertices minus a perfect
matching, the octahedral graph; Euler and face degrees give eight
triangular faces. Every cusp is square, as follows directly without
appealing to a classification or general rigidity theorem.

Send a vertex u to infinity and its unique nonneighbor w to0 in upper
half-space. The four vertical walls through u bound a rectangle
[x_-,x_+] times[y_-,y_+]. Since w is on none of those four faces,
x_-<0<x_+, y_-<0<y_+. The four neighbors of u are its four corners.
The remaining facet across x=x_+ contains w and the two corners on
that side. Right angles put its hemisphere center on x=x_+; passing
through both corners gives center y=(y_-+y_+)/2 and radius
(y_+-y_-)/2. Passing through w=0 then gives

    x_+^2=-y_-*y_+.

The opposite side gives x_-^2=-y_-*y_+, so x_-=-x_+. The horizontal
pair similarly gives y_-=-y_+ and y_+=x_+. Thus the rectangle is square.

At u in R, the four normals through u of the three-facet h_u are the
same four actual facet normals as those of R through u: the two extra
ambient normals n,h_u are an opposite cusp pair. The sum coefficients
of those four normals use the same null vector t_u. Square cusp shape
therefore implies

    alpha_u=beta_u  at every GOOD vertex of R.       (B)

This is the only octahedral metric fact needed below. It is not a claim
that Q is an octahedral-facet garland.

## 3. A normalized ideal quadrilateral identity

For any quadrilateral face of R, all four consecutive t-products are-2
because q=0. In the spinor normalization T(v).T(w)=-2det(v,w)^2 choose
cyclic lifts, preserving the endpoint signs, as

    v1=(1,0), v2=(0,1), v3=(-1,t), v4=(-2/t,1), t>0.

All four consecutive determinants, including det(v4,-v1), are1. The
reviewed side-normal-sum formula gives cyclic corner coefficients

    2/t, t, 2/t, t.

Consequently a corner coefficient B has adjacent coefficients2/B,
and if u,v are the two opposite corners having coefficient B, then

    -(t_u,t_v)=8/B^2.                              (C)

There is no independently chosen rescaling of t_u here; the unit adjacent
determinants fix exactly the scale used in the cusp-defect proof.

## 4. No two actual opposite facets h_u,h_v are orthogonal

This is already known when u is good. Suppose it occurs for a nongood u.
Then -(t_u,t_v)=1, so u,v cannot be adjacent, whose product would be2.
Write A<2<=B for u's two coefficients. The opposite pair of R-facet
normals with sum A*t_u cannot both miss v: otherwise their two products
with t_v would each be at most-1 by the actual-wall gap, forcing A>=2.
Thus v lies in a face of this pair through u.

Distinct nonadjacent vertices share at most one face of a convex
three-polyhedron. Equivalently here, any three of the four cusp facets
through u include an opposite pair, whose hyperplanes meet at infinity
only at u. In the selected face, the intrinsic corner coefficient at u
is the OTHER pair coefficient B. That face cannot be a triangle, whose
coefficient is1. It is a quadrilateral, and v, not adjacent to u, is its
opposite vertex. Formula(C) gives B=2sqrt(2). Its two other vertices
therefore have a small coefficient1/sqrt(2).

At either such other vertex, if it is good, (B) makes its two coefficients
equal1/sqrt(2). If it is nongood, its other coefficient is>=2. In both
cases neither coefficient equals1. No triangular face can meet it.
Thus each of these TWO DISTINCT vertices is incident to four quadrilateral
faces. There are only five quadrilateral faces in R, so their two
four-element incident-facet sets overlap in at least three facets.
This is impossible: among any three facets through an ideal rectangular
cusp, an opposite pair has that cusp as its unique common ideal point.
Two distinct vertices cannot belong to the same three facets.

This contradiction proves h_u.h_v!=0 also when u is nongood. Therefore
ALL eleven opposite facets are pairwise nonorthogonal.

## 5. The universal witness argument now contradicts N0=3

Each h_u has at least six vertices, whereas u and its four neighbors
account for only five. Choose a further vertex w of h_u not adjacent to u.
It is outside R. If w had an R-neighbor v, its outgoing edge wv would
lie in h_v. At w, the incident h_u,h_v are nonorthogonal and hence form
an opposite pair of its rectangular three-box link.

The SAME actual-wall argument as in the reviewed good-cusp lemma now
applies: send w to infinity. R's hemisphere is tangent to those two
parallel vertical planes at the known actual points u,v. Since v is a
box corner, the other tangency point u is also a box corner. Convexity
against every other wall preserves the geodesic uw; the three common
vertical walls cut it out as an actual edge. This contradicts the choice
of w. Thus w belongs to N0.

Goodness was only used in that inherited argument to exclude orthogonal
opposite facets. Section4 has now established the exclusion for ALL h_u,
so every one of the eleven facets supplies an N0 incidence. At an N0
vertex, three distinct incident facets include an orthogonal pair in its
rectangular three-box link. At most two h_u can therefore occur there.
We obtain11<=2N0=6, a contradiction. Equality V(Q)=36 is impossible,
proving(A).

## 6. Conditional numerical consequences

The reviewed Figure7/carrier bridge gives a positive-kappa four-facet Q
in every ideal right-angled five-polytope and a three-facet R of Q with
M>=11 vertices. If M=11, (A) gives V(Q)>=37. If M>=12, the existing
universal +3 bound already gives V(Q)>=39. The separate actual-facet
injection still gives F(Q)>=2M+3>=25. Hence

    F(P5)>=1+37+25=63.

The unchanged reviewed propagation inequality
3V(P5)>=13V(Q)+6F(Q) gives3V(P5)>=631, so

    V(P5)>=211.

Each five-facet of an ideal six-polytope satisfies these bounds; ambient
facet injection gives F6>=1+63+211=275. The unchanged sharp transport
99f4<=99F6+136V6, together with2f4>=63F6, gives

    272V6>=6039F6,
    272*6105=1660560<6039*275=1660725<=1660832=272*6106,
    V6>=6106.

These prospective63/211 and275/6106 counts need full independent review
of the local equality exclusion before promotion. They do not prove an
ideal dimension6 or mixed dimension12 exclusion. No published11-vertex
classification, numerical census, Lean proof or global garland assertion
was used in the new local argument.

# Independent review of the ideal-five facet cusp-defect extension

19 October 2026. **Written mathematical review: PASS within the inputs
listed below.** This reviews the coordinator's proposed extension of the
ideal-four cusp-defect argument to a four-facet of an entirely ideal
right-angled five-polytope. No calculation script, proof process, optimizer,
Lean job or build was run. No existing evidence was edited. I subsequently
read the complete final coordinator source
`results/right_angled_ideal5_facet_cusp_propagation_20261019.md`, SHA256
`e2850f435989c08d2de3f8932596e7c8fd6067fb94263676c44a5855e9b107a1`.
The derivation below reviews that source, including its adopted elementary
18-vertex dependency reduction.

The reviewed geometric statement is the following. If P is an entirely
ideal finite-volume right-angled five-polytope, and Q is any four-facet
with m vertices and f facets, then

    3 V(P) >= 13m+6f.                                  (A)

In particular, the already reviewed carrier with m>=36 and f>=25 gives
V(P)>=206. Together with the existing F(P)>=62 and the unchanged ideal-six
injection and reciprocal transport, this gives the necessary counts
F(P6)>=269 and V(P6)>=5875. These are not dimension exclusions.

## 1. Inputs and what is new

I read the complete reports

* `results/right_angled_outside_cusp_capacity_20261019.md`;
* `results/right_angled_ideal4_cusp_defect_research_20261019.md`;
* `results/review_right_angled_ideal4_cusp_defect_finite57_20261019.md`.

The inherited geometric inputs are acute actual-wall products, the
all-wall positive-definite projection lemma, rectangular cusp links,
face inheritance, and the proved outside-neighbor capacity of two. The
spinor polygon lemma is the one proved in the ideal-four cusp-defect
report, with its fixed-scale normalization and positive wrap signs.

The new argument is not an assumption that all ideal polytopes decompose
into 24-cells. It applies to arbitrary ideal right-angled P and Q.
Furthermore, the proof below does not need a published minimum of
24 vertices for ideal right-angled four-polytopes: the already reviewed
capacity bound supplies the weaker minimum 18, which is ample here.

## 2. The exact outside partition generalizes

Let n be Q's actual unit normal. For each cusp u of Q let h_u be the
unique actual wall opposite n in the four-box cusp of P. The null
vector t_u=n+h_u is future-pointing, is orthogonal to n, and has a
geometrically fixed scale. Different u give different h_u.

There are eight edges leaving Q at each of its vertices. If N_i counts
outside vertices with i neighbors in Q, capacity gives i in {0,1,2}
and

    N_1+2N_2=8m.

An edge uv of Q has three other supporting actual walls. Both h_u and
h_v are orthogonal to those three walls. If they were orthogonal to
one another, the five walls would give an actual finite vertex of P.
Thus, by entire ideality and the actual-wall gap,

    -(t_u,t_v)=1-(h_u,h_v)>=2.

The outgoing two-face of uv is triangular exactly when (h_u,h_v)=-1.
For this direction the third vertex is h_u+h_v, a feasible future null
vector, and the three known side inequalities already enclose the
ideal triangle. Conversely, the outside-capacity box proof puts the
two Q-neighbors of an N_2 vertex at the two endpoints of one box edge.
The three common vertical walls cut out precisely that outgoing
two-face. This is the same bijection as in dimension four, with three
common walls in place of two.

The graph of Q has degree eight, hence 4m edges. If q is the number
of edges with nontriangular outgoing two-face, then

    q=4m-N_2=N_1/2,
    V(P)=5m+N_0+q.                                    (1)

No outside vertex or outgoing edge was omitted.

## 3. Polygon corners and coefficient multiplicities

At a cusp u of Q the six incident intrinsic walls form three opposite
pairs, with sums alpha_j t_u, j=1,2,3. A two-face through u is obtained
by choosing one wall from each of two pairs. Its two intrinsic side
normals are the remaining opposite pair. Thus each alpha_j occurs
exactly four times among the twelve polygon corners of Q at u.

Call u good if all three alpha_j are strictly below 2, and let G count
good u. If S counts polygon corners with coefficient below 2, then

    S<=8m+4G.                                         (2)

Indeed a nongood cusp contributes at most eight small corners, whereas
a good cusp contributes twelve. This is the necessary higher-dimensional
replacement of the two-pair count; using two occurrences here would
be incorrect.

The spinor polygon lemma applies to each two-face of Q because every
adjacent pair of its fixed-scale t_u vectors has product at most -2.
Each polygon has at least two small corners. A triangular polygon
whose three outgoing two-faces are triangular has all coefficients
equal to 1, so contributes an additional small corner.

Every three-facet of Q is an ideal right-angled three-polyhedron and
has at least eight triangular faces. Every two-face of Q has two
three-cofaces. Hence Q has at least 4f triangular two-faces. Each of
the q exceptional edges spoils at most three triangles, because an
edge of a four-polytope has exactly three incident two-faces at its
relative interior. Euler and degree eight give f_2(Q)=3m+f. Consequently

    S>=2f_2(Q)+4f-3q=6m+6f-3q.

The right side is valid even if 4f-3q is negative: the number of extra
small corners is nonnegative and is at least that weaker lower bound.
Together with (2), this gives

    4G+3q>=6f-2m.                                    (3)

## 4. The good-wall argument does not lose a shared-facet case

Suppose h_u is good and orthogonal to h_v. Then (t_u,t_v)=-1, so u,v
are not adjacent in Q by the edge inequality above. They share at
most two facets of Q. To see this without a general simple-polytope
assumption at the cusp, inspect any three common facets among the
six walls at u. If they include an opposite pair, their hyperplanes
are asymptotic at the unique point u and cannot also contain v. If
not, the three are mutually orthogonal. Their intersection with Q
has dimension one by the all-wall projection lemma; since it contains
both u and v, it is their actual edge. This also contradicts
nonadjacency.

There are three opposite pairs at u, and v is incident to at most two
of their six walls. One pair therefore avoids v altogether. If its
sum is alpha_j t_u, each of its wall products with t_v equals its
product with h_v, is nonzero, and is at most -1. Pairing the sum with
t_v forces alpha_j>=2, contrary to goodness. Therefore a good h_u
is orthogonal to no other h_v.

## 5. A good wall has enough N_0 vertices

The facet h_u is an ideal right-angled four-polytope. The old capacity
bound applied to any three-facet R of h_u gives

    V(h_u)>=3V(R)>=18,

since an ideal right-angled three-polyhedron has at least six vertices.
At u the facet h_u has exactly eight neighbors, the eight edges
leaving Q. It has at least nine further vertices w not adjacent to u.
All of them lie outside Q, because Q and h_u meet only at u.

Each such w belongs to N_0. Otherwise let v in Q be a neighbor of w.
The edge wv lies on h_v. At w the walls h_u,h_v have product either
0 or -1; goodness excludes 0, so they are opposite sides of its
four-box cusp. Send w to infinity. The hemisphere of n is tangent
to these parallel walls at the known actual points u and v. Tangency
to both puts its center halfway between them and makes the two
tangency points have the same other coordinates. Since v is the box
corner of an edge from w, u is the opposite box corner and also
defines an edge from w. Convexity preserves that entire edge because
both endpoints are actual vertices. For the dimension count, that geodesic
lies on h_u and the three other common vertical walls: these four independent
walls cut out an edge in dimension five. This spells out the implicit h_u
wall in the coordinator's phrase "three common vertical facets." This
contradicts the choice of w.

Finally an N_0 vertex is incident to at most two good h_u: any three
distinct sides of a rectangular four-box include an orthogonal pair,
and such a pair of good h_u was excluded. Counting the at least nine
choices at each good wall gives

    9G<=2N_0, hence 4G<=8N_0/9<=3N_0.                 (4)

The coordinator's stronger intermediate 15G<=2N_0 also follows if the
published 24-vertex minimum is explicitly supplied. It is unnecessary
for (A), and this review does not rely on it.

## 6. Conclusion and numerical scope

Equations (3),(4) imply 3(N_0+q)>=6f-2m. Substituting (1) proves

    3V(P)>=15m+6f-2m=13m+6f.

For the already reviewed large carrier, m>=36 and f>=25. Thus
13m+6f>=618 and V(P5)>=206. The existing independent facet injection
still gives F(P5)>=62. Applying it in dimension six gives

    F(P6)>=1+62+206=269.

The unchanged reciprocal/incidence argument with at least 62 facets
on every ideal five-facet is V(P6)>=(1485/68)F(P6). Therefore

    1485*269=399465,
    68*5874=399432<399465<=399500=68*5875,

so V(P6)>=5875. These products were checked by hand, not by a new
executed certificate. Their further inputs are the existing large
carrier theorem, facet injection, and 136/99 ideal-six transport.

This review establishes no new mixed or ideal dimension ceiling,
no existence result, no published-priority claim, and no new Lean
coverage. Its universal mathematical claim is precisely (A) under
the foundational geometric inputs listed in Section 1.

## Source clarification addendum

The coordinator subsequently made the supporting-wall count explicit:
the geodesic u-w lies on h_u and three other common vertical facets,
whose four independent walls make it an actual edge in dimension five.
I reread that amended paragraph and checked the source binding. The
final source SHA256 is
`69453207b7ad935c741ed8f31f31bc61f6a8dc12e96675eb2f063aafb5a929bb`.
The earlier `e2850f...` binding above records the full initially reviewed
version. This clarification changes no mathematical claim or review
conclusion; the final source remains a written PASS. No job was executed.

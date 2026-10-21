# A stronger vertex propagation inequality across an ideal four-facet

19 October2026. New root written proof, **awaiting independent review**.
No numerical search or proof process was used in deriving this argument.
The previously reviewed outside-capacity and actual-cusp/acute facts are
explicit inputs. No published four-dimensional classification is needed
for inequality(1). No dimension exclusion or example is claimed.

**Proposed theorem.** Let P be an entirely ideal finite-volume right-angled
five-polytope and Q any of its four-facets. Put m=V(Q), f=F(Q). Then

    3V(P) >= 13m+6f.                                  (1)

Combining this with the reviewed positive-carrier bounds m>=36 and f>=25
would give V(P)>=206. The existing facet bound F(P)>=62 would be unchanged.
The resulting ideal6 implications would be F>=269 and V>=5875. These
applications remain pending independent review of the new geometry below.

## 1. Normalization and the exact outside-vertex defect

Give every ambient wall a unit outward normal. Write n for Q's normal,
and h_u for the actual facet opposite Q in the rectangular cusp link
at a vertex u of Q. Set t_u=n+h_u. This is a future null vector in
the Lorentz model of Q, with its scale fixed by the two unit normals.
As in the [reviewed cusp-defect proof](right_angled_ideal4_cusp_defect_research_20261019.md),
all distinct actual normals have product0 or at most-1.

Every edge uv of Q belongs to three other ambient facets whose normals
are mutually orthogonal and orthogonal to both h_u,h_v. If h_u,h_v
were orthogonal as well, the five actual orthogonal normals would force
a finite vertex of P by acute projection against every other wall.
Thus h_u.h_v<=-1 and

    -(t_u,t_v)=1-(h_u,h_v)>=2.                        (2)

Let N_i count vertices outside Q with i graph-neighbors on Q. The
reviewed outside-cusp capacity gives i in{0,1,2}. There are8 outgoing
edges at every vertex of Q: the respective vertex links have16 and8
vertices. Hence N1+2N2=8m.

An outside vertex with two Q-neighbors gives exactly the triangle
outgoing from the Q-edge joining them. Indeed, send that vertex to
infinity. The capacity proof puts the hemisphere center of Q on one
edge of its four-box link, with the two ideal contact points the
endpoints of that box edge. Three common vertical walls cut out the
actual two-face, and Q cuts out its opposite side. Conversely a
triangular outgoing two-face has one such outside vertex. The edge
determines the triangle uniquely, so this is a bijection. Since Q
has4m edges, the number q of edges with nontriangular outgoing
two-face satisfies

    q=4m-N2=N1/2,
    V(P)=5m+N0+q.                                    (3)

## 2. Small polygon corners and good cusps

The six facet normals of Q through u occur in three opposite pairs.
Their sums are alpha_i*t_u, i=1,2,3, with alpha_i>0. Call u good
when all three alpha_i are<2, and let G count good vertices.

At a corner of any polygonal two-face of Q, the two intrinsic side
normals are one of these opposite pairs. The supporting normals of
the polygon are orthogonal to them, so the intrinsic normals have
the same unit normalization. The [spinor polygon lemma](right_angled_ideal4_cusp_defect_research_20261019.md)
applies by(2): every polygon has at least two corner coefficients<2.
A triangle all of whose outgoing edge-faces are triangular has all
three coefficients equal1, and contributes one further small corner.

Write r for the number of polygonal two-faces of Q and T for their
triangular subset. The ideal4 incidence identities give r=3m+f and
total ridge defect3(m-f). Every nontriangle has defect at least1,
so T>=r-3(m-f)=4f. Each exceptional edge of Q belongs to exactly
three two-faces, since its link is a spherical triangle. At least
T-3q triangles therefore have no exceptional edge. If S counts
small polygon corners, we obtain

    S >= 2r+max(0,T-3q) >= 6m+6f-3q.                 (4)

At a vertex of Q its three-box link has12 edges. A fixed one of
the three opposite-pair coefficients occurs four times among these
polygon corners. A good vertex contributes12 small corners; every
other vertex contributes at most8. Consequently

    S<=8m+4G,
    4G+3q>=6f-2m.                                   (5)

## 3. Good opposite facets force many vertices with no Q-neighbor

A good h_u is orthogonal to no other h_v. For otherwise t_u.t_v=-1,
so(2) shows that u,v are not adjacent in Q. Distinct vertices of Q
cannot share an opposite pair of facets through u: those hyperplanes
are asymptotic at u alone. If they shared three facets, those normals
would therefore belong to three different coordinate directions of
u's three-box link and be mutually orthogonal. Their actual common
face is one-dimensional, by acute projection and independence. Since
its closure contains both u,v, they would be its edge endpoints,
contrary to nonadjacency. Thus they share at most two of the six
facets through u.

At least one opposite pair through u consequently misses v entirely.
Each of those two facet normals has strictly negative product with
t_v, hence product at most-1 by the actual-wall product gap. If their
sum is alpha*t_u, its product with t_v is-alpha<=-2. This forces
alpha>=2, contradicting goodness.

The ambient facet h_u is an entirely ideal right-angled four-polytope
and has at least18 vertices by the earlier outside-capacity theorem:
any three-facet has at least6 vertices by degree4 and Euler, and its
vertex count propagates by a factor3. This elementary lower bound
suffices; the sharper published ideal4 minimum is not needed. The
eight edges of h_u at u are precisely the outgoing Q-edges. It thus
contains at least9 vertices
w other than u and these eight neighbors. None lies on Q, because
the hyperplanes of n and h_u meet only at u.

Each such w has no Q-neighbor. If it had a neighbor v on Q, then
the outgoing edge wv would lie on h_v, so h_u,h_v meet at w. They
are not orthogonal by goodness and must therefore be opposite faces
of its four-box link. Send w to infinity. Q's hemisphere is tangent
to these two parallel vertical planes at u and v. These tangency
points have the same other coordinates. Since v is a box corner,
u is the opposite box corner. The geodesic u-w is feasible by
convexity and lies on h_u and the three other common vertical facets.
These four independent walls cut out an actual edge. This contradicts
the choice of w as a nonneighbor of u.

At a vertex with no Q-neighbor, at most two good h_u can occur.
Any three distinct faces of a rectangular four-box include an
orthogonal pair, which the preceding goodness argument forbids.
Counting the incidences of good h_u with these additional vertices gives

    9G<=2N0.                                        (6)

All walls in this argument are actual ambient facets, and the
edge assertions use feasibility against every wall.

## 4. Propagation and exact prospective constants

From(6), 4G<=(8/9)N0<=3N0. Combining with(5) yields
3(N0+q)>=6f-2m. Equation(3) then proves(1):

    3V(P)=15m+3(N0+q)>=13m+6f.

The existing reviewed carrier has m>=36 and f>=25, so(1) gives
3V(P)>=618, or V(P)>=206. Ambient facet injection still gives
F(P)>=1+36+25=62.

In an entirely ideal six-polytope, a five-facet and its opposite-cusp
and ridge-neighbor facets then give F>=1+62+206=269. The unchanged
universal transport relation99f4<=99F+136V and2f4>=62F imply
68V>=1485F. At F269,1485*269=399465 and

    68*5874=399432 <399465 <=399500=68*5875.

Thus V>=5875 would follow. No new optimizer or census is needed for
these substitutions. The new universal assertion requiring review is(1)
and its actual-wall/count interfaces, not a numerical fit. The ideal
and mixed dimension ceilings remain6 and12.

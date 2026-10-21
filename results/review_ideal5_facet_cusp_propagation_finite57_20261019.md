# Independent review of ideal-five facet cusp propagation

19 October 2026. **Written mathematical review: PASS.** The final reviewed
source is `results/right_angled_ideal5_facet_cusp_propagation_20261019.md`,
SHA256 `69453207b7ad935c741ed8f31f31bc61f6a8dc12e96675eb2f063aafb5a929bb`.
No scientific computation, checker, Lean process or build ran in this
review. No source was edited by this reviewer.

The new universal inequality is

    3 V(P5) >= 13 V(Q4)+6 F(Q4)

for EVERY four-facet Q of an entirely ideal finite-volume right-angled
five-polytope P. The proof checks, using the explicitly inherited
outside-capacity and acute-wall foundations. Combining it with the
reviewed carrier m>=36,f>=25 gives V(P5)>=206; F(P5)>=62 is unchanged.
The next-dimensional conclusions F(P6)>=269,V(P6)>=5875 follow. These
are necessary size bounds, not dimension exclusions or existence claims.

## 1. Actual outgoing two-faces and the exact partition

The normalization t_u=n+h_u is fixed by actual opposite unit walls.
It is future null and lies in the Lorentz space of Q. For an edge uv
of Q, its three other ambient supporting normals are mutually
orthogonal and orthogonal to h_u,h_v. An orthogonal h_u,h_v pair
would make five actual orthogonal walls, giving an actual finite
ambient vertex by the all-wall projection lemma. This contradicts
entire ideality and establishes -t_u.t_v>=2.

The links at u in P and Q have16 and8 vertices, so exactly8 ambient
edges leave Q. These are distinct edges with outside endpoints. An
edge with both endpoints in Q lies in Q, because it is a face.
The inherited outside-vertex capacity therefore gives

    N1+2N2=8m.

The N2 correspondence is valid also in this dimension. At an outside
vertex sent to infinity, the two Q-neighbors are the contact corners
of one edge of a FOUR-box. Their three common vertical wall normals
cut out an actual two-face. The Q hemisphere gives its side between
the two ideal neighbors, and the two remaining coordinate walls give
its other sides. The three known ideal vertices fill that triangle;
convexity against ALL other chamber walls prevents another side
from truncating it. Conversely an outgoing triangular two-face has
one outside vertex with exactly these two Q-neighbors. The Q edge
determines that face uniquely. Thus N2 is exactly the number of
triangular outgoing edge-faces, not merely a lower bound for it.

Since Q has4m edges, q=4m-N2=N1/2 and

    V(P)=5m+N0+q

are exact. No dimension-four adjacency degree was incorrectly reused.

## 2. Polygon coefficients and their multiplicities

At a cusp of Q its SIX facet normals form three opposite pairs. A
polygonal two-face through that cusp fixes two coordinate directions
of the three-box link and leaves the third direction free. Its two
intrinsic side normals are exactly the remaining opposite pair, at
their original unit normalization: the supporting normals are
orthogonal to them. Thus its coefficient gamma at t_u is exactly
one of alpha1,alpha2,alpha3. It is not a projection with an additional
scale factor.

Each of these three coefficients occurs FOUR times, on the four
parallel edges of the three-box. A good cusp consequently supplies12
small polygon corners; a nongood cusp supplies at most8. This gives
S<=8m+4G with the correct factor4.

The spinor polygon lemma applies at the fixed t_u scale because
every Q edge satisfies -t_u.t_v>=2. A triangular polygon all of whose
outgoing edge-faces are triangular has all three side determinants1
and all three gamma coefficients1, just as in the independently
reviewed ideal-four defect proof. Hence it supplies its extra small
corner correctly.

For completeness the ideal-four counts used here follow directly:
Q has4m edges, and Euler gives r=3m+f polygonal two-faces. Its cusp
links have12 edges, so their total corner incidence is12m. Therefore
the total polygon defect is12m-3r=3(m-f). Nontriangles have defect at
least1, giving T>=4f. Each Q edge lies in exactly THREE two-faces:
its transverse spherical link is a triangle. At most3q triangles
are spoiled. Consequently

    S>=2r+max(0,T-3q)>=6m+6f-3q,
    4G+3q>=6f-2m.

The last inequality remains valid if T-3q is negative; its use of the
maximum has the correct direction.

## 3. Good opposite facets and the9G incidence count

If h_u were orthogonal to another h_v, the edge inequality makes u,v
nonadjacent in Q. They cannot share an opposite pair among the six
facets through u: those two hyperplanes are asymptotic only at u.
If they shared three facets, those facets would belong to three
different link directions and be mutually orthogonal. Their actual
common face in Q is one-dimensional, by independence and acute
projection. Its closure contains both u,v, making them edge endpoints.
This contradicts nonadjacency. Thus at most two facets are shared.

One of the three opposite pairs through u misses v completely. Its
two actual normals each have strictly negative product with t_v;
otherwise their hyperplane would contain v. Since they are
orthogonal to n, these products equal their products with h_v and
are at most-1 by the global actual-wall gap. Their positive sum then
forces the associated alpha at least2, contrary to goodness. The
claimed good-h_u nonorthogonality holds globally, not just among
facets already sharing a chosen cusp.

Every h_u has at least18 vertices without using the published
ideal-four minimum: any ideal right-angled three-facet has at least6
vertices by degree4/Euler, and the previously reviewed outside
capacity multiplies this by3. The eight neighbors of u inside h_u
are precisely the outgoing Q-neighbors. Thus there are at least9
additional vertices w not adjacent to u; none is on Q.

If such w had a Q-neighbor v, its outgoing edge would lie on h_v.
At w, h_u and h_v must then be an opposite pair, since goodness
excludes their being orthogonal. With w at infinity, Q's hemisphere
is tangent to the two parallel walls at the already actual points
u,v. The tangency points have the same remaining three coordinates.
As v is a box corner, u is the opposite box corner. The geodesic uw
is contained in every chamber halfspace by convexity. It lies on
h_u AND the three other common vertical walls. These FOUR independent
walls define an actual edge in P5, contradicting the choice of w.
The final source explicitly includes this fourth wall; three alone
would only define a two-face. This was the sole exposition correction
requested during review, and is included in the bound source hash.

Every good h_u therefore contains at least9 N0 vertices. Among any
three distinct facets in a four-box link, two are orthogonal.
The good-h_u exclusion permits at most two good opposite facets at
each N0 vertex. This proves9G<=2N0 even when witnesses coincide or
there are more than9 additional vertices. It counts incidences, not
an unjustified injection of all witnesses.

## 4. The claimed constants and scope

From9G<=2N0 one gets4G<=(8/9)N0<=3N0. Combining with the small-corner
inequality gives3(N0+q)>=6f-2m. Substitution in the exact partition
yields

    3V(P)=15m+3(N0+q)>=13m+6f.

At m>=36,f>=25, its right side is618, giving V(P5)>=206. The separate
facet injection still gives F(P5)>=62. Injection in the next dimension
then gives F(P6)>=1+62+206=269. The unchanged sharp transport and
two-coface incidence give68V>=1485F. Finally

    1485*269=399465,
    68*5874=399432<399465<=399500=68*5875.

Thus V(P6)>=5875 follows. This arithmetic uses the new geometry
universally for every five-facet; it does not extrapolate from a
sample carrier or a finite census.

The core inequality needs no published ideal-four classification
or24-facet minimum. Its numerical application retains the separately
reviewed positive-carrier11-vertex bridge, ideal-four+3 defect and
facet-injection inputs, and the earlier universal136/99 transport.
No hyperbolic geometry has been promoted to formal proof here; no
new Lean coverage or current source compilation is claimed.

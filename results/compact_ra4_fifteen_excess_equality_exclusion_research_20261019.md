# Compact RA4: excluding the fifteen-excess endpoint

19 October2026. **New written candidate, awaiting independent review.**
No scientific calculation, graph census, optimizer, Lean or build ran.
The reviewed W>=15 packet and primary TeX remain unchanged. In particular
the primary source remains pinned at SHA256
`9fe28180a9097d8ee99c1898dc64e0716abf55c9c9f8cb950a6db59576d47738`.

The proposed successor is W=0 or W>=16. This note proves only that
successor claim, conditional on the explicit inherited W>=15 proof and
its metric inputs. It does not establish a compact volume minimum or gap.

## 1. Reviewed interfaces and the exact endpoint reduction

Use the definitions d(F),W,G,H from
`compact_ra4_positive_facet_link_gap_fifteen_20261019.md`, SHA256
`8201b2f5cc7da5988b6d164b258da8e14981cf94c6b5e980a5a87e9deb70dcd7`.
Graph G records nonpentagonal ridges, with weight k-5; H records FULL
adjacency among positive-excess facets. The reviewed source proves:

* sum_F d(F)=2W; every positive facet has d>=2;
* G has no degree-one vertex;
* a degree-two vertex, with equal edge weights w and k=w+5, forces
  W>=5k-15 via its Loebell facet and2k lateral positive neighbors;
* if G has minimum degree three, H has minimum degree five;
* no graph H on at most nine vertices can have minimum degree five,
  no K5 and no induced C4;
* every induced2K2-free subcubic graph on nine vertices has an
  independent set of five;
* a zero-excess facet is a regular right-angled dodecahedron;
* orthogonal actual normals give actual faces, and every compact face
  inherits right angles and the acute-normal property.

The conditional third bullet about H depends only on the minimum degree
of G, not on a separate strict W<15 assumption. Its proof was the
four-neighbor eligibility argument. We reuse it precisely in that form.
Let L be the edge length of the regular right-angled pentagon and
phi=cosh(L)=(1+sqrt(5))/2.

Assume W=15. If G has a degree-two vertex, the reviewed bound forces
k=6 and equality throughout its distinct-neighbor sum. The original
facet, its two base neighbors and its twelve lateral neighbors are
exactly fifteen positive facets, every one of excess2. All their G
vertices therefore have degree two with unit weights, and every facet
is Loebell L6. Each has fourteen positive neighbors, by the same
nonregular-lateral-pentagon argument. Thus H=K15, impossible. This
closes the entire degree-two case without assuming any outside facet
was already zero.

Otherwise G has minimum degree three. Let N,E be its number of vertices
and unweighted edges. Then3N<=2E<=2W=30, so N<=10. Minimum degree
five in H and the reviewed small-graph argument give N>=10. Hence

    N=10, E=15, every G-degree is3, every edge weight is1.       (1)

In particular every positive facet has excess3 and exactly three
nonpentagonal faces, all hexagons. Its other twelve faces are pentagons.

Some H-vertex has degree five. Otherwise H has minimum degree six,
its complement on ten vertices is subcubic and induced2K2-free, and
any nine-vertex induced subgraph has an independent five by the reviewed
lemma. That would give a forbidden K5 in H.

Fix such a degree-five positive facet F. Its five positive neighbor
faces have an adjacency graph J on five vertices, with exactly three
marked vertices corresponding to F's three hexagons. It has no K4 and
no induced C4. All other faces of F border zero-excess facets and are
regular pentagons. We now exclude every possible J, using actual polygon
geometry rather than the weaker abstract ten-vertex graph conditions.

## 2. Three-edge closure and two polygon metric facts

**Three-edge lemma.** A compact right-angled polygon with three
consecutive edges of length L is the regular right-angled pentagon.

Match their inward-oriented chain to three consecutive edges of that
pentagon. Fixed edge lengths and inward right-angle turns determine
the chain. At its two ends the next supporting lines are uniquely
perpendicular through their endpoints. These five lines are therefore
the five lines of the regular pentagon. The two outer ones are
orthogonal. They are ACTUAL wall lines of the original polygon, with
the same inward orientation as in the matched chain. Project a strict
timelike interior point orthogonally off their two outward unit normals.
Acuteness makes every other inequality strict and the result future
timelike. Their intersection is consequently an actual vertex of the
polygon, closing its five-edge boundary. Thus an apparent longer
polygon cannot truncate this vertex by some additional wall. All five
edges are those of the regular pentagon.

**Two-adjacent-edge pentagon lemma.** A compact right-angled pentagon
with two consecutive edges of length L is the regular pentagon. Match
these two edges, fixing the three vertices. The predecessor and successor
support lines are the prescribed perpendiculars at their outer endpoints.
As in the regular model they are ultraparallel. The fifth side must
be their unique common perpendicular. This determines the whole pentagon
and proves the claim without a three-dimensional rigidity theorem.
In particular a right-angled pentagon with three edges of length L is
regular: among three edges of a pentagon two are adjacent.

**Four-fixed-edge hexagon lemma.** A compact right-angled hexagon with
four specified edges of length L and two other side positions must
have those other positions opposite, and both other lengths equal2L.

If the two positions are adjacent or separated by only one edge, the
four specified edges contain three consecutive ones, forbidden by the
three-edge lemma. They are therefore opposite. Write its cyclic side
lengths as a,L,L,d,L,L and put A=cosh(a),D=cosh(d),p=cosh(L)=phi.
The right-angled hexagon cosine relations give

    A=(p^2+D)/(p^2-1),       D=(p^2+A)/(p^2-1).

For clarity these formulas can be derived directly in the Lorentz
plane. If sides b,f neighbor side a and d is opposite a, let x,y be
the unit future vertices at its ends and n_a its outward unit normal.
Then x.y=-cosh(a). The outward normals of sides c,e after traversing
b,f away from a are

    n_c=sinh(b)*x-cosh(b)*n_a,
    n_e=sinh(f)*y-cosh(f)*n_a.

Their product is cosh(b)cosh(f)-cosh(a)sinh(b)sinh(f).
Since d is their perpendicular connecting edge, it equals-cosh(d).
This gives the first cosine relation, and the opposite one is identical.
As p^2-1=p, subtraction forces A=D, and

    A=D=p^2/(p-1)=p^3=cosh(2L).

Both lengths are positive, so a=d=2L.

## 3. Exact marked-vertex restrictions in J

For a hexagonal neighbor A of F, its degree in J is the number of
positive third-facet edge carriers of the hexagon R=F intersect A.
Common neighbors and their adjacencies are actual edges and vertices
of R, by orthogonal projection of three and four normals respectively.
All other edge carriers have zero excess and give length L.
Consequently a marked vertex must satisfy:

1. its J-degree is at least two;
2. at degree two its neighbors are not adjacent, and its corresponding
   two edge carriers are opposite in the hexagon, with shared lengths2L;
3. at degree three its neighbor graph is NOT a three-vertex path.

The first two statements follow from the hexagon lemma (also from the
reviewed four-edge eligibility test for the adjacency clause). For the
third, a path of three carriers gives three consecutive positive edges;
the other three are consecutive zero-carrier edges of length L. The
three-edge lemma excludes this hexagon. A triangle of carriers is already
forbidden by the absence of K4 in J. No restriction beyond these is
assumed at degree four.

## 4. Every five-neighbor graph is covered

If J has a universal vertex v, let B be its graph on the other four.
It is trianglefree (otherwise with v there is K4) and has no induced
C4. Thus B is a forest. A vertex of B-degree zero has J-degree one
and is unmarked. A vertex of B-degree one has two adjacent J-neighbors
and is unmarked. A vertex of B-degree two has three J-neighbors
forming a path centered at v, and is unmarked by restriction3.
Only a B-degree-three vertex can still be marked, and a forest on
four vertices has at most one. Together with v this gives at most
two marks, whereas three are required.

Hence J has maximum degree at most three. Suppose it contains a
triangle T, with outside vertices d,e. If de is absent, an outside
vertex with at most one T-contact has degree at most one; with two
contacts its degree-two neighbors are adjacent; with three contacts
it would give K4. Thus neither outside vertex can be marked. A marked
triangle vertex must have a unique outside contact whose outside end
contacts NO other triangle vertex: otherwise its three neighbors form
a path. There are only two outside vertices, so at most two marks.

If de is present, let D,E be their respective contact sets in T.
They are disjoint: a common triangle neighbor would have degree four,
contrary to our current maximum-degree-three case. They must also
be nested. Indeed x in D\E and y in E\D would give the induced
four-cycle x-d-e-y-x. Thus one is empty; say E is empty. Set D has
size at most two. At size zero there are no marks. At size one the
only possible marks are d and its unique triangle neighbor. At size
two, d is the only possible mark: its two triangle neighbors have
three-neighbor paths and the remaining triangle vertex has two
adjacent neighbors. Vertex e always has degree one. Again there
are fewer than three marks. This exhausts every graph containing
a triangle, not merely one selected triangle's marks.

Finally suppose J is trianglefree. With no induced C4, it is either
C5 or a forest: any cycle on at most five vertices must have length
five, and a chord would create a triangle. A forest on five vertices
with at least three marked vertices must have at least three vertices
of degree at least two. The edge sum at most eight forces it to be
the path P5, with its three internal vertices marked. Thus the only
remaining cases are P5 and C5.

In P5, every marked hexagon has degree two. Its two positive carrier
edges have length2L by the hexagon lemma. An endpoint of the path
represents an unmarked positive pentagonal face; it has one positive
neighbor, a marked hexagon, and four zero-carrier edges of length L.
It is the regular pentagon by the reviewed four-edge lemma, contradicting
its shared edge of length2L.

In C5, again every marked hexagon has two positive carrier edges of
length2L. There are three marked and two unmarked vertices, so a
pentagonal face is adjacent to a marked hexagon. It has two positive
carriers and three zero-carrier edges of length L. Among these three
there are two consecutive edges, forcing this pentagon to be regular
by the two-adjacent-edge lemma. Its shared edge of length2L is again
impossible. This closes both remaining graph cases.

## 5. Conclusion and limit

Every possible five-neighbor facet F has been excluded. But the W=15
reduction forced such a facet. Therefore W=15 is impossible, giving

    W=0 or W>=16.

The weaker abstract ten-vertex graph conditions were not declared
impossible. For example C5[K2] satisfies minimum degree five, no K5
and no induced C4; the present proof uses marked hexagon geometry
and its actual pentagonal consumers to exclude the required cubic G.

No change to the primary manuscript or bounds table is made by this
candidate. Even if independently accepted, the local46/46/96 facet
cutoffs stay the same, and f+W>=120, compact volume minimality and
the double-volume gap remain unproved. All inputs are those of the
reviewed W>=15 theorem plus the explicit elementary polygon arguments
above; there is no new external classification or execution claim.

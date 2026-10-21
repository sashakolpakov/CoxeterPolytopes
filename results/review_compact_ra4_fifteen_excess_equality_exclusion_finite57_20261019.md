# Independent written review: excluding compact RA4 excess fifteen

19 October 2026. **PASS relative to the explicitly reviewed W>=15 and
three-dimensional metric inputs.** I read the complete proposed proof
and independently rederived its polygon and graph arguments. No checker,
matrix program, graph census, optimizer, Lean or build ran. No pinned
source or manuscript was edited.

Reviewed source:

* `results/compact_ra4_fifteen_excess_equality_exclusion_research_20261019.md`,
  SHA256 `237d9ccc7c91062b268124fc78c4fa82dbcf4bbf6b0b6b32c8eae013c713ff87`.

The successor conclusion is precisely `W=0 or W>=16`. This review does
not establish a volume minimum, a second-volume gap, the inequality
`f+W>=120`, an equality classification at16, or a dimension bound. The
primary manuscript remains the separate source bound at9fe28180 and
still states the earlier W>=15 theorem until its owner updates it.

## 1. The complete W15 reduction

The distinction between weighted nonpentagonal adjacency G and full
positive-facet adjacency H is maintained. The old argument that minimum
G-degree three implies minimum H-degree five uses only eligibility in
the four-vertex neighbor graph; its proof does not intrinsically require
strict W<15. Its use at the endpoint is therefore valid.

At a degree-two G vertex, the two weights agree. The old distinct-facet
sum gives `W>=5k-15`. At W15 this forces k6 and equality in every
summand. The original facet, its two base neighbors and its twelve
lateral neighbors exhaust all fifteen positive facets, each of excess2.
Every positive facet then has G-degree two with unit weights, hence the
same Loebell L6 type. Each has fourteen distinct positive neighbors by
the old nonregular lateral-pentagon argument. Thus H would be K15,
contradicting the positive index4. This does not declare an uncounted
outside facet zero without the exact equality sum.

Otherwise `3N<=2E<=30`. The reviewed small-H argument excludes N<=9,
so N10, E15, every G-degree3 and every edge weight1. Thus every positive
facet has excess3, exactly three hexagonal faces, and twelve pentagons.
If H had minimum degree at least six, its complement on ten vertices
would be subcubic and induced-2K2-free. Any induced nine-vertex subgraph
would have an independent five by the old lemma, yielding a forbidden
K5 in H. Hence a positive facet of H-degree exactly five really exists.
Its five positive neighboring faces have exactly the three specified
hexagons. Every other face borders a zero-excess regular dodecahedron.

## 2. Actual three-edge closure and the two metric consumers

For the three-edge lemma, the prescribed inward chain fixes both
adjacent support lines at its ends. Their orientations agree with the
regular model because the chain lies in the same inward halfplanes.
They are actual polygon wall lines and are mutually perpendicular.
Projection of a strict interior point off their two outward normals
stays future timelike and strictly satisfies every other wall inequality.
Their intersection is therefore an actual vertex, not an intersection
outside the polygon that some omitted wall could truncate. The cyclic
boundary closes with these five lines. This is the essential stronger
step beyond the earlier four-edge lemma, and it is justified.

The two-adjacent-L-edge pentagon lemma is also valid. Its two fixed
edges determine their predecessor and successor perpendicular lines.
These are ultraparallel in the matched regular model. The remaining
fifth side must be their unique common perpendicular, determining the
whole pentagon. Three L-edges in a pentagon contain an adjacent pair,
so the stated consequence follows without a three-dimensional
classification input beyond the common value L.

For a hexagon with four L sides, two unspecified side positions which
are not opposite leave three consecutive L sides, already impossible.
For cyclic lengths `a,L,L,d,L,L`, the actual outward-normal computation
gives

    cosh(a) = (cosh(L)^2+cosh(d))/sinh(L)^2,
    cosh(d) = (cosh(L)^2+cosh(a))/sinh(L)^2.

The signs agree with the actual ultraparallel c,e walls connected by
side d. If `p=cosh(L)=phi`, subtracting these equations gives equal
unknown coshes, and either equation then gives
`p^2/(p-1)=p^3=2p^2-1=cosh(2L)`. Both side lengths are positive,
so both equal2L. This is an exact metric identity, with no numerical
angle approximation or unchecked hexagon symmetry assumption.

## 3. Marked-vertex restrictions are actual face restrictions

In the five-neighbor graph J, the degree of a marked vertex A is
exactly the number of positive edge carriers of the actual hexagon
`F intersect A`. Three orthogonal facet normals give actual edges;
four give actual vertices. Thus adjacency of carrier vertices is
adjacency of their polygon edges, with the converse coming from
simple actual incidence.

Degree zero or one is excluded by three consecutive zero-carrier
L-edges. At degree two, the hexagon lemma forces the two positive
carriers to be opposite and their lengths to be2L; adjacency of the
carrier vertices is impossible. At degree three, a three-vertex path
among the carriers gives a chain of three consecutive positive edges
and hence a complementary chain of three L-edges, which the
three-edge lemma excludes. A triangle among them gives K4 in J
already. No further restriction is being silently imposed at degree
four.

## 4. Every five-neighbor graph branch

I checked the graph split independently:

* With a universal vertex, its other four vertices form a triangle-free,
  induced-C4-free graph, hence a forest. Vertices of its degree zero
  or one are ineligible; degree two gives exactly the forbidden
  three-neighbor path centered at the universal vertex. Only a
  degree-three forest vertex can still be marked. A four-vertex
  forest has at most one such vertex, so even with the universal
  vertex there are at most two marks.
* Without a universal vertex the maximum degree is three. With a
  triangle and two nonadjacent outside vertices, neither outside
  vertex is marked. A marked triangle vertex must have an exclusive
  contact to one outside vertex; a shared contact gives its forbidden
  three-neighbor path. Two outside vertices support at most two marks.
* If those outside vertices are adjacent, their triangle-contact sets
  are disjoint by the degree bound and nested by induced-C4 exclusion.
  Thus one is empty. The remaining contact-set sizes0,1,2 give at
  most0,2,1 marks respectively, exactly as the proof states.
* In the remaining triangle-free case, induced-C4 exclusion leaves a
  forest or C5. A forest on five vertices with three vertices of
  degree at least two must be P5, with precisely its three internal
  vertices marked. This follows also from the elementary leaf count;
  disconnected forests cannot supply the three internal vertices.

Thus P5 and C5 are the only unexcluded abstract cases; no computational
graph classification is being assumed. In P5 an endpoint pentagon
has four zero-carrier L-edges and shares a2L edge with its marked
hexagonal neighbor, contradicting regular closure. In C5, some
unmarked pentagon is adjacent to a marked hexagon. Its three
zero-carrier L-edges contain an adjacent pair, so the pentagon is
regular and again cannot have that shared2L edge. These are whole
actual edges of both adjacent faces, not unrelated metric segments.

## Conclusion

No missed equality case, carrier identification, graph branch, sign in
the hexagon formula, or hidden use of an omitted geometric classification
was found. The paper proves the claimed successor `W=0 or W>=16`
relative to the inherited W>=15 theorem, actual-face/acute-wall inputs,
and its explicitly stated three-dimensional metric premises. Its
separation from volume minimization and its unchanged46/46/96 local
cutoffs are correct. This is written mathematical verification only.

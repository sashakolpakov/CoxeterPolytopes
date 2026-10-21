# Root review: exclusion of compact RA4 excess fifteen

19 October 2026. **Full written mathematical review PASS.** I read the
complete source and checked its composition with the reviewed gap-15
interfaces. No graph census, symbolic program, enumeration, Lean batch
or TeX build ran. The primary manuscript at SHA `9fe28180...76d47738`
is unchanged and still states the earlier gap.

Reviewed source:
`compact_ra4_fifteen_excess_equality_exclusion_research_20261019.md`,
SHA256 `237d9ccc7c91062b268124fc78c4fa82dbcf4bbf6b0b6b32c8eae013c713ff87`.
Inherited complete proof:
`compact_ra4_positive_facet_link_gap_fifteen_20261019.md`, SHA256
`8201b2f5cc7da5988b6d164b258da8e14981cf94c6b5e980a5a87e9deb70dcd7`.

## Equality reduction

At W=15 a degree-two G vertex forces equality in the distinct-neighbor
excess sum. All fifteen positive facets then have excess two. The
degree-one prohibition makes every G vertex degree two with unit
weights, so all facets are Loebell L6. Each has fourteen distinct
positive neighbors, giving the impossible K15 in H. No outside-zero
assumption enters that equality argument.

Otherwise `3N<=2E<=30` gives N at most ten. The earlier eligibility
proof gives minimum H degree five whenever minimum G degree is three;
it did not need strict W<15 for any of its geometric steps. The
at-most-nine-vertex graph obstruction thus gives N=10, E=15, cubic G
and unit edge weights. Every positive facet has three hexagons and
twelve pentagons. If minimum H degree were six, a nine-vertex induced
subgraph of its subcubic, induced-2K2-free complement would have an
independent five. Hence some F has exactly five positive neighbors.

## Polygon statements

Three consecutive regular-pentagon lengths determine an oriented
right-angle chain and the two prescribed outer perpendicular support
lines. In the matched regular pentagon those lines are orthogonal.
They are actual walls of the given polygon. The inherited acute
orthogonal-projection lemma makes their intersection an actual vertex,
so an extra wall cannot truncate it. This closes the five-edge boundary.
The source correctly supplies this actual-intersection step; merely
matching an extended chain would not suffice.

For a pentagon, two adjacent fixed edges prescribe the predecessor and
successor support lines, whose unique common perpendicular is the fifth
side. Thus the pentagon is regular. This also applies to any three
specified regular-length edges, since two among them are adjacent.

In a hexagon, four regular-length sides leave two unspecified positions.
Unless opposite, the four fixed positions contain three consecutive
edges, contradicting the first lemma. With cyclic lengths
`a,L,L,d,L,L`, the source derives the right-angled cosine relation
directly from the two outward normals. Writing p=cosh(L), it gives
`A=(p^2+D)/(p^2-1)` and its opposite relation. Since p is the golden
ratio, A=D=p^3=cosh(2L). The two exceptional lengths are exactly 2L.

## Exhaustive five-neighbor link argument

J is the induced graph on the five positive neighbors of F. Three
vertices are marked hexagons. Actual third-/four-wall intersections
identify their neighbor graphs with the corresponding edge carriers.
Zero-excess carriers give edges of length L. A marked vertex has
degree at least two; at degree two its neighbors are nonadjacent and
the two positive edges have length 2L. At degree three, a neighbor path
would give three consecutive positive edges and three consecutive
L edges, which is impossible. A neighbor triangle gives K4 in J.

I checked every case of the finite graph proof:

* A universal vertex leaves a four-vertex forest. Vertices of residual
  degree zero/one/two fail the marked tests. At most one residual
  degree-three vertex and the universal vertex can be marked.
* With maximum degree three and a triangle, outside vertices with no
  mutual edge cannot be marked. A marked triangle vertex needs an
  outside contact unique to it, giving at most two marks.
* If the outside edge exists, its two triangle-contact sets are
  disjoint by maximum degree three. If both were nonempty they would
  form an induced four-cycle. One set is therefore empty; the three
  remaining contact-set sizes give at most two marks as listed.
* A triangle-free graph has no cycle except a possible C5. A forest
  on five vertices with three eligible vertices must be the P5: a
  tree component with three nonleaves requires at least five vertices,
  and its three nonleaves all have degree two. These are the only
  cases left, not an assumed graph census.

In P5 an endpoint pentagon has four L edges and is regular. In C5
a pentagon adjacent to a marked hexagon has three L edges and is
regular. Both contradict the shared edge of length 2L. Thus the
degree-five facet forced by equality cannot exist.

## Conclusion and verification boundary

The proof establishes `W=0 or W>=16`, conditional only on the already
reviewed geometric interfaces and explicitly named metric inputs.
It does not classify arbitrary ten-vertex opposition graphs; in
particular it does not falsely reject the weaker C5[K2] graph model.

This is a written theorem with root review, pending the separate
independent review and successor TeX integration. It is not a geometric
Lean result, a new dimension bound, a proof of the 120-cell minimum,
or a second-volume gap. The local size cutoffs remain 46/46/96.

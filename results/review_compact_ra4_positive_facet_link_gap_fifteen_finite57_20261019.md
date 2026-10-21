# Independent written review: the compact RA4 positive-excess gap

19 October 2026. **PASS within the explicit geometric inputs below.** This is
a full source and mathematical review, not an execution, graph census, kernel
check, or independent classification of compact polyhedra. No pinned source
was changed.

## Source and conclusion

Reviewed in full:

* `results/compact_ra4_positive_facet_link_gap_fifteen_20261019.md`, SHA256
  `8201b2f5cc7da5988b6d164b258da8e14981cf94c6b5e980a5a87e9deb70dcd7`.

The source proves, for a compact right-angled hyperbolic four-polytope,

    W = sum_ridges (number_of_sides - 5) = 0 or W >= 15.

Its zero case is the inherited regular 120-cell recognition. The argument
does not prove a minimum-volume theorem, a second-volume theorem, the
inequality `f + W >= 120`, classification at `W = 15`, or a dimension bound.

The first section and the four-edge polygon closure were also reviewed in
`results/review_compact_ra4_pentagon_propagation_finite57_20261019.md`, SHA256
`d1cde5b9e9f34fd06471c329f85f865c7f4990ed270e7225bcbd263729326bec`.
The new review below checks the stronger full-adjacency argument and every
case of its small-graph proof independently.

## Geometric interfaces and their use

The inherited three-dimensional capacity result, its Loebell equality
recognition, actual-face projection, compact strict normal gap, and
zero-excess recognition remain mathematical premises supplied by the named
reviewed papers. They are not consequences of the graph lemma. In
particular, the source explicitly uses two three-dimensional metric inputs:
zero-excess facets are regular right-angled dodecahedra, and marked Andreev
uniqueness realizes the Loebell cyclic symmetry. This review does not turn
those external inputs into formalized or computation-verified facts.

With these inputs, Section 1 is sound. At a degree-two vertex of the
nonpentagonal-ridge graph `G`, capacity makes the two positive integer edge
weights equal. Equality identifies the entire facet as Loebell type. Its
base is regular, with

    cosh(L_k) = 1 + 2 cos(2 pi/k) > cosh(L_5),  k >= 6.

Consequently every one of its `2k` lateral pentagons is nonregular and
cannot belong to a zero-excess neighboring facet. Convexity makes all its
ridge-neighbors distinct. The total distinct-facet excess is therefore at
least `10k - 30`, giving `W >= 15`. Thus the contradiction hypothesis
`0 < W < 15` really gives ordinary minimum degree three in `G`.

The four-edge closure in Section 2 is a metric polygon argument. An
oriented chain of four right-angled sides of the regular-pentagon side
length is fixed up to isometry. At its endpoints the next two supporting
lines coincide with the fifth pentagon line. A convex polygon has only one
boundary edge on that line, so insertion of further sides is impossible.
An edge carried by a zero-excess third facet has exactly this side length:
it is an entire edge of a regular pentagonal face of that facet.

The passage from actual geometry to the neighbor graph uses the full
adjacency graph `H`, not `G`. Three pairwise adjacent facets have mutually
orthogonal unit normals; projection of an interior point onto their
orthogonal complement produces the actual edge in their intersection.
Four such facets similarly produce an actual finite vertex. Compact
right-angled simplicity makes each edge of a fixed ridge have exactly one
third-facet carrier. Hence, in the graph on `N_H(F)`, the degree of `A`
is exactly the number of positive carriers on `F intersect A`.
Adjacency of two carriers gives adjacent edges of that polygon, not merely
a formal diagram contact. These points justify both uses of the
four-zero-edge prohibition and the eligibility definition.

## The minimum-degree step

I checked the four-vertex graph fact separately. If a degree-three vertex
exists, the other three vertices have zero, one, or two edges, unless a
forbidden `K4` occurs. The first two cases give only the original eligible
vertex; the two-edge path gives exactly its middle and the original
vertex. The other path vertices have adjacent neighbors and are
ineligible. Without degree three, every component is a path or cycle;
triangles have no eligible vertices, induced squares are excluded, and
paths on at most four vertices have at most two eligible internal
vertices. This covers also graphs of fewer than four vertices.

Each positive facet has at least three `G`-neighbors, all eligible. Its
`H`-neighbor graph has no `K4`, since that would give an ambient `K5`,
and it has no induced square. Its order therefore cannot be at most four.
This proves `minimum_degree(H) >= 5` without assuming that all positive
adjacencies themselves carry nonpentagonal ridges.

## The complete graph obstruction

Under `0 < W < 15`, integrality gives

    3N <= 2|E(G)| <= 2W <= 28,

so `N <= 9`. For `C = complement(H)`, the exact transferred conditions are
`Delta(C) <= N-6`, no induced `2K2`, and no independent set of size five.
The cases `N <= 6` are immediate. For `N = 7`, at most one complement edge
is possible, leaving an independent six. For `N = 8`, at most one
component can contain an edge. With maximum degree two it is a path on at
most four vertices or a cycle on at most five; adjoining the isolated
vertices gives the asserted independent sets. These bounds do not use a
computer enumeration.

For `N = 9`, I checked both branches of the stated lemma:

* A complement triangle has at most three outgoing degree slots. Its
  six outside vertices contain an independent set `S` of at least three
  with no triangle contact. If `|S| >= 4`, adjoining a triangle vertex
  works. Otherwise all three slots give a perfect matching to the other
  three outside vertices. An edge from one of these to `S`, together
  with the triangle edge avoiding its unique triangle neighbor, would be
  an induced `2K2`. Thus `S` consists of isolated vertices; an unmatched
  triangle/outside pair completes an independent five.
* In the triangle-free nonbipartite case, a shortest odd cycle is an
  induced `C5`: every induced cycle of length at least six contains an
  induced `2K2`. Each of the four outside vertices has either zero or
  two cycle neighbors. A single neighbor gives the explicitly disjoint
  cycle edge obstruction, and three neighbors contradict
  trianglefreeness. The cycle has at most five outside degree slots, so
  at least two outside vertices have zero contact. Each is isolated:
  the at most two cycle neighbors of a proposed outside partner cannot
  cover all edges of `C5`. Three isolated vertices suffice immediately;
  with exactly two, take another outside vertex and an independent pair
  among the three cycle vertices avoiding its two neighbors. Together
  they give five. The bipartite case has an independent part of size
  at least five directly.

Thus every allowed `N` is excluded. No missing finite graph case or
geometric-to-graph identification was found.

## Status

The gap `W = 0 or W >= 15` passes this independent written review relative
to the explicitly inherited capacity, actual-face, and metric inputs. It
is stronger than the earlier `W >= 9` note. The source's caution at the
endpoint and its separation of this gap from volume minimization are
necessary and correctly retained. No proof program, build, optimizer, or
scientific checker was run for this review.

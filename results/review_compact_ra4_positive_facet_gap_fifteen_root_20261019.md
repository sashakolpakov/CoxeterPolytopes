# Root review: compact RA4 positive-excess gap fifteen

19 October2026. **Complete written mathematical review PASS**, relative
to the named inherited actual-face/capacity lemmas and explicit3D metric
uniqueness/classification inputs. No census, script, build or Lean job
was run for this review. Independent peer review is separate.

Read the complete
`compact_ra4_positive_facet_link_gap_fifteen_20261019.md`, SHA256
`8201b2f5cc7da5988b6d164b258da8e14981cf94c6b5e980a5a87e9deb70dcd7`, all its
small-graph cases, and the earlier Loebell/four-edge propagation sources.
The statement is W=0 or W>=15, not universal minimum volume.

The distinction between G, the weighted nonpentagonal-ridge graph, and
H, FULL adjacency on the positive-excess facets, is correct and essential.
A positive facet has excess at least2; an edge of weight w forces
endpoint excess at least2w. Degree-two vertices therefore have equal
weights and Loebell combinatorics by the own3D closing/equality proof.
The stated Andreev uniqueness clause realizes base rotation, so each
base is a regular right-angled k-gon. Its edge length exceeds that of
a regular pentagon. Thus every lateral pentagon forces a distinct
positive neighbor. Summing those actual facets gives2W>=10k-30,
so a degree-two vertex requires W>=15.

The four-consecutive-regular-pentagon-edge closure proof is valid:
coherent right-angle turns fix the oriented chain; the two remaining
support lines are the same fifth line, and convexity forces closure.
Inside a nonpentagonal ridge, a zero-excess third facet forces an
entire regular-pentagon edge. Common H-neighbors correspond to actual
ridge edges by orthogonal projection; adjacent common neighbors
correspond to adjacent edges, by projection of four orthogonal normals.
No supporting-plane intersection is substituted for an actual face.

These facts give the exact eligibility condition on a positive facet's
neighbor graph. I checked all graphs on at most four vertices in the
written argument: with a degree-three vertex the residual graph has
zero, one or two edges; without one, only paths/triangles occur after
excluding C4. At most two vertices are eligible. Three G-neighbors
therefore force minimum H-degree five under0<W<15.

Then3N<=2E<=2W<=28 gives at most nine positive facets. In the complement
C, maximum degree is at mostN-6; there is no induced2K2 and no independent
five. The proof covers N<=6, seven and eight explicitly. For nine,
I checked every step of the triangle/bipartite/shortest-odd-cycle proof.
In the triangle case, three outside contact slots force either an
immediate independent five or a perfect matching to three outside
vertices, with the other three isolated. In the trianglefree nonbipartite
case the induced cycle must be C5; each outside vertex has zero or two
cycle contacts, at least two have zero, and induced2K2 exclusion makes
them isolated. The residual cases all produce an independent five.
This is a complete own graph proof, not an unverified graph census.

The contradiction establishes W=0 or W>=15. Recognition of W=0 as
the whole120-cell remains the earlier theorem with its explicit3D
metric classification. A distinct polytope at or below that volume
must have W>=15 and every facet at most46 vertices. This still leaves
the global f+W>=120 inequality unproved, as well as the double-volume
gap and any new dimension ceiling. The W=15 endpoint is not excluded.
No geometric Lean or script-based graph verification is claimed.

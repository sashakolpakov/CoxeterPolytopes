# Root review: whole dodecahedral blocks in uniform neighbors

20 October 2026. Full written review of
`compact_ra4_uniform_dodecahedron_neighbor_block_20261020.md`, SHA256
`69d900288cd778e6cda38cb6eedb932ef0583eb363649f60727be41dd98b9da3`.
Paper only; no new script, primary manuscript or Lean result.

The flower hypothesis concerns six ENTIRE regular pentagonal faces.
Their planes and inward polygons agree uniquely with a regular
dodecahedron. In a simple right-angled polyhedron each outer edge
therefore forces its true neighboring support plane. The ten such
outer edges give the five distinct canonical B_i planes, since the
two occurrences of the same plane belong to the same convex face.
Together with the flower this gives eleven actual canonical planes
and fifteen actual canonical vertices.

I checked the orientation step explicitly. The five B_i normals in
the opposite plane's Lorentz space sum to a positive multiple of
its future center. Their homogeneous feasible cone has positive
time except at zero. Its time-one Euclidean section is bounded,
and compactness of the regular pentagon puts that entire section
strictly inside the Klein disk. Thus every nonzero vector in the
cone is future timelike. Projecting any other actual normal into
that space therefore makes its product strictly negative throughout
the opposite pentagon. A zero projection leaves only the two
possible orientations of its unit perpendicular; feasibility of
the already actual flower selects the canonical one. This is not
an unsupported choice of future sheet.

All five opposite vertices are consequently feasible. With the
fifteen flower vertices they prove containment of the ENTIRE
dodecahedron by convex hull. If its opposite plane is an actual
wall, the eleven other actual halfspaces give equality. Otherwise
the cut has a strict interior and exactly its five B_i side walls.
Those walls are perpendicular to the cut, so the far piece is
compact and right-angled. Euler, trivalence and face size at least5
give its vertex count at least20. The seam's five new vertices
occur once in each piece and are not vertices before cutting;
the count V(Q)=20+V(far)-10 is therefore correct. At equality,
the stated20-vertex metric recognition gives the reflected double.
No garland claim for arbitrary larger Q is used.

For the H4 application, at each edge of the shared pentagon the
two adjacent sides in its perpendicular polygon are precisely the
outgoing edges of the marked parent dodecahedron. Uniformity gives
all three consecutive lengths L. The imported actual polygon
closure theorem makes each entire neighboring polygon regular.
This verifies ALL five members of the required flower, rather
than only its first local edge chain. Each non-dodecahedral
neighbor of this UNIFORM facet has d>=5.

Verdict: PASS within the stated marked uniform-facet class. Root
discussed the flower/cone mechanism with the author, so this is not
an independent outside-author audit or external refereeing. The
argument does not give d>=5 for every positive facet of an arbitrary
four-polytope, prove f+W>=120, or identify a whole120-cell inside P.

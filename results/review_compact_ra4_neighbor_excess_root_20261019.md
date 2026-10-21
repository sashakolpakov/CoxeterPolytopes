# Root review: neighboring pentagon excess in compact RA4

19 October2026. **Written mathematical review PASS**, with the explicit
inherited capacity proof and the named external3D metric classification
as specified below. No census, numerical check, build or Lean job ran.

Read the complete source
`compact_ra4_neighbor_excess_and_zero_excess_rigidity_20261019.md` and
reread the original3D closing/equality argument and the scope audit.
The earlier frozen compact-capacity manuscript and proof packet remain
unchanged. This review does not assert an established universal compact
minimum or second-volume theorem.

For a facet F, d(F)=(V(F)-20)/2 is the sum of shared-ridge excesses.
Each adjacent facet across a k-gon has d>=2(k-5) by the already proved
3D capacity. Summing the selected DISTINCT facets against the identity
sum d=2W gives2W>=3d(F). Shared ridges cause no overcount beyond that
identity. Combined with f>=2V+5 this proves

    f+W >= 2V+5+ceil(3(V-20)/4).

At V=48 the lower bound is122>120; at V=98 it is260>256.
Since facet vertex counts are even, every facet has at most46 vertices
at or below the120-cell volume, and at most96 strictly below the double
volume. The endpoints46/96 are not ruled out. This improves the older
48/50/104 necessary local cutoffs, not a dimension ceiling.

The positive-excess graph argument also checks. An edge weight w forces
weighted degree at least2w at either endpoint. A degree-two vertex must
have equal edge weights; equality in the original closing proof gives
Loebell combinatorics with the two nonpentagonal faces disjoint. Thus
its two positive neighbors cannot meet anywhere in the actual polytope,
by acute orthogonal projection. A graph with at most four edges and
minimum degree two is a triangle or a four-cycle. The former contradicts
that disjointness; the latter is induced in the full facet graph and
would give two negative directions in the unit Gram matrix. Hence
W=0 or W>=5. The original source explicitly proves the needed equality
clause, so no unproved extra3D classification is hidden here.

The W=0 recognition separately uses the NAMED external fact that every
compact right-angled3-polyhedron with12 faces is the regular metric
dodecahedron. Under this premise, all facets, not one selected facet,
are regular dodecahedra. A matched whole facet fixes each adjacent
supporting hyperplane by orthogonality, and its matched pentagonal
ridge fixes the adjacent whole regular dodecahedron on the inward side.
Connectivity in both facet graphs makes the matched set all facets of
the regular right-angled120-cell. This is a valid global recognition
argument; it does not propagate an arbitrary lone dodecahedral facet.

A possible smaller or distinct equal-volume polytope must therefore
have W>=5 and all facets with at most46 vertices. Many such positive-
excess configurations remain unresolved. No universal minimum, volume
gap, garland theorem, new example or formal geometric verification is
claimed. These results belong in the existing compact-capacity paper
when its frozen predecessor is revised with appropriate provenance;
they do not need a separate manuscript.

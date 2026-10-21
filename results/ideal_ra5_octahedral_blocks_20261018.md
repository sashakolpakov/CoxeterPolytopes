# Octahedral blocks and a stronger restricted ideal RA5 carrier bound

18 October 2026. This is a written proof with an exact local-data certificate.
It assumes the same explicit three-type facet spectrum as the
[preceding cusp theorem](ideal_ra5_known_carrier_cusp_compatibility_20261018.md):
the regular ideal right-angled24-cell C, its facet double D, and the
particular sixteen-cell polytope T in the
[positive-carrier construction](ideal_ra4_positive_kappa_doubling_20261018.md).
It does not classify arbitrary ideal right-angled four-polytopes.

**Theorem.** In an entirely ideal finite-volume right-angled hyperbolic
five-polytope with that facet spectrum, the T facets split into blocks of
six under adjacency along24-vertex ridges. Each block induces an
octahedral graph. At most two facets can be adjacent to every member of
one block. Consequently there are at least36 T facets and84 D facets.
In particular, the old necessary count assignment with30 T facets and
66 D facets cannot be realized by a polytope.

These are restricted facet-count improvements. No dimension ceiling,
nonexistence of this whole spectrum, or new geometric example is asserted.

## Exact intrinsic data

Use outward unit normals and Lorentz signature(4,1) in T. The inherited
certificate lists all198 walls and258 vertices. Exactly four walls have
24 vertices, with indices5,14,30,37. In that order their Gram matrix is

\[
 \begin{pmatrix}1&0&0&-1\\0&1&-1&0\\0&-1&1&0\\-1&0&0&1\end{pmatrix}.
 \tag{1}
\]

The two opposite-pair sums are positive multiples of the same null vector,
the saved ideal vertex85. Their span therefore has two positive directions
and a one-dimensional radical. Every24-vertex three-ridge has26 polygon
faces: twelve triangles, twelve quadrilaterals and two hexagons. Precisely
the two hexagons are its intersections with other24-vertex ridges of T.

The [new producer](../scripts/fresh_ideal_ra5_octahedral_blocks_20261018.py)
reconstructs every facet-vertex incidence from exact rational coordinates,
checks the complete104 polygon-face records of these four ridges, (1),
and the common null ray. It inherits the actual completeness of T from
the reflection-doubling proof; it does not replace that proof with a new
bounded incidence census. C and D have no24-vertex ridges.

## The graph of24-vertex ridge adjacencies

Let H have one vertex for each T facet, with an edge when their common
three-ridge has24 vertices. Each vertex has degree four, because T has
exactly four such ridges, all necessarily adjacent to T facets.
Different ridges have different neighbors: the intersection of two
convex facets is one face.

Fix one T facet A. Its four neighbors have native normals given by (1),
since the polytope is right-angled. At each of its24-vertex ridges R,
the two hexagonal two-faces are intrinsically distinguished by their
six vertices. The same two hexagons are seen from the T facet on the
other side of R. In either T they are exactly the intersections with
the two other24-vertex ridges meeting R. Thus the two corresponding
neighbor facets are adjacent in H to both ends of the edge through R.
The opposite pairs in (1) have product minus one, so their facets are
not ridge-adjacent even in the full facet graph. Consequently the graph
induced by the four neighbors of A is exactly a four-cycle.

Here a two-face at the intersection of three mutually orthogonal facets
has exactly those three containing facets: its relative interior is
finite and the right-angled local link there is a simplex. Hence the
hexagon correspondence identifies actual neighboring facets; it does
not merely identify unlabeled polygon counts.

A finite simple four-regular graph whose every neighborhood induces a
four-cycle is a disjoint union of octahedral graphs. For a direct proof,
write the neighbors of v in cyclic order as a,b,c,d. The fourth neighbor
w of a is outside this set, since a and c are nonadjacent. The four-cycle
in the neighborhood of a forces w adjacent to b and d and not to v.
The neighbors of b are now v,a,c,w. Their four-cycle forces w adjacent
to c. The same argument gives all edges of the octahedron on
v,a,b,c,d,w; all six degrees are filled, so this is a component.
This proves the claimed partition into six-element blocks.

## At most two common neighbors per block

Work in the ambient Lorentz space of signature(5,1). Choose one normal n
in a block and the four normals r of its neighbors in H. The latter
span a three-dimensional space with two positive directions and null
radical v, by (1). Since they lie in n-perpendicular, adjoining n gives
a four-dimensional space U of signature(3,0,1). The sixth normal m is
orthogonal to all four r. The octahedral incidence inside the neighboring
T facets and (1) give n dot m=-1. Thus m+n is perpendicular to U and
has square zero. The orthogonal complement of U has one positive
direction and null radical v, so m+n is a multiple of v. In particular
the six normals span U.

Any facet adjacent to all six has its unit normal in U-perpendicular.
Choose a unit positive vector e there. Every such normal is of the form
plus or minus e plus a multiple of v. Two with the same sign have product
one, contradicting the nonpositive product of distinct outward normals
of an acute-angled polytope. There are therefore at most two common
neighbors. This applies to all facets, and hence also to T facets.

## A four-color obstruction to30 T facets

Let G be the full ridge-adjacency graph on T facets. Every T has22 ridges
with at least14 vertices. C and D have no such ridges, so G has minimum
degree at least22. We first strengthen the previous positive-index
clique observation: G has no five-clique because the ambient polytope
is entirely ideal.

Indeed, if five outward unit normals n_i are mutually orthogonal, take
a strict timelike interior vector q and set

\[
 p=q-\sum_{i=1}^{5}(q\mathbin\cdot n_i)n_i.
\]

Then p squared equals q squared minus the sum of the five squared
coefficients, so p is timelike and in the same time orientation as q.
It is perpendicular to the five selected normals. For any other outward
normal m, both q dot n_i and n_i dot m are nonpositive, and therefore
p dot m is at most q dot m, which is strictly negative. After
normalization p is an actual finite vertex, contrary to the ideal
hypothesis. This proves the absence of a five-clique without assuming
that pairwise intersections automatically have a common face.

The preceding cusp theorem already proves6 divides z and z at least30,
where z is the number of T facets. Suppose z=30. The published
[Andrasfai--Erdos--Sos theorem, Theorem1.1](https://real.mtak.hu/110445/1/1974-18.pdf)
with r=5 says that a graph without a five-clique and minimum degree
greater than8z/11 is four-colorable. Here22 is greater than240/11.
This graph theorem is an explicit external input, not a new internal
lemma or an asserted Lean proof.

In a proper four-coloring every color class has size at most eight,
because each vertex has at least22 neighbors among the other30 vertices.
There must be four nonempty classes, and their sorted sizes are either
(6,8,8,8) or(7,7,8,8). A vertex in a class of size eight is adjacent to
every vertex outside that class. A vertex in a class of size seven can
miss at most one vertex outside its class.

Consider any one of the five octahedral blocks. A proper coloring of its
three opposite pairs uses either three colors, with each pair
monochromatic, or four colors, with exactly one opposite pair
bichromatic. In the three-color case let the omitted class have size s.
If s=8, all eight vertices are common neighbors of the block. If the
class sizes are(6,8,8,8) and s=6, all six are common neighbors. In the
remaining case s=7 in(7,7,8,8), at most two vertices of the omitted class
can miss a block vertex: only the two block vertices in the other
seven-element class can have missing external edges, and each can miss
at most one. Thus there are at least five common neighbors. Every case
contradicts the upper bound two.

The block therefore uses four colors. Its bichromatic opposite pair is
a nonedge of G, so neither endpoint can be in an eight-element class.
This excludes(6,8,8,8). In(7,7,8,8) the split pair must use the two
seven-element classes. Each other opposite pair is monochromatic in
one of the eight-element classes. Every block thus contains color
counts(1,1,2,2). Five disjoint blocks would give class sizes(5,5,10,10),
contrary to(7,7,8,8). This excludes z=30.

Together with6 divides z, this proves z at least36. The inherited
carrier identity y=3z-24 gives at least84 D facets. No other facet-type
restriction or census is silently added to the three-type hypothesis.

## Reproducibility and verification scope

Run the producer once to generate a new certificate, or use the following
command to replay the saved exact data without modifying it:

```sh
python3 -B scripts/fresh_ideal_ra5_octahedral_blocks_20261018.py --verify
```

Generation and a separate-process replay both passed in about0.07 seconds
per process, in chunks fb0c22 and44d545. Alongside the exact local geometry,
the script exhausts4096 colorings for each possible part-size vector and
checks every finite numerical step in the final coloring argument. It
does not encode the universal Lorentz or local-to-global graph proofs.
Those are the written arguments above. An independent abstract graph
review is recorded in
[the separate review](fresh_ra5_t30_abstract_graph_review_20261018.md).
The [independent geometric review](ideal_ra5_octahedral_blocks_geometric_review_20261018.md)
reads the complete proof and producer, checks the primary coloring theorem,
and verifies all geometric reductions. Its
[standalone checker](../scripts/check_ideal_ra5_octahedral_blocks_saved_20261018.py)
imports no producer, reconstructs all788 relevant wall-pair intersections,
checks polygon rank and all104 polygon-face records, and verifies the source
hash bindings. Its final fresh replay passed in chunk87bd0c. Root separately
read the checker and ran its pre-binding version in c8051f; that version
already checked all geometric records. The three-page manuscript compiled
in eeadd1 with no final-pass warnings; temporary build files were removed.
No complete Lean or external human review is claimed.

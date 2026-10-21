# Adversarial written review: octahedral blocks form a straight path

18 October2026. **PASS**, as a strengthening of the restricted
[Oct/D9 block theorem](ideal_ra4_octahedral_block_decomposition_research_20261018.md).
This review checks the new path argument directly. It assumes the preceding
written decomposition proof, including disjointness of the octahedral cut
sections even at ideal vertices. No program, census, compilation or Lean
was run. Earlier frozen evidence was not changed.

## The root calculation

Use the regular ideal24-cell whose outward unit normals are

    n_r=(1,r),  r in {+/-e_i +/-e_j : 1<=i<j<=4},

in Lorentz signature(-,+,+,+,+). Thus r^2=2 and
n_r dot n_s=-1+r dot s. Distinct D4 roots have products-2,-1,0,1.
Product1 means the two facets meet orthogonally along a triangle;
product0 means they share an ideal cusp. These intersections are actual:
for instance the Klein vertices are the eight vectors+/-e_i and the
sixteen vectors(+/-1,+/-1,+/-1,+/-1)/2. Consequently two cut facets whose
closures are disjoint have root product-1 or-2.

If r dot s=-1, c=r+s is again a D4 root: writing the two signed-coordinate
roots shows their unique common coordinate cancels. It has norm2 and

    c dot r=c dot s=1.

The24-cell facet C with normal n_c is therefore orthogonal to both cut
facets B_r,B_s. In particular C cannot itself be a cut facet, because
the cut sections are closure-disjoint. Its whole octahedron is part of
one original boundary facet F of P.

## The three-piece boundary contradiction

Across B_r lies another block. The regular octahedral facet on the other
side sharing the triangle C cap B_r is the reflection of C across B_r.
Indeed reflection fixes C's supporting wall, since its wall is orthogonal
to B_r; and a regular ideal octahedron is determined by a triangular face
and the side on which it lies. This reflected octahedron is also uncut:
it meets the existing cut B_r, so it cannot be a second internal cut.
It is therefore another whole boundary piece in the same original
supporting wall as F. The same reasoning across B_s gives a third whole
octahedron in F.

These three relative interiors are disjoint. Otherwise two blocks would
share a relatively open three-dimensional boundary set in the same wall.
Their interiors both lie on the P side of that wall, so a small inward
neighborhood there would belong to both block interiors, contrary to the
decomposition into blocks. Equivalently the induced boundary subdivision
is a cellulation. If v_O>0 is the hyperbolic three-volume of a regular
ideal octahedron, this proves Vol_3(F)>=3v_O. But the original hypothesis
allows only O and D9 facets, with volumes v_O and2v_O. This contradiction
excludes root product-1.

It follows that any two cut facets of one block have root product-2,
which means s=-r: they are the unique opposite pair. There cannot be
three pairwise opposite roots. Hence every node of the previously proved
finite connected block tree has degree at most2. The tree is a path.

## Straightness and uniqueness of placement

In a middle block, the outgoing facet is the unique facet opposite the
incoming one. A regular24-cell is uniquely placed once a full octahedral
facet and its ambient side are specified. To justify the absence of
twists: every isometry of an octahedral facet extends to a symmetry of the
regular24-cell. Choose a full flag of the facet; flag transitivity gives
a cell symmetry with the prescribed flag action, and an octahedral
isometry fixing a full flag is the identity. After canceling this facet
action, an ambient isometry fixes its hyperplane pointwise, and is either
the identity or reflection. The specified side decides between them.

Thus, starting from one end block and its first cut, every subsequent
block is forced by reflection across the preceding cut, and every new
cut is the opposite facet in that block. Existing polytopes with a fixed
number of blocks are isometric to this straight garland. This conclusion
classifies the placements that occur under the Oct/D9 hypothesis; no
separate existence assertion for every length is needed for the proof.

The earlier phrase that branching was permitted must consequently be
removed from the theorem and scope statements. The stronger conclusion
is still restricted to the Oct/D9 facet class and gives no universal
24-cell decomposition theorem for arbitrary entirely ideal RA H4
polytopes. This is a written mathematical review, not a priority claim or
end-to-end kernel verification.

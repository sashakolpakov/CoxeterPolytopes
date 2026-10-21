# Dimension restrictions and a changed H22 cusp construction

No new finite-volume example or dimension is established. The literature
excludes the noncompact arithmetic route above dimension 21. The constructive
advance in this checkpoint is an exact H22 chamber that removes one old wall
and cuts an inherited cusp, escaping a restriction of the previous fixed
210-wall core. Its volume is still infinite. A separate H9 construction
closes a compact hexagonal face but also retains infinite volume.

## Literature and search restrictions

Esselmann's maximal reflective integral Lorentzian lattice rank is 22,
meaning hyperbolic dimension 21. The noncompact arithmetic dimension set is
`n<=19` or `n=21`, with Borcherds attaining 21. Vinberg excludes all arithmetic
reflection groups, and all compact Coxeter polytopes, for `n>=30`. The
dimension-21 theorem therefore does not exclude compact arithmetic examples
in dimensions 22–29 over nonrational fields. See
[Esselmann's original article](https://doi.org/10.1006/jnth.1996.0141),
[Vinberg's original article](https://www.mathnet.ru/eng/mmo444), and the explicit
distinction in [Belolipetsky, Section 1](https://arxiv.org/pdf/1506.03111).

These exclusions now appear at the top of README, GOAL and RESEARCH_PLAN,
with the full [dimension audit](arithmetic_dimension_bounds_and_effective_high_dimension_search.md)
and [algorithm audit](vinberg_allcock_search_strategy_20261013.md) linked.
Allcock's published edgewalking algorithm works over integral lattices;
using spherical extensions and edge geometry without a fixed lattice is a
separate construction method requiring metric and global certification.
No noncompact arithmetic H22 host search is warranted.

## A single old-wall deletion permits a cusp cut

The [813-facet construction](fresh_h22_core_replacement_20261013.md) starts
from the 812-wall incidence branch: 210 old vertical walls and 602 parent
caps. It deletes old wall 1 and inserts two new norm-two walls. It does not
retain the 48 subsequent caps of the 860-wall branch.

At cusp orbit 14, the old active diagram has one affine component. Deleting
wall 1 leaves a positive-definite set I of 20 normals. Let delta be its
old cusp ray and p the projection of the old center h perpendicular to I.
The choice of cusp is structural: every irreducible affine component has
a strictly positive linear combination equal to a positive multiple of
delta. Retaining that entire component forces every acute-compatible new
normal to pair nonpositively with delta. To cut a cusp by such normals,
one must delete at least one wall from every affine component. A
single-component cusp permits the smallest possible old-wall deletion.

The exact data give

```text
p² = -2484,   p·delta = -46,
u = -p + 350 delta,   U = u² = 29716,
r+ = (u/sqrt(U),  1/sqrt(2)),
r- = (u/sqrt(U), -1/sqrt(2))
```

in the form `Q direct-sum <2>`. Both new normals have norm two; they are
mutually perpendicular and perpendicular to I. All other new pairings
are strictly dotted. The complete new-pair census is 41 right and 1,582
dotted relations. A strict common timelike center proves that all 813
inequalities are actual facets.

Their old projections have positive norm one and pair positively with delta,
so they cut the old cusp. This was impossible while retaining all 210 old
walls as acute-compatible facets. The canonical new vertex
`(p+296 delta,0)` has norm `-29716`; its 22 incident normals form a
positive-definite Gram matrix. Every other wall inequality is strict there.

A two-cycle generates `Q(sqrt(59432))`, while the retained rational short-wall
block has a negative quadratic value `-118`. Consequently a finite-volume
completion retaining these witnesses would be non-quasi-arithmetic. This
construction therefore respects the literature restriction on H22 arithmetic
examples.

The previous rational strict escape survives, with norm 1586. Thus the
construction has infinite volume, and no global cusp census or completed
edge queue is claimed. The old lower bound of 2,191,374 additional facets
depended on retaining the old core; its premise is absent here.

The [separate radical checker and review](fresh_h22_core_replacement_20261013_review.md)
verify the new norms, all 1,623 new pairings, strict center, actual finite
vertex, cusp cut, field witness and escape. The standalone producer also
reconstructs and checks the complete 812-wall parent's Gram data.

The [supplementary plane analysis](fresh_h22_core_replacement_20261013_planes.md)
archives the preceding 809-facet five-wall-deletion experiment and two
scoped continuous obstructions. These explain why the successful cusp cut
required leaving those restricted projection planes. Their exact replays
pass; they do not exclude general wall replacements.

## The fixed-core nonright search did not produce a wall

The [nonright endpoint census](fresh_h22_nonright_bridge_20261013.md)
exhausts 6,480 explicit spherical profiles on the 88 selected edge supports,
with 217 unbounded dihedral families kept explicit. The exact profile replay
and [source review](fresh_h22_nonright_bridge_20261013_review.md) pass.

A bounded numerical search sampled 9,084 profiles against all 305 selected
anchors with anchor labels 2 through 32. None of 85,889,220 parameter
proposals survived every Coxeter-angle gate. This is a numerical observation,
not an exact exclusion of these profiles or of every nonright bridge.
Further unrestricted work should change the core or the geometric family,
rather than interpret these finite failures as a theorem.

## H9: a compact hexagonal face, with a larger open-edge deficit

The [125-wall construction](fresh_h9_h5_closure_next_20261013.md) adds two
coupled endpoint caps to the 123-wall chamber. All 247 new pairs and the
complete 1,024-subset finite-neighbor affine census pass. A complete
six-facet, six-vertex boundary cycle certifies an actual compact H2 hexagon;
its intrinsic angles are not asserted to be Coxeter angles.

The H5 face now has 26 facets, 29 finite vertices, 46 finite edges and 53
open edges, compared with 46 open edges previously. Strict spacelike escapes
remain in H5 and H9. Closing this polygon is local geometric progress, not
an improvement in the total edge deficit or a compact H9 example.

The [160-digit replay and source review](fresh_h9_h5_closure_next_20261013_review.md)
reconstruct the two new caps using shared routines against rescaled frozen
predecessor enclosures. They pass, with that dependence explicit.

The next construction problem is global: compatible replacements that remove
several escapes while controlling the new edges and cusp links they create.
Both chamber dimensions remain open targets; neither certificate is a
finite-volume example.

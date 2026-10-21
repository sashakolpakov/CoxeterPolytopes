# Independent review: compact dodecahedral cutting and wall radius

20 October 2026. **PASS within the stated hypotheses and imported metric
inputs.** Independent written mathematical and source review only. No proof,
numerical, enumeration, Lean or typesetting job was run; the reviewed files
were not changed.

## Exact sources

* `results/compact_ra4_dodecahedral_double_spectrum_cutting_20261020.md`,
  SHA256 `08129e5aa206016d758114fda97df482c78de5397f7a3631bc542f79b0cbb34b`.
* `results/compact_ra4_farthest_dodecahedron_wall_radius_20261020.md`,
  SHA256 `d082e4a918d4504a59ac200f32bbb141f1f2cb52d701006eaea7ba5b0e5bcc22`.
* The capacity, metric-input and all-dodecahedral continuation passages in
  `manuscripts/compact_ra4_volume_constraints.tex`, SHA256
  `afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f`.

Both new notes were read in full. The conclusions below were derived from
their geometry and algebra, rather than inferred from another review.

## Cutting theorem

The hypothesis is that **every intrinsic facet is the regular dodecahedron
D or its specified reflection double D2**. It is a metric hypothesis, not a
consequence of counting pentagons and hexagons.

The doubled pentagon is intrinsically recovered: each hexagon has exactly
two opposite length-2L sides, their midpoint segment is the old gluing edge,
and these segments give the canonical regular pentagon. At a seam vertex,
the actual parent edge has three incident facets. The shared long-edge
condition forces all three pairwise ridges to be hexagonal. Thus the
pentagons join along their actual geometric edges with no boundary.

The common-span step is sufficient and nondegenerate. At the two endpoints
of a shared seam edge, the third-facet normals are nonadjacent side normals
of a regular right-angled pentagon, hence span a Lorentz two-plane. The two
orthogonal adjacent parent normals complete a Lorentz four-space. This is
both neighboring spaces U_F and U_G; its orthogonal line is spacelike.
Consequently the seam lies in a genuine common hyperplane, not merely an
abstract closed surface.

The section argument closes the main recognition risk: the union of
canonical pentagons is open as well as closed in the boundary of the
compact three-dimensional section. Interior, edge and vertex strata have
exactly the listed one, two and three parent walls; no omitted wall can
pass through them. Connectedness gives the **whole** section. Regular
pentagon continuation identifies it with D. An outside wall cannot meet
the section's interior either, since its nonpositive affine functional
would then vanish identically on the section plane, contradicting strict
inequalities on a pentagon interior.

The cuts therefore introduce D facets, replace exactly twelve D2 facets
by their D halves, and leave every other component strictly on one side.
Later cuts remain disjoint. The all-D final-piece recognition uses unique
orthogonal continuation of whole regular dodecahedra, in both adjacency
graphs, rather than a face-count recognition. It gives N_D2=12(m-1), a
tree of m regular 120-cells, and the unique reflection double when m=2.
The displayed counts check: D facets 94m+26, total facets 106m+14,
vertices 560m+40 and W=30(m-1). In the volume formula, f denotes facets.

## Farthest-wall theorem

The maximizing projection is an actual relative interior point: for
every other actual wall m, both terms in
`<q+Tn,m>=<q,m>+T<n,m>` have the required strict/nonpositive signs.
Each adjacent intrinsic wall height is at most T/sqrt(1+T^2)<1, giving
distinct S2 directions with pair products strictly below 1/2.

I independently expanded the quartic as
`x^4+(13/10)x^3+(3/50)x^2-(8/25)x-2/25` and recovered all five stated
Legendre coefficients. The unnormalized rank-four trace is
`u^2-(7u^2+I)/7+5I/35=0`; contraction gives
`x^4-6x^2/7+3/35`. Thus the positive-kernel sum and N<=14 are valid.
At N=14, every off-diagonal product must be -1 or -2/5; the squared norm
of the sum is then negative. The strict code therefore has N<=13.

The imported polygon capacity bound excludes any hexagonal or larger
face at that size. Euler then forces twelve pentagons. Identification
with regular D uses the explicitly isolated twelve-face metric theorem;
the quartic alone is not a metric classification.

For D, h^2=phi/2 follows from adjacent icosahedral directions. Opposite
walls give maximum intrinsic height at least ah, where a>=1 is the time
coordinate relative to the center. Combining with the projection bound
gives T^2>=phi/(2-phi)=phi^3. Equality forces the projected point to be the
center and **every** adjacent parent height to equal the same T. Applying
the same argument at these maximizing facets propagates all-D status
through the connected facet graph. The continuation theorem then gives
the regular 120-cell; its matched centered foot and normal distance fix q
as its center. Conversely symmetry attains this value there.

## Scope retained

The first theorem proves a volume spectrum and unique first double only
inside the D/D2 facet class. The second theorem is a universal statement
about maximum supporting-wall distance, not minimum wall distance. A
nonmaximizing non-D facet is not thereby redundant. Neither theorem
forces the D/D2 spectrum from a volume bound, supplies an inscribed
120-cell, or proves the unrestricted compact minimum-volume conjecture.
The radius theorem retains the stated external twelve-face metric input;
the cutting theorem uses the given regular D and 120-cell and geometric
continuation. Worldwide novelty and primary-source comparison are outside
this written review.

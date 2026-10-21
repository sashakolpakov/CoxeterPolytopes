# Finite-edge extensions, cusp repair and compact tree incidences

**No new finite-volume Coxeter polytope or dimension is established.** The H20 construction now has 95 actual facets, at least 116 verified finite vertices, and three specified bounded full cusp links. An exact strict spacelike escape still proves infinite volume. The compact H9 search has produced further exact exclusions and unresolved incidence candidates.

This continues the [contact-search checkpoint](e7_contact_search_and_structured_h9_2026_10_14.md), preserving its frozen artifacts. The [literature restrictions](arithmetic_dimension_bounds_and_effective_high_dimension_search.md) remain decisive: noncompact arithmetic H20 and dimensions above 21 are excluded. H20 is an unfilled dimension below the known noncompact maximum 21; compact H9 would improve the known compact maximum eight. The searches here do not impose arithmeticity.

## H20: actual edge walking beyond the completed face

The [two-cap construction](fresh_h20_e7_finite_edgewalk_20261014.md) retains the previous 88 walls and adds two mutually orthogonal norm-two normals with nonright contacts along the shared edge walls. They lie outside the earlier all-right transverse family. They close 31 outgoing edges from the 44 known face vertices and create exactly 40 finite vertices incident to the new caps, including nine on both caps.

The [independent review](fresh_h20_e7_finite_edgewalk_20261014_review.md) reconstructs both formulas with separate sparse radical arithmetic. It checks all 177 new wall pairs, all 90 actual facets, and 84 distinct finite points. Its independent finite-support enumeration uses exhaustive subsets and FLINT principal determinants, reproducing 3,080 candidates with exactly 40 positive-definite supports. There are 145 verified finite edges with a new endpoint. Of the 800 directed germs at the 40 new vertices, 259 have a second finite vertex and 541 do not. Ideal endpoints of those 541 were not all classified.

The new parallel pair `(9,88)` introduces a deficient cusp. Its 20 incident normals have horospherical rank 18, and an explicit tangent vector of square `98/27` proves an unbounded direction. This defect is then repaired in the [five-wall extension](fresh_h20_e7_edge_cusp_repair_20261014.md). The added quotient roots complete the finite A5, A1 and B3 factors and bound the remaining lineality direction. Their shifts are chosen at finite-angle contacts with existing nonincident walls.

The resulting 95-wall chamber preserves all 84 previous vertices and both earlier bounded cusps. The newly repaired cusp has 25 incident facets in six affine components of ranks `8,5,3,1,1,1`, totaling 19. A bounded continuation of the 541 earlier germs finds twelve distinct new finite endpoints. Six of these also have actual edges to the now-bounded cusp, giving finite–finite–proper-ideal paths. The [separate repair review](fresh_h20_e7_edge_cusp_repair_20261014_review.md) verifies the exact five-normal construction, all 460 new pairs, the complete cusp link, the twelve accepted endpoints and the six ideal edges. It does not independently re-enumerate the whole 541-germ sweep or assert a complete 95-wall vertex census.

The [complete edge census at the repaired cusp](fresh_h20_e7_repaired_cusp_walk_20261014.md) subsequently examines all 1,728 supports, one omission per affine component. It finds 26 distinct finite endpoints, including twenty additional vertices, sixteen ideal contact points whose links are unbounded, and 1,686 open opposite ends. The [independent review](fresh_h20_e7_repaired_cusp_walk_20261014_review.md) reconstructs all 25 component projections, checks the finite and ideal endpoint geometry, verifies explicit recessions at all sixteen ideal contacts, and proves 128,136 strict inequalities at the open null ends using outward rational bounds. Thus there are at least 116 verified finite vertices. The complete star at this cusp does not supply a complete global vertex census.

The same rational escape survives both extensions:

```text
E² = 259854917001622621 / 62500000000000000 > 0.
```

It is strict at every wall and has the correct time orientation relative to the common strict timelike center. Thus the chamber contains an open set of ideal directions and has infinite volume. The additional corners and bounded links are starting points for further search, not a finite-volume certificate.

## Compact H9: different cut trees and exact metric obstructions

The [tree-incidence batch](fresh_h9_tree_cube_wedges_20261014.md) generates 240 star and T-tree cut/wedge recipes and reduces them, with literal isomorphism witnesses, to 33 rationally realized nonpath incidence classes. Separate primal checks verify all 14,988 vertices, actual facets, full minimal-nonface data and closed connected edge graphs.

Complete necessary label searches exclude four star classes, accounting for 96 recipes. Six classes, accounting for 52 recipes, have partial-label survivors. The other 23 T classes, accounting for 92 recipes, reached the stated search limit and remain unresolved. These limits are not exclusions, and surviving label assignments are not metric realizations.

The [first star assignment](fresh_h9_star_cube_metric_20261014.md) fails a fully exact rank obstruction. Four positive B2 blocks have total rank eight, leaving at most two residual dimensions in H9. Explicit actual-vertex tests reduce the remaining relevant contacts to 1,372 cases. Every case forces a nonzero three-dimensional residual determinant, hence total rank at least eleven. The [independent review](fresh_h9_star_cube_metric_20261014_review.md) checks every contact-domain guard and determinant.

A larger [cut-contact family](fresh_h9_star_cube_metric_20261014_cut_profiles_review.md), fixing those B2 blocks and the positive leaves' two labels at three, has 6,300 residual cases. Of these, 6,296 fail rank; the four zero determinants force affine A3 diagrams inside actual finite vertices. Its [complete cut-profile extension](fresh_h9_star_cube_metric_20261014_all_cut_profiles.md) also handles own-right contacts, with the necessary angle bounds re-established for those cases. It has exactly eleven ordered profiles per cut and 8,092 relevant residual cases: 8,088 rank contradictions and the four affine-vertex failures. Thus the entire branch with these B2 and positive-leaf labels is excluded, while the cut angles are unrestricted. Changed positive-leaf or B2 labels remain open; the next saved necessary assignment changes one positive-leaf contact from three to four.

The [sixteen-orientation extension](fresh_h9_star_cube_metric_20261014_orientations.md) permits either of the two leaves at each B2 block to be the fixed `(3,3)` leaf. It rebuilds actual-vertex bounds for every choice, without assuming a facet-swap symmetry. The separate review checks all 91,392 residual occurrences, including 64 rank-zero occurrences rejected by finite-vertex minors. The next assignment outside these sixteen families is separately excluded by 28 residual cases. The entire star incidence remains open.

## Next work and verification scope

The H20 census now provides explicit open supports for a directed general spherical-extension search. Prioritize caps that connect existing geometry and cut the surviving escape. Repeated local cusp repairs alone cannot substitute for global finite volume. Any bounded numerical screen must retain its stated scope, including untested large-label families.

The compact search continues through remaining star assignments and T-tree candidates, applying proved metric gates only on their stated domains. Arbitrary large dihedral labels must remain available whenever no actual-vertex bound excludes them.

The [checkpoint audit](edgewalking_cusp_repair_and_tree_incidence_2026_10_14_audit.json) binds the completed files, source manifests and local links while preserving the preceding archive. It is an integrity check; mathematical verification scopes are stated in the individual reports.

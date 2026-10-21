# Complete source census, changed angle field, and global covering constraints

**No new compact or finite-volume Coxeter example or dimension is established.** The new exact constructions still have strict spacelike escapes. This checkpoint completes the known H21 source's corner census, replaces the H22 π/7 branch by a verified multiquadratic π/6 construction, and closes another actual edge in the H9 construction. The global covering analysis explains why local corner and cusp repairs are insufficient.

## Complete H21 source census

The [exact corner walk](fresh_h21_source_corner_walk_20261013.md) exhausts the known 210-wall Borcherds chamber modulo its recorded 241,920-element symmetry group. It reaches 377 finite-corner orbits and 18 cusp orbits, closes all 7,528 directed edge orbits, and leaves an empty queue. Orbit–stabilizer gives 70,235,026 finite vertices, 1,095,990 cusps, and 834,554,238 undirected edges. A separate rational checker verifies every group word, complete active set, edge support and endpoint.

The earlier 15 sampled cusp orbits missed three types: D4^5, A1^3+A17, and D7+A11. The [independent lattice and literature audit](fresh_h22_extension_feasibility_20261013.md) matches the complete 18-type list to the published classification. This is an exhaustive computation of a known chamber, not a new polytope. Its completeness applies to the old H21 source; checking only those representatives does not cover an asymmetric H22 extension's entire boundary.

## One-sided construction and changed angle field

The [one-side reduction](fresh_h22_extension_feasibility_20261013.md) shows that it suffices to seek a coherent upper cap family covering the entire old chamber. Reflecting it automatically gives allowable cross-pairs. Its [862-facet π/7 application](fresh_h22_one_side_20261013.md) retains all preceding 839 facets, two bounded cusp links and a compact quadrilateral, with a reflected second copy. Its exact escape survives.

A finite-volume completion retaining the π/7 pair and old source subgroup would have to be non-pseudo-arithmetic. The [860-facet π/6 replacement](fresh_h22_multiquadratic_20261013.md) avoids that cubic-field obstruction. It changes the cap recipes, discards one provably forbidden angle, and reconstructs the coupled endpoint caps in an explicit multiquadratic field. Four actual finite vertices again bound a compact quadrilateral. A separate exact squarefree-radical implementation verifies all new norms, pairings, vertices and the strict escape.

The π/6 branch still cannot have a quasi-arithmetic finite-volume completion retaining its field witnesses: one Gram cyclic product is `208+64sqrt(3)`, whereas a retained rational principal block remains indefinite under every embedding. Absence of the cubic obstruction does not certify pseudo-arithmeticity or completion. The displayed integer escape has squared norm 1,586 and is strict on all 860 facets.

## Inherited cusps impose a large retained-core lower bound

The [exact inherited-cusp footprint](fresh_h22_inherited_cusp_footprint_20261013.md) establishes that only 303 of the 1,095,990 old cusps have bounded links in the frozen 860-facet chamber. They are the 301 rational parent null directions and two further null directions. The other 44 signed new normals have timelike old projections and meet no inherited cusp.

Any new acute-compatible normal has a future causal old projection. It can be incident to an inherited cusp only when that projection is proportional to the cusp ray; pure vertical normals conflict with retained caps of the same height sign. Each of the 1,095,687 missing cusp links therefore needs two new facets of opposite height signs. A finite-volume completion retaining all 860 actual facets requires at least **2,191,374 additional facets**, or **2,192,234 in total**. The argument permits arbitrary real coordinates and arbitrary angle labels.

This is a conditional lower bound on retaining this core, not a nonexistence theorem or a bound on arbitrary H22 polytopes. Replacing retained walls can evade its assumptions. Repairing all inherited cusps would still leave the whole-chamber covering requirement.

## Exact covering tests and necessary finite-angle bridges

The [covering-cell criterion](fresh_h22_covering_cells_20261013.md) subdivides the old cone by which cap's normalized height is smallest. Coverage is equivalent to checking a quadratic inequality at every extreme ray of every cell. The proof uses the reverse Lorentz triangle inequality and includes null rays and lower-dimensional cells. It gives a finite certificate protocol; the complete H21 subdivision has not been computed.

A cheaper necessary test already detects 19 disconnected components in the 325 retained upper caps' finite-angle contact graph: one component of size 305, one triangle and 17 singletons. Any completion retaining them needs finite-angle bridges. More caps dotted or parallel to every retained upper cap cannot finish the construction.

An exact nine-facet H3 counterexample demonstrates why covering all original finite vertices and bounding all original cusp links is insufficient. Its two ideal links are bounded, yet it has eight spacelike cone vertices and a strict escape. Its complete minimum-cell census finds four uncovered new rays. The [source/proof and replay review](fresh_h22_covering_cells_20261013_review.md) passes and independently reproduces the H22 graph.

The first targeted [finite-contact bridge search](fresh_h22_contact_bridge_20261013.md) now has an exact negative result over all real positions and all integer angle labels within its stated scope. It tests positive-height norm-two caps perpendicular to each of the 88 edge supports at the four verified quadrilateral vertices, seeking finite or parallel contact with the large 305-cap component. Ten edge families fail continuous acuteness outright. For the 301 rational parent anchors, 23,003 cases have empty allowable intervals and 475 force another normalized cosine strictly between 0 and 1/2, where no Coxeter angle exists. A separate exact extension proves all 312 intervals for the four other cap types empty, completing the test against the whole component.

This excludes direct bridges from those right-angle endpoint families, not general bridges. Nonright edge profiles, other vertices, indirect chains of new caps, and wall replacements remain outside the theorem. The [independent source review and replay](fresh_h22_contact_bridge_20261013_review.md) corrected a normalization error in the exploratory filter before the final exact certificate was produced; the corrected equations and saved replay pass.

## Coupled edge closure in H9

The [123-facet H9 construction](fresh_h9_h5_closure_20261013.md) retains all 118 preceding walls and adds five compatible caps. Two impose new π/5 and π/4 couplings. The last closes an actual H5 edge and supplies one finite endpoint. Its H5 face now has 24 facets, 26 finite vertices, 42 finite edges and 46 open rays.

All 600 new wall pairs and all 2,560 subsets of each new cap's finite neighbors pass exact checks, with no new affine diagram. The full 24-choose-5 face census and explicit vertex inequalities verify the reported vertices. Both H9 and the H5 face have strict spacelike escapes. The [scoped review](fresh_h9_h5_closure_20261013_review.md) reproduces the new formulas and gates at 160-digit working precision while preserving the rigorously enclosed predecessor data.

## Algorithm interpretation

The [algorithm audit](vinberg_allcock_search_strategy_20261013.md) remains the governing distinction: Vinberg enumerates mirrors by distance, while Allcock traverses chamber edges. Allcock's published algorithm is itself arithmetic, explicitly for rational lattices over Z. The broader construction uses its geometric idea with additional metric variables and independent whole-chamber verification. The completed source walk illustrates exact edge closure for a fixed known wall set; it does not supply a complete search through all possible new walls.

The research objective remains open. Further additions must address whole-chamber coverage and compatible finite-angle connections, not merely increase the wall count or the number of locally bounded links.

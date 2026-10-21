# Fixed-point symmetry, finite-edge replacements and incidence exhaustion

**No new finite-volume Coxeter polytope or dimension is established.** This checkpoint records a useful full-chamber symmetry in the compact-H9 arithmetic search, a changed-support H20 construction with actual finite edges, a continuous-family obstruction, and completion work on the remaining compact-H9 incidence cases.

## Compact H9: a chamber symmetry beyond the enumerated corner radius

The [cubic-field report](fresh_h9_cubic_corner_20261014.md) validates a lattice involution through a nearby fixed face. The fixed point has exact distance height `2+7θ+3θ²<16`, while the corner it transports lies at height about 1009.44. A complete height-100 Vinberg prefix therefore suffices to recognize the fixed sector and certify the full-chamber symmetry. This avoids an unjustified claim that the distant corner was already known to belong to the full chamber.

The symmetry first gives 22 actual walls. Two specifically selected additional Vinberg levels, up to height 141, give a 23-wall prefix and then 32 walls after applying the involution. An exact complete elliptic-face census finds 24 finite prefix vertices, sixteen certified full-chamber corners, and 200 rank-eight edges. Twenty-eight edges have two finite endpoints. The saved strict spacelike escape proves infinite volume of the 32-wall prefix. Neither reflectivity nor nonreflectivity of the full host is decided.

The [independent arithmetic review](fresh_h9_cubic_corner_20261014_review.md) verifies the involution, fixed point, oriented local sector and all pairs in the first extension. The two new Vinberg levels were fully replayed. A [separate packet and orbit check](fresh_h9_cubic_corner_20261014_h141_review.md) verifies all 496 final pairs and reviews the continuation logic. The final escape signs were checked with two field implementations. No independent exhaustive root census is claimed.

## Nonarithmetic H20: move the endpoint by replacing walls

The [contact-replacement report](fresh_h20_contact_closure_20261014.md) builds two exact 76-facet chambers from the retained sixty-wall H18 lift. The first introduces rational pi/3 anchor contacts and eight finite edges. Its full local walk at sixteen specified vertices still has 304 open outgoing germs, and the previous strict escape survives.

The second removes all eight large orthogonal mates and replaces them with eight norm-two caps using `sqrt(79)`. It has eight actual finite edges joining `D13+D7` corners to `A12+E8` corners. All 572 new pairs and all 76 actual facets are certified exactly. The three previously bounded cusp links persist. Both distinct height-sign images of the earlier escape are cut, but a new strict spacelike witness has square

```text
269344645863312077/40000000000000000 > 0.
```

Thus the final chamber still has infinite volume. This is a new local construction in the search, not a new finite-volume example. Keeping the old sixty facets preserves the conditional non-quasi-arithmetic obstruction: any finite-volume completion retaining them would avoid the prohibited arithmetic H20 route.

## Complete an E8 face using a known H12 chamber

At one `A12+E8` corner, the [exact E8 projection calculation](fresh_h20_contact_closure_20261014_e8_face.json) identifies a face with precisely twelve current walls, forming `A12`. Every other nonincident projected normal is future timelike and imposes no further face boundary. The intrinsic form simplifies to `Cartan(A12) + <-1/13>`, although its embedding in H20 uses `sqrt(79)`.

There is also a useful all-cap restriction. The E8 inverse Gram has positive entries and minimum diagonal two. If an acute Coxeter cap of norm q is not orthogonal to E8, some squared pairing is at least `q/2`. Its projection to the E8 span consequently has square at least q; its residual normal on the H12 face is nonspacelike. Thus every new wall actually bounding that face must be orthogonal to all eight E8 normals. The face needs a finite-volume completion if the ambient chamber is ever to have finite volume; no compactness requirement or global nonexistence follows.

The [separate face review](fresh_h20_contact_closure_20261014_e8_face_review.md) verifies the projection arithmetic, orientation argument, 14-root isometry, and precise scope of these statements. For an arbitrary nonempty acute extension, the orientation proof projects its own strict center; it does not assume that the extension retains the predecessor's center.

The [standalone H12 source](fresh_h12_a12_source_20261014.md) reconstructs the known 14-facet chamber of `I_(12,1)`. Exact enumeration and a complete projective vertex census certify 35 finite vertices, two bounded cusp links and an actual A12 corner. With `T²=-1/13`, its two missing normals are `-omega3+2T` of norm two and `-omega11+3T` of norm one. Both orientations of the A12 chain were transported and tested against every H20 wall.

The norm-one caps pass, but the norm-two caps have forbidden angles with retained walls. In the reversed orientation there is exactly one conflict, with original wall zero. Removing that wall and adjoining the reversed pair gives the [77-wall construction](fresh_h20_contact_closure_20261014_h12_replacement.json). All 2,926 pairs and all 77 actual facets pass exact checks. Its E8 face is isometric to the complete H12 source, and all 35 finite face vertices are actual finite H20 vertices. Its two face cusps initially retain an unbounded E8 direction in their ambient links.

The old strict global escape survives, so the 77-wall chamber remains infinite-volume. Fresh retained pairs verify the same irrational cycle `12+8sqrt(2)` and rational indefinite block of determinant minus 32 after the deletion; any finite-volume completion retaining these witnesses would still be non-quasi-arithmetic. The old three cusp-link claims from the sixty-wall lift are not automatically inherited after deleting a wall.

The [two-cusp completion](fresh_h20_e8_cusp_completion_20261014.md) then adds exactly two walls. For each of the two intrinsic cusp rays nu, take `-omega+8nu`, where omega is the E8 fundamental weight of square two. Each has norm two, one pi/3 contact with E8, and right contacts with its other incident cusp walls. All 155 new pairs pass: 39 right, two pi/3 and 114 dotted. The final chamber has **79 actual facets**, retains the complete H12 face and all 35 certified finite vertices on it, and has bounded full H20 links at these two cusps. The two link decompositions have affine block sizes `4+9+9` and `12+9`; positive null relations and positive principal minors certify the full quotient rank nineteen.

The same strict spacelike escape remains. Thus even this complete face and these two repaired ambient cusps do not establish finite volume. The final 79-wall certificate explicitly proves infinite volume; no further cusp census or exhaustive vertex census is claimed.

The [separate cusp-completion review](fresh_h20_e8_cusp_completion_20261014_review.md) recomputes the two cap formulas, all 155 new pairs, every norm/center/escape sign, both positive affine-link decompositions and the 70 new vertex-preservation inequalities. Actual facet geometry and the inherited complete H12 source remain explicit checked dependencies.

After the fourteen complete H12 face walls are retained, no further spacelike cap orthogonal to E8 can be acute to all of them. Its transported vector would lie in the source's homogeneous cone: the B12 sector and first positive-height wall force nonnegative time, with time zero possible only for the zero vector; the certified time-one slice is contained in the closed future light cone. Positive norm is impossible. Further compatible caps must therefore have a nonzero E8 component. This closes the pure intrinsic-face addition route, not the ambient completion problem.

## Stop the old paired-curve search

The [continuous obstruction](fresh_h20_structured_extension_20261014.md) proves that all admissible complete pairs on the four original fixed-orientation support curves, at every real parameter and under all four height reflections, retain one common strict escape. It permits compatible deletion and replacement of the old optional pairs. All 120 pairs of different curve families stay dotted on their necessary parameter tails.

The proof uses the gap between a right angle and pi/3 to eliminate acute but non-Coxeter candidates; acuteness alone would not suffice. The [separate review](fresh_h20_structured_extension_20261014_review.md) checks the rational frames and symbolic identities and records a fresh primary replay. The theorem does not cover arbitrary unpaired caps, rotations in the support complement, changed supports or replacement of the old sixty walls. The new H20 construction changes support and lies outside it.

## Fifteen-facet incidence search

The [exact completion report](fresh_h9_fifteen_completion_20261014.md), [primary pass](fresh_h9_fifteen_completion_20261014.json) and [separate-arithmetic audit](fresh_h9_fifteen_completion_20261014_independent.json) close all 1,334 cases left open in the [previous checkpoint](finite_ridge_corners_and_fifteen_facet_search_2026_10_14.md). Every case terminates with an exact exclusion: 15,263 states, no node-limited cases, no full category assignments. The separate arithmetic replay agrees on all 1,334 cases. These cases close before metric variables or modular-minor tests are reached.

Combined with the earlier exclusions, all 45,679 distinct fifteen-facet nerves from a codimension-two-through-five cut of the 770 saved Gale parents are excluded. A cut of codimension at least six creates a minimal nonface of that order, impossible for a compact Coxeter polytope because Lannér diagrams have at most five nodes. Thus every standard nontrivial single-face truncation of these specific parents is excluded. This does not exclude all compact-H9 polytopes with fifteen facets, arbitrary larger hyperplane cuts or later surgeries. All finite labels and missing triangles were allowed. The replay uses separate finite-diagram arithmetic while sharing search control and the inherited relation tables; it is not a wholly independent combinatorial enumerator.

A concrete unsearched alternative is direct rank-five Gale data for fifteen facets, with exact positive-circuit realization and a filter excluding the vertex-link structure forced by a standard face truncation. This changes the incidence family rather than revisiting a proved empty catalogue.

## Literature constraints and next work

The [dimension and algorithm check](arithmetic_dimension_bounds_and_effective_high_dimension_search.md) remains in force. Noncompact arithmetic reflection groups occur only in dimensions at most nineteen or exactly twenty-one; all arithmetic and all compact reflection groups have the different upper bound twenty-nine. The arithmetic searches here concern compact H9, while H20 requires a nonarithmetic completion. Allcock's published edgewalking operates on rational Lorentzian lattices; the nonarithmetic corner constructions here require their own metric and global checks.

The next H20 step must target the explicit remaining global escape using compatible caps with a nonzero E8 component. The complete H12 face and the two full cusp links now provide firmer constraints than an isolated vertex. The fixed-orientation paired curves and further purely intrinsic H12 additions are closed. In compact H9, use the certified corners and open edges or seek incidence families outside the completed single-cut catalogue. A bounded prefix, a local endpoint or a cut escape is never a finite-volume certificate.

The [checkpoint auditor](../scripts/audit_fixed_point_and_contact_checkpoint_20261014.py) verifies local links, recursive declared source hashes and nested per-case record hashes. Its [saved archive certificate](fixed_point_symmetry_and_contact_replacements_2026_10_14_audit.json) concerns provenance; the mathematical verification scopes remain those of the individual reports. Earlier checkpoint files and their auditors remain unchanged.

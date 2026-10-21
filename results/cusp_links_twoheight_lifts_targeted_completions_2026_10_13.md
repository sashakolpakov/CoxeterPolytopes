# New-dimension search: cusp links, changed lifts, and targeted completion

**No new compact or finite-volume Coxeter dimension has been established.** This round tests new constructions and resolves a previously unfinished finite cap dictionary. All positive root and angle claims below have exact replay; every tested complete-looking candidate retains an exact obstruction to the required volume or Coxeter angles.

## A changed geometric H9 seed

The [two-height construction](fresh_h8_twoheight_lift_20261013.md) deletes original wall 3 of Bugaenko's H8 chamber before lifting walls 0 and 9 into an additional positive coordinate. Their opposite heights have squares `2φ` and `(2+8φ)/11`. After normalization, edges 0–1 and 0–9 have label three, while the squared dotted pairing 9–10 becomes `3+2φ`. This is outside the earlier obstruction that retained all eleven old directions.

The resulting ten-wall Gram has signature `(9,1)`. Rescaling four normals puts it over `Q(φ)`; its conjugate also has signature `(9,1)`, so this is not an admissible arithmetic host. The first dictionary gives 62 caps. A targeted free-distance norm equation adds 55 distinct mirrors, each cutting a witness that survived all initial caps. Seventy-one tested compatible configurations have 40–51 facets and exact infinite-volume witnesses.

The final, stronger certificate is one strict spacelike escape surviving **all 117 mirrors simultaneously**. Thus this finite dictionary cannot give a finite-volume completion, regardless of subset selection. It does not exclude further mirrors or the unrestricted new seed. A standalone Fraction replay checks the lift, both signatures, all roots, selected angles, mirror distinctions, all volume witnesses and two deliberate corruptions. Another agent reviewed the construction and independently reran the final verifier.

## Resolving the old dictionary and preserving a partial configuration

The [completion search](fresh_h8_cap_completion_20261013.md) changes the selection strategy: multiple exact escape constraints are extracted from candidate affine-section vertices, with complete cutter lists across all 1035 saved cap entries. It then prioritizes coverage of the most restrictive gaps.

Thirteen exact witness clauses reduce compact completion to 147 cutter nodes. All 10,731 pair relations are classified with a complete safe finite-angle list; a 45-state Boolean certificate rules out any compatible selection meeting all clauses. This proves that **no subset of the 1035-entry bank gives a compact Coxeter completion**, even with the polar cap optional. The proof accounts for redundant optional inequalities. It is compact-only: allowing parallel mirrors leaves the finite-volume covering problem satisfiable. One such 22-wall configuration was tested and escaped, without closing that larger problem.

The bottleneck then generates 49 cap entries outside the closed dictionary. Nine are exactly compatible with every wall of a fixed 30-facet configuration. All eight maximal extensions from those nine are tested, with 31–32 facets; each still escapes. This gives coherent new extensions rather than replacing the entire partial configuration. The expanded 1084-entry bank and further caps remain unresolved. Producer-free exact replays cover both the compact obstruction and all eight extensions.

## A changed arithmetic coefficient and saturated spatial searches

The [prime-five continuation](fresh_cubic9_positive_20261013.md) searches `⟨−T⟩+B²+A4+A1`, with `T=2cos(π/9)` and `B=[[2,−T],[-T,2]]`. Its rank-ten determinant squareclass distinguishes it from the preceding coefficient `2T` and from the four earlier `B⁴+⟨1⟩` hosts. This comparison concerns the explicitly searched spaces, not worldwide priority.

Full bounded lower-dimensional faces suggest A4 directions for completion. Starting with the complete 13-wall H9 height-30 prefix, the search frees two independent A4 Cartan covectors and the A1 scalar, retaining the exact mod-five integrality condition. At height 1000 these six saturated slices produce **nine additional mutually compatible walls**. Several remove the earlier exact escape directions. The resulting 22-wall cone nevertheless retains a strict A4/scalar escape.

Independent Fraction implementations repeat all three initial bounded censuses and the targeted slice censuses, then verify every angle and the surviving escape. This does not prove nonreflectivity of the full host. No larger uniform H9 height cutoff was used to obtain the additional walls.

## Entire Coxeter cusp links in H20

The [cusp-adapted sections](fresh_h20_new_route_20261013.md) enforce a full Euclidean Coxeter cusp link before solving global angle conditions. Cutting a bridge in an affine component can split it into two affine Coxeter components; one verified example replaces an affine D7 component with two affine B3 factors. The normal moves in a pencil `e(t)=e₀+t u`, where `u` is the original cusp vector; the entire rank-19 cusp link stays constant.

Eight declared pencils and their target-angle equations have 46 real algebraic roots. Thirty-five give genuine finite-volume H20 sections with a full Coxeter cusp. Eight also have a certified repaired global ridge. **Every one of the 35 has another actual forbidden ridge.** The remaining roots are excluded by exact cusp-cone identities or incorrect angle branches.

Twenty-seven irrational normals are proved to require at least five original walls, putting these sections outside the preceding four-wall work. The proof tests rational normal planes against all spans of at most four walls, using a complete modular quotient-rank filter and exact rational checks of its residual candidates. A standalone Fraction replay checks the complete root lists, cusp links, interiors, ridge witnesses, support bounds and deliberate corruptions. This is a bounded construction result, not a general exclusion of cusp-adapted sections.

The goal remains open. Surviving directions include new caps on the changed two-height seed, further compatible extensions outside the 1035-entry bank, roots outside the tested arithmetic slices, and other cusp parameters or link decompositions. The certificates identify precisely which finite constructions should not be repeated.

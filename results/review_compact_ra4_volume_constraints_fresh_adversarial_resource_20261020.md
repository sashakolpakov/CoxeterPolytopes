# Compact RA4 constraints: fresh adversarial source audit

20 October 2026. **PASS at the stated written-geometric and published-input boundary. No missing implication or countercase found.** This is a direct read of the complete manuscript and supporting proofs, not a deduction from earlier PASS reports. I have previously contributed to compact RA4 supporting research; this review is fresh re-derivation, not a claim of independent discovery or external refereeing. No scientific script, project import, Lean or TeX job was run. Only this review file was written.

Reviewed manuscript: `manuscripts/compact_ra4_volume_constraints.tex`, SHA256 `afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f` (1188 lines, 53508 bytes). Its reproduction map is `docs/compact_ra4_volume_constraints_reproduce_20261019.md`, SHA256 `d4a170d3c7b35a16d5ec1086dcd449fba5ab5d8eb8e83a657e7b44b38bdc26a6`. Existing build/page claims were not re-executed.

## Geometry, capacity and integer accounting

The positive-inverse projection is valid with the displayed signs: from a strict future point it adds nonnegative outward-normal multiples, keeps every unselected inequality strict, lowers the squared norm, and preserves future orientation. It excludes products in `(-1,0)`. At product `-1`, the nonzero feasible null sum contradicts compactness. Mutually orthogonal actual walls therefore give actual faces; finite vertices are simple. These facts justify intrinsic right angles and every later graph-to-face incidence. They are not transfers to arbitrary metric faces of a general Coxeter polytope.

For an outgoing edge, `r_x=sinh(L_x)x-cosh(L_x)a` has the correct orientation. Projection onto `a`-perpendicular proves distinctness and the exact neighbor-incidence equivalence. Any competing actual wall along that perpendicular would have positive product `cosh(L-L_x)` with the first wall. The footprint clique follows from the positive-definite orthogonal complement of the disjoint normal pair.

In the polygon dual-disc proof, the known neighbor links force the indicated outer cycle. No-extra and one-extra-touching cases leave a triangulated polygon with an impermissible degree-two ear. In the equality case, a chord cuts off a disc avoiding the sole interior vertex; an ear away from the chord endpoints again has total degree four. This proves the claimed Löbell lattice rather than assuming it from a census. In dimension four, an extra footprint of size three would already be an outgoing opposite. A neighbor with exactly one extra forces that extra to have singleton footprint. The resulting incidence count gives at least `M` extras and `f>=2V+5`.

Recomputed `e=2v`, `r=v+f`, `sum k=6v`, `W=v-5f`, `d(F)=(V-20)/2=M-12`, and `sum d=2W`. The orbifold face sum and four-dimensional Gauss–Bonnet give the displayed volume formula. Neighbor capacity yields `2W>=3d(F)` by summing distinct **facets**; shared ridges need not be disjoint. Facet `V` is even and `W` integral. At `V=48` the lower bound is `101+21=122`; at `V=98` it is `201+59=260`. Thus the 46/96 cutoffs follow with exactly the stated weak/strict volume hypotheses. The comparison double has `(v,f,W)=(1160,226,30)` and `f+W=256`; no minimum property is used.

## Metric and graph endpoints

At degree two the edge-weight inequalities force equal weights and equality in polygon capacity. Marked uniqueness realizes the Löbell rotation/reflection; its bases are regular. The regular-base formula and pentagon identities give the L6 lateral lengths with coshes `(2,sqrt(2),sqrt(3),sqrt(3),sqrt(2))`. Its unique longest edge identifies the same third parent facet on both sides. Simplicity rules out an additional parent along that edge. In the excess-twenty proof the two rows/bases are distinct selected facet sets: baseline30 plus at least5 from each row gives `2W>=40`, even when some lateral neighbors have larger excess or other positive facets occur elsewhere.

The three-consecutive-`L` argument uses actual orthogonal closing walls and projection, so an additional wall cannot truncate their vertex. The two-adjacent-edge pentagon and four-fixed-edge hexagon arguments are valid; solving the latter relations gives `cosh a=cosh d=phi^3=cosh(2L)`. Common carriers and their adjacencies are actual edges/vertices by the orthogonal-wall lemma. This supplies the marked-degree restrictions used in the full five-neighbor split: universal vertex, triangle with/without its outside edge, forest, and C5. In the last two cases the unmarked pentagons have enough `L`-edges to force regularity, contradicting a shared `2L` edge.

Re-derived the four-vertex eligibility bound, nine-vertex complement lemma, and ten-vertex lemma. In the last lemma the two exceptional nonadjacent degree-four vertices leave common-neighbor size `t=2,3,4`; the stated spare-degree arguments cover all three cases. It is essential that `H` is **full** positive-facet adjacency: absence of an induced square comes from two orthogonal Lorentzian normal planes. It is not inferred for `G` after deleting pentagonal diagonals. For `0<W<=16`, the degree-two result forces minimum G-degree3; the nine-vertex obstruction forces exactly10 positive facets, and the sum `2W<=32` leaves at most two of excess above3. The five-neighbor exclusion then gives precisely the hypotheses of the ten-vertex complement lemma. The final contradiction is valid.

## External inputs and scope

Checked the primary [Alexandrov–Bogachev–Vesnin–Egorov paper, §1.2](https://www.mathnet.ru/php/archive.phtml?jrnid=sm&option_lang=eng&paperid=9740&wshow=paper): the compact20-vertex case is exactly the regular right-angled dodecahedron. Checked [Roeder–Hubbard–Dunbar, Theorem1.4, printed p.832](https://www.numdam.org/article/AIF_2007__57_3_825_0.pdf): compact nonobtuse marked realizations with more than four faces are unique up to isometry. All uses here satisfy that face-count hypothesis; the uniqueness applies to the prescribed combinatorial symmetries. The actual-wall intersection and simplicity scope is consistent with [Potyagailo–Vinberg, pp.63–64 and Proposition4.1](https://ems.press/content/serial-article-files/43032). Standard convex duality, the hyperbolic polygon identities, and orbifold Gauss–Bonnet remain explicit mathematical inputs, not Lean-verified facts.

Zero-excess recognition requires **every** facet to be a regular dodecahedron. Matching a whole facet fixes the perpendicular neighboring hyperplane and its inward side; the whole regular neighboring dodecahedron is then fixed. Connectedness gives both inclusions of facet sets in the comparison120-cell. The proof does not propagate from a lone regular facet. The paper neither assumes nor proves `f+W>=120`, universal120-cell minimality, a second-volume gap, or decomposition into120-cells. Later C26/C28 or uniform-neighbor refinements are not silently imported into this held manuscript.

## Supporting-source bindings

All paths below are in `results/`; each source was read directly. Their historical candidate/status wording is not used as evidence for the audit conclusion.

|Source|SHA256|
|---|---|
|`compact_ra4_perpendicular_facet_capacity_research_20261018.md`|`8e6db87b001b59a16c1b9b55f03fa9373630bcfb20c1697796f4fd88fbaff8ff`|
|`compact_ra4_neighbor_excess_and_zero_excess_rigidity_20261019.md`|`9c670a7b3c25dda48997c3f70a6d1be9be88bcdefb475b7710f5011068aebcbe`|
|`compact_ra4_positive_facet_link_gap_fifteen_20261019.md`|`8201b2f5cc7da5988b6d164b258da8e14981cf94c6b5e980a5a87e9deb70dcd7`|
|`compact_ra4_fifteen_excess_equality_exclusion_research_20261019.md`|`237d9ccc7c91062b268124fc78c4fa82dbcf4bbf6b0b6b32c8eae013c713ff87`|
|`compact_ra4_sixteen_excess_equality_exclusion_research_20261019.md`|`8caf24e085b204375e3a36b108513e55d8fb0cc66559b44f3e3caf1bd17977e8`|
|`compact_ra4_degree_two_excess_twenty_20261019.md`|`05afd260aa2201b8385f2b1aa397e9896bd6d6e590ae2fe5300d5de9a1546a04`|

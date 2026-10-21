# Compact RA4 primary-paper consolidation

19 October2026. Source preparation and static checks only. No scientific
job, Lean, TeX compilation or PDF rendering ran in this lane.

The revised primary source for this paper family is
`manuscripts/compact_ra4_volume_constraints.tex`:

* SHA256 `afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f`;
* 53,508bytes,1188lines,10sections,4inlineTikZ figures;
* static label/reference, citation, environment, brace, dollar and control
  character checks passed;37unique labels;
* at the author's handoff,no PDF had been built from this revision.
  The coordinator subsequently completed its bounded14-page build,
  zero warnings,and all four figure-page visual reviews. See
  `results/compact_ra4_volume_constraints_W17_build_20261019.json` and
  `results/compact_ra4_w17_figure_visual_review_20261019.json`.

This supersedes the old capacity note as the PRIMARY source in the same
paper family. It does not create a second concurrent primary paper. The
old source and PDF remain the frozen predecessor, not an updated PDF for
the successor results. Root owns catalogue/correspondence/status changes
and all execution scheduling.

## Preservation

The predecessor remains byte-identical to its existing freeze
`results/compact_ra4_facet_capacity_20261018_freeze.json`:

| Frozen predecessor | SHA256 |
|---|---|
| `manuscripts/compact_ra4_facet_capacity.tex` | `2c087727ed59b473d9c25154669af55d3a7d6da6af0c47c8b03d71867dce5928` |
| `manuscripts/build/compact_ra4_facet_capacity.pdf` | `957e1f5f710b2d4333e89c432bb9912e6a59ce5835743e57be7c900ed433492c` |

Both hashes were checked by streaming after source preparation. The old
freeze, supporting result sources, reviews and reproduction document were
not edited. The new source was created at a new path, reusing the complete
old geometric/counting proofs and both old TikZ figures.

Before the authorized W>=16 revision, the complete W>=15 primary source
at SHA256 `9fe28180a9097d8ee99c1898dc64e0716abf55c9c9f8cb950a6db59576d47738`
was preserved in
`manuscripts/source_archives/compact_ra4_volume_constraints_W15_20261019.tex.gz`.
The archive is14,344bytes with SHA256
`fda18e557af25e4a3a3d04f1678fb03188d2f0373680884d13f00a71368cb92b`.
The decompressed bytes were verified equal to the full prior source
before editing. No predecessor PDF was replaced.

Before the authorized W>=17 revision, both the complete W>=16 source
and its actually built13-page PDF were archived and their gzip roundtrip
bytes checked exactly:

| Archived W>=16 evidence | Original SHA256 | Gzip archive SHA256 |
|---|---|---|
| `compact_ra4_volume_constraints_W16_20261019.tex.gz` (16,965bytes) | `20e5b4a3aff508520b7db3c4ae1a8904c0aa2d3299b4164a657bb6432433c49e` | `cbfcd5d076ab66a2361fd2ac20e795be0725cc621d34d2bcb33fb5a150a5a54e` |
| `compact_ra4_volume_constraints_W16_20261019.pdf.gz` (389,731bytes) | `a43b94282d047a14286ea9e7deed3e96bd312c2c50131a8cd1f1cd4f1410eb9f` | `1f59c03988a7b1f2ccae83e8268bd137cf44bc4c52bf89d589361153fcfc037b` |

Both archives are under `manuscripts/source_archives/`. The prior PDF
at the normal build path was not overwritten during preparation. It was
subsequently replaced by the reviewed W>=17 build after these archives
were verified. The archived13-page PDF remains W>=16 evidence. Its root build/visual-review
ledgers remain separate unchanged records. The frozen7-page predecessor
and every reviewed result packet also remain unchanged.

## Proof map and exact scope

| Result/interface | Source location/label | Evidence/dependency |
|---|---|---|
| Acute normal gap, actual orthogonal faces, compact simplicity | Section2, `lem:gap` | Full inherited projection/null proof; named acute-chamber background |
| Common-neighbor clique, perpendicular-edge opposite, actual incidences/distinctness | Section2, `lem:footprint`, `lem:opposite` | Complete old proofs retained |
|3D capacity2k+2 and Loebell equality face lattice | Section3, `lem:polygon` and following paragraph | Complete old triangulated-sphere ear proof retained |
|4D capacity f>=2V+5 | Section4, `thm:capacity` | Complete singleton-footprint proof retained |
| W=v-5f, sum d=2W, exact volume identity | Section5, `eq:volume`, `eq:facet-excess` | Euler, corner counts and explicitly named orbifold Gauss--Bonnet |
| Neighboring excess W>=ceil(3(V-20)/4) | Section5, `thm:neighbor-excess` | Distinct-neighbor excess summed in exact2W identity |
| Current local cutoffs46/46/96 | Section5, `cor:local-size` | Exact120/256 thresholds, parity and integer endpoints |
| W=0 identifies whole120-cell | Section6, `prop:zero` | Explicit3D12-face metric classification; full whole-facet propagation and connectedness |
| Degree2positive-facet structure and W>=5k-15 | Section7, `lem:small-degree` | Loebell equality; named compact3D metric uniqueness; regular-base length; distinct positive neighbors |
| Unique longest base edge of each L6 lateral pentagon | Section7, `lem:loebell-six` | Explicit p_i/q_i reflection plus metric uniqueness; own Lorentz derivation of pentagon identity and exact cosh values2,sqrt2,sqrt3 |
| Degree2positive facet forces W>=20 | Section7, `prop:degree-two-twenty` | Full actual-third-parent incidence; baseline30 plus at least5 extra units from each disjoint lateral row/base |
| Three-regular-pentagon-edge closure | Section7, `lem:three-edges` | Oriented chain fixes two orthogonal ACTUAL closing walls; acute projection proves their actual vertex |
| Eligibility and min degreeH>=5 whenever min degreeG>=3 | Section7, `lem:eligibility`, `lem:four-graph`, `cor:min-degree` | Actual edge carriers and four-neighbor graph proof; G/H distinction explicit |
| All<=9positive-facet graph cases | Section8, `lem:nine`, `prop:ten-positive` | Full triangle/bipartite/C5 independent-set argument, isolated from global W chronology |
| Pentagon/hexagon metric consumer lemmas | Section9, `lem:pentagon-two`, `lem:hexagon-four` | Unique common perpendicular; explicit Lorentz derivation of both hexagon equations and exceptional lengths2L |
| Local exclusion of d=3 facet with five positive neighbors | Section9, `lem:five-neighbor-exclusion` | Full universal/triangle/forest/C5 marked-link proof and actual pentagonal consumers; no assumption on excesses of other facets |
| Ten-vertex complement with at most two degree4vertices has independent5 | Section9, `lem:ten-complement` | Delete-a-high-vertex reductions; explicit common-neighbor t=2,3,4 cases and S-edge capacity proof |
| W=0 or W>=17 | `thm:gap` | ONE final proof for every0<W<=16, using degree2>=20 and shared local/graph lemmas |

Supporting packets used:

* Frozen capacity proof/review and old manuscript.
* `results/compact_ra4_neighbor_excess_and_zero_excess_rigidity_20261019.md`,
  SHA256 `9c670a7b3c25dda48997c3f70a6d1be9be88bcdefb475b7710f5011068aebcbe`,
  with `review_compact_ra4_neighbor_excess_root_20261019.md`.
* `results/compact_ra4_positive_facet_link_gap_fifteen_20261019.md`,
  SHA256 `8201b2f5cc7da5988b6d164b258da8e14981cf94c6b5e980a5a87e9deb70dcd7`,
  with `review_compact_ra4_positive_facet_gap_fifteen_root_20261019.md`.
* `results/compact_ra4_fifteen_excess_equality_exclusion_research_20261019.md`,
  SHA256 `237d9ccc7c91062b268124fc78c4fa82dbcf4bbf6b0b6b32c8eae013c713ff87`,
  with root and finite57 full written reviews. These reviewed result
  packets remain unchanged; the present source is their human translation.
* `results/compact_ra4_sixteen_excess_equality_exclusion_research_20261019.md`,
  SHA256 `8caf24e085b204375e3a36b108513e55d8fb0cc66559b44f3e3caf1bd17977e8`,
  with root and finite57 full written PASS reviews.
* `results/compact_ra4_degree_two_excess_twenty_20261019.md`,
  SHA256 `05afd260aa2201b8385f2b1aa397e9896bd6d6e590ae2fe5300d5de9a1546a04`,
  reviewed by root before integration. All the needed metric and counting
  arguments are now in the paper itself.

The intermediate W5/W9 reports are preserved research support, not extra
primary manuscripts or chronology reproduced in the article. The final
argument includes the complete degree-two and polygon lemmas it needs;
it does not rely on an omitted W5 metric Gram computation.

## External inputs and illustrations

The article distinguishes its own capacity/graph arguments from the
following named inputs:

* the acute-chamber wall property and elementary Lorentz/convex-polytope
  facts, cited through Potyagailo--Vinberg and the standard geometry text;
* compact3D12-face metric dodecahedron classification, explicitly cited
  to Alexandrov--Bogachev--Vesnin--Egorov, Section1.2;
* compact nonobtuse3D Andreev uniqueness, explicitly cited to
  Roeder--Hubbard--Dunbar, Theorem1.4;
* reflection-orbifold Euler and four-dimensional Chern--Gauss--Bonnet.

Primary-literature comparison is recorded separately in
`results/compact_ra4_primary_literature_comparison_20261019.md`.
The source now credits Kellerhals--Perren's prior perpendicular-opposite
mechanism and Inoue's prior3D excess identity. The published volume
reference was corrected from Proposition3.5(d) to3.2(d). No worldwide
novelty claim follows from the limited comparison.

There is no four-dimensional minimum/classification input. All figures
are source-native vector diagrams:

1. inherited dual ear/closing argument;
2. inherited singleton-versus-two-incidence capacity count;
3. three-edge pentagon closure by two actual orthogonal walls, and a
   forbidden zero-carrier block;
4. new independent-five constructions for the nine-vertex graph lemma.

Captions distinguish incidence/chain diagrams from metric hyperbolic
models. They illustrate arguments rather than replace proofs.

## Claims deliberately left open

The paper does not prove f+W>=120, universal120-cell volume minimality,
a double-volume gap/equality classification,
a120-cell tessellation or garland theorem, a new example, or a new
dimension ceiling. The new source has no compiled-PDF or end-to-end
geometric Lean claim. Written mathematical reviews, future build ledgers
and formal interfaces are separate verification layers.

Final source clarification: the right panel of Figure4 now includes the
U--V edge joining its two outside nonisolated vertices, so the depicted
subcubic trianglefree component is itself induced2K2-free. The highlighted
independent five is unchanged. Static checks were repeated after this
one-edge illustration correction; no scientific proof or claim changed.

Current revision: the same four illustrations were retained; the polygon
chain panel now shows three fixed edges and its two prescribed closing
walls. No new figure claims a Euclidean metric model of hyperbolic right
angles. The W>=17 gap is integrated with full case arguments, while
the local46/46/96 facet cutoffs remain unchanged. No proof/build job ran.

## Explicit revision map from20e to the current source

* The old degree-two preliminary estimate is retained as structural
  Lemma`small-degree`, followed by the L6 metric lemma and universal
  degree-two>=20 proposition. No W=16-only slack argument is repeated.
* The earlier staged W>=15 proposition became the reusable statement
  that min degreeG>=3 forces at least ten positive facets.
* The full five-neighbor marked-link proof became a LOCAL lemma for an
  excess3 facet; all universal/triangle/forest/C5 cases remain written.
  The earlier global W=15 reduction is not repeated.
* The new ten-vertex complement lemma handles both at most32-unit cases.
  A single final proof contradicts every0<W<=16, rather than separate
  successive W15/W16 endpoint proofs.
* The same four vector figures remain. Following root's visual review
  of the prior PDF, the lower explanatory text in both Figure4 panels
  was shifted downward by0.35cm (about3mm after scaling), increasing
  its separation from the s_i labels. Mathematical diagram edges are
  unchanged, including the U--V correction.
* Literature credit, complete capacity/neighbor-excess/zero-recognition
  proofs, explicit metric inputs, and all remaining-open limitations
  remain present. No source/PDF build or new formal claim is inferred
  from these static preparation checks.

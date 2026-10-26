# Independent source review of the right-angled and volume Sphinx chapters

26 October 2026. Reviewer `/root/finite57_completion_audit`.

Reviewed the complete new chapters:

- `docs/site/right_angled.rst`, 13,725 bytes,
  SHA256 `558c132c19a9af62217afc80f264a810b8cb0d76dff7c357929052ed48f971b3`.
- `docs/site/volume_problems.rst`, 22,972 bytes,
  SHA256 `cd5ac8103a5cb4850a7b168ada9164184f122ab50fc8941a53a0ed4c8517ddf5`.

The comparison sources are the accepted sharp-dimension source `0c29a786…`,
ideal-volume source `e0c2ca24…` and compact-volume source `afe80767…`, together with their
actual proof/reproduction maps and conditional Lean source interfaces.
The three accepted mathematical TeX files retain their original bytes.

## Sharp dimensions

The chapter correctly derives the projection norm, selected-wall product,
future-sheet sign and strict feasibility against every other original wall.
Maximum-height selection and orthogonal ridge incidence have their distinct
roles. The height recurrence, strict finite-start bound and radial-code
product inequality agree with the paper. Actual face inheritance is confined
to the right-angled setting and faces of dimension at least three.

The cubic polynomial, its rank-two and trace-free rank-three contractions,
all four sum-of-squares coefficients and the ratio 27 are correct. The first
inequality gives 27; the strict endpoint requires the separate four-vector
contradiction and gives 26. The universal six-dimensional 27-facet input is
correctly distinguished from the dimension-eight construction used for
sharpness. The entirely ideal argument uses Musin for cardinality, metric
24-cell equality for recognition and actual opposite normal pairs for the
height contradiction. It assumes no uniqueness of a kissing configuration.

The exact workers are described as fixed-coordinate/universal-polynomial
checks, not spherical-code or arbitrary-polytope enumeration. The four Lean
rows retain the real-instance, geometric and contraction parameters of their
sources. Their 12, 10, 8 and 5 declaration counts and the separate supplement
scope agree with the preserved records. The chapter does not enlarge a
historical core's manuscript-source binding.

## Ideal volume and restricted cutting

All ideal incidence, defect, fiber and volume formulas agree with the paper.
The opposite-cusp count retains distinct actual facets and includes both
nontriangular neighbors and nonsquare-cusp opposites. The complete minimum-
defect split yields the stated gap with the divisibility condition retained.
The classification uses the specified small three-dimensional facet cases;
it is not presented as a four-dimensional census.

The equality spectrum is established before the Oct/D hypothesis is used.
The seven-positive-facet branch, Figure 7 obstruction and cuboctahedral metric
identification retain their actual embedding and cusp-shape conditions.
The separate chain theorem requires every facet to be isometric to an
octahedron or its triangular-face double. Its surface counts, Euler bound,
same-component cusp injection, six-cusp graph, cube equations, Lorentz span,
actual section and outside-wall argument agree with the written proof.

Cut termination, direct all-octahedral recognition, opposite intermediate
cuts and extension of gluing symmetries have the correct scopes. The cell
count, volume and f-vector formulas are unchanged. The chapter explicitly
does not classify all ideal four-polytopes as chains or assert existence of
every chain length. It attributes the minimum separately from the internal
second-volume theorem. All fixed enumeration domains and recorded outputs
agree with the source; they are not credited with universal cutting or
recognition. ManuscriptIdeal4 remains a conditional composition interface.

## Compact capacity and excess

The compact normalization of W differs correctly from the ideal one. The
incidence identities, orbifold Euler formula, volume expression, facet
excess and neighboring-excess inequality agree with the source. The
singleton-extra argument is retained in the capacity count. The comparison
thresholds 46 and 96 are necessary local cutoffs, with no assertion of their
realization or further exclusion.

The zero-excess argument requires every facet to be a metric dodecahedron;
it does not propagate from only one dodecahedral facet. In the positive-
excess proof, weighted nonpentagonal graph G and full positive adjacency
graph H remain distinct. The degree-two/Löbell argument, N=10 reduction,
actual five-neighbor metric exclusion and ten-vertex complement lemma are
used in the correct order. No induced-square property of H is silently
assigned to G. The stated `C_5[K_2]` counterexample preserves the reason that
weaker abstract graph conditions alone do not finish the proof.

The result remains `W=0` or `W≥17`, with `f+W≥120` unproved. The chapter claims neither a
120-cell volume minimum, a compact second-volume theorem, a unique double
nor a general 120-cell decomposition. It retains the published metric and
marked-uniqueness inputs and the explicit absence of a complete geometric
Lean verification of this argument.

## Required local presentation repairs

No substantive mathematical or quantified-scope blocker was found. Three
local repairs were reported to the coordinating author before publication:

1. In the height display of `right_angled.rst`, replace the literal
   `}quad(j\geq1)` by `}\quad(j\geq1)`. The missing backslash otherwise
   prints mathematical letters where spacing is intended.
2. Replace “eight actually printed axiom audits” for ManuscriptIdeal4 by
   “eight declarations with actual printed axiom audits” or equivalent.
   The record has one separate actual audit ledger covering eight declarations;
   it is not eight separate audit executions.
3. Replace the sentence fragment beginning “With N positive facets and E
   unweighted edges. Then” by a complete definition of N and E.

These repairs change no proof or execution claim. This review is source-only;
it neither compiles Sphinx nor certifies the rendered figures or pages. No
worker was imported, no mathematical check/test/build/PDF process ran and no
accepted manuscript, ledger or Git state was changed by the reviewer.

## Confirmation of the three local repairs

The author applied precisely the requested corrections. The current chapter
hashes are `1ada5d3946709c1ad01b1d27fe0dbdb2ecd83819860362bbb6d4cc924317dda2`
for `right_angled.rst` and
`d2f3cc9173d953609294e881408322d01386cb0f1599fe8f6abf5c3fc01a446e`
for `volume_problems.rst`. Independently reversing the missing-backslash
repair and the two wording repairs restores both reviewed predecessor hashes
byte for byte. No other source delta is present. The corrected chapters have
no remaining mathematical or quantified-scope blocker from this review;
rendering and deployment remain separate coordinating gates.

## Figure captions and related explanatory chapters

The author then corrected only the two figure alt/caption blocks to describe
the actual SVG flowcharts. The right-angled figure distinguishes the shared
actual projection from subsequent logical arrows. The volume figure
distinguishes its upper proof boxes from the bottom cell-adjacency row and
retains the restricted Oct/D hypothesis. The current hashes are
`3322501fbcf56821bc279bb48f4d914a7050e199c5747a774efa69be0015e1c4`
and `251369fe86d33b2cd20b18d182cb79c11090fd6bb807e0c5cf95b62d60b108b1`.
Independently restoring those two former caption blocks restores the preceding
repaired hashes exactly. No mathematical text changed.

I also read the complete root-authored `formal_logic.rst`, `examples.rst`
and `sources.rst`. The former's Lean version and displayed
`global_dimension_fifteen` theorem agree with the actual toolchain and source.
Its declaration counts, parameter boundaries, distinct compilation/audit
stages and cache interpretation agree with the selected receipts. It does
not infer actual geometry or field/completeness premises from axiom prints.

The examples chapter retains the exact eleven-wall H8 reconstruction and
its separate endpoint criterion, the source-audited K3 status, and the seven
known Enriques types with their published completeness boundary. The P6,10
lineage discussion is confined to its pinned source audit. A finite Gram or
a partially enumerated chamber is not promoted to a new example or complete
root search. The stated fixed counts agree with the saved artifacts read
as metadata; they were not recomputed for this review.

The source chapter agrees with the external-input register and the nine
byte-preserved predecessor provenance READMEs. It retains the distinction
between full lattice and period-preserving groups, covectors and roots,
the OSCAR object ordering and terminal-newline adjustment, and original
third-party licensing/retrieval boundaries. It records local historical
versions without promising unchanged upstream bytes or new redistribution
rights. No external collection was redownloaded by this review.

Across all eight new mathematical/source chapters, exact comparison of
download targets with the scientific bulk manifest found one bulk target:
the hostile Enriques audit JSON in `examples.rst`. Root replaced that
download directive by the literal repository path and its availability
qualification. The other inspected local download targets exist and are
outside that exact bulk manifest. No source-only review authorizes uploading
the bulk payload through Sphinx or any other route.

The two finite-payment/shared-budget SVG sources have the correct logical
scope: one corner rule branches into finite and ideal estimates, and one
original compact norm budget supports all component costs. The bottom notes
retain the independent higher tail and the arithmetic restrictions without
a second budget. Their arrows depict proof dependencies, not realized
polytopes. This is diagram-source review, not a rendered visual inspection.

The expected concise formal README rewrite changes a documentation hash in
the dated right-angled assembly557ad65a. I reported this current-file binding
issue to root; the accepted mathematical sources and executions remain
unchanged. The dated assembly must retain its historical identity, with the
new documentation layer bound explicitly in the successor delivery record.

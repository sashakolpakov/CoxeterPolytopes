# Sharp right-angled paper: completion and readability audit

21 October 2026. Manual full-source read; no mathematical worker, Lean,
TeX, optimizer or literature search was run. No publication source or
dependency was changed.

## Bound sources and verdict

|Input|SHA-256|
|---|---|
|`manuscripts/right_angled_sharp_dimension_bounds_publication.tex`|`0d607b48616429ae8ff7cf93f8cd8ead4564d4a2515f74178848c84fc6f86bab`|
|`docs/right_angled_sharp_dimensions_publication_proof_map_20261021.md`|`cdc764a6e94ce0ee3b24ef33b38152881d54d55d8a71017fbe7d0c97ca27157c`|
|`docs/manuscript_verification.md`|`140417647ef6dee6e31b7faa372c580cf39c6d06f69144945ffd488770293cfd`|

**Complete written argument at the stated geometric and published-input
boundaries.** The current text supports finite-volume right-angled
dimension at most eight, including mixed finite/ideal vertices, and
entirely ideal dimension at most four. The named examples establish
sharpness. I found no unfinished argument, missing case or false
assumed-versus-derived assertion. The two exposition suggestions below
are nonblocking; neither repairs a mathematical gap.

This is a fresh read of all 378 source lines and the current map,
including the projection/height assumption map and the existing
primary-source audits. I previously authored the facet/acuteness audit
and the predecessor final read. This is not external refereeing or a new
independent retrieval of those published sources.

## Mathematical completion

The proof is organized around one actual-facet construction. The
all-wall implication is explicitly attributed to Vinberg (source
93–102), rather than inferred from ridge angles without justification.
The projection proof (115–147) checks future unit norm, incidence with
the chosen wall, and a strict inequality at **every** other actual wall.
Thus the projected point is in the relative interior of the actual
facet. Its intrinsic walls are exactly the orthogonal ridge neighbors;
their unit scaling and height identities are retained. A maximum exists
because there are finitely many facets. Ideal nonsimplicity does not
invalidate this argument: the selected points and intrinsic ridges have
finite relative-interior points.

The recurrence (149–160) retains strictness from the finite initial
maximum. After three projections its bound gives radial products
strictly below one quarter. Distinctness of radial directions is proved
explicitly (196–215), so code cardinality equals the number of actual
facets. The tensor contractions, nonnegative square sum and strict
27-point exclusion are displayed in full (219–281). The contradiction
uses the 27-facet minimum for the actual finite-volume six-face. Face
inheritance also supplies the reduction from any higher dimension to
nine; no unspecified classification of formal Gram configurations is
used.

In the ideal branch (299–344), one projection gives an actual entirely
ideal four-dimensional facet with all heights below one. Musin supplies
only the code cardinality bound; Kolpakov supplies both the opposite
facet bound and **hyperbolic metric** equality. The opposite-normal
calculation then excludes that metric 24-cell at the projected point.
Kissing-configuration uniqueness is never substituted for metric
rigidity. Five-face inheritance covers all higher dimensions. No descent
needs the false assertion that an ideal polygon has right angles.

## Exact dependency and verification boundaries

The general exclusion's external numerical input is the universal
six-dimensional 27-facet minimum. The text accurately cites Dufour's
consequence **via Alexandrov**. The existing primary audit directly
checked Alexandrov's all-vertex-type statement and its derivation;
Dufour's original full text was not obtained. The ideal exclusion adds
Musin's main cardinality theorem and Kolpakov Theorem 4's facet minimum
and isometry equality. Their original sources and the sharp examples
were checked in the existing ideal-input audit. Standard Lorentz/Klein
geometry, all-wall acuteness and face inheritance are geometric inputs,
not outputs of the finite scripts. No unresolved source issue requires
a new literature search for this audit.

The map correctly distinguishes universal exact polynomial expansion
from sampling, and the displayed 24-cell model check from complete
vertex enumeration or arbitrary-polytope classification. Neither worker
is an unpublished finite census on which the dimension theorem depends.
The written proof displays the universal algebra and the equality
calculation itself. Real positivity, actual face incidence, geometric
interpretation and the published theorems remain outside those workers.

The three selected Lean roots check conditional internal logic. The
projection/height root derives its coordinate identities and scalar
descent from primitive coordinate, wall-sign and normalization inputs;
it does not instantiate actual hyperbolic polytopes, face inheritance or
published theorems. The dimension/count root also retains its explicit
geometric interfaces. The map does not call these end-to-end geometric
Lean proofs. The 12+10+8 actual prints are repeated core outputs, not 30
new declarations added to the seven-kernel/62-print publication receipt.

The current reproduction command selects six serial jobs with PDF.
Recorded four- and six-job core runs belong to the preserved predecessor
edition. The current publication successor has the fresh five-page,
zero-warning bounded build recorded separately; the map explicitly says
that its terminology delta did not rerun exact workers or extend Lean
coverage. The unchanged figure's earlier visual record is not presented
as a fresh all-page inspection. No fresh combined successor core PASS is
claimed by this audit.

## Concrete proposed edits, before any source mutation

1. **Optional short geometric explanation, source 99–102.** The present
   inheritance statement is correct and sufficient as a named standard
   foundation, but it carries both dimension reductions. A reader would
   benefit from replacing its first two sentences with:

   > Faces of dimension at least three inherit right angles: every
   > intrinsic ridge has a finite relative-interior point, where the
   > incident walls are orthogonal. In the Klein closure their vertices
   > are vertices of the original polytope, so entire ideality is
   > inherited; triangulation into finite or ideal simplices gives finite
   > volume.

   Retain the following warning about two-faces. This explains the
   already-used foundation; it adds no new hypothesis or result.

2. **Optional notation clarification, source 151–152.** After defining
   `A_0`, write: “Let `A_j` be the maximum height after `j` projections.”
   Then continue with the displayed recurrence. The present meaning is
   recoverable from context and the figure, but not expressly defined.

3. **Map chronology, lines 125–132.** Label the first linked record
   “historical four-job core” and use past tense. Later paragraphs already
   explain its predecessor scope accurately; this small edit would make
   that scope visible immediately under the verification heading.

4. **Reproduction-guide pointer.** In
   `docs/manuscript_verification.md`, change the sharpRA proof-map link
   from the preserved `right_angled_sharp_dimensions_proof_map_20261020.md`
   to the current publication map. Describe the recorded six-job run as
   the predecessor core run and retain the separately recorded successor
   PDF build. This is documentation maintenance, not a request to rerun
   successful jobs or rebind old receipts.

No broader expansion, extra appendix, new manuscript, new finite-data
claim or new formalization is needed to complete this paper. These are
proposals only; no TeX, map, catalogue, registry or guide was edited.

## Research freeze

The arithmetic research packet
`results/arithmetic_dimension19_forest_frontier_and_actual_terminal_exclusions_20261021.md`
is frozen at SHA-256
`c5992464c6d81f38c22b75b8259518f52c69806105e7e364acc325c1147f2211`.
It remains an unreviewed candidate, receives no promotion from this
publication audit, and was not extended. Current accepted bounds and
all six manuscript/PDF snapshots remain unchanged.

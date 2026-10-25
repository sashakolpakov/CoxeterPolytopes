# Six manuscripts: audit boundary and consolidation map

24 October 2026. This is an implementation map for the new audit and
consolidation task. It records inspection of theorem statements, section
structure, proof maps, accepted composition interfaces and execution
receipts. **It is not a fresh mathematical acceptance, a completed audit,
a revised manuscript or a new verification run.** The six delivered
TeX/PDF pairs have not been edited by this preparation. Their byte
bindings and the sources inspected here are recorded in the
[source snapshot](../results/six_manuscript_consolidation_source_snapshot_20261024.json).

The first audit priorities are the general finite-volume bound and the
sharp right-angled bounds. Existing reviews are evidence to inspect, not
premises asserting that a new audit must pass. A counterexample, missing
case, unaffordable simultaneous payment or incorrect external input takes
precedence over an old acceptance label.

## The six results and their revision boundaries

| Existing paper family | Delivered edition | Accepted result to be audited for consolidation | Material outside the revised mathematical argument |
|---|---|---|---|
| [General finite volume](../manuscripts/finite_volume57.tex) | 48 pages; dimension at most 57 | The same theorem, excluding ambient dimension 58 and the separate higher tail | The unfinished ambient-57 campaign toward 56, including literal-heavy and zero-arrival port refinements, even where individual lemmas have been accepted |
| [Sharp right-angled dimensions](../manuscripts/right_angled_sharp_dimension_bounds_publication.tex) | 5 pages; finite volume at most 8, entirely ideal at most 4 | The same two sharp bounds, subject to the fresh geometry and input audits | Older mixed10 and ideal5/6 approaches; no new general-dimensional search |
| [General compact](../manuscripts/compact_farthest_vertex_bounds.tex) | 36 pages; dimension at most 23 | Dimension at most 17, from the accepted complete composition and its dependency closure | The exact17 five-row campaign toward 16, including accepted local row exclusions; spectral and construction research |
| [Compact arithmetic](../manuscripts/compact_arithmetic_dimension_bound.tex) | 21 pages; dimension at most 20 | Dimension at most 15, from the exact17 and exact16 exclusions with a proved global tail | The exact15 ten-group campaign toward 14; parked nine-field19 refinements; unfinished A14/rank-one work |
| [Ideal RA4 volume and structure](../manuscripts/ideal_ra4_volume_and_structure_publication.tex) | 17 pages; published minimum, second volume with unique double, restricted chains | The same coherent sequence, with an independent audit of the structural steps | General garland/decomposition conjectures and unnecessary antiprism extensions |
| [Compact RA4 constraints](../manuscripts/compact_ra4_volume_constraints.tex) | 14 pages; capacity, excess, zero-excess recognition, positive gap 17 | Retain these results; audit the accepted C26/C28 successor before adding its unconditional facet and volume bounds | Conditional inball, central-collar and unresolved separating-homology refinements unless needed for a retained theorem |

The accepted compact17 and arithmetic15 results are not sharpness claims.
There is no new polytope example. General finite-volume, general compact
and right-angled arguments have no arithmetic hypothesis. The
quadratic-form model of an arithmetic reflection lattice is automatic;
“first type” must not be presented as a further restriction.

## Common writing and proof interface

Each paper should start with its geometric idea, then identify the first
place where different cases need different methods. Present each common
lemma once. A chronological sequence of intermediate dimension bounds
would obscure both the method and its dependency closure.

Each proof-relevant computation needs five mathematical statements:
the admissible domain and its geometric coverage; the algorithm; a
termination argument; soundness and completeness of every transition,
pruning rule and quotient; and the exact output used in the theorem.
The [computational proof contract](computational_proof_contract.md)
provides the model. Code and certificates belong in the reproducibility
supplement. A list of filenames or one maximizing witness is not a
proof of exhaustive coverage or of a universal upper bound.

For each retained theorem the dependency map must end at one of four
explicit boundaries: a stated published theorem; a proved internal
geometric/analytic lemma with a TeX proof home; a stated finite computed
proposition with an algorithm and execution evidence; or a conditional
Lean theorem with all ordinary parameters listed. A receipt can bind
those objects but cannot replace the mathematical statement. A prepared
module, a successful build or a source hash adds no mathematical or
execution credit.

## 1. General finite volume: keep the theorem 57

The existing argument already has the right separation. Sections 1--2
define actual face corners and the payment comparison. The local rule
then branches into finite spherical forests and ideal affine links;
the branches reunite in the shape-curvature contradiction. The bounded
and unbounded higher-dimensional tails are separate proof obligations.

The revision should preserve that architecture:

1. State a single simultaneous nonnegative payment rule and every
   priority/transfer operation. Prove triangle payment at least 2 and
   quadrilateral payment at least 1 for actual faces.
2. Define the full finite alphabet, context state, forest/fiber domain
   and recurrence. Prove that every actual rank58 vertex enters it and
   that the computed charge is at most 377.
3. Define the affine-rank57 cusp domain, nonsimple product-link fibers,
   factor envelopes and concentration reduction. Prove actual charge
   at most curvature minus 36.
4. Derive the correction 4205 and conclude the exact dimension58
   contradiction. Keep the bounded tail through 163 and the unbounded
   narrow-facet argument, including the compact-facet alternative.

The exact sources and obligations are in the
[lemma map](finite_volume57_lemma_map.md) and
[reproduction guide](finite_volume57_reproduce.md). The audit must check
the whole-component interpretation, all additional white roots,
simultaneous rather than separately optimized charges, actual receiving
incidences and the absence of unproved Coxeter-face inheritance.

Runtime evidence has a precise limit. All 25 bounded B13/frontier
replacement units actually passed and are integrated in the primary
catalogue. The full selection has 136 scripts, nine Lean roots and one
PDF job. The 111 other command matches remain inherited execution
records. There is no combined 146-job runtime PASS or complete transitive
regeneration. The audit should decide which specific inherited results
need stronger evidence before scheduling a new job; it should not rerun
successful modules simply to obtain a new date.

The selected Lean kernels cover conditional payment composition,
concentration descent, residual Bellman logic, scalar tail arithmetic
and shape-counting steps. Actual geometry, interpretation of the
residual alphabet and contexts, analytic convex-order comparison,
terminal data and geometric tails remain application inputs. Review
the parameters in the [assumption register](manuscript_assumptions.md)
and logic maps; printed axioms do not list them.

The ambient-57 work toward 56 uses different ranks and budgets. Its
accepted weighted ideal calculations, literal-heavy financing and
covered-port lemmas are valuable research, but they do not belong in
the proof of the current 57 theorem merely because they are newer.

## 2. Sharp right-angled dimensions: one projection, two branches

Keep the five-page paper focused. The shared argument projects a finite
interior point to an actual farthest facet, checks all original walls,
and decreases the maximum intrinsic facet height. The universal branch
then uses three projections from a hypothetical nine-polytope to a
six-face. Its tensor certificate gives at most 26 facets, contradicting
the stated universal minimum of 27. The entirely ideal branch uses one
projection from dimension five to four, Musin's bound, the metric
24-cell equality case and the opposite-normal obstruction.

The fresh audit must validate actual-face inheritance, finite volume,
entire ideality, strict height inequalities and the equality step of
the tensor bound. It must verify the external statements at their exact
scope: Dufour's six-dimensional facet minimum, Musin's cardinality
bound, Kolpakov's facet minimum with metric equality, and published
examples for sharpness. Neither a construction with 27 facets nor a
volume minimum substitutes for the required facet minimum.

The [proof map](right_angled_sharp_dimensions_publication_proof_map_20261021.md)
identifies two exact workers and three selected conditional Lean roots.
The polynomial worker checks universal algebra; the normal-model worker
checks a fixed 24-cell model. They do not prove actual face geometry or
the published theorems. The expanded historical core contains six
sequential jobs including a temporary PDF build. Its receipt must not be
silently rebound to an edited source. Successful unchanged Lean modules
need no rerun. Any new correction must have a new source-bound written
review, relevant changed checks and its own build receipt.

## 3. General compact: replace increment chapters by component arguments

Retain the common original-vector geometry in `sec:geometry`, including
positive projection, reciprocal shadow and all-wall feasible points.
Retain ordinary costs and the precise I/J and typed-seven defect
accounting. Consolidate the large-A arguments into one chain-pressure
section and B/D arguments into one prefix/head/fork section. Each local
statement must retain the rank, spectator, comparison factor and
coefficient-price hypotheses actually proved.

The revised final composition should use a finite necessary-cover table
for the dimensions it excludes, with each row pointing to a complete
component theorem. It must include the higher-dimensional tail and not
infer a global bound from nonexistence in one dimension. The accepted
[compact17 composition](../results/compact_dimension17_global_written_acceptance_20261024.json)
has a 470-line final source and 21 direct bound inputs. Its final
whole-D18 argument belongs after the common prefix lemmas: actual
dotted-head point, four complete initial references, individually funded
fork heights, the local actual-A18 lemma and the last near/far
contradiction. Finite and parallel heads and every other fork family
need their complete proof homes, not a statement that earlier receipts
removed them.

The old global23 proof is part of the current dependency history. Do not
delete a proof merely because its theorem is numerically weaker: first
show that its essential component or higher-tail facts have complete
homes in the consolidated argument. Conversely, remove an obsolete
specialized endpoint argument if a retained stronger lemma proves the
same obligation with all hypotheses satisfied. Preserve its source as
history rather than repeating it as another increment section.

The older compact18 kernel has actual compilation and a separate audit
of five declarations under 23 written premises. The two prepared
[compact17 modules](compact_dimension17_conditional_logic_map_20261024.md)
have static reviews but no actual compile or axiom audit at this
snapshot. Their eventual serial compilation does not formalize the
geometry, calculus, individual prices or complete inventories supplied
as theorem parameters.

The exact17 cover and its later five-row refinement concern a hoped-for
16 theorem. They are unnecessary to proving 17 and should remain outside
the revised TeX.

## 4. Compact arithmetic: a possible shorter route to 15

The accepted [global15 proof](../results/review_arithmetic16_whole_d16_manual_dyadic_exclusion_and_global15_root_acceptance_20261023.md)
combines a global16 tail, the complete exact16 frontier and the manual
whole-D16 dyadic exclusion. Its new final source has eight direct
accepted interfaces. The reduction to actual quadratic fields, the
inert/ramified/split cases, full maximal-order strips, field-square
exceptions and separately priced geometric residues must be written as
mathematical arguments. The prepared D16 worker was not executed and is
not a premise of that proof.

A substantial consolidation opportunity is to use the accepted general
compact17 theorem as the global17 tail. Then prove the complete exact17
arithmetic exclusion, giving 16, and the exact16 exclusion, giving 15.
The 396-line [global16 source](../results/arithmetic17_whole_a17_complete_endpoint_census_and_global16_composition_20261022.md)
already separates its global17 input from its 14-group exact17
exhaustion. This replacement could make the old arithmetic20
touching-star/angle-rule and field-by-field global-tail chapters
unnecessary in the primary paper.

**This deletion is a proposal, not an established dependency reduction.**
First inspect the local exact17 and exact16 source closures for any
remaining field, census or angle-rule input; give each retained input
a proof or exact algorithm home. A named internal theorem may be cited
from the general compact companion, with its full statement and proof
location; it must not be relabelled as published or as a Lean axiom.

The proposed chapter order is: common actual-vector framework cited from
the general companion; arithmetic admissibility, actual-field equality
and cyclic integrality; geometric degree/contact reductions; complete
exact17 forest dispatch and its finite-data propositions; complete
exact16 forest dispatch; the quadratic dyadic obstruction; final 15
composition and reproducibility. Finite data belong where their
mathematical proposition is first used. Avoid chapters for each past
bound 20, 19, 18, 17 and 16.

The [dyadic kernel](quadratic_dyadic_triple_kernel_logic_20261023.md)
has an [actual compilation and separate 14-declaration audit](../results/lean_quadratic_dyadic_triple_verification_20261023/manifest.json).
Its preparation map still has a historical source-preparation header;
add current execution status without erasing that boundary. It proves
discrete valuation/norm implications and a conditional numerical
composition. Actual fields, valuation interpretation, maximal-order
strips, square tests, geometry and older finite data remain external.
The older exact17/global16 composition kernel has its own actual scope.
Any new final15 interface must expose the complete dispatch rather than
assume the final impossibility in a differently named field.

The exact15 ten-group research, original15 A14/rank-one work and parked
field-specific19 refinements are not needed for the 15 theorem.

## Section-level acyclicity check for the proposed arithmetic tail

The following is the declared-interface path found by source inspection;
it is not a full transitive mathematical audit.

| Node used by the general17 proof | Exact scope retained | Scope in the same files that is not imported |
|---|---|---|
| Final whole-D18/global17 source `08c3a04e…` | Its final general composition uses G18 and the exact18 geometric frontier | No arithmetic exact17 exclusion or arithmetic15 theorem |
| Global18 source `9a855860…` | General G19, complete F19 and geometric row exclusions | No arithmetic counterpart of a mixed source |
| Mixed global19/arithmetic17 source `0ac80210…` | Sections 2--3: common geometry and exclusion of general dimension20 | Section 4: arithmetic dimension18 exclusion |
| Mixed global20/arithmetic19 source `c6b4873a…` | General dimension21 exhaustion and exclusion in Sections 3--4 | Its arithmetic dimension20 exclusion |
| Mixed general20/arithmetic18 frontier `415e2289…` | General dimension20 row and geometric A19+I subtraction | Arithmetic whole-B18 subtraction |
| General F19 source `c443832c…` | Mixed `b95245cf…` only for its first, nonarithmetic A-rank13 theorem | Arithmetic A15 packing in the same source |
| General A16 consumer `a707755e…` | `arithmetic_dimension18_fresh_forest_frontier` only for its rank-free geometric Section 2 I barrier and I+T2 estimate | Field conditions, the rank18 window and arithmetic forest list |

Thus the inspected composition interfaces do not require
`arithmetic17 => arithmetic15 => general17`. The revised papers must
preserve these section boundaries. Audit the deeper local inputs before
certifying full acyclicity or removing the old global arithmetic tail.

## 5. Ideal RA4: preserve the order of the three main results

The current theorem sequence is already coherent: the attributed
published minimum; the new second-volume gap and equality spectrum;
then the restricted octahedral/doubled-octahedral chain theorem, applied
to the equality spectrum. The cube-normal/cutting construction is
proved once. Do not assume the chain theorem to obtain the facet
spectrum that is needed to apply it.

The [proof map](ideal_ra4_volume_and_structure_proof_map.md) gives the
dependency split. Keep volume and defect identities, opposite-cusp
injection and complete small-facet classifications before equality.
Then prove closed internal surfaces, the six-cusp/eight-normal step,
one actual cutting hyperplane, disjoint sections, all-octahedral
recognition, opposite intermediate cuts and canonical gluing. Retain
figures that explain these geometric steps. The theorem concerns
existing polytopes with the stated intrinsic facet spectrum; it neither
asserts arbitrary garlands nor existence for every chain length.

External inputs are the explicitly stated Kolpakov incidence/volume and
antiprism results, the nine- and ten-vertex three-dimensional cases and
right-angled rigidity. Fixed double/cuboctahedron/root-model checks
support those models; they do not replace universal injection, surfaces,
cutting or gluing. The historical core ledger binds the 16-page
predecessor, while the 17-page publication successor has a separate
source-bound build. The catalogue lists five script selections and two
conditional Lean roots; a new source-binding assessment is needed
before calling any combined publication-successor replay complete.
Unchanged successful logic modules need not be recompiled.

## 6. Compact RA4: a bounded successor, not a minimum-volume paper

Keep the shared perpendicular-edge geometry, polygon capacity,
neighboring excess, volume identity and zero-excess recognition. The
positive-excess graph argument should prove `W=0` or `W>=17` once,
without repeating the historical W15 and W16 stages.

There is a further accepted unconditional result outside the delivered
TeX: the [C28 successor](compact_ra4_c28_successor_reproduce_20261020.md)
proves `f>=65` and `Vol>=49*pi^2/6`. Its proposed integration point is
after the positive-excess gap. It needs a proof that a dodecahedral
facet exists, the C26 and C28 no-regular-pentagon metric lemmas, and
the existing capacity/volume identities. A positive facet adjacent to
a dodecahedron then has at least 30 vertices and excess at least 5;
this is not an excess bound for every positive facet.

The new classification premise must be stated exactly: the published
fullerene count is two at 28 vertices, from Brinkmann--Goedgebeur--McKay,
Table 2, page 15. A primary-source text check confirms that count and
the definition as cubic plane graphs with only pentagonal and hexagonal
faces. [The Generation of Fullerenes](https://arxiv.org/pdf/1207.7010v2)
also gives zero at 22 vertices and one each at 24 and 26. This text
check does not verify the geometric reduction or the metric argument.
The proof must also justify that the remaining
compact RA3 cases reduce to the Loebell type and these two fullerene
types. The worker checks two fixed dual spheres, symmetries,
propagation and exact rational/Q(sqrt5) identities. It does not
enumerate all fullerenes or prove metric uniqueness. The
[27-file immutable successor packet](../results/compact_ra3_c28_successor_20261020_freeze.json)
binds the complete written/root/outside reviews, unchanged C26 inputs,
generation and separate replay. No new execution is supplied here.

The compact RA4 catalogue currently has no Lean root and no numerical
check for its universal geometry. If the C28 result is retained, add
its exact check and algorithm statement to this same paper's catalogue
entry. Prepare a useful conditional logic interface for the excess
identities, capacity consequences, finite graph dispatch and final
volume lower bound; keep all geometric and classification premises
visible. Do not label a successful TeX build as a proof replay.

The regular-inball theorem remains conditional on an actual inball of
the stated radius. A farthest-wall radius does not provide it. The
central-collar and disconnected-dodecahedron/homology investigations
should remain outside the TeX unless selected as an explicit retained
theorem. Neither `f+W>=120`, universal 120-cell minimality, a second
compact volume or a unique double is proved by the accepted successor.

## Catalogue, documentation and release synchronization

Keep the same six primary IDs. An edited source gets its own source/PDF
and mathematical-review binding; preserve the delivered pair and all
old evidence before replacement. Do not retag inherited execution
records as a fresh combined successor run.

The following observed drift needs correction during integration:

* `docs/site/papers.rst` and `scope.rst` still describe later compact18
  and arithmetic17 results; the current accepted targets are 17 and 15.
* The compact23 proof map and compact20 proof map correctly describe
  their printed editions but contain stale descriptions of the latest
  later bounds. Keep printed and current scopes explicit.
* Catalogue bundle notes/titles still say written compact20/PDF23 and
  arithmetic18/PDF20. Primary manuscript entries bind the delivered
  editions correctly; update target metadata without changing historical
  result entries into claims about those PDFs.
* Prepared logic-map headers and later execution manifests disagree in
  apparent current status. Add a clear current-status pointer and retain
  the exact preparation snapshot bound by older receipts.
* The compact RA4 primary entry needs the retained C26/C28 checks and
  conditional logic only after that successor passes the requested audit.

Then synchronize the six TeX/PDF entries, `manuscripts/README.md`, the
result correspondence, current bounds, assumption register, logic maps,
Sphinx downloads and `formal/README.md`. A strict Sphinx build checks
documentation rendering and links without launching mathematical jobs.
Before a push, inspect the actual Git index through the publication
watchdog and exclude dependencies, objects, caches, intermediates and
generated sites. Intentional final PDFs and verification stdout remain
evidence.

All scientific, Lean and TeX jobs are root-owned serial batches under
the existing shared lock and resource gates. Successful unchanged
modules and checks are not repeated merely for documentation. Each
new or substantively changed module needs a bounded compilation and a
separate actual axiom audit. Remove only owned unused build products
after each batch. The bulk research archives remain unuploaded; this
map neither authorizes an upload nor converts restoration into
mathematical verification.

## Pending reproduction-driver successor

Root is preparing a new audited driver and public bounded proof-batch
entry point while preserving the legacy driver bytes that existing
workers bind. This preparation does not change the historical 146-job
plan or its execution scope. The proposed new plan separates each Lean
compilation from its actual axiom audit. If all nine finite57 roots have
the required audit commands, that would make its new full/core plans
155/20 jobs respectively. These are **conditional proposed counts**,
not an executed plan, a runtime PASS or a claim that all those modules
need recompilation. The generated plan and actual source inspection
must determine the successor's eventual counts and selection.

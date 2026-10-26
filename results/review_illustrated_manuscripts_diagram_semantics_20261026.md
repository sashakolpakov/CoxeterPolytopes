# Independent semantic review of the illustrated manuscript sources

26 October 2026. Reviewer `/root/ra_release_completion_audit`; figure and
proposal author `/root/compact_manuscript_completion_audit`.
**PASS for the exact source proposals below, at the underlying manuscripts'
declared mathematical inputs.** Root acceptance, bounded compilation and
rendered inspection remain separate; no illustrated PDF page count is
credited by this review.

## Exact source boundary and byte restoration

| New illustrated source | Reviewed SHA256 | Unchanged base SHA256 restored by removing the marked blocks |
|---|---|---|
| `manuscripts/finite_volume57_illustrated.tex` | `107c31528f80fb1e0e49381546cb439b1b4180c1c6449b9716994a83bcccea18` | `38186aa7b05251cd8a6b781d751d8eadbd70eba9a49f603f0d4ba9c7eb66ee80` |
| `manuscripts/compact_farthest_vertex_bounds_illustrated.tex` | `0f5ae4c89a9848da41582621f2be65558845890972f99889f7312bfb4429958b` | `6222f5c6cee8b30b2c1019d648bc31cde8bfcddf08f031a0355a3e7137c107ac` |
| `manuscripts/compact_arithmetic_dimension_bound_illustrated.tex` | `1eda39b845e0d6353d2447bc288f49aeb5c44523ff1c4bbf5d2097995e50eab5` | `d0f040b253146f318bafd2391c8809c9d995f53fcd920477f78f42038d424295` |

I read all six complete TikZ figure blocks and both complete package
blocks, together with the relevant actual proof paragraphs and referenced
labels. Independently removing the two finite57 blocks, three compact17
blocks and three arithmetic15 blocks restores the corresponding base
bytes exactly. Every individual marked block also matches its recorded
byte count and SHA256. This is a source/byte comparison, not execution of
the mathematical algorithms.

The two package blocks add only TikZ and `arrows.meta`; finite57 already
loads both. Each figure's small font is scoped inside its TikZ picture.
No proof-body font, equation, hypothesis, source interface or existing
figure was changed. Six new figure labels are unique and do not collide
with base labels. Every reference inside the new blocks resolves to a
label in the corresponding base source.

## The six drawings and their logical scope

| Figure block suffix | Source-semantic assessment |
|---|---|
| `FV57_CORNERS` | The triangle receives total payment at least2, the quadrilateral at least1, and the same actual corner ledger is counted as a face sum or a vertex sum. The displayed D5 upper bound is exactly Lemma `lem:payment-comparison`; larger faces are covered by nonnegativity. Dots do not prescribe equal weights. The caption explicitly distinguishes schematic polygons from metric Coxeter realizations and keeps the separate incidence lower bound separate. |
| `FV57_TAIL` | The branches partition n>=58 into exact58, 59--163 and all n>=164. The bounded family appendix and the all-rank narrow-facet/two-budget appendix have separate proof homes in the actual source. Exclusion of58 is not used as a hereditary higher-dimensional theorem. The caption retains geometric and finite-data premises and says the arrows are logical, not passages to Coxeter facets. |
| `C17_BUDGET` | The left panel is an abstract orthogonal decomposition in the original p-perpendicular Euclidean space. The sum of whole-component energies is the norm budget of the same x. The actual feasible future point gives the displayed pressure from Lemma `lem:pointpressure`; the final square uses one funded reference W. The caption retains individual coefficient debits and all cross terms. W is not asserted to be an actual hyperbolic point, and overlapping pressures are not independently optimized and added. |
| `C17_CONJUNCTION` | Exact18, independently proved own-rank19--29 cases and the published compact29 ceiling cover all n>=18. The actual source's complete-conjunction paragraph supplies the separate 24--29,22--23,20,21 and whole A19/B19/D19 exclusions. Necessary forests are not examples. No arithmetic premise or face-heredity reduction is introduced. The original farthest vertex/vector stays fixed inside each case. |
| `A15_FILTER` | One original-vector budget is filtered by arithmetic structure over the actual field K; there is no second budget. The root-line caveat avoids claiming every physical unit normal itself lies in K. The finite-domain/recorded-output lane is tied in the caption to the four computed propositions and their actual-domain inclusions. Algebraic survivors are expressly not realized polytopes. The separate next figure identifies the manual whole-D16 lane. |
| `A15_DISPATCH` | The independent general compact17 companion is a one-way input. Exact17 and exact16 have separate complete forest dispatches. The four numerical-output premises feed their specified domains and consumers; the last whole-D16 step is explicitly manual, splitting odd head norm from the even inert/ramified/split-at2 cases. No new whole-D16 census or arithmetic face inheritance is asserted. This agrees with `app:arithfinalimplication` and the separate whole-D16 reduction, not the D16+I theorem at ambient17. |

The source records used for placement and provenance are
`results/illustrated_manuscript_source_only_proposal_20261026.json` and
`results/illustrated_manuscript_figure_source_notes_20261026.md`. Their
pending-review/build fields describe the authoring stage; this separate
review does not mutate them into execution receipts.

## Findings and acceptance boundary

No mathematical-content or logical-scope correction is required by this
review of the exact source hashes above. Captions consistently identify
schematic, combinatorial or dependency drawings. None supplies a missing
geometric realization, a complete census, a new bound, a proof of an
external literature input, or new numerical/formal verification credit.
The underlying exact-data and geometric premises retain their existing
scope; byte restoration does not prove them anew.

Actual figure size, float placement, label legibility, overlap and final
page count still require root's monitored builds and rendered inspection.
This source-only clearance does not claim those visual checks. The six
canonical editions and three existing base working sources/PDFs are
preserved. I changed no authoritative manuscript and ran no worker import,
enumeration, test, Lean/TeX compilation or PDF/image renderer.
